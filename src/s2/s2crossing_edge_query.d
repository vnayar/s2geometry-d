// Copyright 2013 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2crossing_edge_query;

import s2.r2point;
import s2.r2rect;
import s2.s2cell_id;
import s2.s2edge_clipping;
import s2.s2edge_crosser;
import s2.s2padded_cell;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.shapeutil.count_edges;
import s2.shapeutil.shape_edge;
import s2.shapeutil.shape_edge_id;
import s2.util.container.btree_map;

import std.algorithm : sort, uniq;
import std.array : array;
import std.exception : enforce;
import std.range : empty;

/**
 * A parameter that controls the reporting of edge intersections.
 *
 *  - CrossingType::INTERIOR reports intersections that occur at a point
 *    interior to both edges (i.e., not at a vertex).
 *
 *  - CrossingType::ALL reports all intersections, even those where two edges
 *    intersect only because they share a common vertex.
 */
enum CrossingType { INTERIOR, ALL }

// For small loops it is faster to use brute force.  The threshold below was
// determined using the benchmarks in the unit test.
private enum int MAX_BRUTE_FORCE_EDGES = 27;

/**
 * S2CrossingEdgeQuery is used to find edges or shapes that are crossed by
 * an edge.  Here is an example showing how to index a set of polylines,
 * and then find the polylines that are crossed by a given edge AB:
 *
 * void Test(const vector<S2Polyline*>& polylines,`
 *           const S2Point& a0, const S2Point &a1) {
 *   MutableS2ShapeIndex index;
 *   for (S2Polyline* polyline : polylines) {
 *     index.Add(absl::make_unique<S2Polyline::Shape>(polyline));
 *   }
 *   S2CrossingEdgeQuery query(&index);
 *   for (const auto& edge : query.GetCrossingEdges(a, b, CrossingType::ALL)) {
 *     CHECK_GE(S2::CrossingSign(a0, a1, edge.v0(), edge.v1()), 0);
 *   }
 * }
 *
 * Note that if you need to query many edges, it is more efficient to declare
 * a single S2CrossingEdgeQuery object and reuse it so that temporary storage
 * does not need to be reallocated each time.
 *
 * If you want to find *all* pairs of crossing edges, use
 * s2shapeutil::VisitCrossingEdgePairs() instead.
 */
class S2CrossingEdgeQuery {
public:

  /// Default constructor; requires Init() to be called.
  this() {
    _iter = new S2ShapeIndex.Iterator();
  }

  /// Convenience constructor that calls Init().
  this(S2ShapeIndex index) {
    this();
    initialize(index);
  }

  inout(S2ShapeIndex) index() inout {
    return _index;
  }

  /// REQUIRES: "index" is not modified after this method is called.
  void initialize(S2ShapeIndex index) {
    _index = index;
    _iter.initialize(index);
  }

  /**
   * Returns all edges that intersect the given query edge (a0,a1) and that
   * have the given CrossingType (ALL or INTERIOR).  Edges are sorted and
   * unique.
   */
  ShapeEdge[] getCrossingEdges(in S2Point a0, in S2Point a1, CrossingType type) {
    ShapeEdge[] edges;
    getCrossingEdges(a0, a1, type, edges);
    return edges;
  }

  /**
   * A specialized version of GetCrossingEdges() that only returns the edges
   * that belong to a particular S2Shape.
   */
  ShapeEdge[] getCrossingEdges(
      in S2Point a0, in S2Point a1, in S2Shape shape, CrossingType type) {
    ShapeEdge[] edges;
    getCrossingEdges(a0, a1, shape, type, edges);
    return edges;
  }

  /**
   * These versions can be more efficient when they are called many times,
   * since they do not require allocating a new vector on each call.
   */
  void getCrossingEdges(
      in S2Point a0, in S2Point a1, CrossingType type, ref ShapeEdge[] edges) {
    edges.length = 0;
    getCandidates(a0, a1, _tmpCandidates);
    int min_sign = (type == CrossingType.ALL) ? 0 : 1;
    auto crosser = new S2CopyingEdgeCrosser(a0, a1);
    int shape_id = -1;
    S2Shape shape = null;
    foreach (ShapeEdgeId candidate; _tmpCandidates) {
      if (candidate.shapeId != shape_id) {
        shape_id = candidate.shapeId;
        shape = _index.shape(shape_id);
      }
      int edge_id = candidate.edgeId;
      S2Shape.Edge b = shape.edge(edge_id);
      if (crosser.crossingSign(b.v0, b.v1) >= min_sign) {
        edges ~= new ShapeEdge(shape_id, edge_id, b);
      }
    }
  }

  void getCrossingEdges(
      in S2Point a0, in S2Point a1, in S2Shape shape, CrossingType type, ref ShapeEdge[] edges) {
    edges.length = 0;
    getCandidates(a0, a1, shape, _tmpCandidates);
    int min_sign = (type == CrossingType.ALL) ? 0 : 1;
    auto crosser = new S2CopyingEdgeCrosser(a0, a1);
    foreach (ShapeEdgeId candidate; _tmpCandidates) {
      int edge_id = candidate.edgeId;
      S2Shape.Edge b = shape.edge(edge_id);
      if (crosser.crossingSign(b.v0, b.v1) >= min_sign) {
        edges ~= new ShapeEdge(shape.id(), edge_id, b);
      }
    }
  }

  /////////////////////////// Low-Level Methods ////////////////////////////
  //
  // Most clients will not need the following methods.  They can be slightly
  // more efficient but are harder to use, since they require the client to do
  // all the actual crossing tests.

  // Returns a superset of the edges that intersect a query edge (a0, a1).
  // This method is useful for clients that want to test intersections in some
  // other way, e.g. using S2::EdgeOrVertexCrossing().
  ShapeEdgeId[] getCandidates(in S2Point a0, in S2Point a1) {
    ShapeEdgeId[] edges;
    getCandidates(a0, a1, edges);
    return edges;
  }


  // A specialized version of GetCandidates() that only returns the edges that
  // belong to a particular S2Shape.
  ShapeEdgeId[] getCandidates(in S2Point a0, in S2Point a1, in S2Shape shape) {
    ShapeEdgeId[] edges;
    getCandidates(a0, a1, shape, edges);
    return edges;
  }

  // These versions can be more efficient when they are called many times,
  // since they do not require allocating a new vector on each call.
  void getCandidates(in S2Point a0, in S2Point a1, out ShapeEdgeId[] edges) {
    int num_edges = countEdgesUpTo(_index, MAX_BRUTE_FORCE_EDGES + 1);
    if (num_edges <= MAX_BRUTE_FORCE_EDGES) {
      edges.reserve(num_edges);
    }
    visitRawCandidates(a0, a1, (in ShapeEdgeId id) {
          edges ~= id;
          return true;
        });
    if (edges.length > 1) {
      edges = edges.sort().uniq().array();
    }
  }

  void getCandidates(in S2Point a0, in S2Point a1, in S2Shape shape, out ShapeEdgeId[] edges) {
    int num_edges = shape.numEdges();
    if (num_edges <= MAX_BRUTE_FORCE_EDGES) {
      edges.reserve(num_edges);
    }
    visitRawCandidates(a0, a1, shape, (in ShapeEdgeId id) {
          edges ~= id;
          return true;
        });
    if (edges.length > 1) {
      edges = edges.sort.uniq.array;
    }
  }

  // A function that is called with each candidate intersecting edge.  The
  // function may return false in order to request that the algorithm should
  // be terminated, i.e. no further crossings are needed.
  alias ShapeEdgeIdVisitor = bool delegate(in ShapeEdgeId id);

  // Visits a superset of the edges that intersect the query edge (a0, a1),
  // terminating early if the given ShapeEdgeIdVisitor returns false (in which
  // case this function returns false as well).
  //
  // CAVEAT: Edges may be visited more than once.
  bool visitRawCandidates(in S2Point a0, in S2Point a1, ShapeEdgeIdVisitor visitor) {
    int num_edges = countEdgesUpTo(_index, MAX_BRUTE_FORCE_EDGES + 1);
    if (num_edges <= MAX_BRUTE_FORCE_EDGES) {
      int num_shape_ids = _index.numShapeIds();
      for (int s = 0; s < num_shape_ids; ++s) {
        const(S2Shape) shape = _index.shape(s);
        if (shape is null) continue;
        int num_shape_edges = shape.numEdges();
        for (int e = 0; e < num_shape_edges; ++e) {
          if (!visitor(ShapeEdgeId(s, e))) return false;
        }
      }
      return true;
    }
    return visitCells(a0, a1, (in S2ShapeIndexCell cell) {
          for (int s = 0; s < cell.numClipped(); ++s) {
            const(S2ClippedShape) clipped = cell.clipped(s);
            for (int j = 0; j < clipped.numEdges(); ++j) {
              if (!visitor(ShapeEdgeId(clipped.shapeId(), clipped.edge(j)))) {
                return false;
              }
            }
          }
          return true;
        });
  }

  bool visitRawCandidates(
      in S2Point a0, in S2Point a1, in S2Shape shape, ShapeEdgeIdVisitor visitor) {
    int num_edges = shape.numEdges();
    if (num_edges <= MAX_BRUTE_FORCE_EDGES) {
      for (int e = 0; e < num_edges; ++e) {
        if (!visitor(ShapeEdgeId(shape.id(), e))) return false;
      }
      return true;
    }
    return visitCells(a0, a1, (in S2ShapeIndexCell cell) {
          const(S2ClippedShape)* clipped = cell.findClipped(shape.id());
          if (clipped is null) return true;
          for (int j = 0; j < clipped.numEdges(); ++j) {
            if (!visitor(ShapeEdgeId(shape.id(), clipped.edge(j)))) return false;
          }
          return true;
        });
  }

  // A function that is called with each S2ShapeIndexCell that might contain
  // edges intersecting the given query edge.  The function may return false
  // in order to request that the algorithm should be terminated, i.e. no
  // further crossings are needed.
  alias CellVisitor = bool delegate(in S2ShapeIndexCell);

  // Visits all S2ShapeIndexCells that might contain edges intersecting the
  // given query edge (a0, a1), terminating early if the given CellVisitor
  // returns false (in which case this function returns false as well).
  //
  // NOTE: Each candidate cell is visited exactly once.
  bool visitCells(in S2Point a0, in S2Point a1, CellVisitor visitor) {
    _visitor = visitor;
    FaceSegmentVector segments;
    getFaceSegments(a0, a1, segments);
    foreach (segment; segments) {
      _a0 = segment.a;
      _a1 = segment.b;
      // Optimization: rather than always starting the recursive subdivision at
      // the top level face cell, instead we start at the smallest S2CellId that
      // contains the edge (the "edge root cell").  This typically lets us skip
      // quite a few levels of recursion since most edges are short.
      R2Rect edge_bound = R2Rect.fromPointPair(_a0, _a1);
      auto pcell = new S2PaddedCell(S2CellId.fromFace(segment.face), 0);
      S2CellId edge_root = pcell.shrinkToFit(edge_bound);

      // Now we need to determine how the edge root cell is related to the cells
      // in the spatial index (cell_map_).  There are three cases:
      //
      //  1. edge_root is an index cell or is contained within an index cell.
      //     In this case we only need to look at the contents of that cell.
      //  2. edge_root is subdivided into one or more index cells.  In this case
      //     we recursively subdivide to find the cells intersected by a0a1.
      //  3. edge_root does not intersect any index cells.  In this case there
      //     is nothing to do.
      S2ShapeIndex.CellRelation relation = _iter.locate(edge_root);
      if (relation == S2ShapeIndex.CellRelation.INDEXED) {
        // edge_root is an index cell or is contained by an index cell (case 1).
        enforce(_iter.id().contains(edge_root));
        if (!visitor(_iter.cell())) return false;
      } else if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
        // edge_root is subdivided into one or more index cells (case 2).  We
        // find the cells intersected by a0a1 using recursive subdivision.
        if (!edge_root.isFace()) pcell = new S2PaddedCell(edge_root, 0);
        if (!visitCells(pcell, edge_bound)) return false;
      }
    }
    return true;
  }


  // Visits all S2ShapeIndexCells within "root" that might contain edges
  // intersecting the given query edge (a0, a1), terminating early if the
  // given CellVisitor returns false (in which case this function returns
  // false as well).
  //
  // NOTE: Each candidate cell is visited exactly once.
  bool visitCells(in S2Point a0, in S2Point a1, in S2PaddedCell root, CellVisitor visitor) {
    _visitor = visitor;
    FaceSegmentVector segments;
    getFaceSegments(a0, a1, segments);
    foreach (segment; segments) {
      _a0 = segment.a;
      _a1 = segment.b;

      // Optimization: rather than always starting the recursive subdivision at
      // the top level face cell, instead we start at the smallest S2CellId that
      // contains the edge (the "edge root cell").  This typically lets us skip
      // quite a few levels of recursion since most edges are short.
      R2Rect edge_bound = R2Rect.fromPointPair(_a0, _a1);
      auto pcell = new S2PaddedCell(S2CellId.fromFace(segment.face), 0);
      S2CellId edge_root = pcell.shrinkToFit(edge_bound);

      // Now we need to determine how the edge root cell is related to the cells
      // in the spatial index (cell_map_).  There are three cases:
      //
      //  1. edge_root is an index cell or is contained within an index cell.
      //     In this case we only need to look at the contents of that cell.
      //  2. edge_root is subdivided into one or more index cells.  In this case
      //     we recursively subdivide to find the cells intersected by a0a1.
      //  3. edge_root does not intersect any index cells.  In this case there
      //     is nothing to do.
      S2ShapeIndex.CellRelation relation = _iter.locate(edge_root);
      if (relation == S2ShapeIndex.CellRelation.INDEXED) {
        // edge_root is an index cell or is contained by an index cell (case 1).
        enforce(_iter.id().contains(edge_root));
        if (!visitor(_iter.cell())) return false;
      } else if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
        // edge_root is subdivided into one or more index cells (case 2).  We
        // find the cells intersected by a0a1 using recursive subdivision.
        if (!edge_root.isFace()) pcell = new S2PaddedCell(edge_root, 0);
        if (!visitCells(pcell, edge_bound)) return false;
      }
    }
    return true;
  }


  // Given a query edge AB and a cell "root", returns all S2ShapeIndex cells
  // within "root" that might contain edges intersecting AB.
  void getCells(in S2Point a0, in S2Point a1, in S2PaddedCell root, out const(S2ShapeIndexCell)[] cells) {
    cells.length = 0;
    visitCells(
        a0, a1, root, (in S2ShapeIndexCell cell) {
          cells ~= cell;
          return true;
        });
  }

  ///////////////////// DEPRECATED METHODS //////////////////////////////

  alias EdgeMap = int[][S2Shape];

  deprecated("Use GetCrossingEdges")
  bool getCrossings(
      in S2Point a0, in S2Point a1, in S2Shape shape, CrossingType type, ref int[] edges) {
    edges.length = 0;
    foreach (const(ShapeEdge) edge; getCrossingEdges(a0, a1, shape, type)) {
      edges ~= edge.id().edgeId;
    }
    return !edges.empty();
  }

  deprecated("Use GetCrossingEdges")
  bool getCrossings(in S2Point a0, in S2Point a1, CrossingType type, out EdgeMap edge_map) {
    // Since this API is obsolete, don't worry about reserving vectors, etc.
    foreach (const(ShapeEdge) edge; getCrossingEdges(a0, a1, type)) {
      S2Shape shape = _index.shape(edge.id().shapeId);
      edge_map[shape] ~= edge.id().edgeId;
    }
    return !edge_map.empty();
  }

  deprecated("Use method returning std::vector")
  bool getCandidates(in S2Point a0, in S2Point a1, in S2Shape shape, out int[] edges) {
    foreach (const(ShapeEdgeId) edge; getCandidates(a0, a1, shape)) {
      edges ~= edge.edgeId;
    }
    return !edges.empty();
  }

  deprecated("Use method returning std::vector")
  bool getCandidates(in S2Point a0, in S2Point a1, out EdgeMap edge_map) {
    // Since this API is obsolete, don't worry about reserving vectors, etc.
    foreach (const(ShapeEdgeId) edge; getCandidates(a0, a1)) {
      S2Shape shape = _index.shape(edge.shapeId);
      edge_map[shape] ~= edge.edgeId;
    }
    return !edge_map.empty();
  }

 private:
  // Internal methods are documented with their definitions.

  /**
   * Computes the index cells intersected by the current edge that are
   * descendants of "pcell" and calls visitor_ for each one.
   *
   * WARNING: This function is recursive with a maximum depth of 30.  The frame
   * size is about 2K in versions of GCC prior to 4.7 due to poor overlapping
   * of storage for temporaries.  This is fixed in GCC 4.7, reducing the frame
   * size to about 350 bytes (i.e., worst-case total stack usage of about 10K).
   */
  bool visitCells(S2PaddedCell pcell, in R2Rect edge_bound) {
    _iter.seek(pcell.id().rangeMin());
    if (_iter.done() || _iter.id() > pcell.id().rangeMax()) {
      // The index does not contain "pcell" or any of its descendants.
      return true;
    }
    if (_iter.id() == pcell.id()) {
      return _visitor(_iter.cell());
    }

    // Otherwise, split the edge among the four children of "pcell".
    R2Point center = pcell.middle().lo();
    if (edge_bound[0].hi() < center[0]) {
      // Edge is entirely contained in the two left children.
      return clipVAxis(edge_bound, center[1], 0, pcell);
    } else if (edge_bound[0].lo() >= center[0]) {
      // Edge is entirely contained in the two right children.
      return clipVAxis(edge_bound, center[1], 1, pcell);
    } else {
      R2Rect[2] child_bounds;
      splitUBound(edge_bound, center[0], child_bounds);
      if (edge_bound[1].hi() < center[1]) {
        // Edge is entirely contained in the two lower children.
        return visitCells(new S2PaddedCell(pcell, 0, 0), child_bounds[0])
            && visitCells(new S2PaddedCell(pcell, 1, 0), child_bounds[1]);
      } else if (edge_bound[1].lo() >= center[1]) {
        // Edge is entirely contained in the two upper children.
        return visitCells(new S2PaddedCell(pcell, 0, 1), child_bounds[0])
            && visitCells(new S2PaddedCell(pcell, 1, 1), child_bounds[1]);
      } else {
        // The edge bound spans all four children.  The edge itself intersects
        // at most three children (since no padding is being used).
        return clipVAxis(child_bounds[0], center[1], 0, pcell)
            && clipVAxis(child_bounds[1], center[1], 1, pcell);
      }
    }
  }

  /**
   * Given either the left (i=0) or right (i=1) side of a padded cell "pcell",
   * determine whether the current edge intersects the lower child, upper child,
   * or both children, and call VisitCells() recursively on those children.
   * "center" is the v-coordinate at the center of "pcell".
   */
  bool clipVAxis(in R2Rect edge_bound, double center, int i, S2PaddedCell pcell) {
    if (edge_bound[1].hi() < center) {
      // Edge is entirely contained in the lower child.
      return visitCells(new S2PaddedCell(pcell, i, 0), edge_bound);
    } else if (edge_bound[1].lo() >= center) {
      // Edge is entirely contained in the upper child.
      return visitCells(new S2PaddedCell(pcell, i, 1), edge_bound);
    } else {
      // The edge intersects both children.
      R2Rect[2] child_bounds;
      splitVBound(edge_bound, center, child_bounds);
      return visitCells(new S2PaddedCell(pcell, i, 0), child_bounds[0])
          && visitCells(new S2PaddedCell(pcell, i, 1), child_bounds[1]);
    }
  }

  /**
   * Split the current edge into two child edges at the given u-value "u" and
   * return the bound for each child.
   */
  void splitUBound(in R2Rect edge_bound, double u, ref R2Rect[2] child_bounds) const {
    // See comments in MutableS2ShapeIndex::ClipUBound.
    double v = edge_bound[1].project(interpolateDouble(u, _a0[0], _a1[0], _a0[1], _a1[1]));

    // "diag_" indicates which diagonal of the bounding box is spanned by a0a1:
    // it is 0 if a0a1 has positive slope, and 1 if a0a1 has negative slope.
    int diag = (_a0[0] > _a1[0]) != (_a0[1] > _a1[1]);
    splitBound(edge_bound, 0, u, diag, v, child_bounds);
  }

  /**
   * Split the current edge into two child edges at the given v-value "v" and
   * return the bound for each child.
   */
  void splitVBound(in R2Rect edge_bound, double v, ref R2Rect[2] child_bounds) const {
    double u = edge_bound[0].project(interpolateDouble(v, _a0[1], _a1[1], _a0[0], _a1[0]));
    int diag = (_a0[0] > _a1[0]) != (_a0[1] > _a1[1]);
    splitBound(edge_bound, diag, u, 0, v, child_bounds);
  }

  /**
   * Split the current edge into two child edges at the given point (u,v) and
   * return the bound for each child.  "u_end" and "v_end" indicate which bound
   * endpoints of child 1 will be updated.
   */
  static void splitBound(
      in R2Rect edge_bound, int u_end, double u, int v_end, double v, ref R2Rect[2] child_bounds) {
    child_bounds[0] = edge_bound;
    child_bounds[0][0][1 - u_end] = u;
    child_bounds[0][1][1 - v_end] = v;
    enforce(!child_bounds[0].isEmpty());
    enforce(edge_bound.contains(child_bounds[0]));

    child_bounds[1] = edge_bound;
    child_bounds[1][0][u_end] = u;
    child_bounds[1][1][v_end] = v;
    enforce(!child_bounds[1].isEmpty());
    enforce(edge_bound.contains(child_bounds[1]));
  }

  S2ShapeIndex _index;

  //////////// Temporary storage used while processing a query ///////////

  R2Point _a0, _a1;
  S2ShapeIndex.Iterator _iter;
  CellVisitor _visitor;

  // Avoids repeated allocation when methods are called many times.
  ShapeEdgeId[] _tmpCandidates;
}
