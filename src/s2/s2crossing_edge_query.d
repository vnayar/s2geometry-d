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
import s2.shapeutil.shape_edge_id;

import std.algorithm : sort, uniq;
import std.array : array;

//#include "s2/s2padded_cell.h"
//#include "s2/s2shape_index.h"
//#include "s2/s2shapeutil_shape_edge.h"
//#include "s2/s2shapeutil_shape_edge_id.h"

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType::INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType::ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
enum CrossingType { INTERIOR, ALL }

// For small loops it is faster to use brute force.  The threshold below was
// determined using the benchmarks in the unit test.
private enum int MAX_BRUTE_FORCE_EDGES = 27;

/+
// S2CrossingEdgeQuery is used to find edges or shapes that are crossed by
// an edge.  Here is an example showing how to index a set of polylines,
// and then find the polylines that are crossed by a given edge AB:
//
// void Test(const vector<S2Polyline*>& polylines,`
//           const S2Point& a0, const S2Point &a1) {
//   MutableS2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(absl::make_unique<S2Polyline::Shape>(polyline));
//   }
//   S2CrossingEdgeQuery query(&index);
//   for (const auto& edge : query.GetCrossingEdges(a, b, CrossingType::ALL)) {
//     CHECK_GE(S2::CrossingSign(a0, a1, edge.v0(), edge.v1()), 0);
//   }
// }
//
// Note that if you need to query many edges, it is more efficient to declare
// a single S2CrossingEdgeQuery object and reuse it so that temporary storage
// does not need to be reallocated each time.
//
// If you want to find *all* pairs of crossing edges, use
// s2shapeutil::VisitCrossingEdgePairs() instead.
class S2CrossingEdgeQuery {
public:
  // Convenience constructor that calls Init().
  this(in S2ShapeIndex index) {
    init(index);
  }

  // Default constructor; requires Init() to be called.
  this() { }

  const(S2ShapeIndex) index() const {
    return _index;
  }

  // REQUIRES: "index" is not modified after this method is called.
  void init(in S2ShapeIndex index) {
    _index = index;
    _iter.init(index);
  }

  // Returns all edges that intersect the given query edge (a0,a1) and that
  // have the given CrossingType (ALL or INTERIOR).  Edges are sorted and
  // unique.
  ShapeEdge[] getCrossingEdges(in S2Point a0, in S2Point a1, CrossingType type) {
    ShapeEdge[] edges;
    getCrossingEdges(a0, a1, type, edges);
    return edges;
  }

  // A specialized version of GetCrossingEdges() that only returns the edges
  // that belong to a particular S2Shape.
  ShapeEdge[] getCrossingEdges(
      in S2Point a0, in S2Point a1, in S2Shape shape, CrossingType type) {
    ShapeEdge[] edges;
    getCrossingEdges(a0, a1, shape, type, edges);
    return edges;
  }

  // These versions can be more efficient when they are called many times,
  // since they do not require allocating a new vector on each call.
  void getCrossingEdges(in S2Point a0, in S2Point a1, CrossingType type, ref ShapeEdge[] edges) {
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
    getCandidates(a0, a1, _tmpCandidates);
    int min_sign = (type == CrossingType.ALL) ? 0 : 1;
    auto crosser = new S2CopyingEdgeCrosser(a0, a1);
    int shape_id = -1;
    const(S2Shape)* shape = null;
    for (ShapeEdgeId candidate : _tmpCandidates) {
      if (candidate.shapeId != shape_id) {
        shape_id = candidate.shapeId;
        shape = _index.shape(shapeId);
      }
      int edge_id = candidate.edgeId;
      S2Shape.Edge b = shape.edge(edge_id);
      if (crosser.crossingSign(b.v0, b.v1) >= min_sign) {
        edges ~= new ShapeEdge(shape_id, edge_id, b);
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
  void getCandidates(in S2Point a0, in S2Point a1, ref ShapeEdgeId[] edges) {
    edges.length = 0;
    int num_edges = countEdgesUpTo(_index, MAX_BRUTE_FORCE_EDGES + 1);
    if (num_edges <= MAX_BRUTE_FORCE_EDGES) {
      edges.reserve(num_edges);
    }
    visitRawCandidates(a0, a1, (ShapeEdgeId id) {
          edges ~= id;
          return true;
        });
    if (edges.length > 1) {
      edges = edges.sort().uniq().array();
    }
  }

  void getCandidates(in S2Point a0, in S2Point a1, in S2Shape shape, out ShapeEdgeId[] edges) {
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

  // A function that is called with each candidate intersecting edge.  The
  // function may return false in order to request that the algorithm should
  // be terminated, i.e. no further crossings are needed.
  using ShapeEdgeIdVisitor = bool delegate(in ShapeEdgeId id);

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
          const(S2ClippedShape) clipped = cell.findClipped(shape.id());
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
  using CellVisitor = bool delegate(in S2ShapeIndexCell);

  // Visits all S2ShapeIndexCells that might contain edges intersecting the
  // given query edge (a0, a1), terminating early if the given CellVisitor
  // returns false (in which case this function returns false as well).
  //
  // NOTE: Each candidate cell is visited exactly once.
  bool visitCells(in S2Point a0, in S2Point a1, CellVisitor visitor) {
    _visitor = visitor;
    S2::FaceSegmentVector segments;
    S2::GetFaceSegments(a0, a1, &segments);
    for (const auto& segment : segments) {
      a0_ = segment.a;
      a1_ = segment.b;

      // Optimization: rather than always starting the recursive subdivision at
      // the top level face cell, instead we start at the smallest S2CellId that
      // contains the edge (the "edge root cell").  This typically lets us skip
      // quite a few levels of recursion since most edges are short.
      R2Rect edge_bound = R2Rect::FromPointPair(a0_, a1_);
      S2PaddedCell pcell(S2CellId::FromFace(segment.face), 0);
      S2CellId edge_root = pcell.ShrinkToFit(edge_bound);

      // Now we need to determine how the edge root cell is related to the cells
      // in the spatial index (cell_map_).  There are three cases:
      //
      //  1. edge_root is an index cell or is contained within an index cell.
      //     In this case we only need to look at the contents of that cell.
      //  2. edge_root is subdivided into one or more index cells.  In this case
      //     we recursively subdivide to find the cells intersected by a0a1.
      //  3. edge_root does not intersect any index cells.  In this case there
      //     is nothing to do.
      S2ShapeIndex::CellRelation relation = iter_.Locate(edge_root);
      if (relation == S2ShapeIndex::INDEXED) {
        // edge_root is an index cell or is contained by an index cell (case 1).
        DCHECK(iter_.id().contains(edge_root));
        if (!visitor(iter_.cell())) return false;
      } else if (relation == S2ShapeIndex::SUBDIVIDED) {
        // edge_root is subdivided into one or more index cells (case 2).  We
        // find the cells intersected by a0a1 using recursive subdivision.
        if (!edge_root.is_face()) pcell = S2PaddedCell(edge_root, 0);
        if (!VisitCells(pcell, edge_bound)) return false;
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
  bool VisitCells(const S2Point& a0, const S2Point& a1,
                  const S2PaddedCell& root, const CellVisitor& visitor);


  // Given a query edge AB and a cell "root", returns all S2ShapeIndex cells
  // within "root" that might contain edges intersecting AB.
  void GetCells(const S2Point& a0, const S2Point& a1, const S2PaddedCell& root,
                std::vector<const S2ShapeIndexCell*>* cells);

  ///////////////////// DEPRECATED METHODS //////////////////////////////

  struct CompareBtreeLinearSearch {
    using goog_btree_prefer_linear_node_search = std::true_type;
    bool operator()(const S2Shape* x, const S2Shape* y) const {
      return x->id() < y->id();
    }
  };
  using EdgeMap = gtl::btree_map<const S2Shape*, std::vector<int>,
                                 CompareBtreeLinearSearch>;

  ABSL_DEPRECATED("Use GetCrossingEdges")
  bool GetCrossings(const S2Point& a0, const S2Point& a1,
                    const S2Shape* shape, CrossingType type,
                    std::vector<int>* edges);

  ABSL_DEPRECATED("Use GetCrossingEdges")
  bool GetCrossings(const S2Point& a0, const S2Point& a1, CrossingType type,
                    EdgeMap* edge_map);

  ABSL_DEPRECATED("Use method returning std::vector")
  bool GetCandidates(const S2Point& a0, const S2Point& a1, const S2Shape* shape,
                     std::vector<int>* edges);

  ABSL_DEPRECATED("Use method returning std::vector")
  bool GetCandidates(const S2Point& a0, const S2Point& a1, EdgeMap* edge_map);

 private:
  // Internal methods are documented with their definitions.
  bool VisitCells(const S2PaddedCell& pcell, const R2Rect& edge_bound);
  bool ClipVAxis(const R2Rect& edge_bound, double center, int i,
                 const S2PaddedCell& pcell);
  void SplitUBound(const R2Rect& edge_bound, double u,
                   R2Rect child_bounds[2]) const;
  void SplitVBound(const R2Rect& edge_bound, double v,
                   R2Rect child_bounds[2]) const;
  static void SplitBound(const R2Rect& edge_bound, int u_end, double u,
                         int v_end, double v, R2Rect child_bounds[2]);

  const S2ShapeIndex* index_ = nullptr;

  //////////// Temporary storage used while processing a query ///////////

  R2Point a0_, a1_;
  S2ShapeIndex::Iterator iter_;
  const CellVisitor* visitor_;

  // Avoids repeated allocation when methods are called many times.
  std::vector<s2shapeutil::ShapeEdgeId> tmp_candidates_;
}
+/
