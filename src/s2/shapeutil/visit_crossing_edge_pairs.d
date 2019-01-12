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

module s2.shapeutil.visit_crossing_edge_pairs;

import s2.logger;
import s2.s2cell_id;
import s2.s2crossing_edge_query : CrossingType;
import s2.s2crossing_edge_query;
import s2.s2edge_crosser;
import s2.s2error;
import s2.s2padded_cell;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2wedge_relations : getWedgeRelation, WedgeRelation;
import s2.shapeutil.range_iterator;
import s2.shapeutil.shape_edge;


/**
 * A function that is called with pairs of crossing edges.  The function may
 * return false in order to request that the algorithm should be terminated,
 * i.e. no further crossings are needed.
 *
 * "is_interior" indicates that the crossing is at a point interior to both
 * edges (i.e., not at a vertex).  (The calling function already has this
 * information and it is moderately expensive to recompute.)
 */
alias EdgePairVisitor = bool delegate(in ShapeEdge a, in ShapeEdge b, bool isInterior);

/**
 * Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
 * early if the given EdgePairVisitor function returns false (in which case
 * VisitCrossings returns false as well).  "type" indicates whether all
 * crossings should be visited, or only interior crossings.
 *
 * CAVEAT: Crossings may be visited more than once.
 */
bool visitCrossingEdgePairs(
    S2ShapeIndex index, CrossingType type, EdgePairVisitor visitor) {
  bool need_adjacent = (type == CrossingType.ALL);
  return visitCrossings(index, type, need_adjacent, visitor);
}

/**
 * Like the above, but visits all pairs of crossing edges where one edge comes
 * from each S2ShapeIndex.
 *
 * CAVEAT: Crossings may be visited more than once.
 */
bool visitCrossingEdgePairs(
    S2ShapeIndex a_index, S2ShapeIndex b_index, CrossingType type, EdgePairVisitor visitor) {
  // We look for S2CellId ranges where the indexes of A and B overlap, and
  // then test those edges for crossings.

  // TODO(ericv): Use brute force if the total number of edges is small enough
  // (using a larger threshold if the S2ShapeIndex is not constructed yet).
  auto ai = new RangeIterator(a_index);
  auto bi = new RangeIterator(b_index);
  auto ab = new IndexCrosser(a_index, b_index, type, visitor, false);  // Tests A against B
  auto ba = new IndexCrosser(b_index, a_index, type, visitor, true);   // Tests B against A
  while (!ai.done() || !bi.done()) {
    if (ai.rangeMax() < bi.rangeMin()) {
      // The A and B cells don't overlap, and A precedes B.
      ai.seekTo(bi);
    } else if (bi.rangeMax() < ai.rangeMin()) {
      // The A and B cells don't overlap, and B precedes A.
      bi.seekTo(ai);
    } else {
      // One cell contains the other.  Determine which cell is larger.
      long ab_relation = ai.id().lsb() - bi.id().lsb();
      if (ab_relation > 0) {
        // A's index cell is larger.
        if (!ab.visitCrossings(ai, bi)) return false;
      } else if (ab_relation < 0) {
        // B's index cell is larger.
        if (!ba.visitCrossings(bi, ai)) return false;
      } else {
        // The A and B cells are the same.
        if (ai.cell().numEdges() > 0 && bi.cell().numEdges() > 0) {
          if (!ab.visitCellCellCrossings(ai.cell(), bi.cell())) return false;
        }
        ai.next();
        bi.next();
      }
    }
  }
  return true;
}

// IndexCrosser is a helper class for finding the edge crossings between a
// pair of S2ShapeIndexes.  It is instantiated twice, once for the index pair
// (A,B) and once for the index pair (B,A), in order to be able to test edge
// crossings in the most efficient order.
// TODO: Resume here.
class IndexCrosser {
 public:
  // If "swapped" is true, the loops A and B have been swapped.  This affects
  // how arguments are passed to the given loop relation, since for example
  // A.Contains(B) is not the same as B.Contains(A).
  this(S2ShapeIndex a_index, S2ShapeIndex b_index,
      CrossingType type, EdgePairVisitor visitor, bool swapped) {
    _aIndex = a_index;
    _bIndex = b_index;
    _visitor = visitor;
    _minCrossingSign = type == CrossingType.INTERIOR ? 1 : 0;
    _swapped = swapped;
    _bQuery = new S2CrossingEdgeQuery(_bIndex);
  }

  // Given two iterators positioned such that ai->id().Contains(bi->id()),
  // visits all crossings between edges of A and B that intersect a->id().
  // Terminates early and returns false if visitor_ returns false.
  // Advances both iterators past ai->id().
  bool visitCrossings(RangeIterator ai, RangeIterator bi) {
    logger.logError(!ai.id().contains(bi.id()), "ai must contain bi");
    if (ai.cell().numEdges() == 0) {
      // Skip over the cells of B using binary search.
      bi.seekBeyond(ai);
    } else {
      // If ai->id() intersects many edges of B, then it is faster to use
      // S2CrossingEdgeQuery to narrow down the candidates.  But if it
      // intersects only a few edges, it is faster to check all the crossings
      // directly.  We handle this by advancing "bi" and keeping track of how
      // many edges we would need to test.
      enum int kEdgeQueryMinEdges = 23;
      int b_edges = 0;
      _bCells.length = 0;
      do {
        int cell_edges = bi.cell().numEdges();
        if (cell_edges > 0) {
          b_edges += cell_edges;
          if (b_edges >= kEdgeQueryMinEdges) {
            // There are too many edges, so use an S2CrossingEdgeQuery.
            if (!visitSubcellCrossings(ai.cell(), ai.id())) return false;
            bi.seekBeyond(ai);
            return true;
          }
          _bCells ~= bi.cell();
        }
        bi.next();
      } while (bi.id() <= ai.rangeMax());
      if (_bCells.length != 0) {
        // Test all the edge crossings directly.
        getShapeEdges(_aIndex, ai.cell(), _aShapeEdges);
        getShapeEdges(_bIndex, _bCells, _bShapeEdges);
        if (!visitEdgesEdgesCrossings(_aShapeEdges, _bShapeEdges)) {
          return false;
        }
      }
    }
    ai.next();
    return true;
  }

  // Given two index cells, visits all crossings between edges of those cells.
  // Terminates early and returns false if visitor_ returns false.
  bool visitCellCellCrossings(in S2ShapeIndexCell a_cell, in S2ShapeIndexCell b_cell) {
    // Test all edges of "a_cell" against all edges of "b_cell".
    getShapeEdges(_aIndex, a_cell, _aShapeEdges);
    getShapeEdges(_bIndex, b_cell, _bShapeEdges);
    return visitEdgesEdgesCrossings(_aShapeEdges, _bShapeEdges);
  }

 private:
  bool visitEdgePair(in ShapeEdge a, in ShapeEdge b, bool is_interior) {
    if (_swapped) {
      return _visitor(b, a, is_interior);
    } else {
      return _visitor(a, b, is_interior);
    }
  }

  // Visits all crossings of the current edge with all edges of the given index
  // cell of B.  Terminates early and returns false if visitor_ returns false.
  bool visitEdgeCellCrossings(in ShapeEdge a, in S2ShapeIndexCell b_cell) {
    // Test the current edge of A against all edges of "b_cell".

    // Note that we need to use a new S2EdgeCrosser (or call Init) whenever we
    // replace the contents of b_shape_edges_, since S2EdgeCrosser requires that
    // its S2Point arguments point to values that persist between Init() calls.
    getShapeEdges(_bIndex, b_cell, _bShapeEdges);
    auto crosser = new S2CopyingEdgeCrosser(a.v0(), a.v1());
    foreach (const ShapeEdge b; _bShapeEdges) {
      if (crosser.c() != b.v0()) {
        crosser.restartAt(b.v0());
      }
      int sign = crosser.crossingSign(b.v1());
      if (sign >= _minCrossingSign) {
        if (!visitEdgePair(a, b, sign == 1)) return false;
      }
    }
    return true;
  }

  // Visits all crossings of any edge in "a_cell" with any index cell of B that
  // is a descendant of "b_id".  Terminates early and returns false if
  // visitor_ returns false.
  bool visitSubcellCrossings(in S2ShapeIndexCell a_cell, S2CellId b_id) {
    // Test all edges of "a_cell" against the edges contained in B index cells
    // that are descendants of "b_id".
    getShapeEdges(_aIndex, a_cell, _aShapeEdges);
    auto b_root = new S2PaddedCell(b_id, 0);
    foreach (const ShapeEdge a; _aShapeEdges) {
      // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
      // of B that might contain crossing edges.
      if (!_bQuery.visitCells(
              a.v0(), a.v1(), b_root,
              (const S2ShapeIndexCell cell)  => visitEdgeCellCrossings(a, cell))) {
        return false;
      }
    }
    return true;
  }

  // Visits all crossings of any edge in "a_edges" with any edge in "b_edges".
  bool visitEdgesEdgesCrossings(
      in ShapeEdgeVector a_edges, in ShapeEdgeVector b_edges) {
    // Test all edges of "a_edges" against all edges of "b_edges".
    foreach (const ShapeEdge a; a_edges) {
      auto crosser = new S2EdgeCrosser(a.v0(), a.v1());
      foreach (const ShapeEdge b; b_edges) {
        if (*crosser.c() != b.v0()) {
          crosser.restartAt(b.v0());
        }
        int sign = crosser.crossingSign(b.v1());
        if (sign >= _minCrossingSign) {
          if (!visitEdgePair(a, b, sign == 1)) return false;
        }
      }
    }
    return true;
  }

  S2ShapeIndex _aIndex;
  S2ShapeIndex _bIndex;
  const EdgePairVisitor _visitor;
  const int _minCrossingSign;
  const bool _swapped;

  // Temporary data declared here to avoid repeated memory allocations.
  S2CrossingEdgeQuery _bQuery;
  const(S2ShapeIndexCell)[] _bCells;
  ShapeEdgeVector _aShapeEdges;
  ShapeEdgeVector _bShapeEdges;
}

// Given an S2ShapeIndex containing a single polygonal shape (e.g., an
// S2Polygon or S2Loop), return true if any loop has a self-intersection
// (including duplicate vertices) or crosses any other loop (including vertex
// crossings and duplicate edges) and set "error" to a human-readable error
// message.  Otherwise return false and leave "error" unchanged.
//
// This method is used to implement the FindValidationError methods of S2Loop
// and S2Polygon.
//
// TODO(ericv): Add an option to support S2LaxPolygonShape rules (i.e.,
// duplicate vertices and edges are allowed, but loop crossings are not).
bool findSelfIntersection(S2ShapeIndex index, out S2Error error) {
  if (index.numShapeIds() == 0) return false;
  if (index.numShapeIds() != 1)
    logger.logError("index.numShapeIds()=", index.numShapeIds(), ", expected 1.");
  const(S2Shape) shape = index.shape(0);

  // Visit all crossing pairs except possibly for ones of the form (AB, BC),
  // since such pairs are very common and FindCrossingError() only needs pairs
  // of the form (AB, AC).
  return !visitCrossings(
      index, CrossingType.ALL, false /*need_adjacent*/,
      (in ShapeEdge a, in ShapeEdge b, bool isInterior) =>
          !findCrossingError(shape, a, b, isInterior, error));
}

// Ensure that we don't usually need to allocate memory when collecting the
// edges in an S2ShapeIndex cell (which by default have about 10 edges).
// TODO(vnayar): Investigate if it's worth it to avoid dynamic arrays.
alias ShapeEdgeVector = ShapeEdge[];

// Given a vector of edges within an S2ShapeIndexCell, visit all pairs of
// crossing edges (of the given CrossingType).
private bool visitCrossings(
    in ShapeEdgeVector shape_edges, CrossingType type, bool need_adjacent,
    EdgePairVisitor visitor) {
  const int min_crossing_sign = (type == CrossingType.INTERIOR) ? 1 : 0;
  size_t num_edges = shape_edges.length;
  for (int i = 0; i + 1 < num_edges; ++i) {
    const ShapeEdge a = shape_edges[i];
    int j = i + 1;
    // A common situation is that an edge AB is followed by an edge BC.  We
    // only need to visit such crossings if "need_adjacent" is true (even if
    // AB and BC belong to different edge chains).
    if (!need_adjacent && a.v1() == shape_edges[j].v0()) {
      if (++j >= num_edges) break;
    }
    auto crosser = new S2CopyingEdgeCrosser(a.v0(), a.v1());
    for (; j < num_edges; ++j) {
      const ShapeEdge b = shape_edges[j];
      if (crosser.c() != b.v0()) {
        crosser.restartAt(b.v0());
      }
      int sign = crosser.crossingSign(b.v1());
      if (sign >= min_crossing_sign) {
        if (!visitor(a, b, sign == 1)) return false;
      }
    }
  }
  return true;
}

// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
// early if the given EdgePairVisitor function returns false (in which case
// VisitCrossings returns false as well).  "type" indicates whether all
// crossings should be visited, or only interior crossings.
//
// If "need_adjacent" is false, then edge pairs of the form (AB, BC) may
// optionally be ignored (even if the two edges belong to different edge
// chains).  This option exists for the benefit of FindSelfIntersection(),
// which does not need such edge pairs (see below).
static bool visitCrossings(
    S2ShapeIndex index, CrossingType type, bool need_adjacent, in EdgePairVisitor visitor) {
  // TODO(ericv): Use brute force if the total number of edges is small enough
  // (using a larger threshold if the S2ShapeIndex is not constructed yet).
  ShapeEdgeVector shape_edges;
  for (auto it = new S2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
       !it.done(); it.next()) {
    getShapeEdges(index, it.cell(), shape_edges);
    if (!visitCrossings(shape_edges, type, need_adjacent, visitor)) {
      return false;
    }
  }
  return true;
}

// Returns a vector containing all edges in the given S2ShapeIndexCell.
// (The result is returned as an output parameter so that the same storage can
// be reused, rather than allocating a new temporary vector each time.)
private void getShapeEdges(in S2ShapeIndex index, in S2ShapeIndexCell cell,
    ref ShapeEdgeVector shape_edges) {
  shape_edges.length = 0;
  appendShapeEdges(index, cell, shape_edges);
}

// Returns a vector containing all edges in the given S2ShapeIndexCell vector.
// (The result is returned as an output parameter so that the same storage can
// be reused, rather than allocating a new temporary vector each time.)
private void getShapeEdges(in S2ShapeIndex index, in S2ShapeIndexCell[] cells,
    ref ShapeEdgeVector shape_edges) {
  shape_edges.length = 0;
  for (int c = 0; c < cells.length; ++c) {
    appendShapeEdges(index, cells[c], shape_edges);
  }
}


// Given two loop edges that cross (including at a shared vertex), return true
// if there is a crossing error and set "error" to a human-readable message.
private bool findCrossingError(
    in S2Shape shape, in ShapeEdge a, in ShapeEdge b, bool is_interior, out S2Error error) {
  bool is_polygon = shape.numChains() > 1;
  S2Shape.ChainPosition ap = shape.chainPosition(a.id().edgeId);
  S2Shape.ChainPosition bp = shape.chainPosition(b.id().edgeId);
  if (is_interior) {
    if (ap.chainId != bp.chainId) {
      error.initialize(S2Error.Code.POLYGON_LOOPS_CROSS,
          "Loop %d edge %d crosses loop %d edge %d",
          ap.chainId, ap.offset, bp.chainId, bp.offset);
    } else {
      initLoopError(S2Error.Code.LOOP_SELF_INTERSECTION,
                    "Edge %d crosses edge %d", ap, bp, is_polygon, error);
    }
    return true;
  }
  // Loops are not allowed to have duplicate vertices, and separate loops
  // are not allowed to share edges or cross at vertices.  We only need to
  // check a given vertex once, so we also require that the two edges have
  // the same end vertex.
  if (a.v1() != b.v1()) return false;
  if (ap.chainId == bp.chainId) {
    initLoopError(S2Error.Code.DUPLICATE_VERTICES,
        "Edge %d has duplicate vertex with edge %d",
        ap, bp, is_polygon, error);
    return true;
  }
  int a_len = shape.chain(ap.chainId).length;
  int b_len = shape.chain(bp.chainId).length;
  int a_next = (ap.offset + 1 == a_len) ? 0 : ap.offset + 1;
  int b_next = (bp.offset + 1 == b_len) ? 0 : bp.offset + 1;
  S2Point a2 = shape.chainEdge(ap.chainId, a_next).v1;
  S2Point b2 = shape.chainEdge(bp.chainId, b_next).v1;
  if (a.v0() == b.v0() || a.v0() == b2) {
    // The second edge index is sometimes off by one, hence "near".
    error.initialize(S2Error.Code.POLYGON_LOOPS_SHARE_EDGE,
        "Loop %d edge %d has duplicate near loop %d edge %d",
        ap.chainId, ap.offset, bp.chainId, bp.offset);
    return true;
  }
  // Since S2ShapeIndex loops are oriented such that the polygon interior is
  // always on the left, we need to handle the case where one wedge contains
  // the complement of the other wedge.  This is not specifically detected by
  // GetWedgeRelation, so there are two cases to check for.
  //
  // Note that we don't need to maintain any state regarding loop crossings
  // because duplicate edges are detected and rejected above.
  if (getWedgeRelation(a.v0(), a.v1(), a2, b.v0(), b2) == WedgeRelation.WEDGE_PROPERLY_OVERLAPS
      && getWedgeRelation(a.v0(), a.v1(), a2, b2, b.v0())
          == WedgeRelation.WEDGE_PROPERLY_OVERLAPS) {
    error.initialize(S2Error.Code.POLYGON_LOOPS_CROSS,
        "Loop %d edge %d crosses loop %d edge %d",
        ap.chainId, ap.offset, bp.chainId, bp.offset);
    return true;
  }
  return false;
}

// Appends all edges in the given S2ShapeIndexCell to the given vector.
private void appendShapeEdges(in S2ShapeIndex index,
                             in S2ShapeIndexCell cell,
                             ref ShapeEdgeVector shape_edges) {
  for (int s = 0; s < cell.numClipped(); ++s) {
    const S2ClippedShape clipped = cell.clipped(s);
    const S2Shape shape = index.shape(clipped.shapeId());
    int num_edges = clipped.numEdges();
    for (int i = 0; i < num_edges; ++i) {
      shape_edges ~= new ShapeEdge(shape, clipped.edge(i));
    }
  }
}

// Helper function that formats a loop error message.  If the loop belongs to
// a multi-loop polygon, adds a prefix indicating which loop is affected.
private void initLoopError(S2Error.Code code, string format,
                          S2Shape.ChainPosition ap, S2Shape.ChainPosition bp,
                          bool is_polygon, S2Error error) {
  error.initialize(code, format, ap.offset, bp.offset);
  if (is_polygon) {
    error.initialize(code, "Loop %d: %s", ap.chainId, error.text());
  }
}
