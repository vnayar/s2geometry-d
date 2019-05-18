// Copyright 2017 Google Inc. All Rights Reserved.
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
//

// Original Author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2shape_index_region;

import s2.r2point;
import s2.r2rect;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2contains_point_query;
import s2.s2edge_clipping
    : clipToPaddedFace, intersectsRect, FACE_CLIP_ERROR_UV_COORD, INTERSECTS_RECT_ERROR_UV_DIST;
import s2.s2edge_crosser;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2region;
import s2.s2shape;
import s2.s2shape_index;

import std.exception : enforce;

/**
 * This class wraps an S2ShapeIndex object with the additional methods needed
 * to implement the S2Region API, in order to allow S2RegionCoverer to compute
 * S2CellId coverings of arbitrary collections of geometry.
 *
 * These methods could conceivably be made part of S2ShapeIndex itself, but
 * there are several advantages to having a separate class:
 *
 *  - The class can be templated in order to avoid virtual calls and memory
 *    allocation (for iterators) when the concrete S2ShapeIndex type is known.
 *
 *  - Implementing these methods efficiently requires an S2ShapeIndex iterator,
 *    and this design allows a single iterator to be allocated and reused.
 *
 *  - S2Region::Clone() is not a good fit for the S2ShapeIndex API because
 *    it can't be implemented for some subtypes (e.g., EncodedS2ShapeIndex).
 *
 * Example usage:
 *
 * S2CellUnion GetCovering(const S2ShapeIndex& index) {
 *   S2RegionCoverer coverer;
 *   coverer.mutable_options()->set_max_cells(20);
 *   S2CellUnion covering;
 *   coverer.GetCovering(MakeS2ShapeIndexRegion(&index), &covering);
 *   return covering;
 * }
 *
 * This class is not thread-safe.  To use it in parallel, each thread should
 * construct its own instance (this is not expensive).
 */
final class S2ShapeIndexRegion(IndexT) : S2Region {
public:
  // Rather than calling this constructor, which requires specifying the
  // S2ShapeIndex type explicitly, the preferred idiom is to call
  // MakeS2ShapeIndexRegion(&index) instead.  For example:
  //
  //   coverer.GetCovering(MakeS2ShapeIndexRegion(&index), &covering);
  this(IndexT index) {
    _containsQuery = new S2ContainsPointQuery!IndexT(index);
    _iter = _containsQuery.mutableIter();
  }

  inout(IndexT) index() inout {
    return _containsQuery.index();
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  // Clone() returns a *shallow* copy; it does not make a copy of the
  // underlying S2ShapeIndex.
  override
  S2ShapeIndexRegion!IndexT clone() {
    return new S2ShapeIndexRegion!IndexT(index());
  }

  override
  S2Cap getCapBound() {
    S2CellId[] covering;
    getCellUnionBound(covering);
    return new S2CellUnion(covering).getCapBound();
  }

  override
  S2LatLngRect getRectBound() {
    S2CellId[] covering;
    getCellUnionBound(covering);
    return new S2CellUnion(covering).getRectBound();
  }

  // This method currently returns at most 4 cells, unless the index spans
  // multiple faces in which case it may return up to 6 cells.
  override
  void getCellUnionBound(out S2CellId[] cell_ids) {
    // We find the range of S2Cells spanned by the index and choose a level such
    // that the entire index can be covered with just a few cells.  There are
    // two cases:
    //
    //  - If the index intersects two or more faces, then for each intersected
    //    face we add one cell to the covering.  Rather than adding the entire
    //    face, instead we add the smallest S2Cell that covers the S2ShapeIndex
    //    cells within that face.
    //
    //  - If the index intersects only one face, then we first find the smallest
    //    cell S that contains the index cells (just like the case above).
    //    However rather than using the cell S itself, instead we repeat this
    //    process for each of its child cells.  In other words, for each
    //    child cell C we add the smallest S2Cell C' that covers the index cells
    //    within C.  This extra step is relatively cheap and produces much
    //    tighter coverings when the S2ShapeIndex consists of a small region
    //    near the center of a large S2Cell.
    //
    // The following code uses only a single Iterator object because creating an
    // Iterator may be relatively expensive for some S2ShapeIndex types (e.g.,
    // it may involve memory allocation).
    cell_ids.length = 0;
    cell_ids.reserve(6);

    // Find the last S2CellId in the index.
    _iter.finish();
    if (!_iter.prev()) return;  // Empty index.
    const(S2CellId) last_index_id = _iter.id();
    _iter.begin();
    if (_iter.id() != last_index_id) {
      // The index has at least two cells.  Choose an S2CellId level such that
      // the entire index can be spanned with at most 6 cells (if the index
      // spans multiple faces) or 4 cells (it the index spans a single face).
      int level = _iter.id().getCommonAncestorLevel(last_index_id) + 1;

      // For each cell C at the chosen level, we compute the smallest S2Cell
      // that covers the S2ShapeIndex cells within C.
      const(S2CellId) last_id = last_index_id.parent(level);
      for (auto id = _iter.id().parent(level); id != last_id; id = id.next()) {
        // If the cell C does not contain any index cells, then skip it.
        if (id.rangeMax() < _iter.id()) continue;

        // Find the range of index cells contained by C and then shrink C so
        // that it just covers those cells.
        S2CellId first = _iter.id();
        _iter.seek(id.rangeMax().next());
        _iter.prev();
        coverRange(first, _iter.id(), cell_ids);
        _iter.next();
      }
    }
    coverRange(_iter.id(), last_index_id, cell_ids);
  }


  // Returns true if "target" is contained by any single shape.  If the cell
  // is covered by a union of different shapes then it may return false.
  //
  // The implementation is conservative but not exact; if a shape just barely
  // contains the given cell then it may return false.  The maximum error is
  // less than 10 * DBL_EPSILON radians (or about 15 nanometers).
  override
  bool contains(in S2Cell target) {
    S2ShapeIndex.CellRelation relation = _iter.locate(target.id());

    // If the relation is DISJOINT, then "target" is not contained.  Similarly if
    // the relation is SUBDIVIDED then "target" is not contained, since index
    // cells are subdivided only if they (nearly) intersect too many edges.
    if (relation != S2ShapeIndex.CellRelation.INDEXED) return false;

    // Otherwise, the iterator points to an index cell containing "target".
    // If any shape contains the target cell, we return true.
    enforce(_iter.id().contains(target.id()));
    const(S2ShapeIndexCell) cell = _iter.cell();
    for (int s = 0; s < cell.numClipped(); ++s) {
      const(S2ClippedShape) clipped = cell.clipped(s);
      // The shape contains the target cell iff the shape contains the cell
      // center and none of its edges intersects the (padded) cell interior.
      if (_iter.id() == target.id()) {
        if (clipped.numEdges() == 0 && clipped.containsCenter()) return true;
      } else {
        // It is faster to call AnyEdgeIntersects() before Contains().
        if (index().shape(clipped.shapeId()).hasInterior()
            && !anyEdgeIntersects(clipped, target)
            && _containsQuery.shapeContains(_iter, clipped, target.getCenter())) {
          return true;
        }
      }
    }
    return false;
  }

  // Returns true if any shape intersects "target".
  //
  // The implementation is conservative but not exact; if a shape is just
  // barely disjoint from the given cell then it may return true.  The maximum
  // error is less than 10 * DBL_EPSILON radians (or about 15 nanometers).
  override
  bool mayIntersect(in S2Cell target) {
    S2ShapeIndex.CellRelation relation = _iter.locate(target.id());

    // If "target" does not overlap any index cell, there is no intersection.
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) return false;

    // If "target" is subdivided into one or more index cells, then there is an
    // intersection to within the S2ShapeIndex error bound.
    if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) return true;

    // Otherwise, the iterator points to an index cell containing "target".
    //
    // If "target" is an index cell itself, there is an intersection because index
    // cells are created only if they have at least one edge or they are
    // entirely contained by the loop.
    enforce(_iter.id().contains(target.id()));
    if (_iter.id() == target.id()) return true;

    // Test whether any shape intersects the target cell or contains its center.
    const(S2ShapeIndexCell) cell = _iter.cell();
    for (int s = 0; s < cell.numClipped(); ++s) {
      const(S2ClippedShape) clipped = cell.clipped(s);
      if (anyEdgeIntersects(clipped, target)) return true;
      if (_containsQuery.shapeContains(_iter, clipped, target.getCenter())) {
        return true;
      }
    }
    return false;
  }

  // Returns true if the given point is contained by any two-dimensional shape
  // (i.e., polygon).  Boundaries are treated as being semi-open (i.e., the
  // same rules as S2Polygon).  Zero and one-dimensional shapes are ignored by
  // this method (if you need more flexibility, see S2BooleanOperation).
  override
  bool contains(in S2Point p) {
    if (_iter.locate(p)) {
      const(S2ShapeIndexCell) cell = _iter.cell();
      for (int s = 0; s < cell.numClipped(); ++s) {
        if (_containsQuery.shapeContains(_iter, cell.clipped(s), p)) {
          return true;
        }
      }
    }
    return false;
  }

 private:
  alias Iterator = IndexT.Iterator;

  // Computes the smallest S2Cell that covers the S2Cell range (first, last) and
  // adds this cell to "cell_ids".
  //
  // REQUIRES: "first" and "last" have a common ancestor.
  static void coverRange(S2CellId first, S2CellId last, ref S2CellId[] cell_ids) {
    if (first == last) {
      // The range consists of a single index cell.
      cell_ids ~= first;
    } else {
      // Add the lowest common ancestor of the given range.
      int level = first.getCommonAncestorLevel(last);
      enforce(level >= 0);
      cell_ids ~= first.parent(level);
    }
  }

  // Returns true if the indexed shape "clipped" in the indexed cell "id"
  // contains the point "p".
  //
  // REQUIRES: id.contains(S2CellId(p))
  //bool Contains(S2CellId id, const S2ClippedShape& clipped, const S2Point& p) const;

  // Returns true if any edge of the indexed shape "clipped" intersects the
  // cell "target".  It may also return true if an edge is very close to
  // "target"; the maximum error is less than 10 * DBL_EPSILON radians (about
  // 15 nanometers).
  bool anyEdgeIntersects(in S2ClippedShape clipped, in S2Cell target) const {
    enum double kMaxError = FACE_CLIP_ERROR_UV_COORD + INTERSECTS_RECT_ERROR_UV_DIST;
    const R2Rect bound = target.getBoundUV().expanded(kMaxError);
    const int face = target.face();
    const S2Shape shape = index().shape(clipped.shapeId());
    const int num_edges = clipped.numEdges();
    for (int i = 0; i < num_edges; ++i) {
      const auto edge = shape.edge(clipped.edge(i));
      R2Point p0, p1;
      if (clipToPaddedFace(edge.v0, edge.v1, face, kMaxError, p0, p1)
          && intersectsRect(p0, p1, bound)) {
        return true;
      }
    }
    return false;
  }

  // This class is not thread-safe!
  S2ContainsPointQuery!IndexT _containsQuery;

  // Optimization: rather than declaring our own iterator, instead we reuse
  // the iterator declared by S2ContainsPointQuery.  (This improves benchmark
  // times significantly for classes that create a new S2ShapeIndexRegion
  // object on every call to Contains/MayIntersect(S2Cell).
  Iterator _iter;
}

// Returns an S2ShapeIndexRegion that wraps the given S2ShapeIndex.  Note that
// it is efficient to return S2ShapeIndexRegion objects by value.
S2ShapeIndexRegion!IndexT makeS2ShapeIndexRegion(IndexT)(IndexT index) {
  return new S2ShapeIndexRegion!IndexT(index);
}
