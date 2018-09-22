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
//

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2padded_cell;

import s2.r1interval;
import s2.r2rect;
import s2.s2cell_id;
import s2.s2coords : FaceSiTitoXYZ, IJ_TO_POS, INVERT_MASK, SiTitoST, STtoUV, SWAP_MASK,
  UVtoST, POS_TO_IJ, POS_TO_ORIENTATION, STtoIJ;
import s2.s2point;
import std.algorithm : min, max;
import std.math : log2, floor;

// S2PaddedCell represents an S2Cell whose (u,v)-range has been expanded on
// all sides by a given amount of "padding".  Unlike S2Cell, its methods and
// representation are optimized for clipping edges against S2Cell boundaries
// to determine which cells are intersected by a given set of edges.
//
// This class is intended to be copied by value as desired.
class S2PaddedCell {
public:
  // Construct an S2PaddedCell for the given cell id and padding.
  this(in S2CellId id, double padding) {
    _id = id;
    _padding = padding;

    if (_id.isFace()) {
      // Fast path for constructing a top-level face (the most common case).
      double limit = 1.0 + padding;
      _bound = R2Rect(R1Interval(-limit, limit), R1Interval(-limit, limit));
      _middle = R2Rect(R1Interval(-padding, padding), R1Interval(-padding, padding));
      _ijLo[0] = _ijLo[1] = 0;
      _orientation = _id.face() & 1;
      _level = 0;
    } else {
      int[2] ij;
      id.toFaceIJOrientation(ij[0], ij[1], _orientation);
      _level = id.level();
      _bound = S2CellId.IJLevelToBoundUV(ij, _level).expanded(padding);
      int ij_size = S2CellId.getSizeIJ(_level);
      _ijLo[0] = ij[0] & -ij_size;
      _ijLo[1] = ij[1] & -ij_size;
    }
  }

  // Construct the child of "parent" with the given (i,j) index.  The four
  // child cells have indices of (0,0), (0,1), (1,0), (1,1), where the i and j
  // indices correspond to increasing u- and v-values respectively.
  this(S2PaddedCell parent, int i, int j) {
    _padding = parent._padding;
    _bound = parent._bound;
    _level = parent._level + 1;

    // Compute the position and orientation of the child incrementally from the
    // orientation of the parent.
    int pos = IJ_TO_POS[parent._orientation][2 * i + j];
    _id = parent._id.child(pos);
    int ij_size = S2CellId.getSizeIJ(_level);
    _ijLo[0] = parent._ijLo[0] + i * ij_size;
    _ijLo[1] = parent._ijLo[1] + j * ij_size;
    _orientation = parent._orientation ^ POS_TO_ORIENTATION[pos];
    // For each child, one corner of the bound is taken directly from the parent
    // while the diagonally opposite corner is taken from middle().
    R2Rect middle = parent.middle();
    _bound[0][1 - i] = middle[0][1 - i];
    _bound[1][1 - j] = middle[1][1 - j];
  }

  const(S2CellId) id() const {
    return _id;
  }

  double padding() const {
    return _padding;
  }

  int level() const {
    return _level;
  }

  // Return the bound for this cell (including padding).
  const(R2Rect) bound() const {
    return _bound;
  }

  // Return the "middle" of the padded cell, defined as the rectangle that
  // belongs to all four children.
  //
  // Note that this method is *not* thread-safe, because the return value is
  // computed on demand and cached.  (It is expected that this class will be
  // mainly useful in the context of single-threaded recursive algorithms.)
  const(R2Rect) middle() {
    // We compute this field lazily because it is not needed the majority of the
    // time (i.e., for cells where the recursion terminates).
    if (_middle.isEmpty()) {
      int ij_size = S2CellId.getSizeIJ(_level);
      double u = STtoUV(SiTitoST(2 * _ijLo[0] + ij_size));
      double v = STtoUV(SiTitoST(2 * _ijLo[1] + ij_size));
      _middle = R2Rect(R1Interval(u - _padding, u + _padding),
          R1Interval(v - _padding, v + _padding));
    }
    return _middle;
  }

  // Return the (i,j) coordinates for the child cell at the given traversal
  // position.  The traversal position corresponds to the order in which child
  // cells are visited by the Hilbert curve.
  void getChildIJ(int pos, out int i, out int j) const {
    int ij = POS_TO_IJ[_orientation][pos];
    i = ij >> 1;
    j = ij & 1;
  }


  // Return the smallest cell that contains all descendants of this cell whose
  // bounds intersect "rect".  For algorithms that use recursive subdivision
  // to find the cells that intersect a particular object, this method can be
  // used to skip all the initial subdivision steps where only one child needs
  // to be expanded.
  //
  // Note that this method is not the same as returning the smallest cell that
  // contains the intersection of this cell with "rect".  Because of the
  // padding, even if one child completely contains "rect" it is still
  // possible that a neighboring child also intersects "rect".
  //
  // REQUIRES: bound().Intersects(rect)
  S2CellId shrinkToFit(in R2Rect rect) const
  in {
    assert(bound().intersects(rect));
  } body {
    // Quick rejection test: if "rect" contains the center of this cell along
    // either axis, then no further shrinking is possible.
    int ij_size = S2CellId.getSizeIJ(_level);
    if (_level == 0) {
      // Fast path (most calls to this function start with a face cell).
      if (rect[0].contains(0) || rect[1].contains(0)) return id();
    } else {
      if (rect[0].contains(STtoUV(SiTitoST(2 * _ijLo[0] + ij_size))) ||
          rect[1].contains(STtoUV(SiTitoST(2 * _ijLo[1] + ij_size)))) {
        return id();
      }
    }
    // Otherwise we expand "rect" by the given padding() on all sides and find
    // the range of coordinates that it spans along the i- and j-axes.  We then
    // compute the highest bit position at which the min and max coordinates
    // differ.  This corresponds to the first cell level at which at least two
    // children intersect "rect".

    // Increase the padding to compensate for the error in S2::UVtoST().
    // (The constant below is a provable upper bound on the additional error.)
    R2Rect padded = rect.expanded(padding() + 1.5 * double.epsilon);
    int[2] ij_min;  // Min i- or j- coordinate spanned by "padded"
    int[2] ij_xor;  // XOR of the min and max i- or j-coordinates
    for (int d = 0; d < 2; ++d) {
      ij_min[d] = max(_ijLo[d], STtoIJ(UVtoST(padded[d][0])));
      int ij_max = min(_ijLo[d] + ij_size - 1, STtoIJ(UVtoST(padded[d][1])));
      ij_xor[d] = ij_min[d] ^ ij_max;
    }
    // Compute the highest bit position where the two i- or j-endpoints differ,
    // and then choose the cell level that includes both of these endpoints.  So
    // if both pairs of endpoints are equal we choose kMaxLevel; if they differ
    // only at bit 0, we choose (kMaxLevel - 1), and so on.
    int level_msb = ((ij_xor[0] | ij_xor[1]) << 1) + 1;
    int level = S2CellId.MAX_LEVEL - cast(int) floor(log2(level_msb));
    if (level <= _level) return id();
    return S2CellId.fromFaceIJ(id().face(), ij_min[0], ij_min[1]).parent(level);
  }

  // Return the center of this cell.
  S2Point getCenter() const {
    int ij_size = S2CellId.getSizeIJ(_level);
    uint si = 2 * _ijLo[0] + ij_size;
    uint ti = 2 * _ijLo[1] + ij_size;
    return FaceSiTitoXYZ(_id.face(), si, ti).normalize();
  }

  // Return the vertex where the S2 space-filling curve enters this cell.
  S2Point getEntryVertex() const {
    // The curve enters at the (0,0) vertex unless the axis directions are
    // reversed, in which case it enters at the (1,1) vertex.
    uint i = _ijLo[0];
    uint j = _ijLo[1];
    if (_orientation & INVERT_MASK) {
      int ij_size = S2CellId.getSizeIJ(_level);
      i += ij_size;
      j += ij_size;
    }
    return FaceSiTitoXYZ(_id.face(), 2 * i, 2 * j).normalize();
  }

  // Return the vertex where the S2 space-filling curve exits this cell.
  S2Point getExitVertex() const {
    // The curve exits at the (1,0) vertex unless the axes are swapped or
    // inverted but not both, in which case it exits at the (0,1) vertex.
    uint i = _ijLo[0];
    uint j = _ijLo[1];
    int ij_size = S2CellId.getSizeIJ(_level);
    if (_orientation == 0 || _orientation == SWAP_MASK + INVERT_MASK) {
      i += ij_size;
    } else {
      j += ij_size;
    }
    return FaceSiTitoXYZ(_id.face(), 2 * i, 2 * j).normalize();
  }


 private:
  const(S2CellId) _id;
  double _padding;
  R2Rect _bound;     // Bound in (u,v)-space

  // The rectangle in (u,v)-space that belongs to all four padded children.
  // It is computed on demand by the middle() accessor method.
  R2Rect _middle;

  int[2] _ijLo;     // Minimum (i,j)-coordinates of this cell, before padding
  int _orientation;  // Hilbert curve orientation of this cell (see s2coords.h)
  int _level;        // Level of this cell (see s2coords.h)
}
