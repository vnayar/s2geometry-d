// Copyright 2005 Google Inc. All Rights Reserved.
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

module s2.s2cell_id;

import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s2latlng;
import s2.s2point;
import s2.util.coding.coder;
import s2coords = s2.s2coords;

import algorithm = std.algorithm;
import array = std.array;
import bitop = core.bitop;
import conv = std.conv;
import format = std.format;
import math = std.math;
import range = std.range;

/**
 * An S2CellId is a 64-bit unsigned integer that uniquely identifies a
 * cell in the S2 cell decomposition.  It has the following format:
 *
 *   id = [face][face_pos]
 *
 *   face:     a 3-bit number (range 0..5) encoding the cube face.
 *
 *   face_pos: a 61-bit number encoding the position of the center of this
 *             cell along the Hilbert curve over this face (see the Wiki
 *             pages for details).
 *
 * Sequentially increasing cell ids follow a continuous space-filling curve
 * over the entire sphere.  They have the following properties:
 *
 *  - The id of a cell at level k consists of a 3-bit face number followed
 *    by k bit pairs that recursively select one of the four children of
 *    each cell.  The next bit is always 1, and all other bits are 0.
 *    Therefore, the level of a cell is determined by the position of its
 *    lowest-numbered bit that is turned on (for a cell at level k, this
 *    position is 2 * (MAX_LEVEL - k).)
 *
 *  - The id of a parent cell is at the midpoint of the range of ids spanned
 *    by its children (or by its descendants at any level).
 *
 * Leaf cells are often used to represent points on the unit sphere, and
 * this class provides methods for converting directly between these two
 * representations.  For cells that represent 2D regions rather than
 * discrete point, it is better to use the S2Cell class.
 *
 * This class is intended to be copied by value as desired.  It uses
 * the default copy constructor and assignment operator.
 */
struct S2CellId {
private:
  // The default constructor returns an invalid cell id.
  ulong _id = 0;

public:
  // The extra position bit (61 rather than 60) let us encode each cell as its
  // Hilbert curve position at the cell center (which is halfway along the
  // portion of the Hilbert curve that fills that cell).
  static immutable int FACE_BITS = 3;
  static immutable int NUM_FACES = 6;
  static immutable int MAX_LEVEL = s2coords.MAX_CELL_LEVEL;  // Valid levels: 0..MAX_LEVEL
  static immutable int POS_BITS = 2 * MAX_LEVEL + 1;
  static immutable int MAX_SIZE = 1 << MAX_LEVEL;

  this(ulong id) {
    _id = id;
  }

  // Construct a leaf cell containing the given point "p".  Usually there is
  // is exactly one such cell, but for points along the edge of a cell, any
  // adjacent cell may be (deterministically) chosen.  This is because
  // S2CellIds are considered to be closed sets.  The returned cell will
  // always contain the given point, i.e.
  //
  //   S2Cell(S2CellId(p)).Contains(p)
  //
  // is always true.  The point "p" does not need to be normalized.
  //
  // If instead you want every point to be contained by exactly one S2Cell,
  // you will need to convert the S2CellIds to S2Loops (which implement point
  // containment this way).
  this(in S2Point p) {
    double u, v;
    int face = s2coords.XYZtoFaceUV(p, u, v);
    int i = s2coords.STtoIJ(s2coords.UVtoST(u));
    int j = s2coords.STtoIJ(s2coords.UVtoST(v));
    _id = fromFaceIJ(face, i, j).id();
  }

  this(in S2LatLng ll) {
    this(ll.toS2Point());
  }

  static S2CellId none() {
    return S2CellId();
  }

  // Returns an invalid cell id guaranteed to be larger than any
  // valid cell id.  Useful for creating indexes.
  static S2CellId sentinel() {
    return S2CellId(ulong.max);
  }

  // Return the cell corresponding to a given S2 cube face.
  static S2CellId fromFace(int face) {
    return S2CellId((cast(ulong) face << POS_BITS) + lsbForLevel(0));
  }


  // Return a cell given its face (range 0..5), Hilbert curve position within
  // that face (an unsigned integer with S2CellId::POS_BITS bits), and level
  // (range 0..MAX_LEVEL).  The given position will be modified to correspond
  // to the Hilbert curve position at the center of the returned cell.  This
  // is a static function rather than a constructor in order to indicate what
  // the arguments represent.
  static S2CellId fromFacePosLevel(int face, ulong pos, int level) {
    auto cell = S2CellId((cast(ulong) face << POS_BITS) + (pos | 1));
    return cell.parent(level);
  }


  // Return the direction vector corresponding to the center of the given
  // cell.  The vector returned by ToPointRaw is not necessarily unit length.
  // This method returns the same result as S2Cell::GetCenter().
  //
  // The maximum directional error in ToPoint() (compared to the exact
  // mathematical result) is 1.5 * DBL_EPSILON radians, and the maximum length
  // error is 2 * DBL_EPSILON (the same as Normalize).
  S2Point toS2Point() const {
    return toS2PointRaw().normalize();
  }

  S2Point toS2PointRaw() const {
    int si, ti;
    int face = getCenterSiTi(si, ti);
    return s2coords.FaceSiTitoXYZ(face, si, ti);
  }

  // Return the center of the cell in (s,t) coordinates (see s2coords.h).
  R2Point getCenterST() const {
    int si, ti;
    getCenterSiTi(si, ti);
    return R2Point(s2coords.SiTitoST(si), s2coords.SiTitoST(ti));
  }

  // Return the edge length of this cell in (s,t)-space.
  double getSizeST() const {
    return getSizeST(level());
  }


  // Return the edge length in (s,t)-space of cells at the given level.
  static double getSizeST(int level) {
    return s2coords.IJtoSTMin(getSizeIJ(level));
  }

  // Return the bounds of this cell in (s,t)-space.
  R2Rect getBoundST() const {
    double size = getSizeST();
    return R2Rect.fromCenterSize(getCenterST(), R2Point(size, size));
  }

  // Return the center of the cell in (u,v) coordinates (see s2coords.h).
  // Note that the center of the cell is defined as the point at which it is
  // recursively subdivided into four children; in general, it is not at the
  // midpoint of the (u,v) rectangle covered by the cell.
  R2Point getCenterUV() const {
    int si, ti;
    getCenterSiTi(si, ti);
    return R2Point(s2coords.STtoUV(s2coords.SiTitoST(si)), s2coords.STtoUV(s2coords.SiTitoST(ti)));
  }

  // Return the bounds of this cell in (u,v)-space.
  R2Rect getBoundUV() const {
    int[2] ij;
    toFaceIJOrientation(ij[0], ij[1]);
    return IJLevelToBoundUV(ij, level());
  }

  // Expand a rectangle in (u,v)-space so that it contains all points within
  // the given distance of the boundary, and return the smallest such
  // rectangle.  If the distance is negative, then instead shrink this
  // rectangle so that it excludes all points within the given absolute
  // distance of the boundary.
  //
  // Distances are measured *on the sphere*, not in (u,v)-space.  For example,
  // you can use this method to expand the (u,v)-bound of an S2CellId so that
  // it contains all points within 5km of the original cell.  You can then
  // test whether a point lies within the expanded bounds like this:
  //
  //   R2Point uv;
  //   if (S2::FaceXYZtoUV(face, point, &uv) && bound.Contains(uv)) { ... }
  //
  // Limitations:
  //
  //  - Because the rectangle is drawn on one of the six cube-face planes
  //    (i.e., {x,y,z} = +/-1), it can cover at most one hemisphere.  This
  //    limits the maximum amount that a rectangle can be expanded.  For
  //    example, S2CellId bounds can be expanded safely by at most 45 degrees
  //    (about 5000 km on the Earth's surface).
  //
  //  - The implementation is not exact for negative distances.  The resulting
  //    rectangle will exclude all points within the given distance of the
  //    boundary but may be slightly smaller than necessary.
  static R2Rect expandedByDistanceUV(in R2Rect uv, in S1Angle distance) {
    // Expand each of the four sides of the rectangle just enough to include all
    // points within the given distance of that side.  (The rectangle may be
    // expanded by a different amount in (u,v)-space on each side.)
    double u0 = uv[0][0];
    double u1 = uv[0][1];
    double v0 = uv[1][0];
    double v1 = uv[1][1];
    double max_u = algorithm.max(math.abs(u0), math.abs(u1));
    double max_v = algorithm.max(math.abs(v0), math.abs(v1));
    double sin_dist = sin(distance);
    return R2Rect(
        R1Interval(expandEndpoint(u0, max_v, -sin_dist), expandEndpoint(u1, max_v, sin_dist)),
        R1Interval(expandEndpoint(v0, max_u, -sin_dist), expandEndpoint(v1, max_u, sin_dist)));
  }

  // This is a helper function for expandedByDistanceUV().
  //
  // Given an edge of the form (u,v0)-(u,v1), let max_v = max(abs(v0), abs(v1)).
  // This method returns a new u-coordinate u' such that the distance from the
  // line u=u' to the given edge (u,v0)-(u,v1) is exactly the given distance
  // (which is specified as the sine of the angle corresponding to the distance).
  private static double expandEndpoint(double u, double max_v, double sinDist) {
    // This is based on solving a spherical right triangle, similar to the
    // calculation in S2Cap::GetRectBound.
    double sin_u_shift = sinDist * math.sqrt((1 + u * u + max_v * max_v) / (1 + u * u));
    double cos_u_shift = math.sqrt(1 - sin_u_shift * sin_u_shift);
    // The following is an expansion of tan(atan(u) + asin(sin_u_shift)).
    return (cos_u_shift * u + sin_u_shift) / (cos_u_shift - sin_u_shift * u);
  }

  // Return the (face, si, ti) coordinates of the center of the cell.  Note
  // that although (si,ti) coordinates span the range [0,2**31] in general,
  // the cell center coordinates are always in the range [1,2**31-1] and
  // therefore can be represented using a signed 32-bit integer.
  int getCenterSiTi(out int psi, out int pti) const {
    // First we compute the discrete (i,j) coordinates of a leaf cell contained
    // within the given cell.  Given that cells are represented by the Hilbert
    // curve position corresponding at their center, it turns out that the cell
    // returned by ToFaceIJOrientation is always one of two leaf cells closest
    // to the center of the cell (unless the given cell is a leaf cell itself,
    // in which case there is only one possibility).
    //
    // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
    // jmin) be the coordinates of its lower left-hand corner, the leaf cell
    // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
    // (imin + s/2 - 1, jmin + s/2 - 1).  The first case is the one we want.
    // We can distinguish these two cases by looking at the low bit of "i" or
    // "j".  In the second case the low bit is one, unless s == 2 (i.e. the
    // level just above leaf cells) in which case the low bit is zero.
    //
    // In the code below, the expression ((i ^ (int(_id) >> 2)) & 1) is true
    // if we are in the second case described above.
    int i, j;
    int face = toFaceIJOrientation(i, j);
    int delta = isLeaf() ? 1 : ((i ^ (cast(int) _id >> 2)) & 1) ? 2 : 0;

    // Note that (2 * {i,j} + delta) will never overflow a 32-bit integer.
    psi = 2 * i + delta;
    pti = 2 * j + delta;
    return face;
  }

  // Return the S2LatLng corresponding to the center of the given cell.
  S2LatLng toLatLng() const {
    return S2LatLng(toS2PointRaw());
  }

  // The 64-bit unique identifier for this cell.
  @property
  ulong id() const {
    return _id;
  }

  // Return true if id() represents a valid cell.
  //
  // All methods require isValid() to be true unless otherwise specified
  // (although not all methods enforce this).
  bool isValid() const {
    return (face() < NUM_FACES && (lsb() & 0x1555555555555555uL));
  }

  // Which cube face this cell belongs to, in the range 0..5.
  int face() const {
    return _id >> POS_BITS;
  }

  // The position of the cell center along the Hilbert curve over this face,
  // in the range 0..(2**POS_BITS-1).
  ulong pos() const {
    return _id & (~0uL >> FACE_BITS);
  }

  // Return the subdivision level of the cell (range 0..MAX_LEVEL).
  int level() const
  in {
    // We can't just DCHECK(isValid()) because we want level() to be to be
    // defined for end-iterators, i.e. S2CellId::End(kLevel).  However there is
    // no good way to define S2CellId::None().level(), so we do prohibit that.
    assert(_id != 0);
  } do {
    // A special case for leaf cells is not worthwhile.
    return MAX_LEVEL - (bitop.bsf(_id) >> 1);
  }

  // Return the edge length of this cell in (i,j)-space.
  int getSizeIJ() const {
    return getSizeIJ(level());
  }

  // Like the above, but return the size of cells at the given level.
  static int getSizeIJ(int level) {
    return 1 << (MAX_LEVEL - level);
  }

  // Return true if this is a leaf cell (more efficient than checking
  // whether level() == MAX_LEVEL).
  bool isLeaf() const {
    return cast(int) _id & 1;
  }

  // Return true if this is a top-level face cell (more efficient than
  // checking whether level() == 0).
  bool isFace() const {
    return (_id & (lsbForLevel(0) - 1)) == 0;
  }

  // Return the child position (0..3) of this cell within its parent.
  // REQUIRES: level() >= 1.
  int childPosition() const {
    // No need for a special implementation; the compiler optimizes this well.
    return childPosition(level());
  }

  // Return the child position (0..3) of this cell's ancestor at the given
  // level within its parent.  For example, child_position(1) returns the
  // position of this cell's level-1 ancestor within its top-level face cell.
  // REQUIRES: 1 <= level <= this->level().
  int childPosition(int level) const
  in {
    assert(isValid());
    assert(level >=  1);
    assert(level <= this.level());
  } do {
    return cast(int)(_id >> (2 * (MAX_LEVEL - level) + 1)) & 3;
  }

  // These methods return the range of cell ids that are contained within this
  // cell (including itself).  The range is *inclusive* (i.e. test using >=
  // and <=) and the return values of both methods are valid leaf cell ids.
  // In other words, a.contains(b) if and only if
  //
  //     (b >= a.range_min() && b <= a.range_max())
  //
  // If you want to iterate through all the descendants of this cell at a
  // particular level, use child_begin(level) and child_end(level) instead.
  // Also see maximum_tile(), which can be used to iterate through a range of
  // cells using S2CellIds at different levels that are as large as possible.
  //
  // If you need to convert the range to a semi-open interval [min, limit)
  // (e.g., in order to use a key-value store that only supports semi-open
  // range queries), do not attempt to define "limit" as range_max.next().
  // The problem is that leaf S2CellIds are 2 units apart, so the semi-open
  // interval [min, limit) includes an additional value (range_max.id() + 1)
  // which is happens to be a valid S2CellId about one-third of the time and
  // is *never* contained by this cell.  (It always correpsonds to a cell that
  // is larger than this one.)  You can define "limit" as (range_max.id() + 1)
  // if necessary (which is not always a valid S2CellId but can still be used
  // with FromToken/ToToken), or you can convert range_max() to the key space
  // of your key-value store and define "limit" as Successor(key).
  //
  // Note that Sentinel().range_min() == Sentinel.range_max() == Sentinel().
  S2CellId rangeMin() const {
    return S2CellId(_id - (lsb() - 1));
  }

  S2CellId rangeMax() const {
    return S2CellId(_id + (lsb() - 1));
  }

  // Return true if the given cell is contained within this one.
  bool contains(S2CellId other) const
  in {
    assert(isValid(), "Invalid cell: " ~ toString());
    assert(other.isValid());
  } do {
    return other >= rangeMin() && other <= rangeMax();
  }


  // Return true if the given cell intersects this one.
  bool intersects(S2CellId other) const
  in {
    assert(isValid());
    assert(other.isValid());
  } do {
    return other.rangeMin() <= rangeMax() && other.rangeMax() >= rangeMin();
  }

  // Return the cell at the previous level or at the given level (which must
  // be less than or equal to the current level).
  S2CellId parent() const
  in {
    assert(isValid());
    assert(!isFace());
  } do {
    ulong new_lsb = lsb() << 2;
    return S2CellId((_id & (~new_lsb + 1)) | new_lsb);
  }

  S2CellId parent(int level) const
  in {
    assert(isValid(), format.format("Cannot get parent of invalid cell %s.\n  id=%x, face=%d, lsb=%d", toString(), id(), face(), lsb()));
    assert(level >= 0);
    assert(level <= this.level());
  } do {
    ulong new_lsb = lsbForLevel(level);
    return S2CellId((_id & (~new_lsb + 1)) | new_lsb);
  }

  // Return the immediate child of this cell at the given traversal order
  // position (in the range 0 to 3).  This cell must not be a leaf cell.
  S2CellId child(int position) const
  in {
    assert(isValid());
    assert(!isLeaf());
  } do {
    // To change the level, we need to move the least-significant bit two
    // positions downward.  We do this by subtracting (4 * new_lsb) and adding
    // new_lsb.  Then to advance to the given child cell, we add
    // (2 * position * new_lsb).
    ulong new_lsb = lsb() >> 2;
    return S2CellId(_id + (2 * position + 1 - 4) * new_lsb);
  }

  // Iterator-style methods for traversing the immediate children of a cell or
  // all of the children at a given level (greater than or equal to the current
  // level).  Note that the end value is exclusive, just like standard STL
  // iterators, and may not even be a valid cell id.  You should iterate using
  // code like this:
  //
  //   for(S2CellId c = id.child_begin(); c != id.child_end(); c = c.next())
  //     ...
  //
  // The convention for advancing the iterator is "c = c.next()" rather
  // than "++c" to avoid possible confusion with incrementing the
  // underlying 64-bit cell id.
  S2CellId childBegin() const
  in {
    assert(isValid());
    assert(!isLeaf());
  } do {
    ulong old_lsb = lsb();
    return S2CellId(_id - old_lsb + (old_lsb >> 2));
  }

  S2CellId childBegin(int level) const
  in {
    assert(isValid());
    assert(level >= this.level());
    assert(level <= MAX_LEVEL);
  } do {
    return S2CellId(_id - lsb() + lsbForLevel(level));
  }

  S2CellId childEnd() const
  in {
    assert(isValid());
    assert(!isLeaf());
  } do {
    ulong old_lsb = lsb();
    return S2CellId(_id + old_lsb + (old_lsb >> 2));
  }

  S2CellId childEnd(int level) const
  in {
    assert(isValid());
    assert(level >= this.level());
    assert(level <= MAX_LEVEL);
  } do {
    return S2CellId(_id + lsb() + lsbForLevel(level));
  }

  // Return the next/previous cell at the same level along the Hilbert curve.
  // Works correctly when advancing from one face to the next, but
  // does *not* wrap around from the last face to the first or vice versa.
  S2CellId next() const {
    return S2CellId(_id + (lsb() << 1));
  }

  S2CellId prev() const {
    return S2CellId(_id - (lsb() << 1));
  }

  // This method advances or retreats the indicated number of steps along the
  // Hilbert curve at the current level, and returns the new position.  The
  // position is never advanced past End() or before Begin().
  S2CellId advance(long steps) const {
    if (steps == 0) {
      return this;
    }

    // We clamp the number of steps if necessary to ensure that we do not
    // advance past the End() or before the Begin() of this level.  Note that
    // min_steps and max_steps always fit in a signed 64-bit integer.

    int step_shift = 2 * (MAX_LEVEL - level()) + 1;
    if (steps < 0) {
      long min_steps = -cast(long)(_id >> step_shift);
      if (steps < min_steps) steps = min_steps;
    } else {
      long max_steps = (WRAP_OFFSET + lsb() - _id) >> step_shift;
      if (steps > max_steps) steps = max_steps;
    }
    // If steps is negative, then shifting it left has undefined behavior.
    // Cast to uint64 for a 2's complement answer.
    return S2CellId(_id + (cast(ulong) steps << step_shift));
  }

  // Returns the number of steps that this cell is from Begin(level()). The
  // return value is always non-negative.
  long distanceFromBegin() const {
    const int step_shift = 2 * (MAX_LEVEL - level()) + 1;
    return _id >> step_shift;
  }

  // Like next() and prev(), but these methods wrap around from the last face
  // to the first and vice versa.  They should *not* be used for iteration in
  // conjunction with child_begin(), child_end(), Begin(), or End().  The
  // input must be a valid cell id.
  S2CellId nextWrap() const
  in {
    assert(isValid());
  } do {
    S2CellId n = next();
    if (n._id < WRAP_OFFSET) return n;
    return S2CellId(n._id - WRAP_OFFSET);
  }

  S2CellId prevWrap() const
  in {
    assert(isValid());
  } do {
    S2CellId p = prev();
    if (p._id < WRAP_OFFSET) return p;
    return S2CellId(p._id + WRAP_OFFSET);
  }

  // This method advances or retreats the indicated number of steps along the
  // Hilbert curve at the current level, and returns the new position.  The
  // position wraps between the first and last faces as necessary.  The input
  // must be a valid cell id.
  S2CellId advanceWrap(long steps) const
  in {
    assert(isValid());
  } do {
    if (steps == 0) {
      return this;
    }

    int step_shift = 2 * (MAX_LEVEL - level()) + 1;
    if (steps < 0) {
      long min_steps = -cast(long)(_id >> step_shift);
      if (steps < min_steps) {
        long step_wrap = WRAP_OFFSET >> step_shift;
        steps %= step_wrap;
        if (steps < min_steps) steps += step_wrap;
      }
    } else {
      // Unlike advance(), we don't want to return End(level).
      long max_steps = (WRAP_OFFSET - _id) >> step_shift;
      if (steps > max_steps) {
        long step_wrap = WRAP_OFFSET >> step_shift;
        steps %= step_wrap;
        if (steps > max_steps) steps -= step_wrap;
      }
    }
    return S2CellId(_id + (cast(ulong)(steps) << step_shift));
  }

  // Return the largest cell with the same range_min() and such that
  // range_max() < limit.range_min().  Returns "limit" if no such cell exists.
  // This method can be used to generate a small set of S2CellIds that covers
  // a given range (a "tiling").  This example shows how to generate a tiling
  // for a semi-open range of leaf cells [start, limit):
  //
  //   for (S2CellId id = start.maximum_tile(limit);
  //        id != limit; id = id.next().maximum_tile(limit)) { ... }
  //
  // Note that in general the cells in the tiling will be of different sizes;
  // they gradually get larger (near the middle of the range) and then
  // gradually get smaller (as "limit" is approached).
  S2CellId maximumTile(S2CellId limit) const {
    S2CellId id = this;
    S2CellId start = id.rangeMin();
    if (start >= limit.rangeMin()) {
      return limit;
    }

    if (id.rangeMax() >= limit) {
      // The cell is too large.  Shrink it.  Note that when generating coverings
      // of S2CellId ranges, this loop usually executes only once.  Also because
      // id.range_min() < limit.range_min(), we will always exit the loop by the
      // time we reach a leaf cell.
      do {
        id = id.child(0);
      } while (id.rangeMax() >= limit);
      return id;
    }
    // The cell may be too small.  Grow it if necessary.  Note that generally
    // this loop only iterates once.
    while (!id.isFace()) {
      S2CellId parent = id.parent();
      if (parent.rangeMin() != start || parent.rangeMax() >= limit) {
        break;
      }
      id = parent;
    }
    return id;
  }

  // Returns the level of the lowest common ancestor of this cell and "other",
  // that is, the maximum level such that parent(level) == other.parent(level).
  // Returns -1 if the two cells do not have any common ancestor (i.e., they
  // are from different faces).
  int getCommonAncestorLevel(S2CellId other) const {
    // Basically we find the first bit position at which the two S2CellIds
    // differ and convert that to a level.  The max() below is necessary for the
    // case where one S2CellId is a descendant of the other.
    ulong bits = algorithm.max(id() ^ other.id(), algorithm.max(lsb(), other.lsb()));
    assert(bits != 0);  // Because lsb() is non-zero.

    // Compute the position of the most significant bit, and then map the bit
    // position as follows:
    // {0} -> 30, {1,2} -> 29, {3,4} -> 28, ... , {59,60} -> 0, {61,62,63} -> -1.
    return algorithm.max(60 - bitop.bsr(bits), -1) >> 1;
  }

  // Iterator-style methods for traversing all the cells along the Hilbert
  // curve at a given level (across all 6 faces of the cube).  Note that the
  // end value is exclusive (just like standard STL iterators), and is not a
  // valid cell id.
  static S2CellId begin(int level) {
    return fromFace(0).childBegin(level);
  }

  static S2CellId end(int level) {
    return fromFace(5).childEnd(level);
  }

  // Methods to encode and decode cell ids to compact text strings suitable
  // for display or indexing.  Cells at lower levels (i.e. larger cells) are
  // encoded into fewer characters.  The maximum token length is 16.
  //
  // Tokens preserve ordering, i.e. ToToken(x) < ToToken(y) iff x < y.
  //
  // ToToken() returns a string by value for convenience; the compiler
  // does this without intermediate copying in most cases.
  //
  // These methods guarantee that FromToken(ToToken(x)) == x even when
  // "x" is an invalid cell id.  All tokens are alphanumeric strings.
  // FromToken() returns S2CellId::None() for malformed inputs.
  string toToken() const {
    // Simple implementation: print the id in hex without trailing zeros.
    // Using hex has the advantage that the tokens are case-insensitive, all
    // characters are alphanumeric, no characters require any special escaping
    // in queries for most indexing systems, and it's easy to compare cell
    // tokens against the feature ids of the corresponding features.
    //
    // Using base 64 would produce slightly shorter tokens, but for typical cell
    // sizes used during indexing (up to level 15 or so) the average savings
    // would be less than 2 bytes per cell which doesn't seem worth it.

    // "0" with trailing 0s stripped is the empty string, which is not a
    // reasonable token.  Encode as "X".
    if (_id == 0) {
      return "X";
    }
    const size_t num_zero_digits = bitop.bsf(_id) / 4;
    return hexFormatString(_id >> (4 * num_zero_digits), 16 - num_zero_digits);
  }

  // Print the num_digits low order hex digits.
  private static string hexFormatString(ulong val, size_t num_digits) {
    char[] result = array.replicate([' '], num_digits);
    for (; num_digits--; val >>= 4) {
      result[num_digits] = "0123456789abcdef"[val & 0xF];
    }
    return result.idup;
  }

  static S2CellId fromToken(string token) {
    if (token.length > 16) {
      return S2CellId.none();
    }
    ulong id = 0;
    for (int i = 0, pos = 60; i < token.length; ++i, pos -= 4) {
      ulong d;
      if ('0' <= token[i] && token[i] <= '9') {
        d = token[i] - '0';
      } else if ('a' <= token[i] && token[i] <= 'f') {
        d = token[i] - 'a' + 10;
      } else if ('A' <= token[i] && token[i] <= 'F') {
        d = token[i] - 'A' + 10;
      } else {
        return S2CellId.none();
      }
      id |= d << pos;
    }
    return S2CellId(id);
  }

  /**
   * Use encoder to generate a serialized representation of this cell id.
   * Can also encode an invalid cell.
   */
  void encode(ORangeT)(Encoder!ORangeT encoder) const {
    encoder.ensure(ulong.sizeof);  // A single uint64.
    encoder.put64(_id);
  }

  /// Decodes an S2CellId encoded by Encode(). Returns true on success.
  bool decode(IRangeT)(Decoder!IRangeT decoder) {
    if (decoder.avail() < ulong.sizeof) return false;
    _id = decoder.get64();
    return true;
  }

  // Creates a human readable debug string.  Used for << and available for
  // direct usage as well.  The format is "f/dd..d" where "f" is a digit in
  // the range [0-5] representing the S2CellId face, and "dd..d" is a string
  // of digits in the range [0-3] representing each child's position with
  // respect to its parent.  (Note that the latter string may be empty.)
  //
  // For example "4/" represents S2CellId::FromFace(4), and "3/02" represents
  // S2CellId::FromFace(3).child(0).child(2).
  string toString() const {
    if (!isValid()) {
      return "Invalid: " ~ format.format("%016x", id());
    }
    string s = conv.to!string(face()) ~ "/";
    for (int current_level = 1; current_level <= level(); ++current_level) {
      // Avoid dependencies of SimpleItoA, and slowness of StrAppend &
      // std::to_string.
      s ~= "0123"[childPosition(current_level)];
    }
    return s;
  }

  // Converts a string in the format returned by ToString() to an S2CellId.
  // Returns S2CellId::None() if the string could not be parsed.
  //
  // The method name includes "Debug" in order to avoid possible confusion
  // with FromToken() above.
  static S2CellId fromDebugString(in string str) {
    // This function is reasonably efficient, but is only intended for use in
    // tests.
    int level = cast(int)(str.length - 2);
    if (level < 0 || level > MAX_LEVEL) {
      return S2CellId.none();
    }
    int face = str[0] - '0';
    if (face < 0 || face > 5 || str[1] != '/') {
      return S2CellId.none();
    }
    S2CellId id = S2CellId.fromFace(face);
    for (int i = 2; i < str.length; ++i) {
      int child_pos = str[i] - '0';
      if (child_pos < 0 || child_pos > 3) {
        return S2CellId.none();
      }
      id = id.child(child_pos);
    }
    return id;
  }

  // Return the four cells that are adjacent across the cell's four edges.
  // Neighbors are returned in the order defined by S2Cell::GetEdge.  All
  // neighbors are guaranteed to be distinct.
  void getEdgeNeighbors(out S2CellId[4] neighbors) const {
    int i, j;
    int level = level();
    int size = getSizeIJ(level);
    int face = toFaceIJOrientation(i, j);

    // Edges 0, 1, 2, 3 are in the down, right, up, left directions.
    neighbors[0] = fromFaceIJSame(face, i, j - size, j - size >= 0).parent(level);
    neighbors[1] = fromFaceIJSame(face, i + size, j, i + size < MAX_SIZE).parent(level);
    neighbors[2] = fromFaceIJSame(face, i, j + size, j + size < MAX_SIZE).parent(level);
    neighbors[3] = fromFaceIJSame(face, i - size, j, i - size >= 0).parent(level);
  }

  // Return the neighbors of closest vertex to this cell at the given level,
  // by appending them to "output".  Normally there are four neighbors, but
  // the closest vertex may only have three neighbors if it is one of the 8
  // cube vertices.
  //
  // Requires: level < this->level(), so that we can determine which vertex is
  // closest (in particular, level == MAX_LEVEL is not allowed).
  void appendVertexNeighbors(RangeT)(int level, ref RangeT output) const
  if (range.isOutputRange!(RangeT, S2CellId))
  in {
    // "level" must be strictly less than this cell's level so that we can
    // determine which vertex this cell is closest to.
    assert(level < this.level());
  } do {
    int i, j;
    int face = toFaceIJOrientation(i, j);

    // Determine the i- and j-offsets to the closest neighboring cell in each
    // direction.  This involves looking at the next bit of "i" and "j" to
    // determine which quadrant of this->parent(level) this cell lies in.
    int halfsize = getSizeIJ(level + 1);
    int size = halfsize << 1;
    bool isame, jsame;
    int ioffset, joffset;
    if (i & halfsize) {
      ioffset = size;
      isame = (i + size) < MAX_SIZE;
    } else {
      ioffset = -size;
      isame = (i - size) >= 0;
    }
    if (j & halfsize) {
      joffset = size;
      jsame = (j + size) < MAX_SIZE;
    } else {
      joffset = -size;
      jsame = (j - size) >= 0;
    }

    output ~= parent(level);
    output ~= fromFaceIJSame(face, i + ioffset, j, isame).parent(level);
    output ~= fromFaceIJSame(face, i, j + joffset, jsame).parent(level);
    // If i- and j- edge neighbors are *both* on a different face, then this
    // vertex only has three neighbors (it is one of the 8 cube vertices).
    if (isame || jsame) {
      output ~= fromFaceIJSame(face, i + ioffset, j + joffset, isame && jsame).parent(level);
    }
  }

  // Append all neighbors of this cell at the given level to "output".  Two
  // cells X and Y are neighbors if their boundaries intersect but their
  // interiors do not.  In particular, two cells that intersect at a single
  // point are neighbors.  Note that for cells adjacent to a face vertex, the
  // same neighbor may be appended more than once.
  //
  // REQUIRES: nbr_level >= this->level().
  void appendAllNeighbors(RangeT)(int nbr_level, ref RangeT output) const
  if (range.isOutputRange!(RangeT, S2CellId))
  in {
    assert(nbr_level >= level());
  } do {
    int i, j;
    int face = toFaceIJOrientation(i, j);

    // Find the coordinates of the lower left-hand leaf cell.  We need to
    // normalize (i,j) to a known position within the cell because nbr_level
    // may be larger than this cell's level.
    int size = getSizeIJ();
    i &= -size;
    j &= -size;

    int nbr_size = getSizeIJ(nbr_level);
    assert(nbr_size <= size);

    // We compute the top-bottom, left-right, and diagonal neighbors in one
    // pass.  The loop test is at the end of the loop to avoid 32-bit overflow.
    for (int k = -nbr_size; ; k += nbr_size) {
      bool same_face;
      if (k < 0) {
        same_face = (j + k >= 0);
      } else if (k >= size) {
        same_face = (j + k < MAX_SIZE);
      } else {
        same_face = true;
        // Top and bottom neighbors.
        output ~= fromFaceIJSame(face, i + k, j - nbr_size, j - size >= 0).parent(nbr_level);
        output ~= fromFaceIJSame(face, i + k, j + size, j + size < MAX_SIZE).parent(nbr_level);
      }
      // Left, right, and diagonal neighbors.
      output ~= fromFaceIJSame(face, i - nbr_size, j + k, same_face && i - size >= 0)
          .parent(nbr_level);
      output ~= fromFaceIJSame(face, i + size, j + k, same_face && i + size < MAX_SIZE)
          .parent(nbr_level);
      if (k >= size) {
        break;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  // Low-level methods.

  // Return a leaf cell given its cube face (range 0..5) and
  // i- and j-coordinates (see s2coords.h).
  static S2CellId fromFaceIJ(int face, int i, int j) {
    // Optimization notes:
    //  - Non-overlapping bit fields can be combined with either "+" or "|".
    //    Generally "+" seems to produce better code, but not always.

    // Note that this value gets shifted one bit to the left at the end
    // of the function.
    ulong n = cast(ulong)(face) << (POS_BITS - 1);

    // Alternating faces have opposite Hilbert curve orientations; this
    // is necessary in order for all faces to have a right-handed
    // coordinate system.
    ulong bits = (face & s2coords.SWAP_MASK);

    // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
    // curve position.  The lookup table transforms a 10-bit key of the form
    // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
    // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
    // Hilbert curve orientation respectively.
    int mask;
    static foreach (k; range.iota(7, -1, -1)) {
      mask = (1 << LOOKUP_BITS) - 1;
      bits += ((i >> (k * LOOKUP_BITS)) & mask) << (LOOKUP_BITS + 2);
      bits += ((j >> (k * LOOKUP_BITS)) & mask) << 2;
      bits = LOOKUP_POS[bits];
      n |= (bits >> 2) << (k * 2 * LOOKUP_BITS);
      bits &= (s2coords.SWAP_MASK | s2coords.INVERT_MASK);
    }

    return S2CellId(n * 2 + 1);
  }

  // Return the (face, i, j) coordinates for the leaf cell corresponding to
  // this cell id.  Since cells are represented by the Hilbert curve position
  // at the center of the cell, the returned (i,j) for non-leaf cells will be
  // a leaf cell adjacent to the cell center.  If "orientation" is non-nullptr,
  // also return the Hilbert curve orientation for the current cell.
  int toFaceIJOrientation(out int pi, out int pj) const {
    extractHilbertIJBits(pi, pj);
    return face();
  }

  int toFaceIJOrientation(out int pi, out int pj, out int orientation) const {
    int bits = extractHilbertIJBits(pi, pj);

    // The position of a non-leaf cell at level "n" consists of a prefix of
    // 2*n bits that identifies the cell, followed by a suffix of
    // 2*(kMaxLevel-n)+1 bits of the form 10*.  If n==kMaxLevel, the suffix is
    // just "1" and has no effect.  Otherwise, it consists of "10", followed
    // by (kMaxLevel-n-1) repetitions of "00", followed by "0".  The "10" has
    // no effect, while each occurrence of "00" has the effect of reversing
    // the kSwapMask bit.
    assert(0 == s2coords.POS_TO_ORIENTATION[2]);
    assert(s2coords.SWAP_MASK == s2coords.POS_TO_ORIENTATION[0]);
    if (lsb() & 0x1111111111111110uL) {
      bits ^= s2coords.SWAP_MASK;
    }
    orientation = bits;

    return face();
  }

  private int extractHilbertIJBits(out int pi, out int pj) const {
    int i = 0, j = 0;
    int bits = (face() & s2coords.SWAP_MASK);

    // Each iteration maps 8 bits of the Hilbert curve position into
    // 4 bits of "i" and "j".  The lookup table transforms a key of the
    // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
    // letters [ijpo] represents bits of "i", "j", the Hilbert curve
    // position, and the Hilbert curve orientation respectively.
    //
    // On the first iteration we need to be careful to clear out the bits
    // representing the cube face.
    int nbits;
    static foreach (k; range.iota(7, -1, -1)) {
      nbits = (k == 7) ? (MAX_LEVEL - 7 * LOOKUP_BITS) : LOOKUP_BITS;
      bits += (cast(int)(_id >> (k * 2 * LOOKUP_BITS + 1)) & ((1 << (2 * nbits)) - 1)) << 2;
      bits = LOOKUP_IJ[bits];
      i += (bits >> (LOOKUP_BITS + 2)) << (k * LOOKUP_BITS);
      j += ((bits >> 2) & ((1 << LOOKUP_BITS) - 1)) << (k * LOOKUP_BITS);
      bits &= (s2coords.SWAP_MASK | s2coords.INVERT_MASK);
    }

    pi = i;
    pj = j;
    return bits;
  }

  // Return the lowest-numbered bit that is on for this cell id, which is
  // equal to (uint64(1) << (2 * (MAX_LEVEL - level))).  So for example,
  // a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the
  // first test is more efficient.
  ulong lsb() const {
    return _id & (~_id + 1);
  }

  // Return the lowest-numbered bit that is on for cells at the given level.
  static ulong lsbForLevel(int level) {
    return 1uL << (2 * (MAX_LEVEL - level));
  }

  // Return the bound in (u,v)-space for the cell at the given level containing
  // the leaf cell with the given (i,j)-coordinates.
  static R2Rect IJLevelToBoundUV(int[2] ij, int level) {
    R2Rect bound;
    int cell_size = getSizeIJ(level);
    for (int d = 0; d < 2; ++d) {
      int ij_lo = ij[d] & -cell_size;
      int ij_hi = ij_lo + cell_size;

      bound[d][0] = s2coords.STtoUV(s2coords.IJtoSTMin(ij_lo));
      bound[d][1] = s2coords.STtoUV(s2coords.IJtoSTMin(ij_hi));
    }
    return bound;
  }

  // Support the == and != operators.
  bool opEquals(in S2CellId v) const {
    return id() == v.id();
  }

  // Support <, <=, >, and >= operators.
  int opCmp(in S2CellId v) const {
    if (id() == v.id()) {
      return 0;
    }
    return id() > v.id() ? 1 : -1;
  }

  // The ID serves as a good hash for associative arrays.
  size_t toHash() const @safe pure nothrow {
    return cast(size_t) _id;
  }

 private:
  // This is the offset required to wrap around from the beginning of the
  // Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
  static immutable ulong WRAP_OFFSET = cast(ulong) NUM_FACES << POS_BITS;
  static immutable int LOOKUP_BITS = 4;
  static immutable int LOOKUP_ARRAY_SIZE = 1 << (2 * LOOKUP_BITS + 2);
  static immutable ushort[LOOKUP_ARRAY_SIZE] LOOKUP_POS;
  static immutable ushort[LOOKUP_ARRAY_SIZE] LOOKUP_IJ;

  static this() {
    ushort[LOOKUP_ARRAY_SIZE] lookupPos;
    ushort[LOOKUP_ARRAY_SIZE] lookupIJ;
    initLookupCell(lookupPos, lookupIJ, 0, 0, 0, 0, 0, 0);
    initLookupCell(lookupPos, lookupIJ, 0, 0, 0, s2coords.SWAP_MASK, 0, s2coords.SWAP_MASK);
    initLookupCell(lookupPos, lookupIJ, 0, 0, 0, s2coords.INVERT_MASK, 0, s2coords.INVERT_MASK);
    initLookupCell(lookupPos, lookupIJ, 0, 0, 0, s2coords.SWAP_MASK | s2coords.INVERT_MASK, 0,
        s2coords.SWAP_MASK | s2coords.INVERT_MASK);
    LOOKUP_POS = lookupPos;
    LOOKUP_IJ = lookupIJ;
  }

  // Given a face and a point (i,j) where either i or j is outside the valid
  // range [0..MAX_SIZE-1], this function first determines which neighboring
  // face "contains" (i,j), and then returns the leaf cell on that face which
  // is adjacent to the given face and whose distance from (i,j) is minimal.
  static S2CellId fromFaceIJWrap(int face, int i, int j) {
    // Convert i and j to the coordinates of a leaf cell just beyond the
    // boundary of this face.  This prevents 32-bit overflow in the case
    // of finding the neighbors of a face cell.
    i = algorithm.max(-1, algorithm.min(MAX_SIZE, i));
    j = algorithm.max(-1, algorithm.min(MAX_SIZE, j));

    // We want to wrap these coordinates onto the appropriate adjacent face.
    // The easiest way to do this is to convert the (i,j) coordinates to (x,y,z)
    // (which yields a point outside the normal face boundary), and then call
    // S2::XYZtoFaceUV() to project back onto the correct face.
    //
    // The code below converts (i,j) to (si,ti), and then (si,ti) to (u,v) using
    // the linear projection (u=2*s-1 and v=2*t-1).  (The code further below
    // converts back using the inverse projection, s=0.5*(u+1) and t=0.5*(v+1).
    // Any projection would work here, so we use the simplest.)  We also clamp
    // the (u,v) coordinates so that the point is barely outside the
    // [-1,1]x[-1,1] face rectangle, since otherwise the reprojection step
    // (which divides by the new z coordinate) might change the other
    // coordinates enough so that we end up in the wrong leaf cell.
    static const double SCALE = 1.0 / MAX_SIZE;
    static const double LIMIT = 1.0 + double.epsilon;
    // The arithmetic below is designed to avoid 32-bit integer overflows.
    assert(0 == MAX_SIZE % 2);
    double u = algorithm.max(-LIMIT, algorithm.min(LIMIT, SCALE * (2 * (i - MAX_SIZE / 2) + 1)));
    double v = algorithm.max(-LIMIT, algorithm.min(LIMIT, SCALE * (2 * (j - MAX_SIZE / 2) + 1)));

    // Find the leaf cell coordinates on the adjacent face, and convert
    // them to a cell id at the appropriate level.
    face = s2coords.XYZtoFaceUV(s2coords.FaceUVtoXYZ(face, u, v), u, v);
    return fromFaceIJ(face, s2coords.STtoIJ(0.5*(u+1)), s2coords.STtoIJ(0.5*(v+1)));
  }

  // Inline helper function that calls fromFaceIJ if "same_face" is true,
  // or fromFaceIJWrap if "same_face" is false.
  static S2CellId fromFaceIJSame(int face, int i, int j, bool same_face) {
    if (same_face)
      return S2CellId.fromFaceIJ(face, i, j);
    else
      return S2CellId.fromFaceIJWrap(face, i, j);
  }

  static void initLookupCell(
      ref ushort[LOOKUP_ARRAY_SIZE] lookupPos, ref ushort[LOOKUP_ARRAY_SIZE] lookupIJ,
      int level, int i, int j, int orig_orientation, int pos, int orientation) {
    if (level == LOOKUP_BITS) {
      int ij = (i << LOOKUP_BITS) + j;
      lookupPos[(ij << 2) + orig_orientation] = cast(ushort)((pos << 2) + orientation);
      lookupIJ[(pos << 2) + orig_orientation] = cast(ushort)((ij << 2) + orientation);
    } else {
      level++;
      i <<= 1;
      j <<= 1;
      pos <<= 2;
      const int[4] r = s2coords.POS_TO_IJ[orientation];
      initLookupCell(
          lookupPos, lookupIJ,
          level, i + (r[0] >> 1), j + (r[0] & 1), orig_orientation,
          pos, orientation ^ s2coords.POS_TO_ORIENTATION[0]);
      initLookupCell(
          lookupPos, lookupIJ,
          level, i + (r[1] >> 1), j + (r[1] & 1), orig_orientation,
          pos + 1, orientation ^ s2coords.POS_TO_ORIENTATION[1]);
      initLookupCell(
          lookupPos, lookupIJ,
          level, i + (r[2] >> 1), j + (r[2] & 1), orig_orientation,
          pos + 2, orientation ^ s2coords.POS_TO_ORIENTATION[2]);
      initLookupCell(
          lookupPos, lookupIJ,
          level, i + (r[3] >> 1), j + (r[3] & 1), orig_orientation,
          pos + 3, orientation ^ s2coords.POS_TO_ORIENTATION[3]);
    }
  }
}
