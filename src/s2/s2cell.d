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

module s2.s2cell;

import algorithm = std.algorithm;
import math = std.math;
import s2.r2point;
import s2.r2rect;
import s2.s1chordangle;
import s2.s2cap;
import s2.s2cellid;
import s2.r1interval;
import s2.s1interval;
import s2.s2edgecrosser;
import s2.s2latlng;
import s2.s2latlngrect;
import s2.s2point;
import s2.s2region;
import s2.util.math.vector;
import s2coords = s2.s2coords;
import s2edgedistances = s2.s2edgedistances;
import s2measures = s2.s2measures;
import s2metrics = s2.s2metrics;
import std.exception;

// An S2Cell is an S2Region object that represents a cell.  Unlike S2CellIds,
// it supports efficient containment and intersection tests.  However, it is
// also a more expensive representation (currently 48 bytes rather than 8).

// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator, however it is
// not a "plain old datatype" (POD) because it has virtual functions.
class S2Cell : S2Region {
public:
  // The default constructor is required in order to use freelists.
  // Cells should otherwise always be constructed explicitly.
  this() {}

  // An S2Cell always corresponds to a particular S2CellId.  The other
  // constructors are just convenience methods.
  this(in S2CellId id) {
    _id = id;
    int[2] ij;
    int orientation;
    _face = cast(byte) id.toFaceIJOrientation(ij[0], ij[1], orientation);
    _orientation = cast(byte) orientation;  // Compress int to a byte.
    _level = cast(byte) id.level();
    _uv = S2CellId.IJLevelToBoundUV(ij, _level);
  }

  this(in S2Cell cell) {
    _face = cell._face;
    _level = cell._level;
    _orientation = cell._orientation;
    _id = cell._id;
    _uv = cell._uv;
  }

  // Convenience constructors.  The S2LatLng must be normalized.
  this(in S2Point p) {
    this(S2CellId(p));
  }

  this(in S2LatLng ll) {
    this(S2CellId(ll));
  }

  // Returns the cell corresponding to the given S2 cube face.
  static S2Cell fromFace(int face) {
    return new S2Cell(S2CellId.fromFace(face));
  }

  // Returns a cell given its face (range 0..5), Hilbert curve position within
  // that face (an unsigned integer with S2CellId::kPosBits bits), and level
  // (range 0..kMaxLevel).  The given position will be modified to correspond
  // to the Hilbert curve position at the center of the returned cell.  This
  // is a static function rather than a constructor in order to indicate what
  // the arguments represent.
  static S2Cell fromFacePosLevel(int face, ulong pos, int level) {
    return new S2Cell(S2CellId.fromFacePosLevel(face, pos, level));
  }

  S2CellId id() const {
    return _id;
  }

  int face() const {
    return _face;
  }

  int level() const {
    return _level;
  }

  int orientation() const {
    return _orientation;
  }

  bool isLeaf() const {
    return _level == S2CellId.MAX_LEVEL;
  }

  // These are equivalent to the S2CellId methods, but have a more efficient
  // implementation since the level has been precomputed.
  int getSizeIJ() const {
    return S2CellId.getSizeIJ(level());
  }

  double getSizeST() const {
    return S2CellId.getSizeST(level());
  }

  // Returns the k-th vertex of the cell (k = 0,1,2,3).  Vertices are returned
  // in CCW order (lower left, lower right, upper right, upper left in the UV
  // plane).  The points returned by GetVertexRaw are not normalized.
  // For convenience, the argument is reduced modulo 4 to the range [0..3].
  S2Point getVertex(int k) const {
    return getVertexRaw(k).normalize();
  }

  S2Point getVertexRaw(int k) const {
    return s2coords.FaceUVtoXYZ(_face, _uv.getVertex(k));
  }

  // Returns the inward-facing normal of the great circle passing through the
  // edge from vertex k to vertex k+1 (mod 4).  The normals returned by
  // GetEdgeRaw are not necessarily unit length.  For convenience, the
  // argument is reduced modulo 4 to the range [0..3].
  S2Point getEdge(int k) const {
    return getEdgeRaw(k).normalize();
  }

  S2Point getEdgeRaw(int k) const {
    switch (k & 3) {
      case 0:  return s2coords.GetVNorm(_face, _uv[1][0]);   // Bottom
      case 1:  return s2coords.GetUNorm(_face, _uv[0][1]);   // Right
      case 2:  return -s2coords.GetVNorm(_face, _uv[1][1]);  // Top
      default: return -s2coords.GetUNorm(_face, _uv[0][0]);  // Left
    }
  }

  // If this is not a leaf cell, sets children[0..3] to the four children of
  // this cell (in traversal order) and return true.  Otherwise returns false.
  // This method is equivalent to the following:
  //
  // for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
  //   children[pos] = S2Cell(id);
  //
  // except that it is more than two times faster.
  bool subdivide(out S2Cell[4] children) const {
    // This function is equivalent to just iterating over the child cell ids
    // and calling the S2Cell constructor, but it is about 2.5 times faster.

    if (_id.isLeaf()) {
      return false;
    }

    // Compute the cell midpoint in uv-space.
    R2Point uv_mid = _id.getCenterUV();

    // Create four children with the appropriate bounds.
    S2CellId id = _id.childBegin();
    for (int pos = 0; pos < 4; ++pos, id = id.next()) {
      S2Cell child = new S2Cell();
      child._face = _face;
      child._level = cast(byte)(_level + 1);
      child._orientation = cast(byte)(_orientation ^ s2coords.POS_TO_ORIENTATION[pos]);
      child._id = id;
      // We want to split the cell in half in "u" and "v".  To decide which
      // side to set equal to the midpoint value, we look at cell's (i,j)
      // position within its parent.  The index for "i" is in bit 1 of ij.
      int ij = s2coords.POS_TO_IJ[_orientation][pos];
      int i = ij >> 1;
      int j = ij & 1;
      child._uv[0][i] = _uv[0][i];
      child._uv[0][1-i] = uv_mid[0];
      child._uv[1][j] = _uv[1][j];
      child._uv[1][1-j] = uv_mid[1];
      children[pos] = child;
    }
    return true;
  }

  // Returns the direction vector corresponding to the center in (s,t)-space of
  // the given cell.  This is the point at which the cell is divided into four
  // subcells; it is not necessarily the centroid of the cell in (u,v)-space
  // or (x,y,z)-space.  The point returned by GetCenterRaw is not necessarily
  // unit length.
  S2Point getCenter() const {
    return getCenterRaw().normalize();
  }

  S2Point getCenterRaw() const {
    return _id.toS2PointRaw();
  }

  // Returns the average area for cells at the given level.
  static double averageArea(int level) {
    return s2metrics.AVG_AREA.getValue(level);
  }

  // Returns the average area of cells at this level in steradians.  This is
  // accurate to within a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is
  // extremely cheap to compute.
  double averageArea() const {
    return averageArea(_level);
  }

  // Returns the approximate area of this cell in steradians.  This method is
  // accurate to within 3% percent for all cell sizes and accurate to within
  // 0.1% for cells at level 5 or higher (i.e. squares 350km to a side or
  // smaller on the Earth's surface).  It is moderately cheap to compute.
  double approxArea() const {
    // All cells at the first two levels have the same area.
    if (_level < 2) {
      return averageArea(_level);
    }

    // First, compute the approximate area of the cell when projected
    // perpendicular to its normal.  The cross product of its diagonals gives
    // the normal, and the length of the normal is twice the projected area.
    double flat_area =
        0.5 * (getVertex(2) - getVertex(0)).crossProd(getVertex(3) - getVertex(1)).norm();

    // Now, compensate for the curvature of the cell surface by pretending
    // that the cell is shaped like a spherical cap.  The ratio of the
    // area of a spherical cap to the area of its projected disc turns out
    // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
    // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
    // Here we set Pi*r*r == flat_area to find the equivalent disc.
    return flat_area * 2 / (1 + math.sqrt(1.0 - algorithm.min(math.M_1_PI * flat_area, 1.0)));
  }

  // Returns the area of this cell as accurately as possible.  This method is
  // more expensive but it is accurate to 6 digits of precision even for leaf
  // cells (whose area is approximately 1e-18).
  double exactArea() const {
    // There is a straightforward mathematical formula for the exact surface
    // area (based on 4 calls to asin), but as the cell size gets small this
    // formula has too much cancellation error.  So instead we compute the area
    // as the sum of two triangles (which is very accurate at all cell levels).
    S2Point v0 = getVertex(0);
    S2Point v1 = getVertex(1);
    S2Point v2 = getVertex(2);
    S2Point v3 = getVertex(3);
    return s2measures.area(v0, v1, v2) + s2measures.area(v0, v2, v3);
  }

  // Returns the bounds of this cell in (u,v)-space.
  R2Rect getBoundUV() const {
    return _uv;
  }

  // Returns the distance from the cell to the given point.  Returns zero if
  // the point is inside the cell.
  S1ChordAngle getDistance(in S2Point target) const {
    return getDistanceInternal(target, true /*to_interior*/);
  }

  // Return the distance from the cell boundary to the given point.
  S1ChordAngle getBoundaryDistance(in S2Point target) const {
    return getDistanceInternal(target, false /*to_interior*/);
  }

  // Returns the maximum distance from the cell (including its interior) to the
  // given point.
  S1ChordAngle getMaxDistance(in S2Point target) const {
    // First check the 4 cell vertices.  If all are within the hemisphere
    // centered around target, the max distance will be to one of these vertices.
    S2Point target_uvw = s2coords.FaceXYZtoUVW(_face, target);
    S1ChordAngle max_dist = algorithm.max(
        algorithm.max(vertexChordDist(target_uvw, 0, 0), vertexChordDist(target_uvw, 1, 0)),
        algorithm.max(vertexChordDist(target_uvw, 0, 1), vertexChordDist(target_uvw, 1, 1)));

    if (max_dist <= S1ChordAngle.right()) {
      return max_dist;
    }

    // Otherwise, find the minimum distance d_min to the antipodal point and the
    // maximum distance will be Pi - d_min.
    return S1ChordAngle.straight() - getDistance(-target);
  }


  // Returns the minimum distance from the cell to the given edge AB.  Returns
  // zero if the edge intersects the cell interior.
  S1ChordAngle getDistance(in S2Point a, in S2Point b) const {
    // Possible optimizations:
    //  - Currently the (cell vertex, edge endpoint) distances are computed
    //    twice each, and the length of AB is computed 4 times.
    //  - To fix this, refactor GetDistance(target) so that it skips calculating
    //    the distance to each cell vertex.  Instead, compute the cell vertices
    //    and distances in this function, and add a low-level UpdateMinDistance
    //    that allows the XA, XB, and AB distances to be passed in.
    //  - It might also be more efficient to do all calculations in UVW-space,
    //    since this would involve transforming 2 points rather than 4.

    // First, check the minimum distance to the edge endpoints A and B.
    // (This also detects whether either endpoint is inside the cell.)
    S1ChordAngle min_dist = algorithm.min(getDistance(a), getDistance(b));
    if (min_dist == S1ChordAngle.zero()) {
      return min_dist;
    }

    // Otherwise, check whether the edge crosses the cell boundary.
    // Note that S2EdgeCrosser needs pointers to vertices.
    S2Point[4] v;
    for (int i = 0; i < 4; ++i) {
      v[i] = getVertex(i);
    }
    scope auto crosser = new S2EdgeCrosser(a, b, v[3]);
    for (int i = 0; i < 4; ++i) {
      if (crosser.crossingSign(v[i]) >= 0) {
        return S1ChordAngle.zero();
      }
    }
    // Finally, check whether the minimum distance occurs between a cell vertex
    // and the interior of the edge AB.  (Some of this work is redundant, since
    // it also checks the distance to the endpoints A and B again.)
    //
    // Note that we don't need to check the distance from the interior of AB to
    // the interior of a cell edge, because the only way that this distance can
    // be minimal is if the two edges cross (already checked above).
    for (int i = 0; i < 4; ++i) {
      s2edgedistances.updateMinDistance(v[i], a, b, min_dist);
    }
    return min_dist;
  }


  // Returns the maximum distance from the cell (including its interior) to the
  // given edge AB.
  S1ChordAngle getMaxDistance(in S2Point a, in S2Point b) const {
    // If the maximum distance from both endpoints to the cell is less than Pi/2
    // then the maximum distance from the edge to the cell is the maximum of the
    // two endpoint distances.
    S1ChordAngle max_dist = algorithm.max(getMaxDistance(a), getMaxDistance(b));
    if (max_dist <= S1ChordAngle.right()) {
      return max_dist;
    }

    return S1ChordAngle.straight() - getDistance(-a, -b);
  }

  // Returns the distance from the cell to the given cell.  Returns zero if
  // one cell contains the other.
  S1ChordAngle getDistance(in S2Cell target) const {
    // If the cells intersect, the distance is zero.  We use the (u,v) ranges
    // rather S2CellId::intersects() so that cells that share a partial edge or
    // corner are considered to intersect.
    if (_face == target._face && _uv.intersects(target._uv)) {
      return S1ChordAngle.zero();
    }

    // Otherwise, the minimum distance always occurs between a vertex of one
    // cell and an edge of the other cell (including the edge endpoints).  This
    // represents a total of 32 possible (vertex, edge) pairs.
    //
    // TODO(ericv): This could be optimized to be at least 5x faster by pruning
    // the set of possible closest vertex/edge pairs using the faces and (u,v)
    // ranges of both cells.
    S2Point[4] va, vb;
    for (int i = 0; i < 4; ++i) {
      va[i] = getVertex(i);
      vb[i] = target.getVertex(i);
    }
    S1ChordAngle min_dist = S1ChordAngle.infinity();
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        s2edgedistances.updateMinDistance(va[i], vb[j], vb[(j + 1) & 3], min_dist);
        s2edgedistances.updateMinDistance(vb[i], va[j], va[(j + 1) & 3], min_dist);
      }
    }
    return min_dist;
  }

  // Returns the maximum distance from the cell (including its interior) to the
  // given target cell.
  S1ChordAngle getMaxDistance(in S2Cell target) const {
    // Need to check the antipodal target for intersection with the cell. If it
    // intersects, the distance is S1ChordAngle::Straight().
    if (_face == oppositeFace(target._face) && _uv.intersects(oppositeUV(target._uv))) {
      return S1ChordAngle.straight();
    }

    // Otherwise, the maximum distance always occurs between a vertex of one
    // cell and an edge of the other cell (including the edge endpoints).  This
    // represents a total of 32 possible (vertex, edge) pairs.
    //
    // TODO(user): When the maximum distance is at most Pi/2, the maximum is
    // always attained between a pair of vertices, and this could be made much
    // faster by testing each vertex pair once rather than the current 4 times.
    S2Point[4] va, vb;
    for (int i = 0; i < 4; ++i) {
      va[i] = getVertex(i);
      vb[i] = target.getVertex(i);
    }
    S1ChordAngle max_dist = S1ChordAngle.negative();
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        s2edgedistances.updateMaxDistance(va[i], vb[j], vb[(j + 1) & 3], max_dist);
        s2edgedistances.updateMaxDistance(vb[i], va[j], va[(j + 1) & 3], max_dist);
      }
    }
    return max_dist;
  }


  ////////////////////////////////////////////////////////////////////////
  // Operator implementations:

  override
  bool opEquals(in Object y) {
    S2Cell other = cast(S2Cell) y;
    if (other) {
      return id() == other.id();
    } else {
      return false;
    }
  }

  override
  int opCmp(in Object y) {
    S2Cell other = cast(S2Cell) y;
    enforce(other, "Cannot compare S2Cell with another type.");
    return id > other.id() ? 1
        : id == other.id() ? 0 : -1;
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  override
  S2Region clone() const {
    return new S2Cell(this);
  }

  override
  S2Cap getCapBound() const {
    // Use the cell center in (u,v)-space as the cap axis.  This vector is
    // very close to GetCenter() and faster to compute.  Neither one of these
    // vectors yields the bounding cap with minimal surface area, but they
    // are both pretty close.
    //
    // It's possible to show that the two vertices that are furthest from
    // the (u,v)-origin never determine the maximum cap size (this is a
    // possible future optimization).

    S2Point center = s2coords.FaceUVtoXYZ(_face, _uv.getCenter()).normalize();
    S2Cap cap = S2Cap.fromPoint(center);
    for (int k = 0; k < 4; ++k) {
      cap.addPoint(getVertex(k));
    }
    return cap;
  }

  override
  S2LatLngRect getRectBound() const {
    if (_level > 0) {
      // Except for cells at level 0, the latitude and longitude extremes are
      // attained at the vertices.  Furthermore, the latitude range is
      // determined by one pair of diagonally opposite vertices and the
      // longitude range is determined by the other pair.
      //
      // We first determine which corner (i,j) of the cell has the largest
      // absolute latitude.  To maximize latitude, we want to find the point in
      // the cell that has the largest absolute z-coordinate and the smallest
      // absolute x- and y-coordinates.  To do this we look at each coordinate
      // (u and v), and determine whether we want to minimize or maximize that
      // coordinate based on the axis direction and the cell's (u,v) quadrant.
      double u = _uv[0][0] + _uv[0][1];
      double v = _uv[1][0] + _uv[1][1];
      int i = s2coords.GetUAxis(_face)[2] == 0 ? (u < 0) : (u > 0);
      int j = s2coords.GetVAxis(_face)[2] == 0 ? (v < 0) : (v > 0);
      R1Interval lat = R1Interval.fromPointPair(getLatitude(i, j), getLatitude(1-i, 1-j));
      S1Interval lng = S1Interval.fromPointPair(getLongitude(i, 1-j), getLongitude(1-i, j));

      // We grow the bounds slightly to make sure that the bounding rectangle
      // contains S2LatLng(P) for any point P inside the loop L defined by the
      // four *normalized* vertices.  Note that normalization of a vector can
      // change its direction by up to 0.5 * DBL_EPSILON radians, and it is not
      // enough just to add Normalize() calls to the code above because the
      // latitude/longitude ranges are not necessarily determined by diagonally
      // opposite vertex pairs after normalization.
      //
      // We would like to bound the amount by which the latitude/longitude of a
      // contained point P can exceed the bounds computed above.  In the case of
      // longitude, the normalization error can change the direction of rounding
      // leading to a maximum difference in longitude of 2 * DBL_EPSILON.  In
      // the case of latitude, the normalization error can shift the latitude by
      // up to 0.5 * DBL_EPSILON and the other sources of error can cause the
      // two latitudes to differ by up to another 1.5 * DBL_EPSILON, which also
      // leads to a maximum difference of 2 * DBL_EPSILON.
      return new S2LatLngRect(lat, lng)
          .expanded(S2LatLng.fromRadians(2 * double.epsilon, 2 * double.epsilon))
          .polarClosure();
    }

    // The 4 cells around the equator extend to +/-45 degrees latitude at the
    // midpoints of their top and bottom edges.  The two cells covering the
    // poles extend down to +/-35.26 degrees at their vertices.  The maximum
    // error in this calculation is 0.5 * DBL_EPSILON.
    const double kPoleMinLat = math.asin(math.sqrt(1.0/3)) - 0.5 * double.epsilon;

    // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
    enforce(((_face < 3) ? 1 : -1) == s2coords.GetNorm(_face)[_face % 3]);

    S2LatLngRect bound;
    switch (_face) {
      case 0:
        bound = new S2LatLngRect(
            R1Interval(-math.PI_4, math.PI_4), S1Interval(-math.PI_4, math.PI_4));
        break;
      case 1:
        bound = new S2LatLngRect(
            R1Interval(-math.PI_4, math.PI_4), S1Interval(math.PI_4, 3 * math.PI_4));
        break;
      case 2:
        bound = new S2LatLngRect(R1Interval(kPoleMinLat, math.PI_2), S1Interval.full());
        break;
      case 3:
        bound = new S2LatLngRect(
            R1Interval(-math.PI_4, math.PI_4), S1Interval(3 * math.PI_4, -3 * math.PI_4));
        break;
      case 4:
        bound = new S2LatLngRect(
            R1Interval(-math.PI_4, math.PI_4), S1Interval(-3 * math.PI_4, -math.PI_4));
        break;
      default:
        bound = new S2LatLngRect(R1Interval(-math.PI_2, -kPoleMinLat), S1Interval.full());
        break;
    }
    // Finally, we expand the bound to account for the error when a point P is
    // converted to an S2LatLng to test for containment.  (The bound should be
    // large enough so that it contains the computed S2LatLng of any contained
    // point, not just the infinite-precision version.)  We don't need to expand
    // longitude because longitude is calculated via a single call to atan2(),
    // which is guaranteed to be semi-monotonic.  (In fact the Gnu implementation
    // is also correctly rounded, but we don't even need that here.)
    return bound.expanded(S2LatLng.fromRadians(double.epsilon, 0));
  }

  override
  void getCellUnionBound(out S2CellId[] cellIds) const {
    return getCapBound().getCellUnionBound(cellIds);
  }

  override
  bool contains(in S2Cell cell) const {
    return _id.contains(cell._id);
  }

  override
  bool mayIntersect(in S2Cell cell) const {
    return _id.intersects(cell._id);
  }

  // Returns true if the cell contains the given point "p".  Note that unlike
  // S2Loop/S2Polygon, S2Cells are considered to be closed sets.  This means
  // that points along an S2Cell edge (or at a vertex) belong to the adjacent
  // cell(s) as well.
  //
  // If instead you want every point to be contained by exactly one S2Cell,
  // you will need to convert the S2Cells to S2Loops (which implement point
  // containment this way).
  //
  // The point "p" does not need to be normalized.
  override
  bool contains(in S2Point p) const {
    // We can't just call XYZtoFaceUV, because for points that lie on the
    // boundary between two faces (i.e. u or v is +1/-1) we need to return
    // true for both adjacent cells.
    R2Point uv;
    if (!s2coords.FaceXYZtoUV(_face, p, uv)) return false;

    // Expand the (u,v) bound to ensure that
    //
    //   S2Cell(S2CellId(p)).Contains(p)
    //
    // is always true.  To do this, we need to account for the error when
    // converting from (u,v) coordinates to (s,t) coordinates.  At least in the
    // case of S2_QUADRATIC_PROJECTION, the total error is at most DBL_EPSILON.
    return _uv.expanded(double.epsilon).contains(uv);
  }

  // Appends a serialized representation of the S2Cell to "encoder".
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  //void Encode(Encoder* const encoder) const;

  // Decodes an S2Cell encoded with Encode().  Returns true on success.
  //bool Decode(Decoder* const decoder);

  override
  string toString() const {
    import std.format;
    return format!("[face=%d, level=%d, orientation=%d, id=%s, uv=%s]")(
        _face, _level, _orientation, _id, _uv);
  }

 private:
  // Returns the latitude or longitude of the cell vertex given by (i,j),
  // where "i" and "j" are either 0 or 1.
  double getLatitude(int i, int j) const {
    S2Point p = s2coords.FaceUVtoXYZ(_face, _uv[0][i], _uv[1][j]);
    return S2LatLng.latitude(p).radians();
  }

  double getLongitude(int i, int j) const {
    S2Point p = s2coords.FaceUVtoXYZ(_face, _uv[0][i], _uv[1][j]);
    return S2LatLng.longitude(p).radians();
  }

  S1ChordAngle vertexChordDist(in S2Point p, int i, int j) const {
    S2Point vertex = S2Point(_uv[0][i], _uv[1][j], 1).normalize();
    return S1ChordAngle(p, vertex);
  }

  /**
   * Given a point P and either the lower or upper edge of the S2Cell (specified
   * by setting "v_end" to 0 or 1 respectively), return true if P is closer to
   * the interior of that edge than it is to either endpoint.
   */
  bool UEdgeIsClosest(in S2Point p, int v_end) const {
    double u0 = _uv[0][0], u1 = _uv[0][1], v = _uv[1][v_end];
    // These are the normals to the planes that are perpendicular to the edge
    // and pass through one of its two endpoints.
    auto dir0 = S2Point(v * v + 1, -u0 * v, -u0);
    auto dir1 = S2Point(v * v + 1, -u1 * v, -u1);
    return p.dotProd(dir0) > 0 && p.dotProd(dir1) < 0;
  }

  /**
   * Given a point P and either the left or right edge of the S2Cell (specified
   * by setting "u_end" to 0 or 1 respectively), return true if P is closer to
   * the interior of that edge than it is to either endpoint.
   */
  bool VEdgeIsClosest(in S2Point p, int u_end) const {
    double v0 = _uv[1][0], v1 = _uv[1][1], u = _uv[0][u_end];
    // See comments above.
    auto dir0 = S2Point(-u * v0, u * u + 1, -v0);
    auto dir1 = S2Point(-u * v1, u * u + 1, -v1);
    return p.dotProd(dir0) > 0 && p.dotProd(dir1) < 0;
  }

  // Returns the distance from the given point to the interior of the cell if
  // "to_interior" is true, and to the boundary of the cell otherwise.
  S1ChordAngle getDistanceInternal(in S2Point target_xyz, bool to_interior) const {
    // All calculations are done in the (u,v,w) coordinates of this cell's face.
    S2Point target = s2coords.FaceXYZtoUVW(_face, target_xyz);

    // Compute dot products with all four upward or rightward-facing edge
    // normals.  "dirIJ" is the dot product for the edge corresponding to axis
    // I, endpoint J.  For example, dir01 is the right edge of the S2Cell
    // (corresponding to the upper endpoint of the u-axis).
    double dir00 = target[0] - target[2] * _uv[0][0];
    double dir01 = target[0] - target[2] * _uv[0][1];
    double dir10 = target[1] - target[2] * _uv[1][0];
    double dir11 = target[1] - target[2] * _uv[1][1];
    bool inside = true;
    if (dir00 < 0) {
      inside = false;  // Target is to the left of the cell
      if (VEdgeIsClosest(target, 0)) return edgeDistance(-dir00, _uv[0][0]);
    }
    if (dir01 > 0) {
      inside = false;  // Target is to the right of the cell
      if (VEdgeIsClosest(target, 1)) return edgeDistance(dir01, _uv[0][1]);
    }
    if (dir10 < 0) {
      inside = false;  // Target is below the cell
      if (UEdgeIsClosest(target, 0)) return edgeDistance(-dir10, _uv[1][0]);
    }
    if (dir11 > 0) {
      inside = false;  // Target is above the cell
      if (UEdgeIsClosest(target, 1)) return edgeDistance(dir11, _uv[1][1]);
    }
    if (inside) {
      if (to_interior) return S1ChordAngle.zero();
      // Although you might think of S2Cells as rectangles, they are actually
      // arbitrary quadrilaterals after they are projected onto the sphere.
      // Therefore the simplest approach is just to find the minimum distance to
      // any of the four edges.
      return algorithm.min(
          algorithm.min(
              edgeDistance(-dir00, _uv[0][0]),
              edgeDistance(dir01, _uv[0][1])),
          algorithm.min(
              edgeDistance(-dir10, _uv[1][0]),
              edgeDistance(dir11, _uv[1][1])));
    }
    // Otherwise, the closest point is one of the four cell vertices.  Note that
    // it is *not* trivial to narrow down the candidates based on the edge sign
    // tests above, because (1) the edges don't meet at right angles and (2)
    // there are points on the far side of the sphere that are both above *and*
    // below the cell, etc.
    return algorithm.min(
        algorithm.min(
            vertexChordDist(target, 0, 0),
            vertexChordDist(target, 1, 0)),
        algorithm.min(
            vertexChordDist(target, 0, 1),
            vertexChordDist(target, 1, 1)));
  }

  // Given the dot product of a point P with the normal of a u- or v-edge at the
  // given coordinate value, return the distance from P to that edge.
  static S1ChordAngle edgeDistance(double dirIJ, double uv) {
    // Let P by the target point and let R be the closest point on the given
    // edge AB.  The desired distance PR can be expressed as PR^2 = PQ^2 + QR^2
    // where Q is the point P projected onto the plane through the great circle
    // through AB.  We can compute the distance PQ^2 perpendicular to the plane
    // from "dirIJ" (the dot product of the target point P with the edge
    // normal) and the squared length the edge normal (1 + uv**2).
    double pq2 = (dirIJ * dirIJ) / (1 + uv * uv);

    // We can compute the distance QR as (1 - OQ) where O is the sphere origin,
    // and we can compute OQ^2 = 1 - PQ^2 using the Pythagorean theorem.
    // (This calculation loses accuracy as angle POQ approaches Pi/2.)
    double qr = 1 - math.sqrt(1.0 - pq2);
    return S1ChordAngle.fromLength2(pq2 + qr * qr);
  }

  static int oppositeFace(int face) {
    return face >= 3 ? face - 3 : face + 3;
  }

  // The antipodal UV is the transpose of the original UV, interpreted within
  // the opposite face.
  static R2Rect oppositeUV(in R2Rect uv) {
    return R2Rect(uv[1], uv[0]);
  }

  // This structure occupies 44 bytes plus one pointer for the vtable.
  byte _face;
  byte _level;
  byte _orientation;
  S2CellId _id;
  R2Rect _uv;
}
