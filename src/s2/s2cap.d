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

module s2.s2cap;

import algorithm = std.algorithm;
import format = std.format;
import math = std.math;
import s2.r1interval;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2region;
import s2.util.math.vector;
import s2edge_distances = s2.s2edge_distances;
import s2metrics = s2.s2metrics;
import s2pointutil = s2.s2pointutil;

// S2Cap represents a disc-shaped region defined by a center and radius.
// Technically this shape is called a "spherical cap" (rather than disc)
// because it is not planar; the cap represents a portion of the sphere that
// has been cut off by a plane.  The boundary of the cap is the circle defined
// by the intersection of the sphere and the plane.  For containment purposes,
// the cap is a closed set, i.e. it contains its boundary.
//
// For the most part, you can use a spherical cap wherever you would use a
// disc in planar geometry.  The radius of the cap is measured along the
// surface of the sphere (rather than the straight-line distance through the
// interior).  Thus a cap of radius Pi/2 is a hemisphere, and a cap of radius
// Pi covers the entire sphere.
//
// A cap can also be defined by its center point and height.  The height
// is simply the distance from the center point to the cutoff plane.  There is
// also support for "empty" and "full" caps, which contain no points and all
// points respectively.
//
// This class is intended to be copied by value as desired.  It uses the
// default copy constructor and assignment operator, however it is not a
// "plain old datatype" (POD) because it has virtual functions.
class S2Cap : S2Region {
private:
  // Here are some useful relationships between the cap height (h), the cap
  // radius (r), the maximum chord length from the cap's center (d), and the
  // radius of cap's base (a).
  //
  //     h = 1 - cos(r)
  //       = 2 * sin^2(r/2)
  //   d^2 = 2 * h
  //       = a^2 + h^2

  // Return true if the cap intersects "cell", given that the cap does contain
  // any of the cell vertices (supplied in "vertices", an array of length 4).
  bool intersects(in S2Cell cell, in S2Point[] vertices) const {
    // Return true if this cap intersects any point of 'cell' excluding its
    // vertices (which are assumed to already have been checked).

    // If the cap is a hemisphere or larger, the cell and the complement of the
    // cap are both convex.  Therefore since no vertex of the cell is contained,
    // no other interior point of the cell is contained either.
    if (_radius >= S1ChordAngle.right()) {
      return false;
    }

    // We need to check for empty caps due to the center check just below.
    if (isEmpty()) {
      return false;
    }

    // Optimization: return true if the cell contains the cap center.  (This
    // allows half of the edge checks below to be skipped.)
    if (cell.contains(_center)) {
      return true;
    }

    // At this point we know that the cell does not contain the cap center,
    // and the cap does not contain any cell vertex.  The only way that they
    // can intersect is if the cap intersects the interior of some edge.

    double sin2_angle = _radius.sin2();
    for (int k = 0; k < 4; ++k) {
      S2Point edge = cell.getEdgeRaw(k);
      double dot = _center.dotProd(edge);
      if (dot > 0) {
        // The center is in the interior half-space defined by the edge.  We don't
        // need to consider these edges, since if the cap intersects this edge
        // then it also intersects the edge on the opposite side of the cell
        // (because we know the center is not contained with the cell).
        continue;
      }
      // The Norm2() factor is necessary because "edge" is not normalized.
      if (dot * dot > sin2_angle * edge.norm2()) {
        return false;  // Entire cap is on the exterior side of this edge.
      }
      // Otherwise, the great circle containing this edge intersects
      // the interior of the cap.  We just need to check whether the point
      // of closest approach occurs between the two edge endpoints.
      Vector3_d dir = edge.crossProd(_center);
      if (dir.dotProd(vertices[k]) < 0 && dir.dotProd(vertices[(k + 1) & 3]) > 0)
        return true;
    }
    return false;
  }

  S2Point _center;
  S1ChordAngle _radius;

public:
  // The default constructor returns an empty S2Cap.
  this() {
    _center = S2Point(1, 0, 0);
    _radius = S1ChordAngle.negative();
  }

  this(in S2Cap o) {
    _center = o._center;
    _radius = o._radius;
  }

  // Constructs a cap with the given center and radius.  A negative radius
  // yields an empty cap; a radius of 180 degrees or more yields a full cap
  // (containing the entire sphere).  "center" should be unit length.
  this(in S2Point center, S1Angle radius)
  out {
    assert(isValid());
  } do {
    _center = center;
    // The "min" calculation is necessary to handle S1Angle::Infinity().
    _radius = algorithm.min(radius, S1Angle.fromRadians(M_PI));
  }

  // Constructs a cap where the angle is expressed as an S1ChordAngle.  This
  // constructor is more efficient than the one above.
  this(in S2Point center, S1ChordAngle radius)
  out {
    assert(isValid());
  } do {
    _center = center;
    _radius = radius;
  }

  ~this() {}

  // Convenience function that creates a cap containing a single point.  This
  // method is more efficient that the S2Cap(center, radius) constructor.
  static S2Cap fromPoint(in S2Point center) {
    return new S2Cap(center, S1ChordAngle.zero());
  }

  // Returns a cap with the given center and height (see comments above).  A
  // negative height yields an empty cap; a height of 2 or more yields a full
  // cap.  "center" should be unit length.
  static S2Cap fromCenterHeight(in S2Point center, double height) {
    return new S2Cap(center, S1ChordAngle.fromLength2(2 * height));
  }

  // Return a cap with the given center and surface area.  Note that the area
  // can also be interpreted as the solid angle subtended by the cap (because
  // the sphere has unit radius).  A negative area yields an empty cap; an
  // area of 4*Pi or more yields a full cap.  "center" should be unit length.
  static S2Cap fromCenterArea(in S2Point center, double area) {
    return new S2Cap(center, S1ChordAngle.fromLength2(area / M_PI));
  }

  // Return an empty cap, i.e. a cap that contains no points.
  static S2Cap empty() {
    return new S2Cap();
  }

  // Return a full cap, i.e. a cap that contains all points.
  static S2Cap full() {
    return new S2Cap(S2Point(1, 0, 0), S1ChordAngle.straight());
  }

  // Accessor methods.
  @property
  S2Point center() const {
    return _center;
  }

  @property
  S1ChordAngle radius() const {
    return _radius;
  }

  // Returns the height of the cap, i.e. the distance from the center point to
  // the cutoff plane.
  double height() const {
    return 0.5 * _radius.length2();
  }

  // Return the cap radius as an S1Angle.  (Note that the cap angle is stored
  // internally as an S1ChordAngle, so this method requires a trigonometric
  // operation and may yield a slightly different result than the value passed
  // to the (S2Point, S1Angle) constructor.)
  S1Angle getRadius() const {
    return _radius.toS1Angle();
  }

  // Return the area of the cap.
  double getArea() const {
    return 2 * M_PI * algorithm.max(0.0, height());
  }

  // Return the true centroid of the cap multiplied by its surface area (see
  // s2centroids.h for details on centroids). The result lies on the ray from
  // the origin through the cap's center, but it is not unit length. Note that
  // if you just want the "surface centroid", i.e. the normalized result, then
  // it is much simpler just to call center().
  //
  // The reason for multiplying the result by the cap area is to make it
  // easier to compute the centroid of more complicated shapes.  The centroid
  // of a union of disjoint regions can be computed simply by adding their
  // GetCentroid() results. Caveat: for caps that contain a single point
  // (i.e., zero radius), this method always returns the origin (0, 0, 0).
  // This is because shapes with no area don't affect the centroid of a
  // union whose total area is positive.
  S2Point getCentroid() const {
    // From symmetry, the centroid of the cap must be somewhere on the line
    // from the origin to the center of the cap on the surface of the sphere.
    // When a sphere is divided into slices of constant thickness by a set of
    // parallel planes, all slices have the same surface area. This implies
    // that the radial component of the centroid is simply the midpoint of the
    // range of radial distances spanned by the cap. That is easily computed
    // from the cap height.
    if (isEmpty()) {
      return S2Point();
    }
    double r = 1.0 - 0.5 * height();
    return r * getArea() * _center;
  }


  // We allow negative heights (to represent empty caps) but heights are
  // normalized so that they do not exceed 2.
  bool isValid() const {
    return s2pointutil.isUnitLength(_center) && _radius.length2() <= 4;
  }

  // Return true if the cap is empty, i.e. it contains no points.
  bool isEmpty() const {
    return _radius.isNegative();
  }

  // Return true if the cap is full, i.e. it contains all points.
  bool isFull() const {
    return _radius.length2() == 4;
  }

  // Return the complement of the interior of the cap.  A cap and its
  // complement have the same boundary but do not share any interior points.
  // The complement operator is not a bijection because the complement of a
  // singleton cap (containing a single point) is the same as the complement
  // of an empty cap.
  S2Cap complement() const {
    // The complement of a full cap is an empty cap, not a singleton.
    // Also make sure that the complement of an empty cap is full.
    if (isFull()) {
      return empty();
    }
    if (isEmpty()) {
      return full();
    }
    return new S2Cap(-_center, S1ChordAngle.fromLength2(4 - _radius.length2()));
  }

  // Return true if and only if this cap contains the given other cap
  // (in a set containment sense, e.g. every cap contains the empty cap).
  bool contains(in S2Cap other) const {
    if (isFull() || other.isEmpty()) {
      return true;
    }
    return _radius >= S1ChordAngle(_center, other._center) + other._radius;
  }

  // Return true if and only if this cap intersects the given other cap,
  // i.e. whether they have any points in common.
  bool intersects(in S2Cap other) const {
    if (isEmpty() || other.isEmpty()) {
      return false;
    }
    return _radius + other._radius >= S1ChordAngle(_center, other._center);
  }

  // Return true if and only if the interior of this cap intersects the
  // given other cap.  (This relationship is not symmetric, since only
  // the interior of this cap is used.)
  bool interiorIntersects(in S2Cap other) const {
    // Make sure this cap has an interior and the other cap is non-empty.
    if (_radius.length2() <= 0 || other.isEmpty()) {
      return false;
    }
    return _radius + other._radius > S1ChordAngle(_center, other._center);
  }

  // Return true if and only if the given point is contained in the interior
  // of the cap (i.e. the cap excluding its boundary).  "p" should be be a
  // unit-length vector.
  bool interiorContains(in S2Point p) const
  in {
    assert(s2pointutil.isUnitLength(p));
  } do {
    return isFull() || S1ChordAngle(_center, p) < _radius;
  }

  // Increase the cap height if necessary to include the given point.  If the
  // cap is empty then the center is set to the given point, but otherwise the
  // center is not changed.  "p" should be a unit-length vector.
  void addPoint(in S2Point p)
  in {
    // Compute the squared chord length, then convert it into a height.
    assert(s2pointutil.isUnitLength(p));
  } do {
    if (isEmpty()) {
      _center = p;
      _radius = S1ChordAngle.zero();
    } else {
      // After calling cap.AddPoint(p), cap.Contains(p) must be true.  However
      // we don't need to do anything special to achieve this because Contains()
      // does exactly the same distance calculation that we do here.
      _radius = algorithm.max(_radius, S1ChordAngle(_center, p));
    }
  }

  // Increase the cap height if necessary to include "other".  If the current
  // cap is empty it is set to the given other cap.
  void addCap(in S2Cap other) {
    if (isEmpty()) {
      _center = other._center;
      _radius = other._radius;
    } else {
      // We round up the distance to ensure that the cap is actually contained.
      // TODO(ericv): Do some error analysis in order to guarantee this.
      S1ChordAngle dist = S1ChordAngle(_center, other._center) + other._radius;
      _radius = algorithm.max(_radius, dist.plusError(double.epsilon * dist.length2()));
    }
  }

  // Return a cap that contains all points within a given distance of this
  // cap.  Note that any expansion of the empty cap is still empty.
  S2Cap expanded(S1Angle distance) const
  in {
    assert(distance.radians() >= 0);
  } do {
    if (isEmpty()) {
      return empty();
    }
    return new S2Cap(_center, _radius + S1ChordAngle(distance));
  }

  // Return the smallest cap which encloses this cap and "other".
  S2Cap unite(in S2Cap other) const {
    if (_radius < other._radius) {
      return other.unite(this);
    }
    if (isFull() || other.isEmpty()) {
      return new S2Cap(this);
    }
    // This calculation would be more efficient using S1ChordAngles.
    S1Angle this_radius = getRadius();
    S1Angle other_radius = other.getRadius();
    S1Angle distance = S1Angle(center(), other.center());
    if (this_radius >= distance + other_radius) {
      return new S2Cap(this);
    } else {
      S1Angle result_radius = 0.5 * (distance + this_radius + other_radius);
      S2Point result_center = s2edge_distances.interpolateAtDistance(
          0.5 * (distance - this_radius + other_radius),
          center(),
          other.center());
      return new S2Cap(result_center, result_radius);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.d for details):

  override
  S2Cap clone() const {
    return new S2Cap(this);
  }

  override
  S2Cap getCapBound() const {
    return new S2Cap(this);
  }

  override
  S2LatLngRect getRectBound() const {
    if (isEmpty()) {
      return S2LatLngRect.empty();
    }

    // Convert the center to a (lat,lng) pair, and compute the cap angle.
    S2LatLng center_ll = S2LatLng(_center);
    double cap_angle = getRadius().radians();

    bool all_longitudes = false;
    double[2] lat;
    double[2] lng;
    lng[0] = -M_PI;
    lng[1] = M_PI;

    // Check whether cap includes the south pole.
    lat[0] = center_ll.lat().radians() - cap_angle;
    if (lat[0] <= -math.PI_2) {
      lat[0] = -math.PI_2;
      all_longitudes = true;
    }
    // Check whether cap includes the north pole.
    lat[1] = center_ll.lat().radians() + cap_angle;
    if (lat[1] >= math.PI_2) {
      lat[1] = math.PI_2;
      all_longitudes = true;
    }
    if (!all_longitudes) {
      // Compute the range of longitudes covered by the cap.  We use the law
      // of sines for spherical triangles.  Consider the triangle ABC where
      // A is the north pole, B is the center of the cap, and C is the point
      // of tangency between the cap boundary and a line of longitude.  Then
      // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
      // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
      // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
      // minus the latitude).  This formula also works for negative latitudes.
      //
      // The formula for sin(a) follows from the relationship h = 1 - cos(a).

      double sin_a = _radius.sin();
      double sin_c = math.cos(center_ll.lat().radians());
      if (sin_a <= sin_c) {
        double angle_A = math.asin(sin_a / sin_c);
        lng[0] = math.remainder(center_ll.lng().radians() - angle_A, 2 * M_PI);
        lng[1] = math.remainder(center_ll.lng().radians() + angle_A, 2 * M_PI);
      }
    }
    return new S2LatLngRect(R1Interval(lat[0], lat[1]), S1Interval(lng[0], lng[1]));
  }

  // Computes a covering of the S2Cap.  In general the covering consists of at
  // most 4 cells except for very large caps, which may need up to 6 cells.
  // The output is not sorted.
  override
  void getCellUnionBound(out S2CellId[] cell_ids) {
    // TODO(ericv): The covering could be made quite a bit tighter by mapping
    // the cap to a rectangle in (i,j)-space and finding a covering for that.

    // Find the maximum level such that the cap contains at most one cell vertex
    // and such that S2CellId::AppendVertexNeighbors() can be called.
    int level = s2metrics.MIN_WIDTH.getLevelForMinValue(getRadius().radians()) - 1;

    // If level < 0, then more than three face cells are required.
    if (level < 0) {
      cell_ids.reserve(6);
      for (int face = 0; face < 6; ++face) {
        cell_ids ~= S2CellId.fromFace(face);
      }
    } else {
      // The covering consists of the 4 cells at the given level that share the
      // cell vertex that is closest to the cap center.
      cell_ids.reserve(4);
      S2CellId(_center).appendVertexNeighbors(level, cell_ids);
    }
  }

  override
  bool contains(in S2Cell cell) const  {
    // If the cap does not contain all cell vertices, return false.
    // We check the vertices before taking the Complement() because we can't
    // accurately represent the complement of a very small cap (a height
    // of 2-epsilon is rounded off to 2).
    S2Point[4] vertices;
    for (int k = 0; k < 4; ++k) {
      vertices[k] = cell.getVertex(k);
      if (!contains(vertices[k])) return false;
    }
    // Otherwise, return true if the complement of the cap does not intersect
    // the cell.  (This test is slightly conservative, because technically we
    // want Complement().InteriorIntersects() here.)
    return !complement().intersects(cell, vertices);
  }

  override
  bool mayIntersect(in S2Cell cell) const {
    // If the cap contains any cell vertex, return true.
    S2Point[4] vertices;
    for (int k = 0; k < 4; ++k) {
      vertices[k] = cell.getVertex(k);
      if (contains(vertices[k])) return true;
    }
    return intersects(cell, vertices);
  }


  // The point "p" should be a unit-length vector.
  override
  bool contains(in S2Point p) const
  in {
    assert(s2pointutil.isUnitLength(p));
  } do {
    return S1ChordAngle(_center, p) <= _radius;
  }

  // Appends a serialized representation of the S2Cap to "encoder".
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  // void Encode(Encoder* const encoder) const;

  // Decodes an S2Cap encoded with Encode().  Returns true on success.
  // bool Decode(Decoder* const decoder);

  ///////////////////////////////////////////////////////////////////////
  // The following static methods are convenience functions for assertions
  // and testing purposes only.

  // Return true if two caps are identical.
  override
  bool opEquals(in Object o) const {
    S2Cap other = cast(S2Cap) o;
    if (other !is null) {
      return (_center == other._center && _radius == other._radius)
          || (isEmpty() && other.isEmpty())
          || (isFull() && other.isFull());
    }
    return false;
  }

  // Return true if the cap center and height differ by at most "max_error"
  // from the given cap "other".
  bool approxEquals(in S2Cap other, S1Angle max_error_angle = S1Angle.fromRadians(1e-14)) const {
    const double max_error = max_error_angle.radians();
    const double r2 = _radius.length2();
    const double other_r2 = other._radius.length2();
    return (s2pointutil.approxEquals(_center, other._center, max_error_angle) &&
        math.fabs(r2 - other_r2) <= max_error) ||
        (isEmpty() && other_r2 <= max_error) ||
        (other.isEmpty() && r2 <= max_error) ||
        (isFull() && other_r2 >= 2 - max_error) ||
        (other.isFull() && r2 >= 2 - max_error);
  }

  override
  string toString() const {
    return format.format(
        "[_center=%s, _radius=%s]", center().toString(), getRadius().toString());
  }
}
