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

module s2.s2latlng_rect;

import algorithm = std.algorithm;
import edgedistances = s2.s2edge_distances;
import math = std.math;
import s2.logger;
import s2.r1interval;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2edge_crossings;
import s2.s2latlng;
import s2.s2point;
import s2.s2pointutil;
import s2.s2predicates;
import s2.s2region;
import s2.util.math.vector;
import std.conv;
import std.exception;

/**
 * An S2LatLngRect represents a closed latitude-longitude rectangle.  It is
 * capable of representing the empty and full rectangles as well as single
 * points.  Note that the latitude-longitude space is considered to have a
 * *cylindrical* topology rather than a spherical one, i.e. the poles have
 * multiple lat/lng representations.  An S2LatLngRect may be defined so that
 * includes some representations of a pole but not others.  Use the
 * PolarClosure() method if you want to expand a rectangle so that it contains
 * all possible representations of any contained poles.
 *
 * Because S2LatLngRect uses S1Interval to store the longitude range,
 * longitudes of -180 degrees are treated specially.  Except for empty
 * and full longitude spans, -180 degree longitudes will turn into +180
 * degrees.  This sign flip causes lng_lo() to be greater than lng_hi(),
 * indicating that the rectangle will wrap around through -180 instead of
 * through +179. Thus the math is consistent within the library, but the sign
 * flip can be surprising, especially when working with map projections where
 * -180 and +180 are at opposite ends of the flattened map.  See the comments
 * on S1Interval for more details.
 *
 * This class is intended to be copied by value as desired.  It uses
 * the default copy constructor and assignment operator, however it is
 * not a "plain old datatype" (POD) because it has virtual functions.
 */
class S2LatLngRect : S2Region {
public:
  // Construct a rectangle from minimum and maximum latitudes and longitudes.
  // If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude
  // line. Both points must be normalized, with lo.lat() <= hi.lat().
  // The rectangle contains all the points p such that 'lo' <= p <= 'hi',
  // where '<=' is defined in the obvious way.
  this(in S2LatLng lo, in S2LatLng hi) {
    _lat = R1Interval(lo.lat().radians(), hi.lat().radians());
    _lng = S1Interval(lo.lng().radians(), hi.lng().radians());
    if (!isValid) {
      logger.logError("Invalid rect: ", lo, ", ", hi);
    }
  }

  // Construct a rectangle from latitude and longitude intervals.  The two
  // intervals must either be both empty or both non-empty, and the latitude
  // interval must not extend outside [-90, +90] degrees.
  // Note that both intervals (and hence the rectangle) are closed.
  this(in R1Interval lat, in S1Interval lng) {
    if (!isValid()) logger.logError("Invalid rect: ", lat, ", ", lng);
    _lat = lat;
    _lng = lng;
  }

  this(in S2LatLngRect rect) {
    _lat = rect._lat;
    _lng = rect._lng;
  }

  // The default constructor creates an empty S2LatLngRect.
  this() {
    _lat = R1Interval.empty();
    _lng = S1Interval.empty();
  }

  // Construct a rectangle of the given size centered around the given point.
  // "center" needs to be normalized, but "size" does not.  The latitude
  // interval of the result is clamped to [-90,90] degrees, and the longitude
  // interval of the result is Full() if and only if the longitude size is
  // 360 degrees or more.  Examples of clamping (in degrees):
  //
  //   center=(80,170),  size=(40,60)   -> lat=[60,90],   lng=[140,-160]
  //   center=(10,40),   size=(210,400) -> lat=[-90,90],  lng=[-180,180]
  //   center=(-90,180), size=(20,50)   -> lat=[-90,-80], lng=[155,-155]
  static S2LatLngRect fromCenterSize(in S2LatLng center, in S2LatLng size) {
    return fromPoint(center).expanded(0.5 * size);
  }

  // Construct a rectangle containing a single (normalized) point.
  static S2LatLngRect fromPoint(in S2LatLng p) {
    if (!p.isValid()) {
      logger.logError("Invalid S2LatLng in S2LatLngRect.getDistance: ", p);
    }
    return new S2LatLngRect(p, p);
  }


  // Construct the minimal bounding rectangle containing the two given
  // normalized points.  This is equivalent to starting with an empty
  // rectangle and calling AddPoint() twice.  Note that it is different than
  // the S2LatLngRect(lo, hi) constructor, where the first point is always
  // used as the lower-left corner of the resulting rectangle.
  static S2LatLngRect fromPointPair(in S2LatLng p1, in S2LatLng p2) {
    if (!p1.isValid()) {
      logger.logError("Invalid S2LatLng in S2LatLngRect.fromPointPair: ", p1);
    }
    if (!p2.isValid()) {
      logger.logError("Invalid S2LatLng in S2LatLngRect.fromPointPair: ", p2);
    }

    return new S2LatLngRect(
        R1Interval.fromPointPair(p1.lat().radians(), p2.lat().radians()),
        S1Interval.fromPointPair(p1.lng().radians(), p2.lng().radians()));
  }


  // Accessor methods.
  @property
  S1Angle latLo() const {
    return S1Angle.fromRadians(_lat.lo());
  }
  @property
  S1Angle latHi() const {
    return S1Angle.fromRadians(_lat.hi());
  }
  @property
  S1Angle lngLo() const {
    return S1Angle.fromRadians(_lng.lo());
  }
  @property
  S1Angle lngHi() const {
    return S1Angle.fromRadians(_lng.hi());
  }
  @property
  R1Interval lat() const {
    return _lat;
  }
  @property
  S1Interval lng() const {
    return _lng;
  }

  ref R1Interval mutableLat() { return _lat; }
  ref S1Interval mutableLng() { return _lng; }

  S2LatLng lo() const { return S2LatLng(latLo(), lngLo()); }
  S2LatLng hi() const { return S2LatLng(latHi(), lngHi()); }

  // The canonical empty and full rectangles, as derived from the Empty
  // and Full R1 and S1 Intervals.
  // Empty: lat_lo=1, lat_hi=0, lng_lo=Pi, lng_hi=-Pi (radians)
  static S2LatLngRect empty() {
    return new S2LatLngRect();
  }

  // Full: lat_lo=-Pi/2, lat_hi=Pi/2, lng_lo=-Pi, lng_hi=Pi (radians)
  static S2LatLngRect full() {
    return new S2LatLngRect(fullLat(), fullLng());
  }

  // The full allowable range of latitudes and longitudes.
  static R1Interval fullLat() { return R1Interval(-S1Angle.PI_2, S1Angle.PI_2); }
  static S1Interval fullLng() { return S1Interval.full(); }

  // Returns true if the rectangle is valid, which essentially just means
  // that the latitude bounds do not exceed Pi/2 in absolute value and
  // the longitude bounds do not exceed Pi in absolute value.  Also, if
  // either the latitude or longitude bound is empty then both must be.
  bool isValid() const {
    // The lat/lng ranges must either be both empty or both non-empty.
    return (math.fabs(_lat.lo()) <= S1Angle.PI_2
        && math.fabs(_lat.hi()) <= S1Angle.PI_2
        && _lng.isValid()
        && _lat.isEmpty() == _lng.isEmpty());
  }

  // Returns true if the rectangle is empty, i.e. it contains no points at all.
  bool isEmpty() const {
    return _lat.isEmpty();
  }

  // Returns true if the rectangle is full, i.e. it contains all points.
  bool isFull() const {
    return _lat == fullLat() && _lng.isFull();
  }

  // Returns true if the rectangle is a point, i.e. lo() == hi()
  bool isPoint() const {
    return _lat.lo() == _lat.hi() && _lng.lo() == _lng.hi();
  }

  // Returns true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses
  // the 180 degree longitude line.
  bool isInverted() const { return _lng.isInverted(); }

  // Returns the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order
  // (lower left, lower right, upper right, upper left).  For convenience, the
  // argument is reduced modulo 4 to the range [0..3].
  S2LatLng getVertex(int k) const {
    // Twiddle bits to return the points in CCW order (lower left, lower right,
    // upper right, upper left).
    int i = (k >> 1) & 1;
    return S2LatLng.fromRadians(_lat[i], _lng[i ^ (k & 1)]);
  }

  // Returns the center of the rectangle in latitude-longitude space
  // (in general this is not the center of the region on the sphere).
  S2LatLng getCenter() const {
    return S2LatLng.fromRadians(_lat.getCenter(), _lng.getCenter());
  }

  // Returns the width and height of this rectangle in latitude-longitude
  // space.  Empty rectangles have a negative width and height.
  S2LatLng getSize() const {
    return S2LatLng.fromRadians(_lat.getLength(), _lng.getLength());
  }

  // Returns the surface area of this rectangle on the unit sphere.
  double area() const {
    if (isEmpty()) return 0.0;
    // This is the size difference of the two spherical caps, multiplied by
    // the longitude ratio.
    return lng().getLength() * (latHi().sin() - latLo().sin());
  }

  /**
   * Returns the true centroid of the rectangle multiplied by its surface area
   * (see s2centroids.h for details on centroids).  The result is not unit
   * length, so you may want to normalize it.  Note that in general the
   * centroid is *not* at the center of the rectangle, and in fact it may not
   * even be contained by the rectangle.  (It is the "center of mass" of the
   * rectangle viewed as subset of the unit sphere, i.e. it is the point in
   * space about which this curved shape would rotate.)
   *
   * The reason for multiplying the result by the rectangle area is to make it
   * easier to compute the centroid of more complicated shapes.  The centroid
   * of a union of disjoint regions can be computed simply by adding their
   * GetCentroid() results.
   */
  S2Point getCentroid() const {
    // When a sphere is divided into slices of constant thickness by a set of
    // parallel planes, all slices have the same surface area.  This implies
    // that the z-component of the centroid is simply the midpoint of the
    // z-interval spanned by the S2LatLngRect.
    //
    // Similarly, it is easy to see that the (x,y) of the centroid lies in the
    // plane through the midpoint of the rectangle's longitude interval.  We
    // only need to determine the distance "d" of this point from the z-axis.
    //
    // Let's restrict our attention to a particular z-value.  In this z-plane,
    // the S2LatLngRect is a circular arc.  The centroid of this arc lies on a
    // radial line through the midpoint of the arc, and at a distance from the
    // z-axis of
    //
    //     r * (sin(alpha) / alpha)
    //
    // where r = sqrt(1-z^2) is the radius of the arc, and "alpha" is half of
    // the arc length (i.e., the arc covers longitudes [-alpha, alpha]).
    //
    // To find the centroid distance from the z-axis for the entire rectangle,
    // we just need to integrate over the z-interval.  This gives
    //
    //    d = Integrate[sqrt(1-z^2)*sin(alpha)/alpha, z1..z2] / (z2 - z1)
    //
    // where [z1, z2] is the range of z-values covered by the rectangle.  This
    // simplifies to
    //
    //    d = sin(alpha)/(2*alpha*(z2-z1))*(z2*r2 - z1*r1 + theta2 - theta1)
    //
    // where [theta1, theta2] is the latitude interval, z1=sin(theta1),
    // z2=sin(theta2), r1=cos(theta1), and r2=cos(theta2).
    //
    // Finally, we want to return not the centroid itself, but the centroid
    // scaled by the area of the rectangle.  The area of the rectangle is
    //
    //    A = 2 * alpha * (z2 - z1)
    //
    // which fortunately appears in the denominator of "d".

    if (isEmpty()) return S2Point();
    double z1 = latLo().sin(), z2 = latHi().sin();
    double r1 = latLo().cos(), r2 = latHi().cos();
    double alpha = 0.5 * _lng.getLength();
    double r = math.sin(alpha) * (r2 * z2 - r1 * z1 + _lat.getLength());
    double lng = _lng.getCenter();
    double z = alpha * (z2 + z1) * (z2 - z1);  // scaled by the area
    return S2Point(r * math.cos(lng), r * math.sin(lng), z);
  }

  // More efficient version of Contains() that accepts a S2LatLng rather than
  // an S2Point.  The argument must be normalized.
  bool contains(in S2LatLng ll) const {
    if (!ll.isValid()) logger.logError("Invalid S2LatLng in S2LatLngRect.contains: ", ll);

    return (_lat.contains(ll.lat().radians()) && _lng.contains(ll.lng().radians()));
  }

  // Returns true if and only if the given point is contained in the interior
  // of the region (i.e. the region excluding its boundary).  The point 'p'
  // does not need to be normalized.
  bool interiorContains(in S2Point p) const {
    return interiorContains(S2LatLng(p));
  }

  // More efficient version of InteriorContains() that accepts a S2LatLng
  // rather than an S2Point.  The argument must be normalized.
  bool interiorContains(in S2LatLng ll) const {
    if (!ll.isValid()) logger.logError("Invalid S2LatLng in S2LatLngRect.interiorContains: ", ll);

    return _lat.interiorContains(ll.lat().radians())
        && _lng.interiorContains(ll.lng().radians());
  }

  // Returns true if and only if the rectangle contains the given other
  // rectangle.
  bool contains(in S2LatLngRect other) const {
    return _lat.contains(other._lat) && _lng.contains(other._lng);
  }

  // Returns true if and only if the interior of this rectangle contains all
  // points of the given other rectangle (including its boundary).
  bool interiorContains(in S2LatLngRect other) const {
    return _lat.interiorContains(other._lat) && _lng.interiorContains(other._lng);
  }

  // Returns true if this rectangle and the given other rectangle have any
  // points in common.
  bool intersects(in S2LatLngRect other) const {
    return _lat.intersects(other._lat) && _lng.intersects(other._lng);
  }

  // Returns true if this rectangle intersects the given cell.  (This is an
  // exact test and may be fairly expensive, see also MayIntersect below.)
  bool intersects(in S2Cell cell) const {
    // First we eliminate the cases where one region completely contains the
    // other.  Once these are disposed of, then the regions will intersect
    // if and only if their boundaries intersect.

    if (isEmpty()) return false;
    if (contains(cell.getCenterRaw())) return true;
    if (cell.contains(getCenter().toS2Point())) return true;

    // Quick rejection test (not required for correctness).
    if (!intersects(cell.getRectBound())) return false;

    // Precompute the cell vertices as points and latitude-longitudes.  We also
    // check whether the S2Cell contains any corner of the rectangle, or
    // vice-versa, since the edge-crossing tests only check the edge interiors.

    S2Point[4] cell_v;
    S2LatLng[4] cell_ll;
    for (int i = 0; i < 4; ++i) {
      cell_v[i] = cell.getVertex(i);  // Must be normalized.
      cell_ll[i] = S2LatLng(cell_v[i]);
      if (contains(cell_ll[i])) return true;
      if (cell.contains(getVertex(i).toS2Point())) return true;
    }

    // Now check whether the boundaries intersect.  Unfortunately, a
    // latitude-longitude rectangle does not have straight edges -- two edges
    // are curved, and at least one of them is concave.

    for (int i = 0; i < 4; ++i) {
      S1Interval edge_lng = S1Interval.fromPointPair(
          cell_ll[i].lng().radians(), cell_ll[(i+1) & 3].lng().radians());
      if (!_lng.intersects(edge_lng)) continue;

      const S2Point a = cell_v[i];
      const S2Point b = cell_v[(i+1) & 3];
      if (edge_lng.contains(_lng.lo())) {
        if (intersectsLngEdge(a, b, _lat, _lng.lo())) return true;
      }
      if (edge_lng.contains(_lng.hi())) {
        if (intersectsLngEdge(a, b, _lat, _lng.hi())) return true;
      }
      if (intersectsLatEdge(a, b, _lat.lo(), _lng)) return true;
      if (intersectsLatEdge(a, b, _lat.hi(), _lng)) return true;
    }
    return false;
  }

  // Returns true if and only if the interior of this rectangle intersects
  // any point (including the boundary) of the given other rectangle.
  bool interiorIntersects(in S2LatLngRect other) const {
    return _lat.interiorIntersects(other._lat)
        && _lng.interiorIntersects(other._lng);
  }

  // Returns true if the boundary of this rectangle intersects the given
  // geodesic edge (v0, v1).
  bool boundaryIntersects(in S2Point v0, in S2Point v1) const {
    if (isEmpty()) return false;
    if (!_lng.isFull()) {
      if (intersectsLngEdge(v0, v1, _lat, _lng.lo())) return true;
      if (intersectsLngEdge(v0, v1, _lat, _lng.hi())) return true;
    }
    if (_lat.lo() != -S1Angle.PI_2 && intersectsLatEdge(v0, v1, _lat.lo(), _lng)) {
      return true;
    }
    if (_lat.hi() != S1Angle.PI_2 && intersectsLatEdge(v0, v1, _lat.hi(), _lng)) {
      return true;
    }
    return false;
  }

  /**
   * Increase the size of the bounding rectangle to include the given point.
   * The rectangle is expanded by the minimum amount possible.  The S2LatLng
   * argument must be normalized.
   */
  void addPoint(in S2Point p) {
    addPoint(S2LatLng(p));
  }

  void addPoint(in S2LatLng ll) {
    if (!ll.isValid()) logger.logError("Invalid S2LatLng in S2LatLngRect::AddPoint: ", ll);

    _lat.addPoint(ll.lat().radians());
    _lng.addPoint(ll.lng().radians());
  }

  /**
   * Returns a rectangle that has been expanded by margin.lat() on each side in
   * the latitude direction, and by margin.lng() on each side in the longitude
   * direction.  If either margin is negative, then shrinks the rectangle on
   * the corresponding sides instead.  The resulting rectangle may be empty.
   *
   * As noted above, the latitude-longitude space has the topology of a
   * cylinder.  Longitudes "wrap around" at +/-180 degrees, while latitudes
   * are clamped to range [-90, 90].  This means that any expansion (positive
   * or negative) of the full longitude range remains full (since the
   * "rectangle" is actually a continuous band around the cylinder), while
   * expansion of the full latitude range remains full only if the margin is
   * positive.
   *
   * If either the latitude or longitude interval becomes empty after
   * expansion by a negative margin, the result is empty.
   *
   * Note that if an expanded rectangle contains a pole, it may not contain
   * all possible lat/lng representations of that pole (see header above).
   * Use the PolarClosure() method if you do not want this behavior.
   *
   * If you are trying to grow a rectangle by a certain *distance* on the
   * sphere (e.g. 5km), use the ExpandedByDistance() method instead.
   */
  S2LatLngRect expanded(in S2LatLng margin) const {
    R1Interval lat = _lat.expanded(margin.lat().radians());
    S1Interval lng = _lng.expanded(margin.lng().radians());
    if (lat.isEmpty() || lng.isEmpty()) return empty();
    return new S2LatLngRect(lat.intersection(fullLat()), lng);
  }

  /**
   * If the rectangle does not include either pole, returns it unmodified.
   * Otherwise expands the longitude range to Full() so that the rectangle
   * contains all possible representations of the contained pole(s).
   */
  S2LatLngRect polarClosure() const {
    if (_lat.lo() == -S1Angle.PI_2 || _lat.hi() == S1Angle.PI_2) {
      return new S2LatLngRect(_lat, S1Interval.full());
    }
    return new S2LatLngRect(this);
  }

  /**
   * Returns the smallest rectangle containing the union of this rectangle and
   * the given rectangle.
   */
  S2LatLngRect unite(in S2LatLngRect other) const {
    return new S2LatLngRect(_lat.unite(other._lat), _lng.unite(other._lng));
  }

  /**
   * Returns the smallest rectangle containing the intersection of this
   * rectangle and the given rectangle.  Note that the region of intersection
   * may consist of two disjoint rectangles, in which case a single rectangle
   * spanning both of them is returned.
   */
  S2LatLngRect intersection(in S2LatLngRect other) const {
    R1Interval lat = _lat.intersection(other._lat);
    S1Interval lng = _lng.intersection(other._lng);
    if (lat.isEmpty() || lng.isEmpty()) {
      // The lat/lng ranges must either be both empty or both non-empty.
      return empty();
    }
    return new S2LatLngRect(lat, lng);
  }

  /**
   * Expands this rectangle so that it contains all points within the given
   * distance of the boundary, and return the smallest such rectangle.  If the
   * distance is negative, then instead shrinks this rectangle so that it
   * excludes all points within the given absolute distance of the boundary,
   * and returns the largest such rectangle.
   *
   * Unlike Expanded(), this method treats the rectangle as a set of points on
   * the sphere, and measures distances on the sphere.  For example, you can
   * use this method to find a rectangle that contains all points within 5km
   * of a given rectangle.  Because this method uses the topology of the
   * sphere, note the following:
   *
   *  - The full and empty rectangles have no boundary on the sphere.  Any
   *    expansion (positive or negative) of these rectangles leaves them
   *    unchanged.
   *
   *  - Any rectangle that covers the full longitude range does not have an
   *    east or west boundary, therefore no expansion (positive or negative)
   *    will occur in that direction.
   *
   *  - Any rectangle that covers the full longitude range and also includes
   *    a pole will not be expanded or contracted at that pole, because it
   *    does not have a boundary there.
   *
   *  - If a rectangle is within the given distance of a pole, the result will
   *    include the full longitude range (because all longitudes are present
   *    at the poles).
   *
   * Expansion and contraction are defined such that they are inverses whenver
   * possible, i.e.
   *
   *   rect.ExpandedByDistance(x).ExpandedByDistance(-x) == rect
   *
   * (approximately), so long as the first operation does not cause a
   * rectangle boundary to disappear (i.e., the longitude range newly becomes
   * full or empty, or the latitude range expands to include a pole).
   */
  S2LatLngRect expandedByDistance(S1Angle distance) const {
    if (distance >= S1Angle.zero()) {
      // The most straightforward approach is to build a cap centered on each
      // vertex and take the union of all the bounding rectangles (including the
      // original rectangle; this is necessary for very large rectangles).

      // TODO(ericv): Update this code to use an algorithm like the one below.
      auto radius = S1ChordAngle(distance);
      S2LatLngRect r = new S2LatLngRect(this);
      for (int k = 0; k < 4; ++k) {
        scope S2Cap cap = new S2Cap(getVertex(k).toS2Point(), radius);
        r = r.unite(cap.getRectBound());
      }
      return r;
    } else {
      // Shrink the latitude interval unless the latitude interval contains a pole
      // and the longitude interval is full, in which case the rectangle has no
      // boundary at that pole.
      auto lat_result = R1Interval(
          lat().lo() <= fullLat().lo() && lng().isFull()
              ? fullLat().lo() : lat().lo() - distance.radians(),
          lat().hi() >= fullLat().hi() && lng().isFull()
              ? fullLat().hi() : lat().hi() + distance.radians());
      if (lat_result.isEmpty()) {
        return S2LatLngRect.empty();
      }

      // Maximum absolute value of a latitude in lat_result. At this latitude,
      // the cap occupies the largest longitude interval.
      double max_abs_lat = algorithm.max(-lat_result.lo(), lat_result.hi());

      // Compute the largest longitude interval that the cap occupies. We use the
      // law of sines for spherical triangles. For the details, see the comment in
      // S2Cap::GetRectBound().
      //
      // When sin_a >= sin_c, the cap covers all the latitude.
      double sin_a = math.sin(-distance.radians());
      double sin_c = math.cos(max_abs_lat);
      double max_lng_margin = sin_a < sin_c ? math.asin(sin_a / sin_c) : S1Angle.PI_2;

      S1Interval lng_result = lng().expanded(-max_lng_margin);
      if (lng_result.isEmpty()) {
        return S2LatLngRect.empty();
      }
      return new S2LatLngRect(lat_result, lng_result);
    }
  }

  // Returns the minimum distance (measured along the surface of the sphere) to
  // the given S2LatLngRect. Both S2LatLngRects must be non-empty.
  S1Angle getDistance(in S2LatLngRect other) const
  in {
    assert(!this.isEmpty());
    assert(!other.isEmpty());
  } body {
    const S2LatLngRect a = this;
    const S2LatLngRect b = other;

    // First, handle the trivial cases where the longitude intervals overlap.
    if (a.lng().intersects(b.lng())) {
      if (a.lat().intersects(b.lat()))
        return S1Angle.fromRadians(0);  // Intersection between a and b.

      // We found an overlap in the longitude interval, but not in the latitude
      // interval. This means the shortest path travels along some line of
      // longitude connecting the high-latitude of the lower rect with the
      // low-latitude of the higher rect.
      S1Angle lo, hi;
      if (a.lat().lo() > b.lat().hi()) {
        lo = b.latHi();
        hi = a.latLo();
      } else {
        lo = a.latHi();
        hi = b.latLo();
      }
      return hi - lo;
    }

    // The longitude intervals don't overlap. In this case, the closest points
    // occur somewhere on the pair of longitudinal edges which are nearest in
    // longitude-space.
    S1Angle a_lng, b_lng;
    S1Interval lo_hi = S1Interval.fromPointPair(a.lng().lo(), b.lng().hi());
    S1Interval hi_lo = S1Interval.fromPointPair(a.lng().hi(), b.lng().lo());
    if (lo_hi.getLength() < hi_lo.getLength()) {
      a_lng = a.lngLo();
      b_lng = b.lngHi();
    } else {
      a_lng = a.lngHi();
      b_lng = b.lngLo();
    }

    // The shortest distance between the two longitudinal segments will include at
    // least one segment endpoint. We could probably narrow this down further to a
    // single point-edge distance by comparing the relative latitudes of the
    // endpoints, but for the sake of clarity, we'll do all four point-edge
    // distance tests.
    S2Point a_lo = S2LatLng(a.latLo(), a_lng).toS2Point();
    S2Point a_hi = S2LatLng(a.latHi(), a_lng).toS2Point();
    S2Point b_lo = S2LatLng(b.latLo(), b_lng).toS2Point();
    S2Point b_hi = S2LatLng(b.latHi(), b_lng).toS2Point();
    return algorithm.min(
        edgedistances.getDistance(a_lo, b_lo, b_hi),
        algorithm.min(
            edgedistances.getDistance(a_hi, b_lo, b_hi),
            algorithm.min(
                edgedistances.getDistance(b_lo, a_lo, a_hi),
                edgedistances.getDistance(b_hi, a_lo, a_hi))));
  }

  // Returns the minimum distance (measured along the surface of the sphere)
  // from a given point to the rectangle (both its boundary and its interior).
  // The latlng must be valid.
  S1Angle getDistance(in S2LatLng p) const {
    // The algorithm here is the same as in GetDistance(S2LatLngRect), only
    // with simplified calculations.
    const S2LatLngRect a = this;
    if (a.isEmpty()) logger.logError("Empty S2LatLngRect in S2LatLngRect::GetDistance: ", a);
    if (!p.isValid()) logger.logError("Invalid S2LatLng in S2LatLngRect::GetDistance: ", p);

    if (a.lng().contains(p.lng().radians())) {
      return S1Angle.fromRadians(
          algorithm.max(
              0.0,
              algorithm.max(
                  p.lat().radians() - a.lat().hi(),
                  a.lat().lo() - p.lat().radians())));
    }

    auto interval = S1Interval(a.lng().hi(), a.lng().getComplementCenter());
    double a_lng;
    if (interval.contains(p.lng().radians())) {
      a_lng = a.lng().hi();
    } else {
      a_lng = a.lng().lo();
    }
    S2Point lo = S2LatLng.fromRadians(a.lat().lo(), a_lng).toS2Point();
    S2Point hi = S2LatLng.fromRadians(a.lat().hi(), a_lng).toS2Point();
    return edgedistances.getDistance(p.toS2Point(), lo, hi);
  }

  // Returns the (directed or undirected) Hausdorff distance (measured along the
  // surface of the sphere) to the given S2LatLngRect. The directed Hausdorff
  // distance from rectangle A to rectangle B is given by
  //     h(A, B) = max_{p in A} min_{q in B} d(p, q).
  // The Hausdorff distance between rectangle A and rectangle B is given by
  //     H(A, B) = max{h(A, B), h(B, A)}.
  S1Angle getHausdorffDistance(in S2LatLngRect other) const {
    return algorithm.max(
        getDirectedHausdorffDistance(other),
        other.getDirectedHausdorffDistance(this));
  }

  S1Angle getDirectedHausdorffDistance(in S2LatLngRect other) const {
    if (isEmpty()) {
      return S1Angle.fromRadians(0);
    }
    if (other.isEmpty()) {
      return S1Angle.fromRadians(S1Angle.PI);  // maximum possible distance on S2
    }

    double lng_distance = lng().getDirectedHausdorffDistance(other.lng());
    assert(lng_distance >= 0);
    return getDirectedHausdorffDistance(lng_distance, lat(), other.lat());
  }

  // Returns true if two rectangles contains the same set of points.
  override
  bool opEquals(in Object o) const {
    S2LatLngRect other = cast(S2LatLngRect) o;
    if (other) {
      return lat() == other.lat() && lng() == other.lng();
    } else {
      return false;
    }
  }

  // Returns true if the latitude and longitude intervals of the two rectangles
  // are the same up to the given tolerance (see r1interval.h and s1interval.h
  // for details).
  bool approxEquals(in S2LatLngRect other, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    return _lat.approxEquals(other._lat, max_error.radians())
        && _lng.approxEquals(other._lng, max_error.radians());
  }

  // ApproxEquals() with separate tolerances for latitude and longitude.
  bool approxEquals(in S2LatLngRect other, in S2LatLng max_error) const {
    return _lat.approxEquals(other._lat, max_error.lat().radians())
        && _lng.approxEquals(other._lng, max_error.lng().radians());
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  override
  S2LatLngRect clone() const {
    return new S2LatLngRect(this);
  }

  override
  S2Cap getCapBound() const {
    // We consider two possible bounding caps, one whose axis passes
    // through the center of the lat-long rectangle and one whose axis
    // is the north or south pole.  We return the smaller of the two caps.

    if (isEmpty()) return S2Cap.empty();

    double pole_z, pole_angle;
    if (_lat.lo() + _lat.hi() < 0) {
      // South pole axis yields smaller cap.
      pole_z = -1;
      pole_angle = S1Angle.PI_2 + _lat.hi();
    } else {
      pole_z = 1;
      pole_angle = S1Angle.PI_2 - _lat.lo();
    }
    S2Cap pole_cap = new S2Cap(S2Point(0, 0, pole_z), S1Angle.fromRadians(pole_angle));

    // For bounding rectangles that span 180 degrees or less in longitude, the
    // maximum cap size is achieved at one of the rectangle vertices.  For
    // rectangles that are larger than 180 degrees, we punt and always return a
    // bounding cap centered at one of the two poles.
    double lng_span = _lng.hi() - _lng.lo();
    if (math.remainder(lng_span, 2 * S1Angle.PI) >= 0 && lng_span < 2 * S1Angle.PI) {
      auto mid_cap = new S2Cap(getCenter().toS2Point(), S1Angle.fromRadians(0));
      for (int k = 0; k < 4; ++k) {
        mid_cap.addPoint(getVertex(k).toS2Point());
      }
      if (mid_cap.height() < pole_cap.height())
        return mid_cap;
    }
    return pole_cap;
  }

  override
  S2LatLngRect getRectBound() const {
    return new S2LatLngRect(this);
  }

  override
  void getCellUnionBound(out S2CellId[] cellIds) const {
    return getCapBound().getCellUnionBound(cellIds);
  }

  override
  bool contains(in S2Cell cell) const {
    // A latitude-longitude rectangle contains a cell if and only if it contains
    // the cell's bounding rectangle.  This test is exact from a mathematical
    // point of view, assuming that the bounds returned by S2Cell::GetRectBound()
    // are tight.  However, note that there can be a loss of precision when
    // converting between representations -- for example, if an S2Cell is
    // converted to a polygon, the polygon's bounding rectangle may not contain
    // the cell's bounding rectangle.  This has some slightly unexpected side
    // effects; for instance, if one creates an S2Polygon from an S2Cell, the
    // polygon will contain the cell, but the polygon's bounding box will not.
    return contains(cell.getRectBound());
  }

  // This test is cheap but is NOT exact.  Use Intersects() if you want a more
  // accurate and more expensive test.  Note that when this method is used by
  // an S2RegionCoverer, the accuracy isn't all that important since if a cell
  // may intersect the region then it is subdivided, and the accuracy of this
  // method goes up as the cells get smaller.
  override
  bool mayIntersect(in S2Cell cell) const {
    // This test is cheap but is NOT exact (see s2latlng_rect.h).
    return intersects(cell.getRectBound());
  }

  // The point 'p' does not need to be normalized.
  override
  bool contains(in S2Point p) const  {
    return contains(S2LatLng(p));
  }

  // Appends a serialized representation of the S2LatLngRect to "encoder".
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  //void Encode(Encoder* const encoder) const;

  // Decodes an S2LatLngRect encoded with Encode().  Returns true on success.
  //bool Decode(Decoder* const decoder);

  // Returns true if the edge AB intersects the given edge of constant
  // longitude.
  static bool intersectsLngEdge(
      in S2Point a, in S2Point b, in R1Interval lat, double lng) {
    // Return true if the segment AB intersects the given edge of constant
    // longitude.  The nice thing about edges of constant longitude is that
    // they are straight lines on the sphere (geodesics).

    return crossingSign(
        a, b, S2LatLng.fromRadians(lat.lo(), lng).toS2Point(),
        S2LatLng.fromRadians(lat.hi(), lng).toS2Point()) > 0;
  }

  // Returns true if the edge AB intersects the given edge of constant
  // latitude.  Requires the vectors to have unit length.
  static bool intersectsLatEdge(
      in S2Point a, in S2Point b, double lat, in S1Interval lng)
  in {
    assert(isUnitLength(a));
    assert(isUnitLength(b));
  } body {
    // Return true if the segment AB intersects the given edge of constant
    // latitude.  Unfortunately, lines of constant latitude are curves on
    // the sphere.  They can intersect a straight edge in 0, 1, or 2 points.

    // First, compute the normal to the plane AB that points vaguely north.
    Vector3_d z = robustCrossProd(a, b).normalize();
    if (z[2] < 0) z = -z;

    // Extend this to an orthonormal frame (x,y,z) where x is the direction
    // where the great circle through AB achieves its maximium latitude.
    Vector3_d y = robustCrossProd(z, S2Point(0, 0, 1)).normalize();
    Vector3_d x = y.crossProd(z);
    enforce(isUnitLength(x));
    enforce(x[2] >= 0);

    // Compute the angle "theta" from the x-axis (in the x-y plane defined
    // above) where the great circle intersects the given line of latitude.
    double sin_lat = math.sin(lat);
    if (math.fabs(sin_lat) >= x[2]) {
      return false;  // The great circle does not reach the given latitude.
    }
    enforce(x[2] > 0);
    double cos_theta = sin_lat / x[2];
    double sin_theta = math.sqrt(1 - cos_theta * cos_theta);
    double theta = math.atan2(sin_theta, cos_theta);

    // The candidate intersection points are located +/- theta in the x-y
    // plane.  For an intersection to be valid, we need to check that the
    // intersection point is contained in the interior of the edge AB and
    // also that it is contained within the given longitude interval "lng".

    // Compute the range of theta values spanned by the edge AB.
    S1Interval ab_theta = S1Interval.fromPointPair(
        math.atan2(a.dotProd(y), a.dotProd(x)),
        math.atan2(b.dotProd(y), b.dotProd(x)));

    if (ab_theta.contains(theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      S2Point isect = x * cos_theta + y * sin_theta;
      if (lng.contains(math.atan2(isect[1], isect[0]))) return true;
    }
    if (ab_theta.contains(-theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      S2Point isect = x * cos_theta - y * sin_theta;
      if (lng.contains(math.atan2(isect[1], isect[0]))) return true;
    }
    return false;
  }

  override
  string toString() const {
    return "[Lo" ~ to!string(lo()) ~ ", Hi" ~ to!string(hi()) ~ "]";
  }

 private:
  /**
   * Return the directed Hausdorff distance from one longitudinal edge spanning
   * latitude range 'a_lat' to the other longitudinal edge spanning latitude
   * range 'b_lat', with their longitudinal difference given by 'lng_diff'.
   */
  static S1Angle getDirectedHausdorffDistance(double lng_diff, in R1Interval a, in R1Interval b)
  in {
    assert(lng_diff >= 0);
    assert(lng_diff <= S1Angle.PI);
  } body {
    // By symmetry, we can assume a's longtitude is 0 and b's longtitude is
    // lng_diff. Call b's two endpoints b_lo and b_hi. Let H be the hemisphere
    // containing a and delimited by the longitude line of b. The Voronoi diagram
    // of b on H has three edges (portions of great circles) all orthogonal to b
    // and meeting at b_lo cross b_hi.
    // E1: (b_lo, b_lo cross b_hi)
    // E2: (b_hi, b_lo cross b_hi)
    // E3: (-b_mid, b_lo cross b_hi), where b_mid is the midpoint of b
    //
    // They subdivide H into three Voronoi regions. Depending on how longitude 0
    // (which contains edge a) intersects these regions, we distinguish two cases:
    // Case 1: it intersects three regions. This occurs when lng_diff <= M_PI_2.
    // Case 2: it intersects only two regions. This occurs when lng_diff > M_PI_2.
    //
    // In the first case, the directed Hausdorff distance to edge b can only be
    // realized by the following points on a:
    // A1: two endpoints of a.
    // A2: intersection of a with the equator, if b also intersects the equator.
    //
    // In the second case, the directed Hausdorff distance to edge b can only be
    // realized by the following points on a:
    // B1: two endpoints of a.
    // B2: intersection of a with E3
    // B3: farthest point from b_lo to the interior of D, and farthest point from
    //     b_hi to the interior of U, if any, where D (resp. U) is the portion
    //     of edge a below (resp. above) the intersection point from B2.

    if (lng_diff == 0) {
      return S1Angle.fromRadians(a.getDirectedHausdorffDistance(b));
    }

    // Assumed longtitude of b.
    double b_lng = lng_diff;
    // Two endpoints of b.
    S2Point b_lo = S2LatLng.fromRadians(b.lo(), b_lng).toS2Point();
    S2Point b_hi = S2LatLng.fromRadians(b.hi(), b_lng).toS2Point();

    // Handling of each case outlined at the top of the function starts here.
    // This is initialized a few lines below.
    S1Angle max_distance;

    // Cases A1 and B1.
    S2Point a_lo = S2LatLng.fromRadians(a.lo(), 0).toS2Point();
    S2Point a_hi = S2LatLng.fromRadians(a.hi(), 0).toS2Point();
    max_distance = edgedistances.getDistance(a_lo, b_lo, b_hi);
    max_distance = algorithm.max(max_distance, edgedistances.getDistance(a_hi, b_lo, b_hi));

    if (lng_diff <= S1Angle.PI_2) {
      // Case A2.
      if (a.contains(0) && b.contains(0)) {
        max_distance = algorithm.max(max_distance, S1Angle.fromRadians(lng_diff));
      }
    } else {
      // Case B2.
      const S2Point p = getBisectorIntersection(b, b_lng);
      double p_lat = S2LatLng.latitude(p).radians();
      if (a.contains(p_lat)) {
        max_distance = algorithm.max(max_distance, S1Angle(p, b_lo));
      }

      // Case B3.
      if (p_lat > a.lo()) {
        max_distance = algorithm.max(
            max_distance,
            getInteriorMaxDistance(R1Interval(a.lo(), algorithm.min(p_lat, a.hi())), b_lo));
      }
      if (p_lat < a.hi()) {
        max_distance = algorithm.max(
            max_distance,
            getInteriorMaxDistance(R1Interval(algorithm.max(p_lat, a.lo()), a.hi()), b_hi));
      }
    }

    return max_distance;
  }

  /**
   * Return max distance from a point b to the segment spanning latitude range
   * a_lat on longitude 0, if the max occurs in the interior of a_lat. Otherwise
   * return -1.
   */
  static S1Angle getInteriorMaxDistance(in R1Interval a_lat, in S2Point b) {
    // Longitude 0 is in the y=0 plane. b.x() >= 0 implies that the maximum
    // does not occur in the interior of a_lat.
    if (a_lat.isEmpty() || b.x() >= 0) return S1Angle.fromRadians(-1);

    // Project b to the y=0 plane. The antipodal of the normalized projection is
    // the point at which the maxium distance from b occurs, if it is contained
    // in a_lat.
    S2Point intersection_point = S2Point(-b.x(), 0, -b.z()).normalize();
    if (a_lat.interiorContains(S2LatLng.latitude(intersection_point).radians())) {
      return S1Angle(b, intersection_point);
    } else {
      return S1Angle.fromRadians(-1);
    }
  }


  /**
   * Return the intersection of longitude 0 with the bisector of an edge
   * on longitude 'lng' and spanning latitude range 'lat'.
   */
  static S2Point getBisectorIntersection(in R1Interval lat, double lng) {
    lng = math.fabs(lng);
    double lat_center = lat.getCenter();
    // A vector orthogonal to the bisector of the given longitudinal edge.
    S2LatLng ortho_bisector;
    if (lat_center >= 0) {
      ortho_bisector = S2LatLng.fromRadians(lat_center - S1Angle.PI_2, lng);
    } else {
      ortho_bisector = S2LatLng.fromRadians(-lat_center - S1Angle.PI_2, lng - S1Angle.PI);
    }
    // A vector orthogonal to longitude 0.
    static const S2Point ortho_lng = S2Point(0, -1, 0);
    return robustCrossProd(ortho_lng, ortho_bisector.toS2Point());
  }

  R1Interval _lat;
  S1Interval _lng;
}
