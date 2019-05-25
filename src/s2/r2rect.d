// Copyright 2012 Google Inc. All Rights Reserved.
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

module s2.r2rect;

import s2.r1interval;
import s2.r2point;
import s2.logger;
import conv = std.conv;

/**
 * An R2Rect represents a closed axis-aligned rectangle in the (x,y) plane.
 *
 * This class is intended to be copied by value as desired.  It uses
 * the default copy constructor and assignment operator, however it is
 * not a "plain old datatype" (POD) because it has virtual functions.
 */
struct R2Rect {
private:
  R1Interval[2] _bounds;

public:
  // Construct a rectangle from the given lower-left and upper-right points.
  this(in R2Point lo, in R2Point hi)
  out {
    assert(isValid());
  } do {
    _bounds[0] = R1Interval(lo.x(), hi.x());
    _bounds[1] = R1Interval(lo.y(), hi.y());
  }

  // Construct a rectangle from the given intervals in x and y.  The two
  // intervals must either be both empty or both non-empty.
  this(in R1Interval x, in R1Interval y)
  out {
    if (!isValid()) logger.logError("Invalid: " ~ toString());
  } do {
    _bounds[0] = x;
    _bounds[1] = y;
  }

  // The canonical empty rectangle.  Use isEmpty() to test for empty
  // rectangles, since they have more than one representation.
  static R2Rect empty() {
    return R2Rect(R1Interval.empty(), R1Interval.empty());
  }

  // Construct a rectangle from a center point and size in each dimension.
  // Both components of size should be non-negative, i.e. this method cannot
  // be used to create an empty rectangle.
  static R2Rect fromCenterSize(in R2Point center, in R2Point size) {
    return R2Rect(
        R1Interval(center.x() - 0.5 * size.x(), center.x() + 0.5 * size.x()),
        R1Interval(center.y() - 0.5 * size.y(), center.y() + 0.5 * size.y()));
  }

  // Convenience method to construct a rectangle containing a single point.
  static R2Rect fromPoint(in R2Point p) {
    return R2Rect(p, p);
  }

  // Convenience method to construct the minimal bounding rectangle containing
  // the two given points.  This is equivalent to starting with an empty
  // rectangle and calling AddPoint() twice.  Note that it is different than
  // the R2Rect(lo, hi) constructor, where the first point is always
  // used as the lower-left corner of the resulting rectangle.
  static R2Rect fromPointPair(in R2Point p1, in R2Point p2) {
    return R2Rect(
        R1Interval.fromPointPair(p1.x(), p2.x()),
        R1Interval.fromPointPair(p1.y(), p2.y()));
  }

  // Accessor methods.
  @property
  R1Interval x() const {
    return _bounds[0];
  }

  @property
  R1Interval y() const {
    return _bounds[1];
  }

  @property
  R2Point lo() const {
    return R2Point(x().lo(), y().lo());
  }

  @property
  R2Point hi() const {
    return R2Point(x().hi(), y().hi());
  }

  // Methods that allow the R2Rect to be accessed as a vector.
  ref inout(R1Interval) opIndex(size_t i) inout {
    return _bounds[i];
  }

  void opIndexAssign(in R1Interval value, size_t i)
  in {
    assert(i >= 0 && i < 2);
  } do {
    _bounds[i] = value;
  }

  // Return true if the rectangle is valid, which essentially just means
  // that if the bound for either axis is empty then both must be.
  bool isValid() const {
    // The x/y ranges must either be both empty or both non-empty.
    return x().isEmpty() == y().isEmpty();
  }

  // Return true if the rectangle is empty, i.e. it contains no points at all.
  bool isEmpty() const {
    return x().isEmpty();
  }

  // Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
  // Vertex 0 is in the lower-left corner.  For convenience, the argument is
  // reduced modulo 4 to the range [0..3].
  R2Point getVertex(int k) const {
    // Twiddle bits to return the points in CCW order (lower left, lower right,
    // upper right, upper left).
    int j = (k >> 1) & 1;
    return getVertex(j ^ (k & 1), j);
  }

  // Return the vertex in direction "i" along the x-axis (0=left, 1=right) and
  // direction "j" along the y-axis (0=down, 1=up).  Equivalently, return the
  // vertex constructed by selecting endpoint "i" of the x-interval (0=lo,
  // 1=hi) and vertex "j" of the y-interval.
  R2Point getVertex(int i, int j) const {
    return R2Point(_bounds[0][i], _bounds[1][j]);
  }

  // Return the center of the rectangle in (x,y)-space.
  R2Point getCenter() const {
    return R2Point(x().getCenter(), y().getCenter());
  }

  // Return the width and height of this rectangle in (x,y)-space.  Empty
  // rectangles have a negative width and height.
  R2Point getSize() const {
    return R2Point(x().getLength(), y().getLength());
  }

  // Return true if the rectangle contains the given point.  Note that
  // rectangles are closed regions, i.e. they contain their boundary.
  bool contains(in R2Point p) const {
    return x().contains(p.x()) && y().contains(p.y());
  }

  // Return true if and only if the given point is contained in the interior
  // of the region (i.e. the region excluding its boundary).
  bool interiorContains(in R2Point p) const {
    return x().interiorContains(p.x()) && y().interiorContains(p.y());
  }

  // Return true if and only if the rectangle contains the given other
  // rectangle.
  bool contains(in R2Rect other) const {
    return x().contains(other.x()) && y().contains(other.y());
  }

  // Return true if and only if the interior of this rectangle contains all
  // points of the given other rectangle (including its boundary).
  bool interiorContains(in R2Rect other) const {
    return x().interiorContains(other.x()) && y().interiorContains(other.y());
  }

  // Return true if this rectangle and the given other rectangle have any
  // points in common.
  bool intersects(in R2Rect other) const {
    return x().intersects(other.x()) && y().intersects(other.y());
  }

  // Return true if and only if the interior of this rectangle intersects
  // any point (including the boundary) of the given other rectangle.
  bool InteriorIntersects(in R2Rect other) const {
    return x().interiorIntersects(other.x()) && y().interiorIntersects(other.y());
  }

  // Expand the rectangle to include the given point.  The rectangle is
  // expanded by the minimum amount possible.
  void addPoint(in R2Point p) {
    _bounds[0].addPoint(p[0]);
    _bounds[1].addPoint(p[1]);
  }

  // Expand the rectangle to include the given other rectangle.  This is the
  // same as replacing the rectangle by the union of the two rectangles, but
  // is somewhat more efficient.
  void addRect(in R2Rect other) {
    _bounds[0].addInterval(other[0]);
    _bounds[1].addInterval(other[1]);
  }

  // Return the closest point in the rectangle to the given point "p".
  // The rectangle must be non-empty.
  R2Point Project(in R2Point p) const {
    return R2Point(x().project(p.x()), y().project(p.y()));
  }


  // Return a rectangle that has been expanded on each side in the x-direction
  // by margin.x(), and on each side in the y-direction by margin.y().  If
  // either margin is empty, then shrink the interval on the corresponding
  // sides instead.  The resulting rectangle may be empty.  Any expansion of
  // an empty rectangle remains empty.
  R2Rect expanded(in R2Point margin) const {
    R1Interval xx = x().expanded(margin.x());
    R1Interval yy = y().expanded(margin.y());
    if (xx.isEmpty() || yy.isEmpty()) return empty();
    return R2Rect(xx, yy);
  }

  R2Rect expanded(double margin) const {
    return expanded(R2Point(margin, margin));
  }

  // Return the smallest rectangle containing the union of this rectangle and
  // the given rectangle.
  R2Rect unite(in R2Rect other) const {
    return R2Rect(x().unite(other.x()), y().unite(other.y()));
  }

  // Return the smallest rectangle containing the intersection of this
  // rectangle and the given rectangle.
  R2Rect intersection(in R2Rect other) const {
    R1Interval xx = x().intersection(other.x());
    R1Interval yy = y().intersection(other.y());
    if (xx.isEmpty() || yy.isEmpty()) return empty();
    return R2Rect(xx, yy);
  }

  // Return true if two rectangles contains the same set of points.
  bool opEquals(in R2Rect other) const {
    return x() == other.x() && y() == other.y();
  }

  // Return true if the x- and y-intervals of the two rectangles are the same
  // up to the given tolerance (see r1interval.h for details).
  bool approxEquals(in R2Rect other, double max_error = 1e-15) const {
    return x().approxEquals(other.x(), max_error)
        && y().approxEquals(other.y(), max_error);
  }

  string toString() const {
    return "[Lo" ~ lo().toString() ~ ", Hi" ~ hi().toString() ~ "]";
  }

}
