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

module s2.r1interval;

import s2.util.math.vector;
import algorithm = std.algorithm;
import conv = std.conv;
import math = std.math;

// An R1Interval represents a closed, bounded interval on the real line.
// It is capable of representing the empty interval (containing no points)
// and zero-length intervals (containing a single point).
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
struct R1Interval {
private:
  // The default constructor creates an empty interval.  (Any interval where
  // lo > hi is considered to be empty.)
  //
  // Note: Don't construct an interval using the default constructor and
  // set_lo()/set_hi(), since this technique doesn't work with S1Interval and
  // is bad programming style anyways.  If you need to set both endpoints, use
  // the constructor below:
  //
  //   lat_bounds_ = R1Interval(lat_lo, lat_hi);
  Vector2_d _bounds = Vector2_d([1, 0]);

public:
  // Constructor.  If lo > hi, the interval is empty.
  this(double lo, double hi) {
    _bounds = Vector2_d([lo, hi]);
  }

  // Returns an empty interval.
  static R1Interval empty() {
    return R1Interval();
  }

  // Convenience method to construct an interval containing a single point.
  static R1Interval fromPoint(double p) {
    return R1Interval(p, p);
  }

  // Convenience method to construct the minimal interval containing
  // the two given points.  This is equivalent to starting with an empty
  // interval and calling AddPoint() twice, but it is more efficient.
  static R1Interval fromPointPair(double p1, double p2) {
    if (p1 <= p2) {
      return R1Interval(p1, p2);
    } else {
      return R1Interval(p2, p1);
    }
  }

  // Accessors methods.
  @property
  double lo() const {
    return _bounds[0];
  }

  @property
  double hi() const {
    return _bounds[1];
  }

  // Methods to modify one endpoint of an existing R1Interval.  Do not use
  // these methods if you want to replace both endpoints of the interval; use
  // a constructor instead.  For example:
  //
  //   *lat_bounds = R1Interval(lat_lo, lat_hi);
  void setLo(double p) {
    _bounds[0] = p;
  }

  void setHi(double p) {
    _bounds[1] = p;
  }

  // Methods that allow the R1Interval to be accessed as a vector.  (The
  // recommended style is to use lo() and hi() whenever possible, but these
  // methods are useful when the endpoint to be selected is not constant.)
  ref inout(double) opIndex(size_t i) inout {
    return _bounds[i];
  }

  void opIndexAssign(in double v, size_t i) {
    _bounds[i] = v;
  }

  Vector2_d bounds() const {
    return _bounds;
  }

  ref Vector2_d mutableBounds() {
    return _bounds;
  }

  // Return true if the interval is empty, i.e. it contains no points.
  bool isEmpty() const {
    return lo() > hi();
  }

  // Return the center of the interval.  For empty intervals,
  // the result is arbitrary.
  double getCenter() const {
    return 0.5 * (lo() + hi());
  }

  // Return the length of the interval.  The length of an empty interval
  // is negative.
  double getLength() const {
    return hi() - lo();
  }

  bool contains(double p) const {
    return p >= lo() && p <= hi();
  }

  bool interiorContains(double p) const {
    return p > lo() && p < hi();
  }

  // Return true if this interval contains the interval 'y'.
  bool contains(in R1Interval y) const {
    if (y.isEmpty()) {
      return true;
    }
    return y.lo() >= lo() && y.hi() <= hi();
  }

  // Return true if the interior of this interval contains the entire
  // interval 'y' (including its boundary).
  bool interiorContains(in R1Interval y) const {
    if (y.isEmpty()) {
      return true;
    }
    return y.lo() > lo() && y.hi() < hi();
  }

  // Return true if this interval intersects the given interval,
  // i.e. if they have any points in common.
  bool intersects(in R1Interval y) const {
    if (lo() <= y.lo()) {
      return y.lo() <= hi() && y.lo() <= y.hi();
    } else {
      return lo() <= y.hi() && lo() <= hi();
    }
  }

  // Return true if the interior of this interval intersects
  // any point of the given interval (including its boundary).
  bool interiorIntersects(in R1Interval y) const {
    return y.lo() < hi() && lo() < y.hi() && lo() < hi() && y.lo() <= y.hi();
  }

  // Return the Hausdorff distance to the given interval 'y'. For two
  // R1Intervals x and y, this distance is defined as
  //     h(x, y) = max_{p in x} min_{q in y} d(p, q).
  double getDirectedHausdorffDistance(in R1Interval y) const {
    if (isEmpty()) {
      return 0.0;
    } else if (y.isEmpty()) {
      return double.max;
    } else {
      return algorithm.max(0.0, algorithm.max(hi() - y.hi(), y.lo() - lo()));
    }
  }

  // Expand the interval so that it contains the given point "p".
  void addPoint(double p) {
    if (isEmpty()) {
      setLo(p);
      setHi(p);
    } else if (p < lo()) {
      setLo(p);
    } else if (p > hi()) {
      setHi(p);
    }
  }

  // Expand the interval so that it contains the given interval "y".
  void addInterval(in R1Interval y) {
    if (y.isEmpty()) {
      return;
    } else if (isEmpty()) {
      this = y;
      return;
    }
    if (y.lo() < lo()) {
      setLo(y.lo());
    }
    if (y.hi() > hi()) {
      setHi(y.hi());
    }
  }

  // Return the closest point in the interval to the given point "p".
  // The interval must be non-empty.
  double project(double p) const
  in {
    assert(!isEmpty());
  } body {
    return algorithm.max(lo(), algorithm.min(hi(), p));
  }

  // Return an interval that has been expanded on each side by the given
  // distance "margin".  If "margin" is negative, then shrink the interval on
  // each side by "margin" instead.  The resulting interval may be empty.  Any
  // expansion of an empty interval remains empty.
  R1Interval expanded(double margin) const {
    if (isEmpty()) {
      return this;
    }
    return R1Interval(lo() - margin, hi() + margin);
  }

  // Return the smallest interval that contains this interval and the
  // given interval "y".
  R1Interval unite(in R1Interval y) const {
    if (isEmpty()) {
      return y;
    } else if (y.isEmpty()) {
      return this;
    } else {
      return R1Interval(algorithm.min(lo(), y.lo()), algorithm.max(hi(), y.hi()));
    }
  }

  // Return the intersection of this interval with the given interval.
  // Empty intervals do not need to be special-cased.
  R1Interval intersection(in R1Interval y) const {
    return R1Interval(algorithm.max(lo(), y.lo()), algorithm.min(hi(), y.hi()));
  }

  // Support the == and != operators. Return true if two intervals contain the same set of points.
  bool opEquals(in R1Interval y) const {
    return (lo() == y.lo() && hi() == y.hi()) || (isEmpty() && y.isEmpty());
  }

  // Return true if this interval can be transformed into the given interval
  // by moving each endpoint by at most "max_error".  The empty interval is
  // considered to be positioned arbitrarily on the real line, thus any
  // interval with (length <= 2*max_error) matches the empty interval.
  bool approxEquals(in R1Interval y, double max_error = 1e-15) const {
    if (isEmpty()) {
      return y.getLength() <= 2 * max_error;
    }
    if (y.isEmpty()) {
      return getLength() <= 2 * max_error;
    }
    return (math.fabs(y.lo() - lo()) <= max_error
        && math.fabs(y.hi() - hi()) <= max_error);
  }

  string toString() const {
    return "[" ~ conv.to!string(lo()) ~ ", " ~ conv.to!string(hi()) ~ "]";
  }

 private:
  Vector2_d bounds_;
}
