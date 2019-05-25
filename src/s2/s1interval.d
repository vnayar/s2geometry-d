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

module s2.s1interval;

import s2.util.math.vector;
import algorithm = std.algorithm;
import conv = std.conv;
import math = std.math;
import s2.util.math.s2const;

// An S1Interval represents a closed interval on a unit circle (also known
// as a 1-dimensional sphere).  It is capable of representing the empty
// interval (containing no points), the full interval (containing all
// points), and zero-length intervals (containing a single point).
//
// Points are represented by the angle they make with the positive x-axis in
// the range [-Pi, Pi].  An interval is represented by its lower and upper
// bounds (both inclusive, since the interval is closed).  The lower bound may
// be greater than the upper bound, in which case the interval is "inverted"
// (i.e. it passes through the point (-1, 0)).
//
// Note that the point (-1, 0) has two valid representations, Pi and -Pi.
// The normalized representation of this point internally is Pi, so that
// endpoints of normal intervals are in the range (-Pi, Pi].  However, we
// take advantage of the point -Pi to construct two special intervals:
// the Full() interval is [-Pi, Pi], and the Empty() interval is [Pi, -Pi].
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
struct S1Interval {
private:
  enum ArgsChecked { ARGS_CHECKED }

  // Internal constructor that assumes that both arguments are in the
  // correct range, i.e. normalization from -Pi to Pi is already done.
  this(double lo, double hi, ArgsChecked dummy)
  out {
    assert(isValid());
  } do {
    _bounds = Vector2_d(lo, hi);
  }

  // Return true if the interval (which is closed) contains the point 'p'.
  // Skips the normalization of 'p' from -Pi to Pi.
  bool fastContains(double p) const {
    if (isInverted()) {
      return (p >= lo() || p <= hi()) && !isEmpty();
    } else {
      return p >= lo() && p <= hi();
    }
  }

  static double positiveDistance(double a, double b) {
    // Compute the distance from "a" to "b" in the range [0, 2*Pi).
    // This is equivalent to (remainder(b - a - M_PI, 2 * M_PI) + M_PI),
    // except that it is more numerically stable (it does not lose
    // precision for very small positive distances).
    double d = b - a;
    if (d >= 0) {
      return d;
    }
    // We want to ensure that if b == Pi and a == (-Pi + eps),
    // the return result is approximately 2*Pi and not zero.
    return (b + M_PI) - (a - M_PI);
  }

  Vector2_d _bounds = Vector2_d(M_PI, -M_PI);

public:
  // Constructor.  Both endpoints must be in the range -Pi to Pi inclusive.
  // The value -Pi is converted internally to Pi except for the Full()
  // and Empty() intervals.
  this(double lo, double hi)
  out {
    assert(isValid());
  } do {
    _bounds = Vector2_d(lo, hi);
    if (lo == -M_PI && hi != M_PI) {
      setLo(M_PI);
    }
    if (hi == -M_PI && lo != M_PI) {
      setHi(M_PI);
    }
  }

  // Returns the empty interval.
  static S1Interval empty() {
    return S1Interval();
  }

  // Returns the full interval.
  static S1Interval full() {
    return S1Interval(-M_PI, M_PI, ArgsChecked.ARGS_CHECKED);
  }

  // Convenience method to construct an interval containing a single point.
  static S1Interval fromPoint(double p) {
    if (p == -M_PI) {
      p = M_PI;
    }
    return S1Interval(p, p, ArgsChecked.ARGS_CHECKED);
  }

  // Convenience method to construct the minimal interval containing
  // the two given points.  This is equivalent to starting with an empty
  // interval and calling AddPoint() twice, but it is more efficient.
  static S1Interval fromPointPair(double p1, double p2)
  in {
    assert(math.fabs(p1) <= M_PI);
    assert(math.fabs(p2) <= M_PI);
  } do {
    if (p1 == -M_PI) {
      p1 = M_PI;
    }
    if (p2 == -M_PI) {
      p2 = M_PI;
    }
    if (positiveDistance(p1, p2) <= M_PI) {
      return S1Interval(p1, p2, ArgsChecked.ARGS_CHECKED);
    } else {
      return S1Interval(p2, p1, ArgsChecked.ARGS_CHECKED);
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

  // Methods that allow the S1Interval to be accessed as a vector.  (The
  // recommended style is to use lo() and hi() whenever possible, but these
  // methods are useful when the endpoint to be selected is not constant.)
  //
  // Only const versions of these methods are provided, since S1Interval
  // has invariants that must be maintained after each update.
  double opIndex(size_t i) const {
    return _bounds[i];
  }

  @property
  Vector2_d bounds() const {
    return _bounds;
  }

  // An interval is valid if neither bound exceeds Pi in absolute value,
  // and the value -Pi appears only in the Empty() and Full() intervals.
  bool isValid() const {
    return (math.fabs(lo()) <= M_PI && math.fabs(hi()) <= M_PI
        && !(lo() == -M_PI && hi() != M_PI)
        && !(hi() == -M_PI && lo() != M_PI));
  }

  // Return true if the interval contains all points on the unit circle.
  bool isFull() const {
    return lo() == -M_PI && hi() == M_PI;
  }

  // Return true if the interval is empty, i.e. it contains no points.
  bool isEmpty() const {
    return lo() == M_PI && hi() == -M_PI;
  }

  // Return true if lo() > hi().  (This is true for empty intervals.)
  bool isInverted() const {
    return lo() > hi();
  }

  // Return the midpoint of the interval.  For full and empty intervals,
  // the result is arbitrary.
  double getCenter() const {
    double center = 0.5 * (lo() + hi());
    if (!isInverted()) {
      return center;
    }
    // Return the center in the range (-Pi, Pi].
    return (center <= 0) ? (center + M_PI) : (center - M_PI);
  }

  // Return the length of the interval.  The length of an empty interval
  // is negative.
  double getLength() const {
    double length = hi() - lo();
    if (length >= 0) {
      return length;
    }
    length += 2 * M_PI;
    // Empty intervals have a negative length.
    return (length > 0) ? length : -1;
  }

  // Return the complement of the interior of the interval.  An interval and
  // its complement have the same boundary but do not share any interior
  // values.  The complement operator is not a bijection, since the complement
  // of a singleton interval (containing a single value) is the same as the
  // complement of an empty interval.
  S1Interval complement() const {
    if (lo() == hi()) {
      return full();   // Singleton.
    }
    return S1Interval(hi(), lo(), ArgsChecked.ARGS_CHECKED);  // Handles empty and full.
  }

  // Return the midpoint of the complement of the interval. For full and empty
  // intervals, the result is arbitrary. For a singleton interval (containing a
  // single point), the result is its antipodal point on S1.
  double getComplementCenter() const {
    if (lo() != hi()) {
      return complement().getCenter();
    } else {  // Singleton.
      return (hi() <= 0) ? (hi() + M_PI) : (hi() - M_PI);
    }
  }

  // Return true if the interval (which is closed) contains the point 'p'.
  bool contains(double p) const
  in {
    // Works for empty, full, and singleton intervals.
    assert(math.fabs(p) <= M_PI);
  } do {
    if (p == -M_PI) {
      p = M_PI;
    }
    return fastContains(p);
  }

  // Return true if the interior of the interval contains the point 'p'.
  bool interiorContains(double p) const
  in {
    // Works for empty, full, and singleton intervals.
    assert(math.fabs(p) <= M_PI);
  } do {
    if (p == -M_PI) {
      p = M_PI;
    }

    if (isInverted()) {
      return p > lo() || p < hi();
    } else {
      return (p > lo() && p < hi()) || isFull();
    }
  }

  // Return true if the interval contains the given interval 'y'.
  // Works for empty, full, and singleton intervals.
  bool contains(in S1Interval y) const {
    // It might be helpful to compare the structure of these tests to
    // the simpler Contains(double) method above.
    if (isInverted()) {
      if (y.isInverted()) {
        return y.lo() >= lo() && y.hi() <= hi();
      }
      return (y.lo() >= lo() || y.hi() <= hi()) && !isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return y.lo() >= lo() && y.hi() <= hi();
    }
  }

  // Returns true if the interior of this interval contains the entire
  // interval 'y'.  Note that x.InteriorContains(x) is true only when
  // x is the empty or full interval, and x.InteriorContains(S1Interval(p,p))
  // is equivalent to x.InteriorContains(p).
  bool interiorContains(in S1Interval y) const {
    if (isInverted()) {
      if (!y.isInverted()) {
        return y.lo() > lo() || y.hi() < hi();
      }
      return (y.lo() > lo() && y.hi() < hi()) || y.isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return (y.lo() > lo() && y.hi() < hi()) || isFull();
    }
  }

  // Return true if the two intervals contain any points in common.
  // Note that the point +/-Pi has two representations, so the intervals
  // [-Pi,-3] and [2,Pi] intersect, for example.
  bool intersects(in S1Interval y) const {
    if (isEmpty() || y.isEmpty()) {
      return false;
    }
    if (isInverted()) {
      // Every non-empty inverted interval contains Pi.
      return y.isInverted() || y.lo() <= hi() || y.hi() >= lo();
    } else {
      if (y.isInverted()) {
        return y.lo() <= hi() || y.hi() >= lo();
      }
      return y.lo() <= hi() && y.hi() >= lo();
    }
  }

  // Return true if the interior of this interval contains any point of the
  // interval 'y' (including its boundary).  Works for empty, full, and
  // singleton intervals.
  bool interiorIntersects(in S1Interval y) const {
    if (isEmpty() || y.isEmpty() || lo() == hi()) {
      return false;
    }
    if (isInverted()) {
      return y.isInverted() || y.lo() < hi() || y.hi() > lo();
    } else {
      if (y.isInverted()) {
        return y.lo() < hi() || y.hi() > lo();
      }
      return (y.lo() < hi() && y.hi() > lo()) || isFull();
    }
  }

  // Return the Hausdorff distance to the given interval 'y'. For two
  // S1Intervals x and y, this distance is defined by
  //     h(x, y) = max_{p in x} min_{q in y} d(p, q),
  // where d(.,.) is measured along S1.
  double getDirectedHausdorffDistance(in S1Interval y) const
  out (distance) {
    assert(distance >= 0);
  } do {
    if (y.contains(this)) {
      return 0.0;  // this includes the case this is empty
    }
    if (y.isEmpty()) {
      return M_PI;  // maximum possible distance on S1
    }

    double y_complement_center = y.getComplementCenter();
    if (contains(y_complement_center)) {
      return positiveDistance(y.hi(), y_complement_center);
    } else {
      // The Hausdorff distance is realized by either two hi() endpoints or two
      // lo() endpoints, whichever is farther apart.
      double hi_hi = S1Interval(y.hi(), y_complement_center).contains(hi()) ?
          positiveDistance(y.hi(), hi()) : 0;
      double lo_lo = S1Interval(y_complement_center, y.lo()).contains(lo()) ?
          positiveDistance(lo(), y.lo()) : 0;
      return algorithm.max(hi_hi, lo_lo);
    }
  }

  // Expand the interval by the minimum amount necessary so that it
  // contains the given point "p" (an angle in the range [-Pi, Pi]).
  void addPoint(double p)
  in {
    assert(math.fabs(p) <= M_PI);
  } do {
    if (p == -M_PI) {
      p = M_PI;
    }

    if (fastContains(p)) {
      return;
    }
    if (isEmpty()) {
      setHi(p);
      setLo(p);
    } else {
      // Compute distance from p to each endpoint.
      double dlo = positiveDistance(p, lo());
      double dhi = positiveDistance(hi(), p);
      if (dlo < dhi) {
        setLo(p);
      } else {
        setHi(p);
      }
      // Adding a point can never turn a non-full interval into a full one.
    }
  }

  // Return the closest point in the interval to the given point "p".
  // The interval must be non-empty.
  double project(double p) const
  in {
    assert(!isEmpty());
    assert(math.fabs(p) <= M_PI);
  } do {
    if (p == -M_PI) {
      p = M_PI;
    }
    if (fastContains(p)) {
      return p;
    }
    // Compute distance from p to each endpoint.
    double dlo = positiveDistance(p, lo());
    double dhi = positiveDistance(hi(), p);
    return (dlo < dhi) ? lo() : hi();
  }

  // Return an interval that has been expanded on each side by the given
  // distance "margin".  If "margin" is negative, then shrink the interval on
  // each side by "margin" instead.  The resulting interval may be empty or
  // full.  Any expansion (positive or negative) of a full interval remains
  // full, and any expansion of an empty interval remains empty.
  S1Interval expanded(double margin) const {
    if (margin >= 0) {
      if (isEmpty()) {
        return this;
      }
      // Check whether this interval will be full after expansion, allowing
      // for a 1-bit rounding error when computing each endpoint.
      if (getLength() + 2 * margin + 2 * double.epsilon >= 2 * M_PI) {
        return full();
      }
    } else {
      if (isFull()) {
        return this;
      }
      // Check whether this interval will be empty after expansion, allowing
      // for a 1-bit rounding error when computing each endpoint.
      if (getLength() + 2 * margin - 2 * double.epsilon <= 0) {
        return empty();
      }
    }
    S1Interval result = S1Interval(math.remainder(lo() - margin, 2 * M_PI),
        math.remainder(hi() + margin, 2 * M_PI));
    if (result.lo() <= -M_PI) {
      result.setLo(M_PI);
    }
    return result;
  }

  // Return the smallest interval that contains this interval and the
  // given interval "y".
  S1Interval unite(in S1Interval y) const {
    // The y.is_full() case is handled correctly in all cases by the code
    // below, but can follow three separate code paths depending on whether
    // this interval is inverted, is non-inverted but contains Pi, or neither.

    if (y.isEmpty()) {
      return this;
    }
    if (fastContains(y.lo())) {
      if (fastContains(y.hi())) {
        // Either this interval contains y, or the union of the two
        // intervals is the Full() interval.
        if (contains(y)) {
          return this;  // is_full() code path
        }
        return full();
      }
      return S1Interval(lo(), y.hi(), ArgsChecked.ARGS_CHECKED);
    }
    if (fastContains(y.hi())) {
      return S1Interval(y.lo(), hi(), ArgsChecked.ARGS_CHECKED);
    }

    // This interval contains neither endpoint of y.  This means that either y
    // contains all of this interval, or the two intervals are disjoint.
    if (isEmpty() || y.fastContains(lo())) {
      return y;
    }

    // Check which pair of endpoints are closer together.
    double dlo = positiveDistance(y.hi(), lo());
    double dhi = positiveDistance(hi(), y.lo());
    if (dlo < dhi) {
      return S1Interval(y.lo(), hi(), ArgsChecked.ARGS_CHECKED);
    } else {
      return S1Interval(lo(), y.hi(), ArgsChecked.ARGS_CHECKED);
    }
  }

  // Return the smallest interval that contains the intersection of this
  // interval with "y".  Note that the region of intersection may
  // consist of two disjoint intervals.
  S1Interval intersection(in S1Interval y) const {
    // The y.is_full() case is handled correctly in all cases by the code
    // below, but can follow three separate code paths depending on whether
    // this interval is inverted, is non-inverted but contains Pi, or neither.

    if (y.isEmpty()) {
      return empty();
    }
    if (fastContains(y.lo())) {
      if (fastContains(y.hi())) {
        // Either this interval contains y, or the region of intersection
        // consists of two disjoint subintervals.  In either case, we want
        // to return the shorter of the two original intervals.
        if (y.getLength() < getLength()) {
          return y;  // is_full() code path
        }
        return this;
      }
      return S1Interval(y.lo(), hi(), ArgsChecked.ARGS_CHECKED);
    }
    if (fastContains(y.hi())) {
      return S1Interval(lo(), y.hi(), ArgsChecked.ARGS_CHECKED);
    }

    // This interval contains neither endpoint of y.  This means that either y
    // contains all of this interval, or the two intervals are disjoint.

    if (y.fastContains(lo())) {
      return this;  // is_empty() okay here
    }
    assert(!intersects(y));
    return empty();

  }

  // Return true if two intervals contains the same set of points.
  bool opEquals(in S1Interval y) const {
    return lo() == y.lo() && hi() == y.hi();
  }


  // Return true if this interval can be transformed into the given interval by
  // moving each endpoint by at most "max_error" (and without the endpoints
  // crossing, which would invert the interval).  Empty and full intervals are
  // considered to start at an arbitrary point on the unit circle, thus any
  // interval with (length <= 2*max_error) matches the empty interval, and any
  // interval with (length >= 2*Pi - 2*max_error) matches the full interval.
  bool approxEquals(in S1Interval y, double max_error = 1e-15) const {
    // Full and empty intervals require special cases because the "endpoints"
    // are considered to be positioned arbitrarily.
    if (isEmpty()) {
      return y.getLength() <= 2 * max_error;
    }
    if (y.isEmpty()) {
      return getLength() <= 2 * max_error;
    }
    if (isFull()) {
      return y.getLength() >= 2 * (M_PI - max_error);
    }
    if (y.isFull()) {
      return getLength() >= 2 * (M_PI - max_error);
    }

    // The purpose of the last test below is to verify that moving the endpoints
    // does not invert the interval, e.g. [-1e20, 1e20] vs. [1e20, -1e20].
    return (math.fabs(math.remainder(y.lo() - lo(), 2 * M_PI)) <= max_error
        && math.fabs(math.remainder(y.hi() - hi(), 2 * M_PI)) <= max_error
        && math.fabs(getLength() - y.getLength()) <= 2 * max_error);
  }

  // Low-level methods to modify one endpoint of an existing S1Interval.
  // These methods should really be private because setting just one endpoint
  // can violate the invariants maintained by S1Interval.  In particular:
  //
  //  - It is not valid to call these methods on an Empty() or Full()
  //    interval, since these intervals do not have any endpoints.
  //
  //  - It is not allowed to set an endpoint to -Pi.  (When these methods are
  //    used internally, values of -Pi have already been normalized to Pi.)
  //
  // The preferred way to modify both endpoints of an interval is to use a
  // constructor, e.g. lng = S1Interval(lng_lo, lng_hi).
  void setLo(double p)
  out {
    assert(isValid());
  } do {
    _bounds[0] = p;
  }

  void setHi(double p)
  out {
    assert(isValid());
  } do {
    _bounds[1] = p;
  }

  string toString() const {
    return "[" ~ conv.to!string(lo()) ~ ", " ~ conv.to!string(hi()) ~ "]";
  }
}
