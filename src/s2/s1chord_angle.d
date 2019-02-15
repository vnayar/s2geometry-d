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

module s2.s1chord_angle;

import s2.s1angle;
import s2.s2point;
import s2pointutil = s2.s2pointutil;
import algorithm = std.algorithm;
import math = std.math;

/**
 * S1ChordAngle represents the angle subtended by a chord (i.e., the straight
 * line segment connecting two points on the sphere).  Its representation
 * makes it very efficient for computing and comparing distances, but unlike
 * S1Angle it is only capable of representing angles between 0 and Pi radians.
 * Generally, S1ChordAngle should only be used in loops where many angles need
 * to be calculated and compared.  Otherwise it is simpler to use S1Angle.
 *
 * S1ChordAngle also loses some accuracy as the angle approaches Pi radians.
 * Specifically, the representation of (Pi - x) radians has an error of about
 * (1e-15 / x), with a maximum error of about 2e-8 radians (about 13cm on the
 * Earth's surface).  For comparison, for angles up to 90 degrees (10000km)
 * the worst-case representation error is about 2e-16 radians (1 nanonmeter),
 * which is about the same as S1Angle.
 *
 * This class is intended to be copied by value as desired.  It uses
 * the default copy constructor and assignment operator.
 */
struct S1ChordAngle {
private:
  static immutable double MAX_LENGTH2 = 4.0;

  // S1ChordAngles are represented by the squared chord length, which can
  // range from 0 to 4.  Infinity() uses an infinite squared length.
  this(double length2)
  in {
    assert(isValid());
  } body {
    _length2 = length2;
  }

  double _length2 = 0;

public:
  /// Construct the S1ChordAngle corresponding to the distance between the two
  /// given points.  The points must be unit length.
  this(in S2Point x, in S2Point y)
  in {
    assert(s2pointutil.isUnitLength(x));
    assert(s2pointutil.isUnitLength(y));
  } out {
    isValid();
  } body {
    // The distance may slightly exceed MAX_LENGTH2 due to roundoff errors.
    // The maximum error in the result is 2 * DBL_EPSILON * length2_.
    _length2 = algorithm.min(MAX_LENGTH2, (x - y).norm2());
  }

  this(in S1ChordAngle chordAngle) {
    _length2 = chordAngle.length2();
  }

  /**
   * Conversion from an S1Angle.  Angles outside the range [0, Pi] are handled
   * as follows: Infinity() is mapped to Infinity(), negative angles are
   * mapped to Negative(), and finite angles larger than Pi are mapped to
   * Straight().
   *
   * Note that this operation is relatively expensive and should be avoided.
   * To use S1ChordAngle effectively, you should structure your code so that
   * input arguments are converted to S1ChordAngles at the beginning of your
   * algorithm, and results are converted back to S1Angles only at the end.
   */
  this(in S1Angle angle)
  out {
    assert(isValid(), toString() ~ " is invalid.");
  } body {
    if (angle.radians() < 0) {
      this = negative();
    } else if (angle == S1Angle.infinity()) {
      this = infinity();
    } else {
      // The chord length is 2 * sin(angle / 2).
          algorithm.min(math.PI, angle.radians()));
      double length = 2 * math.sin(0.5 * algorithm.min(
              cast(double) math.PI,
              math.isFinite(angle.radians()) ? angle.radians() : double.max));
      _length2 = length * length;
    }
  }


  /// Return the zero chord angle.
  static S1ChordAngle zero() {
    return S1ChordAngle(0);
  }

  /// Return a chord angle of 90 degrees (a "right angle").
  static S1ChordAngle right() {
    return S1ChordAngle(2);
  }

  /// Return a chord angle of 180 degrees (a "straight angle").  This is the
  /// maximum finite chord angle.
  static S1ChordAngle straight() {
    return S1ChordAngle(4);
  }

  /**
   * Return a chord angle larger than any finite chord angle.  The only valid
   * operations on Infinity() are comparisons, S1Angle conversions, and
   * Successor() / Predecessor().
   */
  static S1ChordAngle infinity() {
    return S1ChordAngle(double.infinity);
  }

  /**
   * Return a chord angle smaller than Zero().  The only valid operations on
   * Negative() are comparisons, S1Angle conversions, and Successor() /
   * Predecessor().
   */
  static S1ChordAngle negative() {
    return S1ChordAngle(-1);
  }

  /// Convenience methods implemented by converting from an S1Angle.
  static S1ChordAngle fromRadians(double radians) {
    return S1ChordAngle(S1Angle.fromRadians(radians));
  }

  static S1ChordAngle fromDegrees(double degrees) {
    return S1ChordAngle(S1Angle.fromDegrees(degrees));
  }

  S1ChordAngle fromE5(int e5) {
    return S1ChordAngle(S1Angle.fromE5(e5));
  }

  S1ChordAngle fromE6(int e6) {
    return S1ChordAngle(S1Angle.fromE6(e6));
  }

  S1ChordAngle fromE7(int e7) {
    return S1ChordAngle(S1Angle.fromE7(e7));
  }

  /**
   * Construct an S1ChordAngle that is an upper bound on the given S1Angle,
   * i.e. such that FastUpperBoundFrom(x).ToAngle() >= x.  Unlike the S1Angle
   * constructor above, this method is very fast, and the bound is accurate to
   * within 1% for distances up to about 3100km on the Earth's surface.
   */
  static S1ChordAngle fastUpperBoundFrom(in S1Angle angle) {
    // This method uses the distance along the surface of the sphere as an upper
    // bound on the distance through the sphere's interior.
    return S1ChordAngle.fromLength2(angle.radians() * angle.radians());
  }

  /**
   * Construct an S1ChordAngle from the squared chord length.  Note that the
   * argument is automatically clamped to a maximum of 4.0 to handle possible
   * roundoff errors.  The argument must be non-negative.
   */
  static S1ChordAngle fromLength2(double length2) {
    return S1ChordAngle(algorithm.min(4.0, length2));
  }

  /**
   * Converts to an S1Angle.
   *
   * Infinity() is converted to S1Angle::Infinity(), and Negative() is
   * converted to an unspecified negative S1Angle.
   *
   * Note that the conversion uses trigonometric functions and therefore
   * should be avoided in inner loops.
   */
  S1Angle toS1Angle() const {
    if (isNegative()) {
      return S1Angle.fromRadians(-1);
    }
    if (isInfinity()) {
      return S1Angle.infinity();
    }
    return S1Angle.fromRadians(2 * math.asin(0.5 * math.sqrt(_length2)));
  }

  /**
   * Convenience methods implemented by calling ToAngle() first.  Note that
   * because of the S1Angle conversion these methods are relatively expensive
   * (despite their lowercase names), so the results should be cached if they
   * are needed inside loops.
   */
  double radians() const {
    return toS1Angle().radians();
  }

  double degrees() const {
    return toS1Angle().degrees();
  }

  int e5() const {
    return toS1Angle().e5();
  }

  int e6() const {
    return toS1Angle().e6();
  }

  int e7() const {
    return toS1Angle().e7();
  }

  // All operators and functions are declared here so that we can put them all
  // in one place.  (The compound assignment operators must be put here.)

  // Comparison operators.
  bool opEquals(in S1ChordAngle x) const {
    return length2() == x.length2();
  }

  int opCmp(in S1ChordAngle x) const {
    if (length2() == x.length2()) {
      return 0;
    }
    return length2() < x.length2() ? -1 : 1;
  }

  /// Comparison predicates.
  bool isZero() const {
    return _length2 == 0;
  }

  bool isNegative() const {
    // TODO(ericv): Consider stricter check here -- only allow Negative().
    return _length2 < 0;
  }

  bool isInfinity() const {
    return _length2 == double.infinity;
  }

  /// Negative or infinity.
  bool isSpecial() const {
    return isNegative() || isInfinity();
  }

  /**
   * Only addition and subtraction of S1ChordAngles is supported.  These
   * methods add or subtract the corresponding S1Angles, and clamp the result
   * to the range [0, Pi].  Both arguments must be non-negative and
   * non-infinite.
   *
   * REQUIRES: !a.is_special() && !b.is_special()
   */
  S1ChordAngle opBinary(string op)(S1ChordAngle b) const
  if (op == "+")
  in {
    assert(!isSpecial());
    assert(!b.isSpecial());
  } body {
    // Note that this method is much more efficient than converting the chord
    // angles to S1Angles and adding those.  It requires only one square root
    // plus a few additions and multiplications.

    // Optimization for the common case where "b" is an error tolerance
    // parameter that happens to be set to zero.
    double a2 = length2();
    double b2 = b.length2();
    if (b2 == 0) {
      return S1ChordAngle(this);
    }

    // Clamp the angle sum to at most 180 degrees.
    if (a2 + b2 >= MAX_LENGTH2) {
      return S1ChordAngle.straight();
    }

    // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
    // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
    // Then the formula below can be derived from c = 2 * sin(A+B) and the
    // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
    //                 cos(X) = sqrt(1 - sin^2(X)) .
    double x = a2 * (1 - 0.25 * b2);  // is_valid() => non-negative
    double y = b2 * (1 - 0.25 * a2);  // is_valid() => non-negative
    return S1ChordAngle(algorithm.min(MAX_LENGTH2, x + y + 2 * math.sqrt(x * y)));
  }

  S1ChordAngle opBinary(string op)(S1ChordAngle b) const
  if (op == "-")
  in {
    assert(!isSpecial());
    assert(!b.isSpecial());
  } body {
    // See comments in opBinary!"+"().
    double a2 = length2(), b2 = b.length2();
    if (b2 == 0) {
      return S1ChordAngle(this);
    }
    if (a2 <= b2) {
      return S1ChordAngle.zero();
    }
    double x = a2 * (1 - 0.25 * b2);
    double y = b2 * (1 - 0.25 * a2);
    return S1ChordAngle(algorithm.max(0.0, x + y - 2 * math.sqrt(x * y)));
  }

  S1ChordAngle opOpAssign(string op)(S1ChordAngle a)
  if (op == "+" || op == "-") {
    return this = mixin("this " ~ op ~ "a");
  }

  /// Trigonmetric functions.  It is more accurate and efficient to call these
  /// rather than first converting to an S1Angle.
  double sin() const {
    return math.sqrt(sin2());
  }

  double cos() const
  in {
    assert(!isSpecial());
  } body {
    // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
    return 1 - 0.5 * length2();
  }

  double tan() const {
    return sin() / cos();
  }

  /// Returns sin(a)**2, but computed more efficiently.
  double sin2() const
  in {
    assert(!isSpecial());
  } body {
    // Let "a" be the (non-squared) chord length, and let A be the corresponding
    // half-angle (a = 2*sin(A)).  The formula below can be derived from:
    //   sin(2*A) = 2 * sin(A) * cos(A)
    //   cos^2(A) = 1 - sin^2(A)
    // This is much faster than converting to an angle and computing its sine.
    return length2() * (1 - 0.25 * length2());
  }


  /// The squared length of the chord.  (Most clients will not need this.)
  @property
  double length2() const {
    return _length2;
  }

  /**
   * Returns the smallest representable S1ChordAngle larger than this object.
   * This can be used to convert a "<" comparison to a "<=" comparison.  For
   * example:
   *
   *   S2ClosestEdgeQuery query(...);
   *   S1ChordAngle limit = ...;
   *   if (query.IsDistanceLess(target, limit.Successor())) {
   *     // Distance to "target" is less than or equal to "limit".
   *   }
   *
   * Note the following special cases:
   *   Negative().Successor() == Zero()
   *   Straight().Successor() == Infinity()
   *   Infinity().Successor() == Infinity()
   */
  S1ChordAngle successor() const {
    if (_length2 >= MAX_LENGTH2) {
      return infinity();
    }
    if (_length2 < 0.0) {
      return zero();
    }
    return S1ChordAngle(math.nextafter(_length2, 10.0));
  }

  /**
   * Like Successor(), but returns the largest representable S1ChordAngle less
   * than this object.
   *
   * Note the following special cases:
   *   Infinity().Predecessor() == Straight()
   *   Zero().Predecessor() == Negative()
   *   Negative().Predecessor() == Negative()
   */
  S1ChordAngle predecessor() const {
    if (_length2 <= 0.0) {
      return negative();
    }
    if (_length2 > MAX_LENGTH2) {
      return straight();
    }
    return S1ChordAngle(math.nextafter(_length2, -10.0));
  }


  /**
   * Returns a new S1ChordAngle that has been adjusted by the given error
   * bound (which can be positive or negative).  "error" should be the value
   * returned by one of the error bound methods below.  For example:
   *    S1ChordAngle a(x, y);
   *    S1ChordAngle a1 = a.PlusError(a.GetS2PointConstructorMaxError());
   */
  S1ChordAngle plusError(double error) const {
    // If angle is Negative() or Infinity(), don't change it.
    // Otherwise clamp it to the valid range.
    if (isSpecial()) {
      return this;
    }
    return S1ChordAngle(algorithm.max(0.0, algorithm.min(MAX_LENGTH2, _length2 + error)));
  }

  /**
   * Return the maximum error in length2() for the S1ChordAngle(x, y)
   * constructor, assuming that "x" and "y" are normalized to within the
   * bounds guaranteed by S2Point::Normalize().  (The error is defined with
   * respect to the true distance after the points are projected to lie
   * exactly on the sphere.)
   */
  double getS2PointConstructorMaxError() const {
    // There is a relative error of 2.5 * DBL_EPSILON when computing the squared
    // distance, plus a relative error of 2 * DBL_EPSILON and an absolute error
    // of (16 * DBL_EPSILON**2) because the lengths of the input points may
    // differ from 1 by up to (2 * DBL_EPSILON) each.  (This is the maximum
    // length error in S2Point::Normalize.)
    return 4.5 * double.epsilon * _length2 + 16 * double.epsilon * double.epsilon;
  }

  /// Return the maximum error in length2() for the S1Angle constructor.
  double getS1AngleConstructorMaxError() const {
    // Assuming that an accurate math library is being used, the sin() call and
    // the multiply each have a relative error of 0.5 * DBL_EPSILON.
    return double.epsilon * _length2;
  }

  /// Return true if the internal representation is valid.  Negative() and
  /// Infinity() are both considered valid.
  bool isValid() const {
    return (_length2 >= 0 && _length2 <= MAX_LENGTH2) || isSpecial();
  }

  string toString() const {
    return toS1Angle().toString();
  }
}
