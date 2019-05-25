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

module s2.s1angle;

import s2.s2point;
import s2.util.math.mathutil;
import s2.util.math.vector;
import s2.s2latlng;
import math = std.math;
import format = std.format;

/**
 * This class represents a one-dimensional angle (as opposed to a
 * two-dimensional solid angle).  It has methods for converting angles to
 * or from radians, degrees, and the E5/E6/E7 representations (i.e. degrees
 * multiplied by 1e5/1e6/1e7 and rounded to the nearest integer).
 *
 * The internal representation is a double-precision value in radians, so
 * conversion to and from radians is exact.  Conversions between E5, E6, E7,
 * and Degrees are not always exact; for example, Degrees(3.1) is different
 * from E6(3100000) or E7(310000000).  However, the following properties are
 * guaranteed for any integer "n", provided that "n" is in the input range of
 * both functions:
 *
 *     Degrees(n) == E6(1000000 * n)
 *     Degrees(n) == E7(10000000 * n)
 *          E6(n) == E7(10 * n)
 *
 * The corresponding properties are *not* true for E5, so if you use E5 then
 * don't test for exact equality when comparing to other formats such as
 * Degrees or E7.
 *
 * The following conversions between degrees and radians are exact:
 *
 *          Degrees(180) == Radians(M_PI)
 *       Degrees(45 * k) == Radians(k * M_PI / 4)  for k == 0..8
 *
 * These identities also hold when the arguments are scaled up or down by any
 * power of 2.  Some similar identities are also true, for example,
 * Degrees(60) == Radians(M_PI / 3), but be aware that this type of identity
 * does not hold in general.  For example, Degrees(3) != Radians(M_PI / 60).
 *
 * Similarly, the conversion to radians means that Angle::Degrees(x).degrees()
 * does not always equal "x".  For example,
 *
 *         S1Angle::Degrees(45 * k).degrees() == 45 * k      for k == 0..8
 *   but       S1Angle::Degrees(60).degrees() != 60.
 *
 * This means that when testing for equality, you should allow for numerical
 * errors (EXPECT_DOUBLE_EQ) or convert to discrete E5/E6/E7 values first.
 *
 * CAVEAT: All of the above properties depend on "double" being the usual
 * 64-bit IEEE 754 type (which is true on almost all modern platforms).
 *
 * This class is intended to be copied by value as desired.  It uses
 * the default copy constructor and assignment operator.
 */
struct S1Angle {
private:
  /// The default constructor yields a zero angle.  This is useful for STL
  /// containers and class methods with output arguments.
  double _radians = 0;

  this(double radians) {
    _radians = radians;
  }

public:
  /// These methods construct S1Angle objects from their measure in radians or degrees.
  static S1Angle fromRadians(double radians) {
    return S1Angle(radians);
  }

  static S1Angle fromDegrees(double degrees) {
    return S1Angle((M_PI / 180) * degrees);
  }

  static S1Angle fromE5(int e5) {
    return S1Angle.fromDegrees(1e-5 * e5);
  }

  static S1Angle fromE6(int e6) {
    return S1Angle.fromDegrees(1e-6 * e6);
  }

  static S1Angle fromE7(int e7) {
    return S1Angle.fromDegrees(1e-7 * e7);
  }

  // Convenience functions -- to use when args have been fixed32s in protos.
  //
  // The arguments are cast into int, so very large unsigned values
  // are treated as negative numbers.
  static S1Angle fromUnsignedE6(uint e6) {
    return S1Angle.fromE6(cast(int) e6);
  }

  static S1Angle fromUnsignedE7(uint e7) {
    return S1Angle.fromE7(cast(int) e7);
  }

  // Return an angle larger than any finite angle.
  static S1Angle infinity() {
    return S1Angle(double.infinity);
  }

  // A explicit shorthand for the default constructor.
  static S1Angle zero() {
    return S1Angle(0);
  }

  // Return the angle between two points, which is also equal to the distance
  // between these points on the unit sphere.  The points do not need to be
  // normalized.  This function has a maximum error of 3.25 * DBL_EPSILON (or
  // 2.5 * DBL_EPSILON for angles up to 1 radian).
  this(in S2Point x, in S2Point y) {
    _radians = x.angle(y);
  }

  // Like the constructor above, but return the angle (i.e., distance) between
  // two S2LatLng points.  This function has about 15 digits of accuracy for
  // small distances but only about 8 digits of accuracy as the distance
  // approaches 180 degrees (i.e., nearly-antipodal points).
  this(in S2LatLng x, in S2LatLng y) {
    _radians = x.getDistance(y).radians();
  }

  double radians() const {
    return _radians;
  }

  double degrees() const {
    return (180 / M_PI) * _radians;
  }

  // Note that the E5, E6, and E7 conversion involve two multiplications rather
  // than one.  This is mainly for backwards compatibility (changing this would
  // break many tests), but it does have the nice side effect that conversions
  // between Degrees, E6, and E7 are exact when the arguments are integers.

  int e5() const {
    return cast(int) math.lround(1e5 * degrees());
  }

  int e6() const {
    return cast(int) math.lround(1e6 * degrees());
  }

  int e7() const {
    return cast(int) math.lround(1e7 * degrees());
  }

  // Return the absolute value of an angle.
  S1Angle abs() const {
    return S1Angle(math.fabs(_radians));
  }


  // Comparison operators.
  bool opEquals(in S1Angle y) const {
    return radians() == y.radians();
  }

  int opCmp(in S1Angle y) const {
    if (radians() == y.radians()) {
      return 0;
    }
    return radians() > y.radians() ? 1 : -1;
  }

  // Implement negation.
  S1Angle opUnary(string op)() const
  if (op == "-") {
    return S1Angle(-radians());
  }

  // Simple arithmetic operators for manipulating S1Angles.
  S1Angle opOpAssign(string op)(in S1Angle v) {
    mixin("_radians " ~ op ~ "= v.radians();");
    return this;
  }

  S1Angle opOpAssign(string op)(in double v) {
    mixin("_radians " ~ op ~ "= v;");
    return this;
  }

  S1Angle opBinary(string op)(in S1Angle v) const {
    return mixin("S1Angle(radians() " ~ op ~ " v.radians())");
  }

  S1Angle opBinary(string op)(in double v) const {
    return mixin("S1Angle(radians() " ~ op ~ " v)");
  }

  S1Angle opBinaryRight(string op)(double v) const {
    return mixin("S1Angle(v " ~ op ~ " radians())");
  }

  // Return the angle normalized to the range (-180, 180] degrees.
  S1Angle normalized() const {
    auto a = S1Angle(_radians);
    a.normalize();
    return a;
  }

  // Normalize this angle to the range (-180, 180] degrees.
  void normalize() {
    _radians = math.remainder(_radians, 2.0 * M_PI);
    if (_radians <= -M_PI) {
      _radians = M_PI;
    }
  }

  string toString() const {
    return format.format("%.7f", degrees());
  }
}

// Trigonmetric functions (not necessary but slightly more convenient).
double sin(in S1Angle a) {
  return math.sin(a.radians());
}

double cos(in S1Angle a) {
  return math.cos(a.radians());
}

double tan(in S1Angle a) {
  return math.tan(a.radians());
}
