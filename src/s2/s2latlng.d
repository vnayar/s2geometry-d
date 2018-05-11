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

module s2.s2latlng;

import s2.r2point;
import s2.s2point;
import s2.s1angle;
import s2.util.math.vector;
import algorithm = std.algorithm;
import conv = std.conv;
import format = std.format;
import math = std.math;

// This class represents a point on the unit sphere as a pair
// of latitude-longitude coordinates.  Like the rest of the "geometry"
// package, the intent is to represent spherical geometry as a mathematical
// abstraction, so functions that are specifically related to the Earth's
// geometry (e.g. easting/northing conversions) should be put elsewhere.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
struct S2LatLng {
private:
  // The default constructor sets the latitude and longitude to zero.  This is
  // mainly useful when declaring arrays, STL containers, etc.
  R2Point _coords;

  // Internal constructor.
  this(in R2Point coords) {
    _coords = coords;
  }

  // This is internal to avoid ambiguity about which units are expected.
  this(double latRadians, double lngRadians) {
    _coords = R2Point([latRadians, lngRadians]);
  }


public:
  // Constructor.  The latitude and longitude are allowed to be outside
  // the is_valid() range.  However, note that most methods that accept
  // S2LatLngs expect them to be normalized (see normalized() below).
  this(in S1Angle lat, in S1Angle lng) {
    _coords = R2Point([lat.radians(), lng.radians()]);
  }

  // Convert a direction vector (not necessarily unit length) to an S2LatLng.
  this(in S2Point p)
  out {
    // The latitude and longitude are already normalized.
    assert(isValid(), "Invalid S2LatLng in constructor: " ~ this.toString());
  } body {
    _coords = R2Point([latitude(p).radians(), longitude(p).radians]);
  }

  // Returns an S2LatLng for which isValid() will return false.
  static S2LatLng invalid() {
    // These coordinates are outside the bounds allowed by is_valid().
    return S2LatLng(math.PI, 2 * math.PI);
  }

  // Convenience functions -- shorter than calling S1Angle::Radians(), etc.
  static S2LatLng fromRadians(double latRadians, double lngRadians) {
    return S2LatLng(latRadians, lngRadians);
  }

  static S2LatLng fromDegrees(double latDegrees, double lngDegrees) {
    return S2LatLng(S1Angle.fromDegrees(latDegrees), S1Angle.fromDegrees(lngDegrees));
  }

  static S2LatLng fromE5(int latE5, int lngE5) {
    return S2LatLng(S1Angle.fromE5(latE5), S1Angle.fromE5(lngE5));
  }

  static S2LatLng fromE6(int latE6, int lngE6) {
    return S2LatLng(S1Angle.fromE6(latE6), S1Angle.fromE6(lngE6));
  }

  static S2LatLng fromE7(int latE7, int lngE7) {
    return S2LatLng(S1Angle.fromE7(latE7), S1Angle.fromE7(lngE7));
  }

  // Convenience functions -- to use when args have been fixed32s in protos.
  //
  // The arguments are cast into int, so very large unsigned values
  // are treated as negative numbers.
  static S2LatLng fromUnsignedE6(uint latE6, uint lngE6) {
    return S2LatLng(S1Angle.fromUnsignedE6(latE6), S1Angle.fromUnsignedE6(lngE6));
  }

  static S2LatLng fromUnsignedE7(uint latE7, uint lngE7) {
    return S2LatLng(S1Angle.fromUnsignedE7(latE7), S1Angle.fromUnsignedE7(lngE7));
  }

  // Methods to compute the latitude and longitude of a point separately.
  static S1Angle latitude(in S2Point p) {
    // We use atan2 rather than asin because the input vector is not necessarily
    // unit length, and atan2 is much more accurate than asin near the poles.
    return S1Angle.fromRadians(math.atan2(p[2], math.sqrt(p[0]*p[0] + p[1]*p[1])));
  }

  static S1Angle longitude(in S2Point p) {
    // Note that atan2(0, 0) is defined to be zero.
    return S1Angle.fromRadians(math.atan2(p[1], p[0]));
  }

  // Accessor methods.
  S1Angle lat() const {
    return S1Angle.fromRadians(_coords[0]);
  }

  S1Angle lng() const {
    return S1Angle.fromRadians(_coords[1]);
  }

  ref R2Point coords() {
    return _coords;
  }

  // Return true if the latitude is between -90 and 90 degrees inclusive
  // and the longitude is between -180 and 180 degrees inclusive.
  bool isValid() const {
    return math.fabs(lat().radians()) <= math.PI_2 && math.fabs(lng().radians()) <= math.PI;
  }

  // Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts
  // a multiple of 360 degrees to the longitude if necessary to reduce it to
  // the range [-180, 180].
  S2LatLng normalized() const {
    // remainder(x, 2 * M_PI) reduces its argument to the range [-M_PI, M_PI]
    // inclusive, which is what we want here.
    return S2LatLng(
        algorithm.max(-math.PI_2, algorithm.min(math.PI_2, lat().radians())),
        math.remainder(lng().radians(), 2 * math.PI));
  }

  // Converts a normalized S2LatLng to the equivalent unit-length vector.
  // The maximum error in the result is 1.5 * DBL_EPSILON.  (This does not
  // include the error of converting degrees, E5, E6, or E7 to radians.)
  S2Point toS2Point() const
  in {
    assert(isValid(), "Invalid S2LatLng in S2LatLng.toPoint() : " ~ this.toString());
  } body {
    double phi = lat().radians();
    double theta = lng().radians();
    double cosphi = math.cos(phi);
    return S2Point([math.cos(theta) * cosphi, math.sin(theta) * cosphi, math.sin(phi)]);
  }

  // Returns the distance (measured along the surface of the sphere) to the
  // given S2LatLng, implemented using the Haversine formula.  This is
  // equivalent to
  //
  //   S1Angle(ToPoint(), o.ToPoint())
  //
  // except that this function is slightly faster, and is also somewhat less
  // accurate for distances approaching 180 degrees (see s1angle.h for
  // details).  Both S2LatLngs must be normalized.
  S1Angle getDistance(in S2LatLng o) const
  in {
    assert(isValid(), "Invalid S2LatLng in S2LatLng.getDistance: " ~ toString());
    assert(o.isValid(), "Invalid S2LatLng in S2LatLng.getDistance: " ~ o.toString());
  } body {
    // This implements the Haversine formula, which is numerically stable for
    // small distances but only gets about 8 digits of precision for very large
    // distances (e.g. antipodal points).  Note that 8 digits is still accurate
    // to within about 10cm for a sphere the size of the Earth.
    //
    // This could be fixed with another sin() and cos() below, but at that point
    // you might as well just convert both arguments to S2Points and compute the
    // distance that way (which gives about 15 digits of accuracy for all
    // distances).

    double lat1 = lat().radians();
    double lat2 = o.lat().radians();
    double lng1 = lng().radians();
    double lng2 = o.lng().radians();
    double dlat = math.sin(0.5 * (lat2 - lat1));
    double dlng = math.sin(0.5 * (lng2 - lng1));
    double x = dlat * dlat + dlng * dlng * math.cos(lat1) * math.cos(lat2);
    return S1Angle.fromRadians(2 * math.asin(math.sqrt(algorithm.min(1.0, x))));
  }

  // Simple arithmetic operations for manipulating latitude-longitude pairs.
  // The results are not normalized (see Normalized()).
  S2LatLng opBinary(string op)(in S2LatLng v) {
    static if (op == "+" || op == "-") {
      return mixin("S2LatLng(coords() " ~ op ~ " v.coords())");
    } else {
      static assert(false, "Unsupporeted operator: " ~ op);
    }
  }

  S2LatLng opBinary(string op)(double m) const
  if (op == "*" || op == "/") {
    return mixin("S2LatLng(_coords " ~ op ~ " m)");
  }

  S2LatLng opBinaryRight(string op)(double m) const
  if (op == "*" || op == "/") {
    return mixin("S2LatLng(m " ~ op ~ " _coords)");
  }

  bool opEquals(in S2LatLng o) const {
    return _coords == o._coords;
  }

  bool opCmp(in S2LatLng o) const {
    if (_coords == o._coords) {
      return 0;
    }
    return _coords > o._coords;
  }

  bool approxEquals(in S2LatLng o, S1Angle maxError = S1Angle.fromRadians(1e-15)) const {
    return _coords.aequal(o._coords, maxError.radians());
  }

  string toString() const {
    return "[" ~ conv.to!string(lat()) ~ ", " ~ conv.to!string(lng()) ~ "]";
  }

  // Exports the latitude and longitude in degrees, separated by a comma.
  // e.g. "94.518000,150.300000"
  string toStringInDegrees() const {
    auto pt = normalized();
    return format.format("%f,%f", pt.lat().degrees(), pt.lng().degrees());
  }
}
