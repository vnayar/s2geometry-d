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

module s2.s2latlng_rect_test;

// Most of the S2LatLngRect methods have trivial implementations that
// use the R1Interval and S1Interval classes, so most of the testing
// is done in those unit tests.

// #include "s2/util/coding/coder.h"
import algorithm = std.algorithm;
import fluent.asserts;
import math = std.math;
import s2.logger;
import s2.r1interval;
import s2.s1angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2edge_distances;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2pointutil;
import s2.s2testing;
import s2.util.math.vector;
import s2.s2predicates : sign;

enum double DOUBLE_ERR = 0.0001;

private S2LatLngRect rectFromDegrees(double lat_lo, double lng_lo, double lat_hi, double lng_hi) {
  // Convenience method to construct a rectangle.  This method is
  // intentionally *not* in the S2LatLngRect interface because the
  // argument order is ambiguous, but hopefully it's not too confusing
  // within the context of this unit test.

  return new S2LatLngRect(
      S2LatLng.fromDegrees(lat_lo, lng_lo).normalized(),
      S2LatLng.fromDegrees(lat_hi, lng_hi).normalized());
}

@("S2LatLngRect.EmptyAndFull")
unittest {
  // Test basic properties of empty and full rectangles.
  S2LatLngRect empty = S2LatLngRect.empty();
  S2LatLngRect full = S2LatLngRect.full();
  Assert.equal(empty.isValid(), true);
  Assert.equal(empty.isEmpty(), true);
  Assert.equal(empty.isPoint(), false);
  Assert.equal(full.isValid(), true);
  Assert.equal(full.isFull(), true);
  Assert.equal(full.isPoint(), false);
  // Check that the default S2LatLngRect is identical to Empty().
  S2LatLngRect default_empty = new S2LatLngRect();
  Assert.equal(default_empty.isValid(), true);
  Assert.equal(default_empty.isEmpty(), true);
  Assert.equal(empty.lat().bounds(), default_empty.lat().bounds());
  Assert.equal(empty.lng().bounds(), default_empty.lng().bounds());
}

@("S2LatLngRect.Accessors")
unittest {
  // Check various accessor methods.
  S2LatLngRect d1 = rectFromDegrees(-90, 0, -45, 180);
  Assert.approximately(d1.latLo().degrees(), -90, DOUBLE_ERR);
  Assert.approximately(d1.latHi().degrees(), -45, DOUBLE_ERR);
  Assert.approximately(d1.lngLo().degrees(), 0, DOUBLE_ERR);
  Assert.approximately(d1.lngHi().degrees(), 180, DOUBLE_ERR);
  Assert.equal(d1.lat(), R1Interval(-math.PI_2, -math.PI_4));
  Assert.equal(d1.lng(), S1Interval(0, math.PI));
}

@("S2LatLngRect.ApproxEquals")
unittest {
  // S1Interval and R1Interval have additional testing.

  Assert.equal(S2LatLngRect.empty().approxEquals(rectFromDegrees(1, 5, 1, 5)), true);
  Assert.equal(rectFromDegrees(1, 5, 1, 5).approxEquals(S2LatLngRect.empty()), true);
  Assert.equal(rectFromDegrees(1, 5, 1, 5).approxEquals(rectFromDegrees(2, 7, 2, 7)), false);

  // Test the max_error (double) parameter.
  Assert.equal(
      rectFromDegrees(10, 10, 20, 20)
          .approxEquals(rectFromDegrees(11, 11, 19, 19), S1Angle.fromDegrees(1.001)),
      true);
  Assert.equal(
      rectFromDegrees(10, 10, 20, 20)
          .approxEquals(rectFromDegrees(11, 11, 19, 19), S1Angle.fromDegrees(0.999)),
      false);

  // Test the max_error (S2LatLng) parameter.
  Assert.equal(
      rectFromDegrees(0, 10, 20, 30)
          .approxEquals(rectFromDegrees(-1, 8, 21, 32), S2LatLng.fromDegrees(1.001, 2.001)),
      true);
  Assert.equal(
      rectFromDegrees(0, 10, 20, 30)
          .approxEquals(rectFromDegrees(-1, 8, 21, 32), S2LatLng.fromDegrees(0.999, 1.999)),
      false);
}

@("S2LatLngRect.FromCenterSize")
unittest {
  Assert.equal(
      S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(80, 170), S2LatLng.fromDegrees(40, 60))
          .approxEquals(rectFromDegrees(60, 140, 90, -160)),
      true);
  Assert.equal(
      S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(10, 40), S2LatLng.fromDegrees(210, 400))
          .isFull(),
      true);
  Assert.equal(
      S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(-90, 180), S2LatLng.fromDegrees(20, 50))
          .approxEquals(rectFromDegrees(-90, 155, -80, -155)),
      true);
}

@("S2LatLngRect.FromPoint")
unittest {
  S2LatLng p = S2LatLng.fromDegrees(23, 47);
  Assert.equal(S2LatLngRect.fromPoint(p), new S2LatLngRect(p, p));
  Assert.equal(S2LatLngRect.fromPoint(p).isPoint(), true);
}

@("S2LatLngRect.FromPointPair")
unittest {
  Assert.equal(
      S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(-35, -140), S2LatLng.fromDegrees(15, 155)),
      rectFromDegrees(-35, 155, 15, -140));
  Assert.equal(
      S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(25, -70), S2LatLng.fromDegrees(-90, 80)),
      rectFromDegrees(-90, -70, 25, 80));
}

@("S2LatLngRect.GetCenterSize")
unittest {
  S2LatLngRect r1 = new S2LatLngRect(R1Interval(0, math.PI_2), S1Interval(-math.PI, 0));
  Assert.equal(r1.getCenter(), S2LatLng.fromRadians(math.PI_4, -math.PI_2));
  Assert.equal(r1.getSize(), S2LatLng.fromRadians(math.PI_2, math.PI));
  Assert.lessThan(S2LatLngRect.empty().getSize().lat().radians(), 0);
  Assert.lessThan(S2LatLngRect.empty().getSize().lng().radians(), 0);
}

@("S2LatLngRect.GetVertex")
unittest {
  S2LatLngRect r1 = new S2LatLngRect(R1Interval(0, math.PI_2), S1Interval(-math.PI, 0));
  Assert.equal(r1.getVertex(0), S2LatLng.fromRadians(0, math.PI));
  Assert.equal(r1.getVertex(1), S2LatLng.fromRadians(0, 0));
  Assert.equal(r1.getVertex(2), S2LatLng.fromRadians(math.PI_2, 0));
  Assert.equal(r1.getVertex(3), S2LatLng.fromRadians(math.PI_2, math.PI));

  // Make sure that GetVertex() returns vertices in CCW order.
  for (int i = 0; i < 4; ++i) {
    double lat = math.PI_4 * (i - 2);
    double lng = math.PI_2 * (i - 2) + 0.2;
    S2LatLngRect r = new S2LatLngRect(
        R1Interval(lat, lat + math.PI_4),
        S1Interval(
            math.remainder(lng, 2 * math.PI), math.remainder(lng + math.PI_2, 2 * math.PI)));
    for (int k = 0; k < 4; ++k) {
      Assert.equal(
          sign(
              r.getVertex(k - 1).toS2Point(),
              r.getVertex(k).toS2Point(),
              r.getVertex(k + 1).toS2Point()),
          1);
    }
  }
}

@("S2LatLngRect.Contains")
unittest {
  // Contains(S2LatLng), InteriorContains(S2LatLng), Contains()
  S2LatLng eq_m180 = S2LatLng.fromRadians(0, -math.PI);
  S2LatLng north_pole = S2LatLng.fromRadians(math.PI_2, 0);
  S2LatLngRect r1 = new S2LatLngRect(eq_m180, north_pole);

  Assert.equal(r1.contains(S2LatLng.fromDegrees(30, -45)), true);
  Assert.equal(r1.interiorContains(S2LatLng.fromDegrees(30, -45)), true);
  Assert.equal(r1.contains(S2LatLng.fromDegrees(30, 45)), false);
  Assert.equal(r1.interiorContains(S2LatLng.fromDegrees(30, 45)), false);
  Assert.equal(r1.contains(eq_m180), true);
  Assert.equal(r1.interiorContains(eq_m180), false);
  Assert.equal(r1.contains(north_pole), true);
  Assert.equal(r1.interiorContains(north_pole), false);
  Assert.equal(r1.contains(S2Point(0.5, -0.3, 0.1)), true);
  Assert.equal(r1.contains(S2Point(0.5, 0.2, 0.1)), false);
}

static void testIntervalOps(
    in S2LatLngRect x, in S2LatLngRect y,
    string expected_relation,
    in S2LatLngRect expected_union,
    in S2LatLngRect expected_intersection) {
  // Test all of the interval operations on the given pair of intervals.
  // "expected_relation" is a sequence of "T" and "F" characters corresponding
  // to the expected results of Contains(), InteriorContains(), Intersects(),
  // and InteriorIntersects() respectively.

  Assert.equal(x.contains(y), expected_relation[0] == 'T');
  Assert.equal(x.interiorContains(y), expected_relation[1] == 'T');
  Assert.equal(x.intersects(y), expected_relation[2] == 'T');
  Assert.equal(x.interiorIntersects(y), expected_relation[3] == 'T');

  Assert.equal(x.contains(y), x.unite(y) == x);
  Assert.equal(x.intersects(y), !x.intersection(y).isEmpty());

  Assert.equal(x.unite(y), expected_union);
  Assert.equal(x.intersection(y), expected_intersection);

  if (y.getSize() == S2LatLng.fromRadians(0, 0)) {
    S2LatLngRect r = new S2LatLngRect(x);
    r.addPoint(y.lo());
    Assert.equal(r, expected_union);
  }
}

@("S2LatLngRect.IntervalOps")
unittest {
  // Contains(S2LatLngRect), InteriorContains(S2LatLngRect),
  // Intersects(), InteriorIntersects(), Union(), Intersection().
  //
  // Much more testing of these methods is done in s1interval_test
  // and r1interval_test.

  // Rectangle "r1" covers one-quarter of the sphere.
  S2LatLngRect r1 = rectFromDegrees(0, -180, 90, 0);

  // Test operations where one rectangle consists of a single point.
  S2LatLngRect r1_mid = rectFromDegrees(45, -90, 45, -90);
  testIntervalOps(r1, r1_mid, "TTTT", r1, r1_mid);

  S2LatLngRect req_m180 = rectFromDegrees(0, -180, 0, -180);
  testIntervalOps(r1, req_m180, "TFTF", r1, req_m180);

  S2LatLngRect rnorth_pole = rectFromDegrees(90, 0, 90, 0);
  testIntervalOps(r1, rnorth_pole, "TFTF", r1, rnorth_pole);

  testIntervalOps(r1, rectFromDegrees(-10, -1, 1, 20), "FFTT",
                  rectFromDegrees(-10, 180, 90, 20),
                  rectFromDegrees(0, -1, 1, 0));
  testIntervalOps(r1, rectFromDegrees(-10, -1, 0, 20), "FFTF",
                  rectFromDegrees(-10, 180, 90, 20),
                  rectFromDegrees(0, -1, 0, 0));
  testIntervalOps(r1, rectFromDegrees(-10, 0, 1, 20), "FFTF",
                  rectFromDegrees(-10, 180, 90, 20),
                  rectFromDegrees(0, 0, 1, 0));

  testIntervalOps(rectFromDegrees(-15, -160, -15, -150),
                  rectFromDegrees(20, 145, 25, 155), "FFFF",
                  rectFromDegrees(-15, 145, 25, -150),
                  S2LatLngRect.empty());
  testIntervalOps(rectFromDegrees(70, -10, 90, -140),
                  rectFromDegrees(60, 175, 80, 5), "FFTT",
                  rectFromDegrees(60, -180, 90, 180),
                  rectFromDegrees(70, 175, 80, 5));

  // Check that the intersection of two rectangles that overlap in latitude
  // but not longitude is valid, and vice versa.
  testIntervalOps(rectFromDegrees(12, 30, 60, 60),
                  rectFromDegrees(0, 0, 30, 18), "FFFF",
                  rectFromDegrees(0, 0, 60, 60), S2LatLngRect.empty());
  testIntervalOps(rectFromDegrees(0, 0, 18, 42),
                  rectFromDegrees(30, 12, 42, 60), "FFFF",
                  rectFromDegrees(0, 0, 42, 60), S2LatLngRect.empty());
}

@("BoundaryIntersects.EmptyRectangle")
unittest {
  S2LatLngRect rect = S2LatLngRect.empty();
  S2Point lo = S2Point(rect.lo().toS2Point()), hi = S2Point(rect.hi().toS2Point());
  Assert.equal(rect.boundaryIntersects(lo, lo), false);
  Assert.equal(rect.boundaryIntersects(lo, hi), false);
}

S2Point makeS2Point(double latDeg, double lngDeg) {
  return S2LatLng.fromDegrees(latDeg, lngDeg).toS2Point();
}

@("BoundaryIntersects.FullRectangle")
unittest {
  S2LatLngRect rect = S2LatLngRect.full();
  S2Point lo = S2Point(rect.lo().toS2Point()), hi = S2Point(rect.hi().toS2Point());
  Assert.equal(rect.boundaryIntersects(lo, lo), false);
  Assert.equal(rect.boundaryIntersects(lo, hi), false);
}

@("BoundaryIntersects.SphericalLune")
unittest {
  // This rectangle only has two non-degenerate sides.
  S2LatLngRect rect = rectFromDegrees(-90, 100, 90, 120);
  Assert.equal(rect.boundaryIntersects(makeS2Point(60, 60), makeS2Point(90, 60)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-60, 110), makeS2Point(60, 110)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-60, 95), makeS2Point(60, 110)), true);
  Assert.equal(rect.boundaryIntersects(makeS2Point(60, 115), makeS2Point(80, 125)), true);
}

@("BoundaryIntersects.NorthHemisphere")
unittest {
  // This rectangle only has only one non-degenerate side.
  S2LatLngRect rect = rectFromDegrees(0, -180, 90, 180);
  Assert.equal(rect.boundaryIntersects(makeS2Point(60, -180), makeS2Point(90, -180)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(60, -170), makeS2Point(60, 170)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-10, -180), makeS2Point(10, -180)), true);
}

@("BoundaryIntersects.SouthHemisphere")
unittest {
  // This rectangle only has only one non-degenerate side.
  S2LatLngRect rect = rectFromDegrees(-90, -180, 0, 180);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-90, -180), makeS2Point(-60, -180)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-60, -170), makeS2Point(-60, 170)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(-10, -180), makeS2Point(10, -180)), true);
}

@("BoundaryIntersects.RectCrossingAntiMeridian")
unittest {
  S2LatLngRect rect = rectFromDegrees(20, 170, 40, -170);
  Assert.equal(rect.contains(makeS2Point(30, 180)), true);

  // Check that crossings of all four sides are detected.
  Assert.equal(rect.boundaryIntersects(makeS2Point(25, 160), makeS2Point(25, 180)), true);
  Assert.equal(rect.boundaryIntersects(makeS2Point(25, -160), makeS2Point(25, -180)), true);
  Assert.equal(rect.boundaryIntersects(makeS2Point(15, 175), makeS2Point(30, 175)), true);
  Assert.equal(rect.boundaryIntersects(makeS2Point(45, 175), makeS2Point(30, 175)), true);

  // Check that the edges on the opposite side of the sphere but at the same
  // latitude do not intersect the rectangle boundary.
  Assert.equal(rect.boundaryIntersects(makeS2Point(25, -20), makeS2Point(25, 0)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(25, 20), makeS2Point(25, 0)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(15, -5), makeS2Point(30, -5)), false);
  Assert.equal(rect.boundaryIntersects(makeS2Point(45, -5), makeS2Point(30, -5)), false);
}

@("S2LatLngRect.AddPoint")
unittest {
  S2LatLngRect p = S2LatLngRect.empty();
  p.addPoint(S2LatLng.fromDegrees(0, 0));
  Assert.equal(true, p.isPoint());
  p.addPoint(S2LatLng.fromRadians(0, -math.PI_2));
  Assert.equal(false, p.isPoint());
  p.addPoint(S2LatLng.fromRadians(math.PI_4, -math.PI));
  p.addPoint(S2Point(0, 0, 1));
  Assert.equal(p, rectFromDegrees(0, -180, 90, 0));
}

@("S2LatLngRect.Expanded")
unittest {
  Assert.equal(
      rectFromDegrees(70, 150, 80, 170).expanded(S2LatLng.fromDegrees(20, 30))
          .approxEquals(rectFromDegrees(50, 120, 90, -160)),
      true);
  Assert.equal(
      S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(20, 30)).isEmpty(), true);
  Assert.equal(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(500, 500)).isFull(), true);
  Assert.equal(
      rectFromDegrees(-90, 170, 10, 20).expanded(S2LatLng.fromDegrees(30, 80))
          .approxEquals(rectFromDegrees(-90, -180, 40, 180)),
      true);

  // Negative margins.
  Assert.equal(
      rectFromDegrees(10, -50, 60, 70).expanded(S2LatLng.fromDegrees(-10, -10))
          .approxEquals(rectFromDegrees(20, -40, 50, 60)),
      true);
  Assert.equal(
      rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(-10, -10))
          .approxEquals(rectFromDegrees(-10, -180, 10, 180)),
      true);
  Assert.equal(
      rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(-30, -30)).isEmpty(),
      true);
  Assert.equal(
      rectFromDegrees(-90, 10, 90, 11).expanded(S2LatLng.fromDegrees(-10, -10)).isEmpty(), true);
  Assert.equal(
      rectFromDegrees(-90, 10, 90, 100).expanded(S2LatLng.fromDegrees(-10, -10))
          .approxEquals(rectFromDegrees(-80, 20, 80, 90)),
      true);
  Assert.equal(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(-50, -500)).isEmpty(), true);
  Assert.equal(
      S2LatLngRect.full().expanded(S2LatLng.fromDegrees(-50, -50))
          .approxEquals(rectFromDegrees(-40, -180, 40, 180)),
      true);

  // Mixed margins.
  Assert.equal(
      rectFromDegrees(10, -50, 60, 70).expanded(S2LatLng.fromDegrees(-10, 30))
          .approxEquals(rectFromDegrees(20, -80, 50, 100)),
      true);
  Assert.equal(
      rectFromDegrees(-20, -180, 20, 180).expanded(S2LatLng.fromDegrees(10, -500))
          .approxEquals(rectFromDegrees(-30, -180, 30, 180)),
      true);
  Assert.equal(
      rectFromDegrees(-90, -180, 80, 180).expanded(S2LatLng.fromDegrees(-30, 500))
          .approxEquals(rectFromDegrees(-60, -180, 50, 180)),
      true);
  Assert.equal(
      rectFromDegrees(-80, -100, 80, 150).expanded(S2LatLng.fromDegrees(30, -50))
          .approxEquals(rectFromDegrees(-90, -50, 90, 100)),
      true);
  Assert.equal(
      rectFromDegrees(0, -180, 50, 180).expanded(S2LatLng.fromDegrees(-30, 500)).isEmpty(),
      true);
  Assert.equal(
      rectFromDegrees(-80, 10, 70, 20).expanded(S2LatLng.fromDegrees(30, -200)).isEmpty(), true);
  Assert.equal(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(100, -100)).isEmpty(), true);
  Assert.equal(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(100, -100)).isFull(), true);
}

@("S2LatLngRect.PolarClosure")
unittest {
  Assert.equal(rectFromDegrees(-89, 0, 89, 1), rectFromDegrees(-89, 0, 89, 1).polarClosure());
  Assert.equal(
      rectFromDegrees(-90, -30, -45, 100).polarClosure(),
      rectFromDegrees(-90, -180, -45, 180));
  Assert.equal(
      rectFromDegrees(89, 145, 90, 146).polarClosure(),
      rectFromDegrees(89, -180, 90, 180));
  Assert.equal(rectFromDegrees(-90, -145, 90, -144).polarClosure(), S2LatLngRect.full());
}

@("ExpandedByDistance.PositiveDistance")
unittest {
  Assert.equal(
      rectFromDegrees(0, 170, 0, -170).expandedByDistance(S1Angle.fromDegrees(15))
          .approxEquals(rectFromDegrees(-15, 155, 15, -155)),
      true);
  Assert.equal(
      rectFromDegrees(60, 150, 80, 10).expandedByDistance(S1Angle.fromDegrees(15))
          .approxEquals(rectFromDegrees(45, -180, 90, 180)),
      true);
}

@("ExpandedByDistance.NegativeDistanceNorthEast")
unittest {
  S2LatLngRect in_rect = rectFromDegrees(0.0, 0.0, 30.0, 90.0);
  S1Angle distance = S1Angle.fromDegrees(5.0);

  S2LatLngRect out_rect = in_rect.expandedByDistance(distance).expandedByDistance(-distance);

  Assert.equal(out_rect.approxEquals(in_rect), true);
}

@("ExpandedByDistance.NegativeDistanceSouthWest")
unittest {
  S2LatLngRect in_rect = rectFromDegrees(-30.0, -90.0, 0.0, 0.0);
  S1Angle distance = S1Angle.fromDegrees(5.0);

  S2LatLngRect out_rect =
      in_rect.expandedByDistance(distance).expandedByDistance(-distance);

  Assert.equal(out_rect.approxEquals(in_rect), true);
}

@("ExpandedByDistance.NegativeDistanceLatWithNorthPole")
unittest {
  S2LatLngRect rect = rectFromDegrees(0.0, -90.0, 90.0, 180.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.approxEquals(rectFromDegrees(5.0, 0.0, 85.0, 90.0)), true);
}

@("ExpandedByDistance.NegativeDistanceLatWithNorthPoleAndLngFull")
unittest {
  S2LatLngRect rect = rectFromDegrees(0.0, -180.0, 90.0, 180.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.approxEquals(rectFromDegrees(5.0, -180.0, 90.0, 180.0)), true);
}

@("ExpandedByDistance.NegativeDistanceLatWithSouthPole")
unittest {
  S2LatLngRect rect = rectFromDegrees(-90.0, -90.0, 0.0, 180.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.approxEquals(rectFromDegrees(-85.0, 0.0, -5.0, 90.0)), true);
}

@("ExpandedByDistance.NegativeDistanceLatWithSouthPoleAndLngFull")
unittest {
  S2LatLngRect rect = rectFromDegrees(-90.0, -180.0, 0.0, 180.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.approxEquals(rectFromDegrees(-90.0, -180.0, -5.0, 180.0)), true);
}

@("ExpandedByDistance.NegativeDistanceLngFull")
unittest {
  S2LatLngRect rect = rectFromDegrees(0.0, -180.0, 30.0, 180.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.approxEquals(rectFromDegrees(5.0, -180.0, 25.0, 180.0)), true);
}

@("ExpandedByDistance.NegativeDistanceLatResultEmpty")
unittest {
  S2LatLngRect rect = rectFromDegrees(0.0, 0.0, 9.9, 90.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  Assert.equal(rect.isEmpty(), true);
}

@("ExpandedByDistance.NegativeDistanceLngResultEmpty")
unittest {
  S2LatLngRect rect = rectFromDegrees(0.0, 0.0, 30.0, 11.0)
      .expandedByDistance(-S1Angle.fromDegrees(5.0));

  // The cap center is at latitude 30 - 5 = 25 degrees. The length of the
  // latitude 25 degree line is 0.906 times the length of the equator. Thus the
  // cap whose radius is 5 degrees covers the rectangle whose latitude interval
  // is 11 degrees.
  Assert.equal(rect.isEmpty(), true);
}

@("S2LatLngRect.GetCapBound")
unittest {
  // Bounding cap at center is smaller:
  Assert.equal(
      rectFromDegrees(-45, -45, 45, 45).getCapBound()
          .approxEquals(S2Cap.fromCenterHeight(S2Point(1, 0, 0), 0.5)),
      true);

  // Bounding cap at north pole is smaller:
  Assert.equal(
      rectFromDegrees(88, -80, 89, 80).getCapBound()
          .approxEquals(new S2Cap(S2Point(0, 0, 1), S1Angle.fromDegrees(2))),
      true);

  // Longitude span > 180 degrees:
  Assert.equal(
      rectFromDegrees(-30, -150, -10, 50).getCapBound()
          .approxEquals(new S2Cap(S2Point(0, 0, -1), S1Angle.fromDegrees(80))),
      true);
}

static void testCellOps(in S2LatLngRect r, in S2Cell cell, int level) {
  // Test the relationship between the given rectangle and cell:
  // 0 == no intersection, 1 == MayIntersect, 2 == Intersects,
  // 3 == Vertex Containment, 4 == Contains

  bool vertex_contained = false;
  foreach (i; 0 .. 4) {
    if (r.contains(cell.getVertexRaw(i)) ||
        (!r.isEmpty() && cell.contains(r.getVertex(i).toS2Point())))
      vertex_contained = true;
  }
  Assert.equal(r.mayIntersect(cell), level >= 1);
  Assert.equal(r.intersects(cell), level >= 2);
  Assert.equal(vertex_contained, level >= 3);
  Assert.equal(r.contains(cell), level >= 4);
}

@("S2LatLngRect.CellOps")
unittest {
  // contains(S2Cell), mayIntersect(S2Cell), intersects(S2Cell)

  // Special cases.
  testCellOps(S2LatLngRect.empty(), S2Cell.fromFacePosLevel(3, 0, 0), 0);
  testCellOps(S2LatLngRect.full(), S2Cell.fromFacePosLevel(2, 0, 0), 4);
  testCellOps(S2LatLngRect.full(), S2Cell.fromFacePosLevel(5, 0, 25), 4);

  // This rectangle includes the first quadrant of face 0.  It's expanded
  // slightly because cell bounding rectangles are slightly conservative.
  S2LatLngRect r4 = rectFromDegrees(-45.1, -45.1, 0.1, 0.1);
  testCellOps(r4, S2Cell.fromFacePosLevel(0, 0, 0), 3);
  testCellOps(r4, S2Cell.fromFacePosLevel(0, 0, 1), 4);
  testCellOps(r4, S2Cell.fromFacePosLevel(1, 0, 1), 0);

  // This rectangle intersects the first quadrant of face 0.
  S2LatLngRect r5 = rectFromDegrees(-10, -45, 10, 0);
  testCellOps(r5, S2Cell.fromFacePosLevel(0, 0, 0), 3);
  testCellOps(r5, S2Cell.fromFacePosLevel(0, 0, 1), 3);
  testCellOps(r5, S2Cell.fromFacePosLevel(1, 0, 1), 0);

  // Rectangle consisting of a single point.
  testCellOps(rectFromDegrees(4, 4, 4, 4), S2Cell.fromFace(0), 3);

  // Rectangles that intersect the bounding rectangle of a face
  // but not the face itself.
  testCellOps(rectFromDegrees(41, -87, 42, -79), S2Cell.fromFace(2), 1);
  testCellOps(rectFromDegrees(-41, 160, -40, -160), S2Cell.fromFace(5), 1);

  // This is the leaf cell at the top right hand corner of face 0.
  // It has two angles of 60 degrees and two of 120 degrees.
  S2Cell cell0tr = new S2Cell(S2Point(1 + 1e-12, 1, 1));
  S2LatLngRect bound0tr = cell0tr.getRectBound();
  S2LatLng v0 = S2LatLng(cell0tr.getVertexRaw(0));
  testCellOps(rectFromDegrees(v0.lat().degrees() - 1e-8,
                              v0.lng().degrees() - 1e-8,
                              v0.lat().degrees() - 2e-10,
                              v0.lng().degrees() + 1e-10),
              cell0tr, 1);

  // Rectangles that intersect a face but where no vertex of one region
  // is contained by the other region.  The first one passes through
  // a corner of one of the face cells.
  testCellOps(rectFromDegrees(-37, -70, -36, -20), S2Cell.fromFace(5), 2);

  // These two intersect like a diamond and a square.
  S2Cell cell202 = S2Cell.fromFacePosLevel(2, 0, 2);
  S2LatLngRect bound202 = cell202.getRectBound();
  testCellOps(rectFromDegrees(bound202.lo().lat().degrees() + 3,
                              bound202.lo().lng().degrees() + 3,
                              bound202.hi().lat().degrees() - 3,
                              bound202.hi().lng().degrees() - 3),
              cell202, 2);
}

/+
// TODO: Enable after encode and decode are added.
@("S2LatLngRect.EncodeDecode")
unittest {
  S2LatLngRect r = rectFromDegrees(-20, -80, 10, 20);
  Encoder encoder;
  r.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2LatLngRect decoded_rect = S2LatLngRect.empty();
  Assert.equal(true, decoded_rect.Decode(&decoder));
  Assert.equal(r, decoded_rect);
}
+/

@("S2LatLngRect.Area")
unittest {
  Assert.equal(S2LatLngRect.empty().area(), 0.0);
  Assert.approximately(S2LatLngRect.full().area(), 4 * math.PI, DOUBLE_ERR);
  Assert.approximately(rectFromDegrees(0, 0, 90, 90).area(), math.PI_2, DOUBLE_ERR);
}

// Recursively verify that when a rectangle is split into two pieces, the
// centroids of the children sum to give the centroid of the parent.
static void testCentroidSplitting(in S2LatLngRect r, int splits_left) {
  S2LatLngRect child0, child1;
  if (S2Testing.rnd.oneIn(2)) {
    double lat = S2Testing.rnd.uniformDouble(r.lat().lo(), r.lat().hi());
    child0 = new S2LatLngRect(R1Interval(r.lat().lo(), lat), r.lng());
    child1 = new S2LatLngRect(R1Interval(lat, r.lat().hi()), r.lng());
  } else {
    if (r.lng().lo() > r.lng().hi()) logger.logError("lng.lo() should be <= lng.hi()");
    double lng = S2Testing.rnd.uniformDouble(r.lng().lo(), r.lng().hi());
    child0 = new S2LatLngRect(r.lat(), S1Interval(r.lng().lo(), lng));
    child1 = new S2LatLngRect(r.lat(), S1Interval(lng, r.lng().hi()));
  }
  Assert.notGreaterThan(
      (r.getCentroid() - child0.getCentroid() - child1.getCentroid()).norm(),
      2e-15);
  if (splits_left > 0) {
    testCentroidSplitting(child0, splits_left - 1);
    testCentroidSplitting(child1, splits_left - 1);
  }
}

@("S2LatLngRect.GetCentroid")
unittest {
  Random rnd = S2Testing.rnd;

  // Empty and full rectangles.
  Assert.equal(S2LatLngRect.empty().getCentroid(), S2Point());
  Assert.notGreaterThan(S2LatLngRect.full().getCentroid().norm(), 1e-15);

  // Rectangles that cover the full longitude range.
  for (int i = 0; i < 100; ++i) {
    double lat1 = rnd.uniformDouble(-math.PI_2, math.PI_2);
    double lat2 = rnd.uniformDouble(-math.PI_2, math.PI_2);
    S2LatLngRect r = new S2LatLngRect(R1Interval.fromPointPair(lat1, lat2), S1Interval.full());
    S2Point centroid = r.getCentroid();
    Assert.approximately(0.5 * (math.sin(lat1) + math.sin(lat2)) * r.area(), centroid.z(), 1e-14);
    Assert.notGreaterThan(Vector2_d(centroid.x(), centroid.y()).norm(), 1e-15);
  }

  // Rectangles that cover the full latitude range.
  for (int i = 0; i < 100; ++i) {
    double lng1 = rnd.uniformDouble(-math.PI, math.PI);
    double lng2 = rnd.uniformDouble(-math.PI, math.PI);
    S2LatLngRect r =
        new S2LatLngRect(S2LatLngRect.fullLat(), S1Interval.fromPointPair(lng1, lng2));
    S2Point centroid = r.getCentroid();
    Assert.notGreaterThan(math.fabs(centroid.z()), 1e-15);
    Assert.approximately(r.lng().getCenter(), S2LatLng(centroid).lng().radians(), 1e-15);
    double alpha = 0.5 * r.lng().getLength();
    Assert.approximately(
        0.25 * math.PI * math.sin(alpha) / alpha * r.area(),
        Vector2_d(centroid.x(), centroid.y()).norm(),
        1e-15);
  }

  // Finally, verify that when a rectangle is recursively split into pieces,
  // the centroids of the pieces add to give the centroid of their parent.
  // To make the code simpler we avoid rectangles that cross the 180 degree
  // line of longitude.
  testCentroidSplitting(
      new S2LatLngRect(S2LatLngRect.fullLat(), S1Interval(-3.14, 3.14)), 10 /*splits_left*/);
}

// Returns the minimum distance from X to the latitude line segment defined by
// the given latitude and longitude interval.
S1Angle getDistance(in S2LatLng x, in S1Angle lat, in S1Interval interval) {
  Assert.equal(x.isValid(), true);
  Assert.equal(interval.isValid(), true);

  // Is X inside the longitude interval?
  if (interval.contains(x.lng().radians()))
    return (x.lat() - lat).abs();

  // Return the distance to the closer endpoint.
  return algorithm.min(
      x.getDistance(S2LatLng(lat, S1Angle.fromRadians(interval.lo()))),
      x.getDistance(S2LatLng(lat, S1Angle.fromRadians(interval.hi()))));
}

static S1Angle bruteForceDistance(in S2LatLngRect a, in S2LatLngRect b) {
  if (a.intersects(b))
    return S1Angle.fromRadians(0);

  // Compare every point in 'a' against every latitude edge and longitude edge
  // in 'b', and vice-versa, for a total of 16 point-vs-latitude-edge tests and
  // 16 point-vs-longitude-edge tests.
  S2LatLng[4] pnt_a, pnt_b;
  pnt_a[0] = S2LatLng(a.latLo(), a.lngLo());
  pnt_a[1] = S2LatLng(a.latLo(), a.lngHi());
  pnt_a[2] = S2LatLng(a.latHi(), a.lngHi());
  pnt_a[3] = S2LatLng(a.latHi(), a.lngLo());
  pnt_b[0] = S2LatLng(b.latLo(), b.lngLo());
  pnt_b[1] = S2LatLng(b.latLo(), b.lngHi());
  pnt_b[2] = S2LatLng(b.latHi(), b.lngHi());
  pnt_b[3] = S2LatLng(b.latHi(), b.lngLo());

  // Make arrays containing the lo/hi latitudes and the lo/hi longitude edges.
  S1Angle[2] lat_a = [ a.latLo(), a.latHi() ];
  S1Angle[2] lat_b = [ b.latLo(), b.latHi() ];
  S2Point[2][2] lng_edge_a = [ [ pnt_a[0].toS2Point(), pnt_a[3].toS2Point() ],
                               [ pnt_a[1].toS2Point(), pnt_a[2].toS2Point() ] ];
  S2Point[2][2] lng_edge_b = [ [ pnt_b[0].toS2Point(), pnt_b[3].toS2Point() ],
                               [ pnt_b[1].toS2Point(), pnt_b[2].toS2Point() ] ];

  S1Angle min_distance = S1Angle.fromDegrees(180.0);
  foreach (i; 0 .. 4) {
    // For each point in a and b.
    S2LatLng current_a = pnt_a[i];
    S2LatLng current_b = pnt_b[i];

    foreach (j; 0 .. 2) {
      // Get distances to latitude and longitude edges.
      S1Angle a_to_lat = getDistance(current_a, lat_b[j], b.lng());
      S1Angle b_to_lat = getDistance(current_b, lat_a[j], a.lng());
      S1Angle a_to_lng = s2.s2edge_distances.getDistance(
          current_a.toS2Point(), lng_edge_b[j][0], lng_edge_b[j][1]);
      S1Angle b_to_lng = s2.s2edge_distances.getDistance(
          current_b.toS2Point(), lng_edge_a[j][0], lng_edge_a[j][1]);

      min_distance = algorithm.min(
          min_distance,
          algorithm.min(a_to_lat, algorithm.min(b_to_lat, algorithm.min(a_to_lng, b_to_lng))));
    }
  }
  return min_distance;
}

static S1Angle bruteForceRectPointDistance(in S2LatLngRect a, in S2LatLng b) {
  if (a.contains(b)) {
    return S1Angle.fromRadians(0);
  }

  S1Angle b_to_lo_lat = getDistance(b, a.latLo(), a.lng());
  S1Angle b_to_hi_lat = getDistance(b, a.latHi(), a.lng());
  S1Angle b_to_lo_lng = s2.s2edge_distances.getDistance(
      b.toS2Point(),
      S2LatLng(a.latLo(), a.lngLo()).toS2Point(),
      S2LatLng(a.latHi(), a.lngLo()).toS2Point());
  S1Angle b_to_hi_lng = s2.s2edge_distances.getDistance(
      b.toS2Point(),
      S2LatLng(a.latLo(), a.lngHi()).toS2Point(),
      S2LatLng(a.latHi(), a.lngHi()).toS2Point());
  return algorithm.min(
      b_to_lo_lat, algorithm.min(b_to_hi_lat, algorithm.min(b_to_lo_lng, b_to_hi_lng)));
}

// This method verifies a.GetDistance(b) by comparing its result against a
// brute-force implementation. The correctness of the brute-force version is
// much easier to verify by inspection.
static void verifyGetDistance(in S2LatLngRect a, in S2LatLngRect b) {
  S1Angle distance1 = bruteForceDistance(a, b);
  S1Angle distance2 = a.getDistance(b);
  Assert.approximately(distance1.radians() - distance2.radians(), 0, 1e-10);
}

static S2LatLngRect pointrectFromDegrees(double lat, double lng) {
  return S2LatLngRect.fromPoint(
      S2LatLng.fromDegrees(lat, lng).normalized());
}

// This method verifies a.GetDistance(b), where b is a S2LatLng, by comparing
// its result against a.GetDistance(c), c being the point rectangle created
// from b.
static void verifyGetRectPointDistance(in S2LatLngRect a, in S2LatLng p) {
  S1Angle distance1 = bruteForceRectPointDistance(a, p.normalized());
  S1Angle distance2 = a.getDistance(p.normalized());
  Assert.approximately(math.fabs(distance1.radians() - distance2.radians()), 0, 1e-10);
}

@("S2LatLngRect.GetDistanceOverlapping")
unittest {
  // Check pairs of rectangles that overlap: (should all return 0):
  S2LatLngRect a = rectFromDegrees(0, 0, 2, 2);
  S2LatLngRect b = pointrectFromDegrees(0, 0);
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(a));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(b));
  Assert.equal(S1Angle.fromRadians(0), b.getDistance(b));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(S2LatLng.fromDegrees(0, 0)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(0, 1, 2, 3)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(0, 2, 2, 4)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(1, 0, 3, 2)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(2, 0, 4, 2)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(1, 1, 3, 3)));
  Assert.equal(S1Angle.fromRadians(0), a.getDistance(rectFromDegrees(2, 2, 4, 4)));
}

@("S2LatLngRect.GetDistanceRectVsPoint")
unittest {
  // Rect that spans 180.
  S2LatLngRect a = rectFromDegrees(-1, -1, 2, 1);
  verifyGetDistance(a, pointrectFromDegrees(-2, -1));
  verifyGetDistance(a, pointrectFromDegrees(1, 2));

  verifyGetDistance(pointrectFromDegrees(-2, -1), a);
  verifyGetDistance(pointrectFromDegrees(1, 2), a);

  verifyGetRectPointDistance(a, S2LatLng.fromDegrees(-2, -1));
  verifyGetRectPointDistance(a, S2LatLng.fromDegrees(1, 2));

  // Tests near the north pole.
  S2LatLngRect b = rectFromDegrees(86, 0, 88, 2);
  verifyGetDistance(b, pointrectFromDegrees(87, 3));
  verifyGetDistance(b, pointrectFromDegrees(87, -1));
  verifyGetDistance(b, pointrectFromDegrees(89, 1));
  verifyGetDistance(b, pointrectFromDegrees(89, 181));
  verifyGetDistance(b, pointrectFromDegrees(85, 1));
  verifyGetDistance(b, pointrectFromDegrees(85, 181));
  verifyGetDistance(b, pointrectFromDegrees(90, 0));

  verifyGetDistance(pointrectFromDegrees(87, 3), b);
  verifyGetDistance(pointrectFromDegrees(87, -1), b);
  verifyGetDistance(pointrectFromDegrees(89, 1), b);
  verifyGetDistance(pointrectFromDegrees(89, 181), b);
  verifyGetDistance(pointrectFromDegrees(85, 1), b);
  verifyGetDistance(pointrectFromDegrees(85, 181), b);
  verifyGetDistance(pointrectFromDegrees(90, 0), b);

  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, 3));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, -1));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 1));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 181));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 1));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 181));
  verifyGetRectPointDistance(b, S2LatLng.fromDegrees(90, 0));

  // Rect that touches the north pole.
  S2LatLngRect c = rectFromDegrees(88, 0, 90, 2);
  verifyGetDistance(c, pointrectFromDegrees(89, 3));
  verifyGetDistance(c, pointrectFromDegrees(89, 90));
  verifyGetDistance(c, pointrectFromDegrees(89, 181));
  verifyGetDistance(pointrectFromDegrees(89, 3), c);
  verifyGetDistance(pointrectFromDegrees(89, 90), c);
  verifyGetDistance(pointrectFromDegrees(89, 181), c);
}

@("S2LatLngRect.GetDistanceRectVsRect")
unittest {
  // Rect that spans 180.
  S2LatLngRect a = rectFromDegrees(-1, -1, 2, 1);
  verifyGetDistance(a, rectFromDegrees(0, 2, 1, 3));
  verifyGetDistance(a, rectFromDegrees(-2, -3, -1, -2));

  // Tests near the south pole.
  S2LatLngRect b = rectFromDegrees(-87, 0, -85, 3);
  verifyGetDistance(b, rectFromDegrees(-89, 1, -88, 2));
  verifyGetDistance(b, rectFromDegrees(-84, 1, -83, 2));
  verifyGetDistance(b, rectFromDegrees(-88, 90, -86, 91));
  verifyGetDistance(b, rectFromDegrees(-84, -91, -83, -90));
  verifyGetDistance(b, rectFromDegrees(-90, 181, -89, 182));
  verifyGetDistance(b, rectFromDegrees(-84, 181, -83, 182));
}

@("S2LatLngRect.GetDistanceRandomPairs")
unittest {
  // Test random pairs.
  foreach (i; 0 .. 10000) {
    S2LatLngRect a = S2LatLngRect.fromPointPair(
        S2LatLng(S2Testing.randomPoint()), S2LatLng(S2Testing.randomPoint()));
    S2LatLngRect b = S2LatLngRect.fromPointPair(
        S2LatLng(S2Testing.randomPoint()), S2LatLng(S2Testing.randomPoint()));
    verifyGetDistance(a, b);

    S2LatLng c = S2LatLng(S2Testing.randomPoint());
    verifyGetRectPointDistance(a, c);
    verifyGetRectPointDistance(b, c);
  }
}

// This function assumes that GetDirectedHausdorffDistance() always returns
// a distance from some point in a to b. So the function mainly tests whether
// the returned distance is large enough, and only does a weak test on whether
// it is small enough.
static void verifyGetDirectedHausdorffDistance(in S2LatLngRect a, in S2LatLngRect b) {
  S1Angle hausdorff_distance = a.getDirectedHausdorffDistance(b);

  static const double kResolution = 0.1;
  // Record the max sample distance as well as the sample point realizing the
  // max for easier debugging.
  S1Angle max_distance;
  double lat_max, lng_max;

  int sample_size_on_lat = cast(int)(a.lat().getLength() / kResolution) + 1;
  int sample_size_on_lng = cast(int)(a.lng().getLength() / kResolution) + 1;
  double delta_on_lat = a.lat().getLength() / sample_size_on_lat;
  double delta_on_lng = a.lng().getLength() / sample_size_on_lng;

  double lng = a.lng().lo();
  for (int i = 0; i <= sample_size_on_lng; ++i, lng += delta_on_lng) {
    double lat = a.lat().lo();
    for (int j = 0; j <= sample_size_on_lat; ++j, lat += delta_on_lat) {
      S2LatLng latlng = S2LatLng.fromRadians(lat, lng).normalized();
      S1Angle distance_to_b = b.getDistance(latlng);

      if (distance_to_b >= max_distance) {
        max_distance = distance_to_b;
        lat_max = lat;
        lng_max = lng;
      }
    }
  }

  Assert.notGreaterThan(max_distance.radians(), hausdorff_distance.radians() + 1e-10);
  Assert.notLessThan(max_distance.radians(), hausdorff_distance.radians() - kResolution);
}

@("S2LatLngRect.GetDirectedHausdorffDistanceRandomPairs")
unittest {
  // Test random pairs.
  const int kIters = 1000;
  foreach (i; 0 .. kIters) {
    S2LatLngRect a = S2LatLngRect.fromPointPair(
        S2LatLng(S2Testing.randomPoint()), S2LatLng(S2Testing.randomPoint()));
    S2LatLngRect b = S2LatLngRect.fromPointPair(
        S2LatLng(S2Testing.randomPoint()), S2LatLng(S2Testing.randomPoint()));
    // a and b are *minimum* bounding rectangles of two random points, in
    // particular, their Voronoi diagrams are always of the same topology. We
    // take the "complements" of a and b for more thorough testing.
    S2LatLngRect a2 = new S2LatLngRect(a.lat(), a.lng().complement());
    S2LatLngRect b2 = new S2LatLngRect(b.lat(), b.lng().complement());

    // Note that "a" and "b" come from the same distribution, so there is no
    // need to test pairs such as (b, a), (b, a2), etc.
    verifyGetDirectedHausdorffDistance(a, b);
    verifyGetDirectedHausdorffDistance(a, b2);
    verifyGetDirectedHausdorffDistance(a2, b);
    verifyGetDirectedHausdorffDistance(a2, b2);
  }
}

@("S2LatLngRect.GetDirectedHausdorffDistanceContained")
unittest {
  // Caller rect is contained in callee rect. Should return 0.
  S2LatLngRect a = rectFromDegrees(-10, 20, -5, 90);
  Assert.equal(
      S1Angle.fromRadians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-10, 20, -5, 90)));
  Assert.equal(
      S1Angle.fromRadians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-10, 19, -5, 91)));
  Assert.equal(
      S1Angle.fromRadians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-11, 20, -4, 90)));
  Assert.equal(
      S1Angle.fromRadians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-11, 19, -4, 91)));
}

@("S2LatLngRect.GetDirectHausdorffDistancePointToRect")
unittest {
  // The Hausdorff distance from a point to a rect should be the same as its
  // distance to the rect.
  S2LatLngRect a1 = pointrectFromDegrees(5, 8);
  S2LatLngRect a2 = pointrectFromDegrees(90, 10);  // north pole

  S2LatLngRect b = rectFromDegrees(-85, -50, -80, 10);
  Assert.approximately(
      a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians(), DOUBLE_ERR);
  Assert.approximately(
      a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians(), DOUBLE_ERR);

  b = rectFromDegrees(4, -10, 80, 10);
  Assert.approximately(
      a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians(), DOUBLE_ERR);
  Assert.approximately(
      a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians(), DOUBLE_ERR);

  b = rectFromDegrees(70, 170, 80, -170);
  Assert.approximately(
      a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians(), DOUBLE_ERR);
  Assert.approximately(
      a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians(), DOUBLE_ERR);
}

@("S2LatLngRect.getDirectedHausdorffDistanceRectToPoint")
unittest {
  S2LatLngRect a = rectFromDegrees(1, -8, 10, 20);
  verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(5, 8));
  verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(-6, -100));
  // south pole
  verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(-90, -20));
  // north pole
  verifyGetDirectedHausdorffDistance(a, pointrectFromDegrees(90, 0));
}

@("S2LatLngRect.getDirectedHausdorffDistanceRectToRectNearPole")
unittest {
  // Tests near south pole.
  S2LatLngRect a = rectFromDegrees(-87, 0, -85, 3);
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-89, 1, -88, 2));
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 1, -83, 2));
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-88, 90, -86, 91));
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, -91, -83, -90));
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-90, 181, -89, 182));
  verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 181, -83, 182));
}

@("S2LatLngRect.getDirectedHausdorffDistanceRectToRectDegenerateCases")
unittest {
  // Rectangles that contain poles.
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(0, 10, 90, 20), rectFromDegrees(-4, -10, 4, 0));
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(-4, -10, 4, 0), rectFromDegrees(0, 10, 90, 20));

  // Two rectangles share same or complement longitudinal intervals.
  S2LatLngRect a = rectFromDegrees(-50, -10, 50, 10);
  S2LatLngRect b = rectFromDegrees(30, -10, 60, 10);
  verifyGetDirectedHausdorffDistance(a, b);
  S2LatLngRect c = new S2LatLngRect(a.lat(), a.lng().complement());
  verifyGetDirectedHausdorffDistance(c, b);

  // rectangle a touches b_opposite_lng.
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(10, 170, 30, 180), rectFromDegrees(-50, -10, 50, 10));
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(10, -180, 30, -170), rectFromDegrees(-50, -10, 50, 10));

  // rectangle b's Voronoi diagram is degenerate (lng interval spans 180
  // degrees), and a touches the degenerate Voronoi vertex.
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(-30, 170, 30, 180), rectFromDegrees(-10, -90, 10, 90));
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(-30, -180, 30, -170), rectFromDegrees(-10, -90, 10, 90));

  // rectangle a touches a voronoi vertex of rectangle b.
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(-20, 105, 20, 110), rectFromDegrees(-30, 5, 30, 15));
  verifyGetDirectedHausdorffDistance(
      rectFromDegrees(-20, 95, 20, 105), rectFromDegrees(-30, 5, 30, 15));
}
