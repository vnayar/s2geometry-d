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
//

// Original Author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2latlng_rect_bounder_test;

import s2.r1interval;
import s2.s1angle;
import s2.s1interval;
import s2.s2edge_distances : interpolateAtDistance;
import s2.s2edge_distances;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2point;
import s2.s2pointutil : isUnitLength, robustCrossProd;
import s2.s2testing;
import s2.logger;

import fluent.asserts;

import std.math;

S2LatLngRect getEdgeBound(in S2Point a, in S2Point b) {
  auto bounder = new S2LatLngRectBounder();
  bounder.addPoint(a);
  bounder.addPoint(b);
  return bounder.getBound();
}

S2LatLngRect getEdgeBound(
    double x1, double y1, double z1,
    double x2, double y2, double z2) {
  return getEdgeBound(S2Point(x1, y1, z1).normalize(), S2Point(x2, y2, z2).normalize());
}

const S2LatLng kRectError = S2LatLngRectBounder.maxErrorForTests();

@("RectBounder.MaxLatitudeSimple")
unittest {
  // Check cases where the min/max latitude is attained at a vertex.
  const double kCubeLat = asin(1 / sqrt(3.0));  // 35.26 degrees
  Assert.equal(
      getEdgeBound(1,1,1, 1,-1,-1).approxEquals(
          new S2LatLngRect(R1Interval(-kCubeLat, kCubeLat), S1Interval(-PI_4, PI_4)), kRectError),
      true);
  Assert.equal(
      getEdgeBound(1,-1,1, 1,1,-1).approxEquals(
          new S2LatLngRect(R1Interval(-kCubeLat, kCubeLat), S1Interval(-PI_4, PI_4)), kRectError),
      true);

  // Check cases where the min/max latitude occurs in the edge interior.
  // These tests expect the result to be pretty close to the middle of the
  // allowable error range (i.e., by adding 0.5 * kRectError).

  // Max latitude, CW edge
  Assert.approximately(
      PI_4 + 0.5 * kRectError.lat().radians(),
      getEdgeBound(1,1,1, 1,-1,1).lat().hi(),
      DOUBLE_ERR);
  // Max latitude, CCW edge
  Assert.approximately(
      PI_4 + 0.5 * kRectError.lat().radians(),
      getEdgeBound(1,-1,1, 1,1,1).lat().hi(),
      DOUBLE_ERR);
  // Min latitude, CW edge
  Assert.approximately(
      -PI_4 - 0.5 * kRectError.lat().radians(),
      getEdgeBound(1,-1,-1, -1,-1,-1).lat().lo(),
      DOUBLE_ERR);
  // Min latitude, CCW edge
  Assert.approximately(
      -PI_4 - 0.5 * kRectError.lat().radians(),
      getEdgeBound(-1,1,-1, -1,-1,-1).lat().lo(),
      DOUBLE_ERR);

  // Check cases where the edge passes through one of the poles.
  Assert.approximately(PI_2, getEdgeBound(.3,.4,1, -.3,-.4,1).lat().hi(), DOUBLE_ERR);
  Assert.approximately(-PI_2, getEdgeBound(.3,.4,-1, -.3,-.4,-1).lat().lo(), DOUBLE_ERR);
}

@("RectBounder.MaxLatitudeRandom")
unittest {
  // Check that the maximum latitude of edges is computed accurately to within
  // 3 * DBL_EPSILON (the expected maximum error).  We concentrate on maximum
  // latitudes near the equator and north pole since these are the extremes.

  auto rnd = S2Testing.rnd;
  const int kIters = 100;
  for (int iter = 0; iter < kIters; ++iter) {
    // Construct a right-handed coordinate frame (U,V,W) such that U points
    // slightly above the equator, V points at the equator, and W is slightly
    // offset from the north pole.
    S2Point u = S2Testing.randomPoint();
    u[2] = double.epsilon * 1e-6 * pow(1e12, rnd.randDouble());  // log is uniform
    u = u.normalize();
    S2Point v = robustCrossProd(S2Point(0, 0, 1), u).normalize();
    S2Point w = robustCrossProd(u, v).normalize();

    // Construct a line segment AB that passes through U, and check that the
    // maximum latitude of this segment matches the latitude of U.
    S2Point a = (u - rnd.randDouble() * v).normalize();
    S2Point b = (u + rnd.randDouble() * v).normalize();
    S2LatLngRect ab_bound = getEdgeBound(a, b);
    Assert.approximately(
        S2LatLng.latitude(u).radians(), ab_bound.lat().hi(), kRectError.lat().radians());

    // Construct a line segment CD that passes through W, and check that the
    // maximum latitude of this segment matches the latitude of W.
    S2Point c = (w - rnd.randDouble() * v).normalize();
    S2Point d = (w + rnd.randDouble() * v).normalize();
    S2LatLngRect cd_bound = getEdgeBound(c, d);
    Assert.approximately(
        S2LatLng.latitude(w).radians(), cd_bound.lat().hi(), kRectError.lat().radians());
  }
}

S2Point perturbATowardsB(in S2Point a, in S2Point b) {
  auto rnd = S2Testing.rnd;
  double choice = rnd.randDouble();
  if (choice < 0.1) {
    return a;
  }
  if (choice < 0.3) {
    // Return a point that is exactly proportional to A and that still
    // satisfies S2::IsUnitLength().
    for (;;) {
      S2Point b2 = (2 - a.norm() + 5*(rnd.randDouble()-0.5) * double.epsilon) * a;
      if (b2 != a && isUnitLength(b2))
        return b2;
    }
  }
  if (choice < 0.5) {
    // Return a point such that the distance squared to A will underflow.
    return interpolateAtDistance(S1Angle.fromRadians(1e-300), a, b);
  }
  // Otherwise return a point whose distance from A is near DBL_EPSILON such
  // that the log of the pdf is uniformly distributed.
  double distance = double.epsilon * 1e-5 * pow(1e6, rnd.randDouble());
  return interpolateAtDistance(S1Angle.fromRadians(distance), a, b);
}

S2Point randomPole() {
  return S2Point(0, 0, S2Testing.rnd.oneIn(2) ? 1 : -1);
}

S2Point pointNearPole() {
  return perturbATowardsB(randomPole(), S2Testing.randomPoint());
}

S2Point pointNearEquator() {
  return perturbATowardsB(
      S2Point(S2Testing.rnd.randDouble(), S2Testing.rnd.randDouble(), 0).normalize(), randomPole());
}

@("RectBounder.NearlyIdenticalOrAntipodalPoints")
unittest {
  // Test pairs of points that are either:
  //  - identical
  //  - nearly or exactly proportional, e.g. (1,0,0) vs. (1+2e-16, 0, 0)
  //  - very close to each other
  // Furthermore we want to test cases where the two points are:
  //  - on a nearly-polar great circle
  //  - on a nearly-equatorial great circle
  //  - near the poles, but on any great circle
  //  - near the equator, but on any great circle
  //  - positioned arbitrarily
  // Also test the corresponding situations for antipodal points, i.e. by
  // negating one of the points so that they are almost 180 degrees apart.

  auto rnd = S2Testing.rnd;
  const int kIters = 10000;
  for (int iter = 0; iter < kIters; ++iter) {
    logger.logTrace("Iteration ", iter);
    S2Point a, b;
    final switch (rnd.uniform(5)) {
      case 0:
        // Two nearby points on a nearly-polar great circle.
        a = S2Testing.randomPoint();
        b = perturbATowardsB(a, pointNearPole());
        break;
      case 1:
        // Two nearby points on a nearly-equatorial great circle.
        a = pointNearEquator();
        b = perturbATowardsB(a, pointNearEquator());
        break;
      case 2:
        // Two nearby points near a pole, but on any great circle.
        a = pointNearPole();
        b = perturbATowardsB(a, S2Testing.randomPoint());
        break;
      case 3:
        // Two nearby points near the equator, but on any great circle.
        a = pointNearEquator();
        b = perturbATowardsB(a, S2Testing.randomPoint());
        break;
      case 4:
        // Two nearby points anywhere on the sphere.
        a = S2Testing.randomPoint();
        b = perturbATowardsB(a, S2Testing.randomPoint());
        break;
    }
    // The two points are chosen to be so close to each other that the min/max
    // latitudes are nearly always achieved at the edge endpoints.  The only
    // thing we need to watch out for is that the latitude error bound is
    // slightly larger if the min/max latitude occurs in the edge interior.
    S2LatLngRect expected_bound = S2LatLngRect.fromPointPair(S2LatLng(a), S2LatLng(b));
    S2LatLngRect bound = getEdgeBound(a, b);
    Assert.equal(bound.contains(expected_bound), true);
    Assert.equal(expected_bound.expanded(kRectError).polarClosure().contains(bound), true);

    // If the two points are close enough and one point is negated (antipodal
    // points), the bound should be the entire sphere.
    if ((a - b).crossProd(a + b).norm() <= 6.110 * double.epsilon) {
      Assert.equal(S2LatLngRect.full(), getEdgeBound(a, -b));
    }
  }
}

S2LatLngRect getSubregionBound(double x_lat, double x_lng, double y_lat, double y_lng) {
  S2LatLngRect input = S2LatLngRect.fromPointPair(
      S2LatLng.fromRadians(x_lat, x_lng),
      S2LatLng.fromRadians(y_lat, y_lng));
  S2LatLngRect output = S2LatLngRectBounder.expandForSubregions(input);

  // Test that the bound is actually expanded.
  Assert.equal(output.contains(input), true);
  if (input.lat() == S2LatLngRect.fullLat()) {
    Assert.equal(input.lat().contains(output.lat()), false);
  }
  return output;
}

@("RectBounder.ExpandForSubregions")
unittest {
  // First we check the various situations where the bound contains
  // nearly-antipodal points.  The tests are organized into pairs where the
  // two bounds are similar except that the first bound meets the
  // nearly-antipodal criteria while the second does not.

  // Cases where the bound does not straddle the equator (but almost does),
  // and spans nearly 180 degrees in longitude.
  Assert.equal(getSubregionBound(3e-16, 0, 1e-14, M_PI).isFull(), true);
  Assert.equal(getSubregionBound(9e-16, 0, 1e-14, M_PI).isFull(), false);
  Assert.equal(getSubregionBound(1e-16, 7e-16, 1e-14, M_PI).isFull(), true);
  Assert.equal(getSubregionBound(3e-16, 14e-16, 1e-14, M_PI).isFull(), false);
  Assert.equal(getSubregionBound(1e-100, 14e-16, 1e-14, M_PI).isFull(), true);
  Assert.equal(getSubregionBound(1e-100, 22e-16, 1e-14, M_PI).isFull(), false);

  // Cases where the bound spans at most 90 degrees in longitude, and almost
  // 180 degrees in latitude.  Note that DBL_EPSILON is about 2.22e-16, which
  // implies that the double-precision value just below Pi/2 can be written as
  // (M_M_PI_2 - 2e-16).
  Assert.equal(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 7e-16, 0).isFull(), true);
  Assert.equal(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 30e-16, 0).isFull(), false);
  Assert.equal(getSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 2e-16, 1e-7).isFull(), true);
  Assert.equal(getSubregionBound(-M_PI_2 + 30e-16, 0, M_PI_2, 1e-7).isFull(), false);
  Assert.equal(getSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 4e-16, M_PI_2).isFull(), true);
  Assert.equal(getSubregionBound(-M_PI_2, 0, M_PI_2 - 30e-16, M_PI_2).isFull(), false);

  // Cases where the bound straddles the equator and spans more than 90
  // degrees in longitude.  These are the cases where the critical distance is
  // between a corner of the bound and the opposite longitudinal edge.  Unlike
  // the cases above, here the bound may contain nearly-antipodal points (to
  // within 3.055 * DBL_EPSILON) even though the latitude and longitude ranges
  // are both significantly less than (Pi - 3.055 * DBL_EPSILON).
  Assert.equal(getSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-8, M_PI - 1e-7).isFull(), true);
  Assert.equal(getSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-7, M_PI - 1e-7).isFull(), false);
  Assert.equal(getSubregionBound(-M_PI_2 + 1e-12, -M_PI + 1e-4, M_PI_2, 0).isFull(), true);
  Assert.equal(getSubregionBound(-M_PI_2 + 1e-11, -M_PI + 1e-4, M_PI_2, 0).isFull(), true);

  // Now we test cases where the bound does not contain nearly-antipodal
  // points, but it does contain points that are approximately 180 degrees
  // apart in latitude.
  Assert.equal(
      getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 2e-16).approxEquals(
          new S2LatLngRect(R1Interval(1.5, 1.5), S1Interval.full()), kRectError),
      true);
  Assert.equal(
      getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 7e-16).approxEquals(
          new S2LatLngRect(R1Interval(1.5, 1.5), S1Interval(-M_PI_2, M_PI_2 - 7e-16)), kRectError),
      true);

  // Test the full and empty bounds.
  Assert.equal(S2LatLngRectBounder.expandForSubregions(S2LatLngRect.full()).isFull(), true);
  Assert.equal(S2LatLngRectBounder.expandForSubregions(S2LatLngRect.empty()).isEmpty(), true);

  // Check for cases where the bound is expanded to include one of the poles.
  Assert.equal(
      getSubregionBound(-M_PI_2 + 1e-15, 0, -M_PI_2 + 1e-15, 0).approxEquals(
          new S2LatLngRect(R1Interval(-M_PI_2, -M_PI_2 + 1e-15), S1Interval.full()), kRectError),
      true);
  Assert.equal(
      getSubregionBound(M_PI_2 - 1e-15, 0, M_PI_2 - 1e-15, 0).approxEquals(
          new S2LatLngRect(R1Interval(M_PI_2 - 1e-15, M_PI_2), S1Interval.full()), kRectError),
      true);
}
