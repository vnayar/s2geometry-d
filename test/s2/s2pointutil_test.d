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

module s2.s2pointutil_test;

import fluent.asserts;
import s2.s2pointutil;
import s2.util.math.matrix3x3;
import s2.s2point;
import s2.s2edge_distances;
import s2.s2testing;
import s2.s1angle;
import s2.s2measures;
import math = std.math;

immutable double DBL_ERR = 0.0001;

@("Frames")
unittest {
  Matrix3x3_d m;
  S2Point z = S2Point(0.2, 0.5, -3.3).normalize();
  getFrame(z, m);
  Assert.equal(approxEquals(m.col(2), z), true);
  Assert.equal(isUnitLength(m.col(0)), true);
  Assert.equal(isUnitLength(m.col(1)), true);
  Assert.approximately(m.det(), 1.0, DBL_ERR);

  Assert.equal(approxEquals(toFrame(m, m.col(0)), S2Point(1, 0, 0)), true);
  Assert.equal(approxEquals(toFrame(m, m.col(1)), S2Point(0, 1, 0)), true);
  Assert.equal(approxEquals(toFrame(m, m.col(2)), S2Point(0, 0, 1)), true);

  Assert.equal(approxEquals(fromFrame(m, S2Point(1, 0, 0)), m.col(0)), true);
  Assert.equal(approxEquals(fromFrame(m, S2Point(0, 1, 0)), m.col(1)), true);
  Assert.equal(approxEquals(fromFrame(m, S2Point(0, 0, 1)), m.col(2)), true);
}

private void testRotate(in S2Point p, in S2Point axis, in S1Angle angle) {
  S2Point result = rotate(p, axis, angle);

  // "result" should be unit length.
  Assert.equal(isUnitLength(result), true);

  // "result" and "p" should be the same distance from "axis".
  double kMaxPositionError = 1e-15;
  Assert.notGreaterThan((S1Angle(result, axis) - S1Angle(p, axis)).abs().radians(),
            kMaxPositionError);

  // Check that the rotation angle is correct.  We allow a fixed error in the
  // *position* of the result, so we need to convert this into a rotation
  // angle.  The allowable error can be very large as "p" approaches "axis".
  double axis_distance = p.crossProd(axis).norm();
  double max_rotation_error;
  if (axis_distance < kMaxPositionError) {
    max_rotation_error = 2 * math.PI;
  } else {
    max_rotation_error = math.asin(kMaxPositionError / axis_distance);
  }
  double actual_rotation = turnAngle(p, axis, result) + math.PI;
  double rotation_error = math.remainder(angle.radians() - actual_rotation, 2 * math.PI);
  Assert.notGreaterThan(rotation_error, max_rotation_error);
}

@("Rotate")
unittest {
  foreach (int iter; 0 .. 1000) {
    S2Point axis = S2Testing.randomPoint();
    S2Point target = S2Testing.randomPoint();
    // Choose a distance whose logarithm is uniformly distributed.
    double distance = math.PI * math.pow(1e-15, S2Testing.rnd.randDouble());
    // Sometimes choose points near the far side of the axis.
    if (S2Testing.rnd.oneIn(5)) {
      distance = math.PI - distance;
    }
    S2Point p = interpolateAtDistance(S1Angle.fromRadians(distance), axis, target);
    // Choose the rotation angle.
    double angle = 2 * math.PI * math.pow(1e-15, S2Testing.rnd.randDouble());
    if (S2Testing.rnd.oneIn(3)) angle = -angle;
    if (S2Testing.rnd.oneIn(10)) angle = 0;
    testRotate(p, axis, S1Angle.fromRadians(angle));
  }
}

/*
// Given a point P, return the minimum level at which an edge of some S2Cell
// parent of P is nearly collinear with S2::Origin().  This is the minimum
// level for which Sign() may need to resort to expensive calculations in
// order to determine which side of an edge the origin lies on.
static int GetMinExpensiveLevel(const S2Point& p) {
  S2CellId id(p);
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2Cell cell(id.parent(level));
    for (int k = 0; k < 4; ++k) {
      S2Point a = cell.GetVertex(k);
      S2Point b = cell.GetVertex(k + 1);
      if (s2pred::TriageSign(a, b, S2::Origin(), a.CrossProd(b)) == 0) {
        return level;
      }
    }
  }
  return S2CellId::kMaxLevel + 1;
}

TEST(S2, OriginTest) {
  // To minimize the number of expensive Sign() calculations,
  // S2::Origin() should not be nearly collinear with any commonly used edges.
  // Two important categories of such edges are:
  //
  //  - edges along a line of longitude (reasonably common geographically)
  //  - S2Cell edges (used extensively when computing S2Cell coverings)
  //
  // This implies that the origin:
  //
  //  - should not be too close to either pole (since all lines of longitude
  //    converge at the poles)
  //  - should not be colinear with edges of any S2Cell except for very small
  //    ones (which are used less frequently)
  //
  // The point chosen below is about 66km from the north pole towards the East
  // Siberian Sea.  The purpose of the STtoUV(2/3) calculation is to keep the
  // origin as far away as possible from the longitudinal edges of large
  // S2Cells.  (The line of longitude through the chosen point is always 1/3
  // or 2/3 of the way across any S2Cell with longitudinal edges that it
  // passes through.)

  EXPECT_EQ(S2Point(-0.01, 0.01 * S2::STtoUV(2./3), 1).Normalize(),
            S2::Origin());

  // Check that the origin is not too close to either pole.  (We don't use
  // S2Earth because we don't want to depend on that package.)
  double distance_km = acos(S2::Origin().z()) * S2Testing::kEarthRadiusKm;
  EXPECT_GE(distance_km, 50.0);
  LOG(INFO) << "\nS2::Origin() coordinates: " << S2LatLng(S2::Origin())
            << ", distance from pole: " << distance_km << " km";

  // Check that S2::Origin() is not collinear with the edges of any large
  // S2Cell.  We do this is two parts.  For S2Cells that belong to either
  // polar face, we simply need to check that S2::Origin() is not nearly
  // collinear with any edge of any cell that contains it (except for small
  // cells < 3 meters across).
  EXPECT_GE(GetMinExpensiveLevel(S2::Origin()), 22);

  // For S2Cells that belong to the four non-polar faces, only longitudinal
  // edges can possibly be colinear with S2::Origin().  We check these edges
  // by projecting S2::Origin() onto the equator, and then testing all S2Cells
  // that contain this point to make sure that none of their edges are nearly
  // colinear with S2::Origin() (except for small cells < 3 meters across).
  S2Point equator_point(S2::Origin().x(), S2::Origin().y(), 0);
  EXPECT_GE(GetMinExpensiveLevel(equator_point), 22);
}
*/
