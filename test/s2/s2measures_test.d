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

module s2.s2measures_test;

import fluent.asserts;
import s2.s2measures;
import s2.s2latlng;
import s2.s2testing;
import s2.s2point;
import algorithm = std.algorithm;
import math = std.math;
import std.stdio;

immutable double DBL_ERR = 0.0001;

@("S2.AngleMethods")
unittest {
  S2Point pz = S2Point(0, 0, 1);
  S2Point p000 = S2Point(1, 0, 0);
  S2Point p045 = S2Point(1, 1, 0).normalize();
  S2Point p090 = S2Point(0, 1, 0);
  S2Point p180 = S2Point(-1, 0, 0);

  Assert.approximately(angle(p000, pz, p045), math.PI_4, DBL_ERR);
  Assert.approximately(turnAngle(p000, pz, p045), -3 * math.PI_4, DBL_ERR);

  Assert.approximately(angle(p045, pz, p180), 3 * math.PI_4, DBL_ERR);
  Assert.approximately(turnAngle(p045, pz, p180), -math.PI_4, DBL_ERR);

  Assert.approximately(angle(p000, pz, p180), math.PI, DBL_ERR);
  Assert.approximately(turnAngle(p000, pz, p180), 0, DBL_ERR);

  Assert.approximately(angle(pz, p000, p045), math.PI_2, DBL_ERR);
  Assert.approximately(turnAngle(pz, p000, p045), math.PI_2, DBL_ERR);

  Assert.approximately(angle(pz, p000, pz), 0, DBL_ERR);
  Assert.approximately(math.fabs(turnAngle(pz, p000, pz)), math.PI, DBL_ERR);
}

@("S2.AreaMethods")
unittest {
  S2Point pz = S2Point(0, 0, 1);
  S2Point p000 = S2Point(1, 0, 0);
  S2Point p045 = S2Point(1, 1, 0).normalize();
  S2Point p090 = S2Point(0, 1, 0);
  S2Point p180 = S2Point(-1, 0, 0);

  Assert.approximately(area(p000, p090, pz), math.PI_2, DBL_ERR);
  Assert.approximately(area(p045, pz, p180), 3 * math.PI_4, DBL_ERR);

  // Make sure that Area() has good *relative* accuracy even for
  // very small areas.
  static const double eps = 1e-10;
  S2Point pepsx = S2Point(eps, 0, 1).normalize();
  S2Point pepsy = S2Point(0, eps, 1).normalize();
  double expected1 = 0.5 * eps * eps;
  Assert.approximately(area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

  // Make sure that it can handle degenerate triangles.
  S2Point pr = S2Point(0.257, -0.5723, 0.112).normalize();
  S2Point pq = S2Point(-0.747, 0.401, 0.2235).normalize();
  Assert.equal(area(pr, pr, pr), 0.0);
  // The following test is not exact due to rounding error.
  Assert.approximately(area(pr, pq, pr), 0, 1e-15);
  Assert.equal(area(p000, p045, p090), 0.0);

  double max_girard = 0;
  for (int i = 0; i < 10000; ++i) {
    S2Point p0 = S2Testing.randomPoint();
    S2Point d1 = S2Testing.randomPoint();
    S2Point d2 = S2Testing.randomPoint();
    S2Point p1 = (p0 + 1e-15 * d1).normalize();
    S2Point p2 = (p0 + 1e-15 * d2).normalize();
    // The actual displacement can be as much as 1.2e-15 due to roundoff.
    // This yields a maximum triangle area of about 0.7e-30.
    Assert.notGreaterThan(area(p0, p1, p2), 0.7e-30);
    max_girard = algorithm.max(max_girard, girardArea(p0, p1, p2));
  }
  // This check only passes if GirardArea() uses RobustCrossProd().
  writeln("Worst case Girard for triangle area 1e-30: ", max_girard);
  Assert.notGreaterThan(max_girard, 1e-14);

  // Try a very long and skinny triangle.
  S2Point p045eps = S2Point(1, 1, eps).normalize();
  double expected2 = 5.8578643762690495119753e-11;  // Mathematica.
  Assert.approximately(area(p000, p045eps, p090), expected2, 1e-9 * expected2);

  // Triangles with near-180 degree edges that sum to a quarter-sphere.
  static const double eps2 = 1e-14;
  S2Point p000eps2 = S2Point(1, 0.1*eps2, eps2).normalize();
  double quarter_area1 = area(p000eps2, p000, p045)
      + area(p000eps2, p045, p180)
      + area(p000eps2, p180, pz)
      + area(p000eps2, pz, p000);
  Assert.approximately(quarter_area1, math.PI, DBL_ERR);

  // Four other triangles that sum to a quarter-sphere.
  S2Point p045eps2 = S2Point(1, 1, eps2).normalize();
  double quarter_area2 = area(p045eps2, p000, p045)
      + area(p045eps2, p045, p180)
      + area(p045eps2, p180, pz)
      + area(p045eps2, pz, p000);
  Assert.approximately(quarter_area2, math.PI, DBL_ERR);

  // Compute the area of a hemisphere using four triangles with one near-180
  // degree edge and one near-degenerate edge.
  for (int i = 0; i < 100; ++i) {
    double lng = 2 * math.PI * S2Testing.rnd.randDouble();
    S2Point p0 = S2LatLng.fromRadians(1e-20, lng).normalized().toS2Point();
    S2Point p1 = S2LatLng.fromRadians(0, lng).normalized().toS2Point();
    double p2_lng = lng + S2Testing.rnd.randDouble();
    S2Point p2 = S2LatLng.fromRadians(0, p2_lng).normalized().toS2Point();
    S2Point p3 = S2LatLng.fromRadians(0, lng + math.PI).normalized().toS2Point();
    S2Point p4 = S2LatLng.fromRadians(0, lng + 5.0).normalized().toS2Point();
    double area = area(p0, p1, p2) + area(p0, p2, p3)
        + area(p0, p3, p4) + area(p0, p4, p1);
    Assert.approximately(area, 2 * math.PI, 2e-15);
  }
}
