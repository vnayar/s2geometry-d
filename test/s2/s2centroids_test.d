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

module s2.s2centroids_test;

import s2.s2centroids;
import s2.s2point;
import s2.s2testing;
import s2.util.math.vector;

import fluent.asserts;
import std.math;

@("TrueCentroid") unittest {
  // Test TrueCentroid() with very small triangles.  This test assumes that
  // the triangle is small enough so that it is nearly planar.
  for (int i = 0; i < 100; ++i) {
    Vector3_d p, x, y;
    S2Testing.getRandomFrame(p, x, y);
    double d = 1e-4 * pow(1e-4, S2Testing.rnd.randDouble());
    S2Point p0 = (p - d * x).normalize();
    S2Point p1 = (p + d * x).normalize();
    S2Point p2 = (p + 3 * d * y).normalize();
    S2Point centroid = trueCentroid(p0, p1, p2).normalize();

    // The centroid of a planar triangle is at the intersection of its
    // medians, which is two-thirds of the way along each median.
    S2Point expected_centroid = (p + d * y).normalize();
    Assert.notGreaterThan(centroid.angle(expected_centroid), 2e-8);
  }
}
