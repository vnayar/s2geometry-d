// Copyright 2016 Google Inc. All Rights Reserved.
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

module s2.s2polyline_simplifier_test;

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2edge_distances;
import s2.s2point;
import s2.s2pointutil;
import s2.s2polyline_simplifier;
import s2.s2testing;
import s2.s2text_format;

import fluent.asserts;

void checkSimplify(string src, string dst,
    string target, string avoid,
    in bool[] disc_on_left,
    double radius_degrees, bool expected_result) {
  auto radius = S1ChordAngle(S1Angle.fromDegrees(radius_degrees));
  auto s = new S2PolylineSimplifier();
  s.initialize(makePoint(src));
  foreach (S2Point p; parsePoints(target)) {
    s.targetDisc(p, radius);
  }
  int i = 0;
  foreach (S2Point p; parsePoints(avoid)) {
    s.avoidDisc(p, radius, disc_on_left[i++]);
  }
  Assert.equal(expected_result, s.extend(makePoint(dst)),
      "\nsrc = " ~ src ~ "\ndst = " ~ dst
      ~ "\ntarget = " ~ target ~ "\navoid = " ~ avoid);
}

@("S2PolylineSimplifier.Reuse") unittest {
  // Check that Init() can be called more than once.
  auto s = new S2PolylineSimplifier();
  auto radius = S1ChordAngle(S1Angle.fromDegrees(10.0));
  s.initialize(S2Point(1, 0, 0));
  Assert.equal(s.targetDisc(S2Point(1, 1, 0).normalize(), radius), true);
  Assert.equal(s.targetDisc(S2Point(1, 1, 0.1).normalize(), radius), true);
  Assert.equal(s.extend(S2Point(1, 1, 0.4).normalize()), false);

  // s.Init(S2Point(0, 1, 0));
  Assert.equal(s.targetDisc(S2Point(1, 1, 0.3).normalize(), radius), true);
  Assert.equal(s.targetDisc(S2Point(1, 1, 0.2).normalize(), radius), true);
  Assert.equal(s.extend(S2Point(1, 1, 0).normalize()), false);
}

@("S2PolylineSimplifier.NoConstraints") unittest {
  // No constraints, dst == src.
  checkSimplify("0:1", "0:1", "", "", [], 0, true);

  // No constraints, dst != src.
  checkSimplify("0:1", "1:0", "", "", [], 0, true);

  // No constraints, (src, dst) longer than 90 degrees (not supported).
  checkSimplify("0:0", "0:91", "", "", [], 0, false);
}

@("S2PolylineSimplifier.TargetOnePoint") unittest {
  // Three points on a straight line.  In theory zero tolerance should work,
  // but in practice there are floating point errors.
  checkSimplify("0:0", "0:2", "0:1", "", [], 1e-10, true);

  // Three points where the middle point is too far away.
  checkSimplify("0:0", "0:2", "1:1", "", [], 0.9, false);

  // A target disc that contains the source vertex.
  checkSimplify("0:0", "0:2", "0:0.1", "", [], 1.0, true);

  // A target disc that contains the destination vertex.
  checkSimplify("0:0", "0:2", "0:2.1", "", [], 1.0, true);
}

@("S2PolylineSimplifier.AvoidOnePoint") unittest {
  // Three points on a straight line, attempting to avoid the middle point.
  checkSimplify("0:0", "0:2", "", "0:1", [true], 1e-10, false);

  // Three points where the middle point can be successfully avoided.
  checkSimplify("0:0", "0:2", "", "1:1", [true], 0.9, true);

  // Three points where the middle point is on the left, but where the client
  // requires the point to be on the right of the edge.
  checkSimplify("0:0", "0:2", "", "1:1", [false], 1e-10, false);
}

@("S2PolylineSimplifier.TargetAndAvoid") unittest {
  // Target several points that are separated from the proposed edge by about
  // 0.7 degrees, and avoid several points that are separated from the
  // proposed edge by about 1.4 degrees.
  checkSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                "4:2, 7:5, 7:9", [true, true, false], 1.0, true);

  // The same example, but one point to be targeted is 1.4 degrees away.
  checkSimplify("0:0", "10:10", "2:3, 4:6, 7:8",
                "4:2, 7:5, 7:9", [true, true, false], 1.0, false);

  // The same example, but one point to be avoided is 0.7 degrees away.
  checkSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                "4:2, 6:5, 7:9", [true, true, false], 1.0, false);
}

@("S2PolylineSimplifier.Precision") unittest {
  // This is a rough upper bound on both the error in constructing the disc
  // locations (i.e., S2::InterpolateAtDistance, etc), and also on the
  // padding that S2PolylineSimplifier uses to ensure that its results are
  // conservative (i.e., the error calculated by GetSemiwidth).
  const S1Angle kMaxError = S1Angle.fromRadians(25 * double.epsilon);

  // We repeatedly generate a random edge.  We then target several discs that
  // barely overlap the edge, and avoid several discs that barely miss the
  // edge.  About half the time, we choose one disc and make it slightly too
  // large or too small so that targeting fails.
  const int kIters = 1000;  // Passes with 1 million iterations.
  auto simplifier = new S2PolylineSimplifier();
  for (int iter = 0; iter < kIters; ++iter) {
    S2Testing.rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point src = S2Testing.randomPoint();
    simplifier.initialize(src);
    S2Point dst = interpolateAtDistance(
        S1Angle.fromRadians(S2Testing.rnd.randDouble()),
        src, S2Testing.randomPoint());
    S2Point n = robustCrossProd(src, dst).normalize();

    // If bad_disc >= 0, then we make targeting fail for that disc.
    const int kNumDiscs = 5;
    int bad_disc = S2Testing.rnd.uniform(2 * kNumDiscs) - kNumDiscs;
    for (int i = 0; i < kNumDiscs; ++i) {
      double f = S2Testing.rnd.randDouble();
      S2Point a = ((1 - f) * src + f * dst).normalize();
      S1Angle r = S1Angle.fromRadians(S2Testing.rnd.randDouble());
      bool on_left = S2Testing.rnd.oneIn(2);
      S2Point x = interpolateAtDistance(r, a, on_left ? n : -n);
      // We grow the radius slightly if we want to target the disc and shrink
      // it otherwise, *unless* we want targeting to fail for this disc, in
      // which case these actions are reversed.
      bool avoid = S2Testing.rnd.oneIn(2);
      bool grow_radius = (avoid == (i == bad_disc));
      auto radius = S1ChordAngle(grow_radius ? r + kMaxError : r - kMaxError);
      if (avoid) {
        simplifier.avoidDisc(x, radius, on_left);
      } else {
        simplifier.targetDisc(x, radius);
      }
    }
    // The result is true iff all the disc constraints were satisfiable.
    Assert.equal(bad_disc < 0, simplifier.extend(dst));
  }
}
