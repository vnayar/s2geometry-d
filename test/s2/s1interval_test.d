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

module s2.s1interval_test;

import s2.s1interval;

import math = std.math;
import s2.util.math.s2const;
import s2.util.math.vector;
import std.conv;

import fluent.asserts;

// Create some standard intervals to use in the tests.  These include the
// empty and full intervals, intervals containing a single point, and
// intervals spanning one or more "quadrants" which are numbered as follows:
//    quad1 == [0, Pi/2]
//    quad2 == [Pi/2, Pi]
//    quad3 == [-Pi, -Pi/2]
//    quad4 == [-Pi/2, 0]
const S1Interval empty = S1Interval.empty();
const S1Interval full = S1Interval.full();
// Single-point intervals:
const S1Interval zero = S1Interval(0, 0);
const S1Interval pi2 = S1Interval(M_PI_2, M_PI_2);
const S1Interval pi = S1Interval(M_PI, M_PI);
const S1Interval mipi = S1Interval(-M_PI, -M_PI);  // Same as "pi" after normalization.
const S1Interval mipi2 = S1Interval(-M_PI_2, -M_PI_2);
// Single quadrants:
const S1Interval quad1 = S1Interval(0, M_PI_2);
const S1Interval quad2 = S1Interval(M_PI_2, -M_PI);
const S1Interval quad3 = S1Interval(M_PI, -M_PI_2);
const S1Interval quad4 = S1Interval(-M_PI_2, 0);
// Quadrant pairs:
const S1Interval quad12 = S1Interval(0, -M_PI);
const S1Interval quad23 = S1Interval(M_PI_2, -M_PI_2);
const S1Interval quad34 = S1Interval(-M_PI, 0);
const S1Interval quad41 = S1Interval(-M_PI_2, M_PI_2);
// Quadrant triples:
const S1Interval quad123 = S1Interval(0, -M_PI_2);
const S1Interval quad234 = S1Interval(M_PI_2, 0);
const S1Interval quad341 = S1Interval(M_PI, M_PI_2);
const S1Interval quad412 = S1Interval(-M_PI_2, -M_PI);
// Small intervals around the midpoints between quadrants, such that
// the center of each interval is offset slightly CCW from the midpoint.
const S1Interval mid12 = S1Interval(M_PI_2 - 0.01, M_PI_2 + 0.02);
const S1Interval mid23 = S1Interval(M_PI - 0.01, -M_PI + 0.02);
const S1Interval mid34 = S1Interval(-M_PI_2 - 0.01, -M_PI_2 + 0.02);
const S1Interval mid41 = S1Interval(-0.01, 0.02);

const double DOUBLE_ERR = 0.0001;
const float FLOAT_ERR = 0.0001;

@("ConstructorsAndAccessors")
unittest {
  // Spot-check the constructors and accessors.
  Assert.equal(quad12.lo(), 0);
  Assert.equal(quad12.hi(), M_PI);
  Assert.equal(quad34[0], M_PI);
  Assert.equal(quad34[1], 0);
  Assert.equal(quad34.bounds(), Vector2_d(M_PI, 0));
  Assert.equal(pi.lo(), M_PI);
  Assert.equal(pi.hi(), M_PI);

  // Check that [-Pi, -Pi] is normalized to [Pi, Pi].
  Assert.equal(mipi.lo(), M_PI);
  Assert.equal(mipi.hi(), M_PI);
  Assert.equal(quad23.lo(), M_PI_2);
  Assert.equal(quad23.hi(), -M_PI_2);

  // Check that the default S1Interval is identical to Empty().
  S1Interval default_empty;
  Assert.equal(default_empty.isValid(), true);
  Assert.equal(default_empty.isEmpty(), true);
  Assert.equal(empty.lo(), default_empty.lo());
  Assert.equal(empty.hi(), default_empty.hi());
}

@("SimplePredicates")
unittest {
  // is_valid(), is_empty(), is_full(), is_inverted()
  Assert.equal(zero.isValid() && !zero.isEmpty() && !zero.isFull(), true);
  Assert.equal(empty.isValid() && empty.isEmpty() && !empty.isFull(), true);
  Assert.equal(empty.isInverted(), true);
  Assert.equal(full.isValid() && !full.isEmpty() && full.isFull(), true);
  Assert.equal(!quad12.isEmpty() && !quad12.isFull() && !quad12.isInverted(), true);
  Assert.equal(!quad23.isEmpty() && !quad23.isFull() && quad23.isInverted(), true);
  Assert.equal(pi.isValid() && !pi.isEmpty() && !pi.isInverted(), true);
  Assert.equal(mipi.isValid() && !mipi.isEmpty() && !mipi.isInverted(), true);
}

@("AlmostEmptyOrFull")
unittest {
  // Test that rounding errors don't cause intervals that are almost empty or
  // full to be considered empty or full.  The following value is the greatest
  // representable value less than Pi.
  const double kAlmostPi = M_PI - 2 * double.epsilon;
  Assert.equal(S1Interval(-kAlmostPi, M_PI).isFull(), false);
  Assert.equal(S1Interval(-M_PI, kAlmostPi).isFull(), false);
  Assert.equal(S1Interval(M_PI, -kAlmostPi).isEmpty(), false);
  Assert.equal(S1Interval(kAlmostPi, -M_PI).isEmpty(), false);
}

@("GetCenter")
unittest {
  Assert.equal(quad12.getCenter(), M_PI_2);
  Assert.approximately(S1Interval(3.1, 2.9).getCenter(), 3.0 - M_PI, DOUBLE_ERR);
  Assert.approximately(S1Interval(-2.9, -3.1).getCenter(), M_PI - 3.0, DOUBLE_ERR);
  Assert.approximately(S1Interval(2.1, -2.1).getCenter(), M_PI, DOUBLE_ERR);
  Assert.equal(pi.getCenter(), M_PI);
  Assert.equal(mipi.getCenter(), M_PI);
  Assert.equal(math.fabs(quad23.getCenter()), M_PI);
  Assert.approximately(quad123.getCenter(), 0.75 * M_PI, DOUBLE_ERR);
}

@("GetLength")
unittest {
  Assert.equal(quad12.getLength(), M_PI);
  Assert.equal(pi.getLength(), 0);
  Assert.equal(mipi.getLength(), 0);
  Assert.approximately(quad123.getLength(), 1.5 * M_PI, DOUBLE_ERR);
  Assert.equal(math.fabs(quad23.getLength()), M_PI);
  Assert.equal(full.getLength(), 2 * M_PI);
  Assert.lessThan(empty.getLength(), 0);
}

@("Complement")
unittest {
  Assert.equal(empty.complement().isFull(), true);
  Assert.equal(full.complement().isEmpty(), true);
  Assert.equal(pi.complement().isFull(), true);
  Assert.equal(mipi.complement().isFull(), true);
  Assert.equal(zero.complement().isFull(), true);
  Assert.equal(quad12.complement().approxEquals(quad34), true);
  Assert.equal(quad34.complement().approxEquals(quad12), true);
  Assert.equal(quad123.complement().approxEquals(quad4), true);
}

@("Contains")
unittest {
  // Contains(double), InteriorContains(double)
  Assert.equal(
      !empty.contains(0) && !empty.contains(M_PI) && !empty.contains(-M_PI), true);
  Assert.equal(!empty.interiorContains(M_PI) && !empty.interiorContains(-M_PI), true);
  Assert.equal(full.contains(0) && full.contains(M_PI) && full.contains(-M_PI), true);
  Assert.equal(full.interiorContains(M_PI) && full.interiorContains(-M_PI), true);
  Assert.equal(
      quad12.contains(0) && quad12.contains(M_PI) && quad12.contains(-M_PI), true);
  Assert.equal(quad12.interiorContains(M_PI_2) && !quad12.interiorContains(0), true);
  Assert.equal(
      !quad12.interiorContains(M_PI) && !quad12.interiorContains(-M_PI), true);
  Assert.equal(quad23.contains(M_PI_2) && quad23.contains(-M_PI_2), true);
  Assert.equal(quad23.contains(M_PI) && quad23.contains(-M_PI), true);
  Assert.equal(!quad23.contains(0), true);
  Assert.equal(
      !quad23.interiorContains(M_PI_2) && !quad23.interiorContains(-M_PI_2), true);
  Assert.equal(quad23.interiorContains(M_PI) && quad23.interiorContains(-M_PI), true);
  Assert.equal(!quad23.interiorContains(0), true);
  Assert.equal(pi.contains(M_PI) && pi.contains(-M_PI) && !pi.contains(0), true);
  Assert.equal(!pi.interiorContains(M_PI) && !pi.interiorContains(-M_PI), true);
  Assert.equal(mipi.contains(M_PI) && mipi.contains(-M_PI) && !mipi.contains(0), true);
  Assert.equal(!mipi.interiorContains(M_PI) && !mipi.interiorContains(-M_PI), true);
  Assert.equal(zero.contains(0) && !zero.interiorContains(0), true);
}

private void testIntervalOps(
    in S1Interval x, in S1Interval y,
    in string expected_relation,
    in S1Interval expected_union,
    in S1Interval expected_intersection) {
  // Test all of the interval operations on the given pair of intervals.
  // "expected_relation" is a sequence of "T" and "F" characters corresponding
  // to the expected results of Contains(), InteriorContains(), Intersects(),
  // and InteriorIntersects() respectively.

  Assert.equal(x.contains(y), expected_relation[0] == 'T');
  Assert.equal(x.interiorContains(y), expected_relation[1] == 'T');
  Assert.equal(x.intersects(y), expected_relation[2] == 'T');
  Assert.equal(x.interiorIntersects(y), expected_relation[3] == 'T');

  // bounds() returns a const reference to a member variable, so we need to
  // make a copy when invoking it on a temporary object.
  Assert.equal(Vector2_d(x.unite(y).bounds()), expected_union.bounds());
  Assert.equal(Vector2_d(x.intersection(y).bounds()), expected_intersection.bounds());

  Assert.equal(x.contains(y), x.unite(y) == x);
  Assert.equal(x.intersects(y), !x.intersection(y).isEmpty());

  if (y.lo() == y.hi()) {
    S1Interval r = x;
    r.addPoint(y.lo());
    Assert.equal(r.bounds(), expected_union.bounds());
  }
}

@("IntervalOps")
unittest {
  // Contains(S1Interval), InteriorContains(S1Interval),
  // Intersects(), InteriorIntersects(), Union(), Intersection()
  testIntervalOps(empty, empty, "TTFF", empty, empty);
  testIntervalOps(empty, full, "FFFF", full, empty);
  testIntervalOps(empty, zero, "FFFF", zero, empty);
  testIntervalOps(empty, pi, "FFFF", pi, empty);
  testIntervalOps(empty, mipi, "FFFF", mipi, empty);

  testIntervalOps(full, empty, "TTFF", full, empty);
  testIntervalOps(full, full, "TTTT", full, full);
  testIntervalOps(full, zero, "TTTT", full, zero);
  testIntervalOps(full, pi, "TTTT", full, pi);
  testIntervalOps(full, mipi, "TTTT", full, mipi);
  testIntervalOps(full, quad12, "TTTT", full, quad12);
  testIntervalOps(full, quad23, "TTTT", full, quad23);

  testIntervalOps(zero, empty, "TTFF", zero, empty);
  testIntervalOps(zero, full, "FFTF", full, zero);
  testIntervalOps(zero, zero, "TFTF", zero, zero);
  testIntervalOps(zero, pi, "FFFF", S1Interval(0, M_PI), empty);
  testIntervalOps(zero, pi2, "FFFF", quad1, empty);
  testIntervalOps(zero, mipi, "FFFF", quad12, empty);
  testIntervalOps(zero, mipi2, "FFFF", quad4, empty);
  testIntervalOps(zero, quad12, "FFTF", quad12, zero);
  testIntervalOps(zero, quad23, "FFFF", quad123, empty);

  testIntervalOps(pi2, empty, "TTFF", pi2, empty);
  testIntervalOps(pi2, full, "FFTF", full, pi2);
  testIntervalOps(pi2, zero, "FFFF", quad1, empty);
  testIntervalOps(pi2, pi, "FFFF", S1Interval(M_PI_2, M_PI), empty);
  testIntervalOps(pi2, pi2, "TFTF", pi2, pi2);
  testIntervalOps(pi2, mipi, "FFFF", quad2, empty);
  testIntervalOps(pi2, mipi2, "FFFF", quad23, empty);
  testIntervalOps(pi2, quad12, "FFTF", quad12, pi2);
  testIntervalOps(pi2, quad23, "FFTF", quad23, pi2);

  testIntervalOps(pi, empty, "TTFF", pi, empty);
  testIntervalOps(pi, full, "FFTF", full, pi);
  testIntervalOps(pi, zero, "FFFF", S1Interval(M_PI, 0), empty);
  testIntervalOps(pi, pi, "TFTF", pi, pi);
  testIntervalOps(pi, pi2, "FFFF", S1Interval(M_PI_2, M_PI), empty);
  testIntervalOps(pi, mipi, "TFTF", pi, pi);
  testIntervalOps(pi, mipi2, "FFFF", quad3, empty);
  testIntervalOps(pi, quad12, "FFTF", S1Interval(0, M_PI), pi);
  testIntervalOps(pi, quad23, "FFTF", quad23, pi);

  testIntervalOps(mipi, empty, "TTFF", mipi, empty);
  testIntervalOps(mipi, full, "FFTF", full, mipi);
  testIntervalOps(mipi, zero, "FFFF", quad34, empty);
  testIntervalOps(mipi, pi, "TFTF", mipi, mipi);
  testIntervalOps(mipi, pi2, "FFFF", quad2, empty);
  testIntervalOps(mipi, mipi, "TFTF", mipi, mipi);
  testIntervalOps(mipi, mipi2, "FFFF", S1Interval(-M_PI, -M_PI_2), empty);
  testIntervalOps(mipi, quad12, "FFTF", quad12, mipi);
  testIntervalOps(mipi, quad23, "FFTF", quad23, mipi);

  testIntervalOps(quad12, empty, "TTFF", quad12, empty);
  testIntervalOps(quad12, full, "FFTT", full, quad12);
  testIntervalOps(quad12, zero, "TFTF", quad12, zero);
  testIntervalOps(quad12, pi, "TFTF", quad12, pi);
  testIntervalOps(quad12, mipi, "TFTF", quad12, mipi);
  testIntervalOps(quad12, quad12, "TFTT", quad12, quad12);
  testIntervalOps(quad12, quad23, "FFTT", quad123, quad2);
  testIntervalOps(quad12, quad34, "FFTF", full, quad12);

  testIntervalOps(quad23, empty, "TTFF", quad23, empty);
  testIntervalOps(quad23, full, "FFTT", full, quad23);
  testIntervalOps(quad23, zero, "FFFF", quad234, empty);
  testIntervalOps(quad23, pi, "TTTT", quad23, pi);
  testIntervalOps(quad23, mipi, "TTTT", quad23, mipi);
  testIntervalOps(quad23, quad12, "FFTT", quad123, quad2);
  testIntervalOps(quad23, quad23, "TFTT", quad23, quad23);
  testIntervalOps(quad23, quad34, "FFTT", quad234, S1Interval(-M_PI, -M_PI_2));

  testIntervalOps(quad1, quad23, "FFTF", quad123, S1Interval(M_PI_2, M_PI_2));
  testIntervalOps(quad2, quad3, "FFTF", quad23, mipi);
  testIntervalOps(quad3, quad2, "FFTF", quad23, pi);
  testIntervalOps(quad2, pi, "TFTF", quad2, pi);
  testIntervalOps(quad2, mipi, "TFTF", quad2, mipi);
  testIntervalOps(quad3, pi, "TFTF", quad3, pi);
  testIntervalOps(quad3, mipi, "TFTF", quad3, mipi);

  testIntervalOps(quad12, mid12, "TTTT", quad12, mid12);
  testIntervalOps(mid12, quad12, "FFTT", quad12, mid12);

  S1Interval quad12eps = S1Interval(quad12.lo(), mid23.hi());
  S1Interval quad2hi = S1Interval(mid23.lo(), quad12.hi());
  testIntervalOps(quad12, mid23, "FFTT", quad12eps, quad2hi);
  testIntervalOps(mid23, quad12, "FFTT", quad12eps, quad2hi);

  // This test checks that the union of two disjoint intervals is the smallest
  // interval that contains both of them.  Note that the center of "mid34"
  // slightly CCW of -Pi/2 so that there is no ambiguity about the result.
  S1Interval quad412eps = S1Interval(mid34.lo(), quad12.hi());
  testIntervalOps(quad12, mid34, "FFFF", quad412eps, empty);
  testIntervalOps(mid34, quad12, "FFFF", quad412eps, empty);

  S1Interval quadeps12 = S1Interval(mid41.lo(), quad12.hi());
  S1Interval quad1lo = S1Interval(quad12.lo(), mid41.hi());
  testIntervalOps(quad12, mid41, "FFTT", quadeps12, quad1lo);
  testIntervalOps(mid41, quad12, "FFTT", quadeps12, quad1lo);

  S1Interval quad2lo = S1Interval(quad23.lo(), mid12.hi());
  S1Interval quad3hi = S1Interval(mid34.lo(), quad23.hi());
  S1Interval quadeps23 = S1Interval(mid12.lo(), quad23.hi());
  S1Interval quad23eps = S1Interval(quad23.lo(), mid34.hi());
  S1Interval quadeps123 = S1Interval(mid41.lo(), quad23.hi());
  testIntervalOps(quad23, mid12, "FFTT", quadeps23, quad2lo);
  testIntervalOps(mid12, quad23, "FFTT", quadeps23, quad2lo);
  testIntervalOps(quad23, mid23, "TTTT", quad23, mid23);
  testIntervalOps(mid23, quad23, "FFTT", quad23, mid23);
  testIntervalOps(quad23, mid34, "FFTT", quad23eps, quad3hi);
  testIntervalOps(mid34, quad23, "FFTT", quad23eps, quad3hi);
  testIntervalOps(quad23, mid41, "FFFF", quadeps123, empty);
  testIntervalOps(mid41, quad23, "FFFF", quadeps123, empty);
}

@("AddPoint")
unittest {
  S1Interval r = empty; r.addPoint(0);
  Assert.equal(r, zero);
  r = empty; r.addPoint(M_PI);
  Assert.equal(r, pi);
  r = empty; r.addPoint(-M_PI);
  Assert.equal(r, mipi);
  r = empty; r.addPoint(M_PI); r.addPoint(-M_PI);
  Assert.equal(r, pi);
  r = empty; r.addPoint(-M_PI); r.addPoint(M_PI);
  Assert.equal(r, mipi);
  r = empty; r.addPoint(mid12.lo()); r.addPoint(mid12.hi());
  Assert.equal(r, mid12);
  r = empty; r.addPoint(mid23.lo()); r.addPoint(mid23.hi());
  Assert.equal(r, mid23);
  r = quad1; r.addPoint(-0.9*M_PI); r.addPoint(-M_PI_2);
  Assert.equal(r, quad123);
  r = full; r.addPoint(0);
  Assert.equal(r.isFull(), true);
  r = full; r.addPoint(M_PI);
  Assert.equal(r.isFull(), true);
  r = full; r.addPoint(-M_PI);
  Assert.equal(r.isFull(), true);
}

@("Project")
unittest {
  S1Interval r = S1Interval(-M_PI, -M_PI);
  Assert.equal(M_PI, r.project(-M_PI));
  Assert.equal(M_PI, r.project(0));
  r = S1Interval(0, M_PI);
  Assert.equal(0.1, r.project(0.1));
  Assert.equal(0.0, r.project(-M_PI_2 + 1e-15));
  Assert.equal(M_PI, r.project(-M_PI_2 - 1e-15));
  r = S1Interval(M_PI - 0.1, -M_PI + 0.1);
  Assert.equal(M_PI, r.project(M_PI));
  Assert.equal(M_PI - 0.1, r.project(1e-15));
  Assert.equal(-M_PI + 0.1, r.project(-1e-15));
  Assert.equal(0.0, S1Interval.full().project(0));
  Assert.equal(M_PI, S1Interval.full().project(M_PI));
  Assert.equal(M_PI, S1Interval.full().project(-M_PI));
}

@("FromPointPair")
unittest {
  Assert.equal(S1Interval.fromPointPair(-M_PI, M_PI), pi);
  Assert.equal(S1Interval.fromPointPair(M_PI, -M_PI), pi);
  Assert.equal(S1Interval.fromPointPair(mid34.hi(), mid34.lo()), mid34);
  Assert.equal(S1Interval.fromPointPair(mid23.lo(), mid23.hi()), mid23);
}

@("Expanded")
unittest {
  Assert.equal(empty.expanded(1), empty);
  Assert.equal(full.expanded(1), full);
  Assert.equal(zero.expanded(1), S1Interval(-1, 1));
  Assert.equal(mipi.expanded(0.01), S1Interval(M_PI - 0.01, -M_PI + 0.01));
  Assert.equal(pi.expanded(27), full);
  Assert.equal(pi.expanded(M_PI_2).approxEquals(quad23), true);
  Assert.equal(pi2.expanded(M_PI_2), quad12);
  Assert.equal(mipi2.expanded(M_PI_2), quad34);

  Assert.equal(empty.expanded(-1), empty);
  Assert.equal(full.expanded(-1), full);
  Assert.equal(quad123.expanded(-27), empty);
  Assert.equal(quad234.expanded(-27), empty);
  Assert.equal(quad123.expanded(-M_PI_2), quad2);
  Assert.equal(quad341.expanded(-M_PI_2).approxEquals(quad4), true);
  Assert.equal(quad412.expanded(-M_PI_2), quad1);
}

@("ApproxEquals")
unittest {
  // Choose two values kLo and kHi such that it's okay to shift an endpoint by
  // kLo (i.e., the resulting interval is equivalent) but not by kHi.
  static const double kLo = 3 * double.epsilon;  // < max_error default
  static const double kHi = 6 * double.epsilon;  // > max_error default

  // Empty intervals.
  Assert.equal(empty.approxEquals(empty), true);
  Assert.equal(zero.approxEquals(empty) && empty.approxEquals(zero), true);
  Assert.equal(pi.approxEquals(empty) && empty.approxEquals(pi), true);
  Assert.equal(mipi.approxEquals(empty) && empty.approxEquals(mipi), true);
  Assert.equal(empty.approxEquals(full), false);
  Assert.equal(empty.approxEquals(S1Interval(1, 1 + 2*kLo)), true);
  Assert.equal(empty.approxEquals(S1Interval(1, 1 + 2*kHi)), false);
  Assert.equal(S1Interval(M_PI - kLo, -M_PI + kLo).approxEquals(empty), true);

  // Full intervals.
  Assert.equal(full.approxEquals(full), true);
  Assert.equal(full.approxEquals(empty), false);
  Assert.equal(full.approxEquals(zero), false);
  Assert.equal(full.approxEquals(pi), false);
  Assert.equal(full.approxEquals(S1Interval(kLo, -kLo)), true);
  Assert.equal(full.approxEquals(S1Interval(2*kHi, 0)), false);
  Assert.equal(S1Interval(-M_PI + kLo, M_PI - kLo).approxEquals(full), true);
  Assert.equal(S1Interval(-M_PI, M_PI - 2*kHi).approxEquals(full), false);

  // Singleton intervals.
  Assert.equal(pi.approxEquals(pi) && mipi.approxEquals(pi), true);
  Assert.equal(pi.approxEquals(S1Interval(M_PI - kLo, M_PI - kLo)), true);
  Assert.equal(pi.approxEquals(S1Interval(M_PI - kHi, M_PI - kHi)), false);
  Assert.equal(pi.approxEquals(S1Interval(M_PI - kLo, -M_PI + kLo)), true);
  Assert.equal(pi.approxEquals(S1Interval(M_PI - kHi, -M_PI)), false);
  Assert.equal(zero.approxEquals(pi), false);
  Assert.equal(pi.unite(mid12).unite(zero).approxEquals(quad12), true);
  Assert.equal(quad2.intersection(quad3).approxEquals(pi), true);
  Assert.equal(quad3.intersection(quad2).approxEquals(pi), true);

  // Intervals whose corresponding endpoints are nearly the same but where the
  // endpoints are in opposite order (i.e., inverted intervals).
  Assert.equal(S1Interval(0, kLo).approxEquals(S1Interval(kLo, 0)), false);
  Assert.equal(S1Interval(M_PI - 0.5 * kLo, -M_PI + 0.5 * kLo)
      .approxEquals(S1Interval(-M_PI + 0.5 * kLo, M_PI - 0.5 * kLo)), false);

  // Other intervals.
  Assert.equal(S1Interval(1 - kLo, 2 + kLo).approxEquals(S1Interval(1, 2)), true);
  Assert.equal(S1Interval(1 + kLo, 2 - kLo).approxEquals(S1Interval(1, 2)), true);
  Assert.equal(S1Interval(2 - kLo, 1 + kLo).approxEquals(S1Interval(2, 1)), true);
  Assert.equal(S1Interval(2 + kLo, 1 - kLo).approxEquals(S1Interval(2, 1)), true);
  Assert.equal(S1Interval(1 - kHi, 2 + kLo).approxEquals(S1Interval(1, 2)), false);
  Assert.equal(S1Interval(1 + kHi, 2 - kLo).approxEquals(S1Interval(1, 2)), false);
  Assert.equal(S1Interval(2 - kHi, 1 + kLo).approxEquals(S1Interval(2, 1)), false);
  Assert.equal(S1Interval(2 + kHi, 1 - kLo).approxEquals(S1Interval(2, 1)), false);
  Assert.equal(S1Interval(1 - kLo, 2 + kHi).approxEquals(S1Interval(1, 2)), false);
  Assert.equal(S1Interval(1 + kLo, 2 - kHi).approxEquals(S1Interval(1, 2)), false);
  Assert.equal(S1Interval(2 - kLo, 1 + kHi).approxEquals(S1Interval(2, 1)), false);
  Assert.equal(S1Interval(2 + kLo, 1 - kHi).approxEquals(S1Interval(2, 1)), false);
}

@("GetDirectedHausdorffDistance")
unittest {
  Assert.approximately(0.0, empty.getDirectedHausdorffDistance(empty), FLOAT_ERR);
  Assert.approximately(0.0, empty.getDirectedHausdorffDistance(mid12), FLOAT_ERR);
  Assert.approximately(M_PI, mid12.getDirectedHausdorffDistance(empty), FLOAT_ERR);

  Assert.equal(0.0, quad12.getDirectedHausdorffDistance(quad123));
  S1Interval interval = S1Interval(3.0, -3.0);  // an interval whose complement center is 0.
  Assert.approximately(
      3.0, S1Interval(-0.1, 0.2).getDirectedHausdorffDistance(interval), FLOAT_ERR);
  Assert.approximately(
      3.0 - 0.1, S1Interval(0.1, 0.2).getDirectedHausdorffDistance(interval), FLOAT_ERR);
  Assert.approximately(
      3.0 - 0.1, S1Interval(-0.2, -0.1).getDirectedHausdorffDistance(interval), FLOAT_ERR);
}
