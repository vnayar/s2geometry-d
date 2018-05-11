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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2predicates_test;

import fluent.asserts;
import s2pred = s2.s2predicates;
import s2.s2point;
import s2.s1angle;
import s2.s1chordangle;
import s2.s2edgedistances;
import s2.s2pointutil;
import s2.s2testing;
import s2.util.math.exactfloat;
import s2.util.math.vector;
import std.exception;
import algorithm = std.algorithm;
import array = std.array;
import format = std.format;
import math = std.math;
import range = std.range;
import traits = std.traits;

import std.stdio;

// Number of iterations for precision consistency tests.
immutable int CONSISTENCY_ITERS = 5000;

@("rounding_epsilon")
unittest {
  // Check that rounding_epsilon<T>() returns the expected value for "float"
  // and "double".  We explicitly do not test "long double" since if this type
  // is implemented using double-double arithmetic then the numeric_limits
  // epsilon() value is completely unrelated to the maximum rounding error.
  Assert.equal(0.5 * float.epsilon, s2pred.roundingEpsilon!float());
  Assert.equal(0.5 * double.epsilon, s2pred.roundingEpsilon!double());
}

@("sign - ColinnearPoints")
unittest {
  // The following points happen to be *exactly collinear* along a line that it
  // approximate tangent to the surface of the unit sphere.  In fact, C is the
  // exact midpoint of the line segment AB.  All of these points are close
  // enough to unit length to satisfy S2::IsUnitLength().
  S2Point a = S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
  S2Point b = S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
  S2Point c = S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
  Assert.equal(c - a, b - c);
  Assert.notEqual(0, s2pred.sign(a, b, c));
  Assert.equal(s2pred.sign(a, b, c), s2pred.sign(b, c, a));
  Assert.equal(s2pred.sign(a, b, c), -s2pred.sign(c, b, a));

  // The points "x1" and "x2" are exactly proportional, i.e. they both lie
  // on a common line through the origin.  Both points are considered to be
  // normalized, and in fact they both satisfy (x == x.Normalize()).
  // Therefore the triangle (x1, x2, -x1) consists of three distinct points
  // that all lie on a common line through the origin.
  S2Point x1 = S2Point(0.99999999999999989, 1.4901161193847655e-08, 0);
  S2Point x2 = S2Point(1, 1.4901161193847656e-08, 0);
  Assert.equal(x1, x1.normalize());
  Assert.equal(x2, x2.normalize());
  Assert.notEqual(0, s2pred.sign(x1, x2, -x1));
  Assert.equal(s2pred.sign(x1, x2, -x1), s2pred.sign(x2, -x1, x1));
  Assert.equal(s2pred.sign(x1, x2, -x1), -s2pred.sign(-x1, x2, x1));

  // Here are two more points that are distinct, exactly proportional, and
  // that satisfy (x == x.Normalize()).
  S2Point x3 = S2Point(1, 1, 1).normalize();
  S2Point x4 = 0.99999999999999989 * x3;
  Assert.equal(x3, x3.normalize());
  Assert.equal(x4, x4.normalize());
  Assert.notEqual(x3, x4);
  Assert.notEqual(0, s2pred.sign(x3, x4, -x3));

  // TODO(vnayar): Why is it important for normalize() to not be idempotent?
  //
  // The following two points demonstrate that Normalize() is not idempotent,
  // i.e. y0.Normalize() != y0.Normalize().Normalize().  Both points satisfy
  // S2::IsNormalized(), though, and the two points are exactly proportional.
  // S2Point y0 = S2Point(1, 1, 0);
  // S2Point y1 = y0.normalize();
  // S2Point y2 = y1.normalize();
  // Assert.notEqual(y1, y2);
  // Assert.equal(y2, y2.normalize());
  // Assert.notEqual(0, sign(y1, y2, -y1));
  // Assert.equal(sign(y1, y2, -y1), sign(y2, -y1, y1));
  // Assert.equal(sign(y1, y2, -y1), -sign(-y1, y2, y1));
}

// This test repeatedly constructs some number of points that are on or nearly
// on a given great circle.  Then it chooses one of these points as the
// "origin" and sorts the other points in CCW order around it.  Of course,
// since the origin is on the same great circle as the points being sorted,
// nearly all of these tests are degenerate.  It then does various consistency
// checks to verify that the points are indeed sorted in CCW order.
//
// It is easier to think about what this test is doing if you imagine that the
// points are in general position rather than on a great circle.
class SignTest {
protected:
  // The following method is used to sort a collection of points in CCW order
  // around a given origin.  It returns true if A comes before B in the CCW
  // ordering (starting at an arbitrary fixed direction).
  static class LessCCW {
  public:
    this(in S2Point origin, in S2Point start) {
      _origin = origin;
      _start = start;
    }
    bool compare()(in S2Point a, in S2Point b) const {
      // OrderedCCW() acts like "<=", so we need to invert the comparison.
      return !s2pred.orderedCCW(_start, b, a, _origin);
    }
  private:
    const S2Point _origin;
    const S2Point _start;
  }

  // Given a set of points with no duplicates, first remove "origin" from
  // "points" (if it exists) and then sort the remaining points in CCW order
  // around "origin" putting the result in "sorted".
  static void sortCCW(in S2Point[] points, in S2Point origin, out S2Point[] sorted) {
    // Make a copy of the points with "origin" removed.
    array.appender(&sorted) ~=
        range.retro(algorithm.filterBidirectional!(a => a != origin)(points));

    // Sort the points CCW around the origin starting at (*sorted)[0].
    LessCCW less = new LessCCW(origin, sorted[0]);
    algorithm.sort!((a, b) => less.compare(a, b))(sorted);
  }

  // Given a set of points sorted circularly CCW around "origin", and the
  // index "start" of a point A, count the number of CCW triangles OAB over
  // all sorted points B not equal to A.  Also check that the results of the
  // CCW tests are consistent with the hypothesis that the points are sorted.
  static int countCCW(in S2Point[] sorted, in S2Point origin, size_t start) {
    int num_ccw = 0;
    int last_sign = 1;
    const size_t n = sorted.length;
    for (int j = 1; j < n; ++j) {
      int sign = s2pred.sign(origin, sorted[start], sorted[(start + j) % n]);
      Assert.notEqual(0, sign);
      if (sign > 0) ++num_ccw;

      // Since the points are sorted around the origin, we expect to see a
      // (possibly empty) sequence of CCW triangles followed by a (possibly
      // empty) sequence of CW triangles.
      Assert.equal(false, sign > 0 && last_sign < 0);
      last_sign = sign;
    }
    return num_ccw;
  }

  // Test exhaustively whether the points in "sorted" are sorted circularly
  // CCW around "origin".
  static void testCCW(in S2Point[] sorted, in S2Point origin) {
    const int n = cast(int) sorted.length;
    int total_num_ccw = 0;
    int last_num_ccw = countCCW(sorted, origin, n - 1);
    for (int start = 0; start < n; ++start) {
      int num_ccw = countCCW(sorted, origin, start);
      // Each iteration we increase the start index by 1, therefore the number
      // of CCW triangles should decrease by at most 1.
      Assert.notLessThan(num_ccw, last_num_ccw - 1);
      total_num_ccw += num_ccw;
      last_num_ccw = num_ccw;
    }
    // We have tested all triangles of the form OAB.  Exactly half of these
    // should be CCW.
    Assert.equal(n * (n-1) / 2, total_num_ccw);
  }

  static void addNormalized(S2Point a, ref S2Point[] points) {
    points ~= a.normalize();
  }

  // Add two points A1 and A2 that are slightly offset from A along the
  // tangent toward B, and such that A, A1, and A2 are exactly collinear
  // (i.e. even with infinite-precision arithmetic).
  static void addTangentPoints(in S2Point a, in S2Point b, ref S2Point[] points) {
    Vector3_d dir = robustCrossProd(a, b).crossProd(a).normalize();
    if (dir == S2Point(0, 0, 0)) return;
    for (;;) {
      S2Point delta = 1e-15 * S2Testing.rnd.randDouble() * dir;
      if ((a + delta) != a && (a + delta) - a == a - (a - delta) &&
          isUnitLength(a + delta) && isUnitLength(a - delta)) {
        points ~= a + delta;
        points ~= a - delta;
        return;
      }
    }
  }

  // Add zero or more (but usually one) point that is likely to trigger
  // Sign() degeneracies among the given points.
  static void addDegeneracy(ref S2Point[] points) {
    Random rnd = S2Testing.rnd;
    S2Point a = points[rnd.uniform(cast(int) points.length)];
    S2Point b = points[rnd.uniform(cast(int) points.length)];
    int coord = rnd.uniform(3);
    final switch (rnd.uniform(8)) {
      case 0:
        // Add a random point (not uniformly distributed) along the great
        // circle AB.
        addNormalized(rnd.uniformDouble(-1, 1) * a + rnd.uniformDouble(-1, 1) * b, points);
        break;
      case 1:
        // Perturb one coordinate by the minimum amount possible.
        a.data[coord] = math.nextafter(a[coord], rnd.oneIn(2) ? 2 : -2);
        addNormalized(a, points);
        break;
      case 2:
        // Perturb one coordinate by up to 1e-15.
        a.data[coord] += 1e-15 * rnd.uniformDouble(-1, 1);
        addNormalized(a, points);
        break;
      case 3:
        // Scale a point just enough so that it is different while still being
        // considered normalized.
        a *= rnd.oneIn(2) ? (1 + 2e-16) : (1 - 1e-16);
        if (isUnitLength(a)) points ~= a;
        break;
      case 4: {
        // Add the intersection point of AB with X=0, Y=0, or Z=0.
        S2Point dir = S2Point(0, 0, 0);
        dir[coord] = rnd.oneIn(2) ? 1 : -1;
        Vector3_d norm = robustCrossProd(a, b).normalize();
        if (norm.norm2() > 0) {
          addNormalized(robustCrossProd(dir, norm), points);
        }
        break;
      }
      case 5:
        // Add two closely spaced points along the tangent at A to the great
        // circle through AB.
        addTangentPoints(a, b, points);
        break;
      case 6:
        // Add two closely spaced points along the tangent at A to the great
        // circle through A and the X-axis.
        addTangentPoints(a, S2Point(1, 0, 0), points);
        break;
      case 7:
        // Add the negative of a point.
        points ~= -a;
        break;
    }
  }

  // Sort the points around the given origin, and then do some consistency
  // checks to verify that they are actually sorted.
  static void sortAndTest(in S2Point[] points, in S2Point origin) {
    S2Point[] sorted;
    sortCCW(points, origin, sorted);
    testCCW(sorted, origin);
  }

  // Construct approximately "n" points near the great circle through A and B,
  // then sort them and test whether they are sorted.
  static void testGreatCircle(S2Point a, S2Point b, int n) {
    a = a.normalize();
    b = b.normalize();
    S2Point[] points;
    points ~= a;
    points ~= b;
    while (points.length < n) {
      addDegeneracy(points);
    }
    // Remove any (0, 0, 0) points that were accidentically created, then sort
    // the points and remove duplicates.
    points = array.array(
        algorithm.uniq(
            algorithm.sort(
                array.array(
                    algorithm.filter!(point => point != S2Point(0, 0, 0))(points)))));

    // enforce(points.length >= n / 2);

    sortAndTest(points, a);
    sortAndTest(points, b);
    for (int k = 0; k < points.length; ++k) {
      sortAndTest(points, points[k]);
    }
  }
}

@("SignTest.StressTest")
unittest {
  // The run time of this test is *cubic* in the parameter below.
  static const int kNumPointsPerCircle = 20;
  SignTest signTest = new SignTest();
  // This test is randomized, so it is beneficial to run it several times.
  for (int iter = 0; iter < 3; ++iter) {
    // The most difficult great circles are the ones in the X-Y, Y-Z, and X-Z
    // planes, for two reasons.  First, when one or more coordinates are close
    // to zero then the perturbations can be much smaller, since floating
    // point numbers are spaced much more closely together near zero.  (This
    // tests the handling of things like underflow.)  The second reason is
    // that most of the cases of SymbolicallyPerturbedSign() can only be
    // reached when one or more input point coordinates are zero.
    signTest.testGreatCircle(S2Point(1, 0, 0), S2Point(0, 1, 0), kNumPointsPerCircle);
    signTest.testGreatCircle(S2Point(1, 0, 0), S2Point(0, 0, 1), kNumPointsPerCircle);
    signTest.testGreatCircle(S2Point(0, -1, 0), S2Point(0, 0, 1), kNumPointsPerCircle);

    // This tests a great circle where at least some points have X, Y, and Z
    // coordinates with exactly the same mantissa.  One useful property of
    // such points is that when they are scaled (e.g. multiplying by 1+eps),
    // all such points are exactly collinear with the origin.
    signTest.testGreatCircle(S2Point(1 << 25, 1, -8), S2Point(-4, -(1 << 20), 1),
                    kNumPointsPerCircle);
  }
}

class StableSignTest {
protected:
  // Estimate the probability that S2::StableSign() will not be able to compute
  // the determinant sign of a triangle A, B, C consisting of three points
  // that are as collinear as possible and spaced the given distance apart.
  double getFailureRate(double km) {
    const int kIters = 1000;
    int failure_count = 0;
    double m = math.tan(S2Testing.kmToAngle(km).radians());
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a, x, y;
      S2Testing.getRandomFrame(a, x, y);
      S2Point b = (a - m * x).normalize();
      S2Point c = (a + m * x).normalize();
      int sign = s2pred.stableSign(a, b, c);
      if (sign != 0) {
        Assert.equal(s2pred.exactSign(a, b, c, true), sign);
      } else {
        ++failure_count;
      }
    }
    double rate = cast(double) failure_count / kIters;
    writeln("StableSign failure rate for ", km, " km = ", rate);
    return rate;
  }
}

@("StableSignTest.FailureRate")
unittest {
  // Verify that StableSign() is able to handle most cases where the three
  // points are as collinear as possible.  (For reference, TriageSign() fails
  // virtually 100% of the time on this test.)
  //
  // Note that the failure rate *decreases* as the points get closer together,
  // and the decrease is approximately linear.  For example, the failure rate
  // is 0.4% for collinear points spaced 1km apart, but only 0.0004% for
  // collinear points spaced 1 meter apart.
  StableSignTest stableSignTest = new StableSignTest();
  Assert.lessThan(stableSignTest.getFailureRate(1.0), 0.01);  //  1km spacing: <  1% (actual 0.4%)
  Assert.lessThan(stableSignTest.getFailureRate(10.0), 0.1);  // 10km spacing: < 10% (actual 4%)
}

// Given 3 points A, B, C that are exactly coplanar with the origin and where
// A < B < C in lexicographic order, verify that ABC is counterclockwise (if
// expected == 1) or clockwise (if expected == -1) using ExpensiveSign().
//
// This method is intended specifically for checking the cases where
// symbolic perturbations are needed to break ties.
static void checkSymbolicSign(int expected, in S2Point a, in S2Point b, in S2Point c) {
  Assert.lessThan(a, b);
  Assert.lessThan(b, c);
  Assert.equal(a.dotProd(b.crossProd(c)), 0.0);

  // Use ASSERT rather than EXPECT to suppress spurious error messages.
  Assert.equal(expected, s2pred.expensiveSign(a, b, c));
  Assert.equal(expected, s2pred.expensiveSign(b, c, a));
  Assert.equal(expected, s2pred.expensiveSign(c, a, b));
  Assert.equal(-expected, s2pred.expensiveSign(c, b, a));
  Assert.equal(-expected, s2pred.expensiveSign(b, a, c));
  Assert.equal(-expected, s2pred.expensiveSign(a, c, b));
}

@("Sign.SymbolicPerturbationCodeCoverage")
unittest {
  // The purpose of this test is simply to get code coverage of
  // SymbolicallyPerturbedSign().  Let M_1, M_2, ... be the sequence of
  // submatrices whose determinant sign is tested by that function.  Then the
  // i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
  //
  //    det(M) = 0
  //    det(M_j) = 0 for j < i
  //    det(M_i) != 0
  //    A < B < C in lexicographic order.
  //
  // I checked that reversing the sign of any of the "return" statements in
  // SymbolicallyPerturbedSign() will cause this test to fail.

  // det(M_1) = b0*c1 - b1*c0
  checkSymbolicSign(1, S2Point(-3, -1, 0), S2Point(-2, 1, 0), S2Point(1, -2, 0));

  // det(M_2) = b2*c0 - b0*c2
  checkSymbolicSign(1, S2Point(-6, 3, 3), S2Point(-4, 2, -1), S2Point(-2, 1, 4));

  // det(M_3) = b1*c2 - b2*c1
  checkSymbolicSign(1, S2Point(0, -1, -1), S2Point(0, 1, -2), S2Point(0, 2, 1));
  // From this point onward, B or C must be zero, or B is proportional to C.

  // det(M_4) = c0*a1 - c1*a0
  checkSymbolicSign(1, S2Point(-1, 2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8));

  // det(M_5) = c0
  checkSymbolicSign(1, S2Point(-4, -2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8));

  // det(M_6) = -c1
  checkSymbolicSign(1, S2Point(0, -5, 7), S2Point(0, -4, 8), S2Point(0, -2, 4));

  // det(M_7) = c2*a0 - c0*a2
  checkSymbolicSign(1, S2Point(-5, -2, 7), S2Point(0, 0, -2), S2Point(0, 0, -1));

  // det(M_8) = c2
  checkSymbolicSign(1, S2Point(0, -2, 7), S2Point(0, 0, 1), S2Point(0, 0, 2));
  // From this point onward, C must be zero.

  // det(M_9) = a0*b1 - a1*b0
  checkSymbolicSign(1, S2Point(-3, 1, 7), S2Point(-1, -4, 1), S2Point(0, 0, 0));

  // det(M_10) = -b0
  checkSymbolicSign(1, S2Point(-6, -4, 7), S2Point(-3, -2, 1), S2Point(0, 0, 0));

  // det(M_11) = b1
  checkSymbolicSign(-1, S2Point(0, -4, 7), S2Point(0, -2, 1), S2Point(0, 0, 0));

  // det(M_12) = a0
  checkSymbolicSign(-1, S2Point(-1, -4, 5), S2Point(0, 0, -3), S2Point(0, 0, 0));

  // det(M_13) = 1
  checkSymbolicSign(1, S2Point(0, -4, 5), S2Point(0, 0, -5), S2Point(0, 0, 0));
}

enum Precision { DOUBLE, REAL, EXACT, SYMBOLIC }

// A helper class that keeps track of how often each precision was used and
// generates a string for logging purposes.
struct PrecisionStats {
public:
  void tally(Precision precision) {
    ++_counts[precision];
  }

  string toString() {
    string result;
    int total = 0;
    foreach (i, precision; traits.EnumMembers!Precision) {
      result ~= format.format("%s=%6d, ", precision, _counts[i]);
      total += _counts[i];
    }
    result ~= format.format("total=%6d", total);
    return result;
  }

 private:
  int[traits.EnumMembers!(Precision).length] _counts;
};

// Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
static S2Point choosePoint() {
  S2Point x = S2Testing.randomPoint();
  for (int i = 0; i < 3; ++i) {
    if (S2Testing.rnd.oneIn(3)) {
      x.data[i] *= math.pow(1e-50, S2Testing.rnd.randDouble());
    }
  }
  return x.normalize();
}

// The following helper classes allow us to test the various distance
// calculation methods using a common test framework.
class Sin2Distances {
public:
  static int triage(T)(in Vector!(T, 3) x, in Vector!(T, 3) a, in Vector!(T, 3) b) {
    return s2pred.triageCompareSin2Distances(x, a, b);
  }
}

class CosDistances {
public:
  static int triage(T)(in Vector!(T, 3) x, in Vector!(T, 3) a, in Vector!(T, 3) b) {
    return s2pred.triageCompareCosDistances(x, a, b);
  }
}

// Compares distances greater than 90 degrees using sin^2(distance).
class MinusSin2Distances {
public:
  static int triage(T)(in Vector!(T, 3) x, in Vector!(T, 3) a, in Vector!(T, 3) b) {
    return -s2pred.triageCompareSin2Distances(-x, a, b);
  }
}

// Verifies that CompareDistances(x, a, b) == expected_sign, and furthermore
// checks that the minimum required precision is "expected_prec" when the
// distance calculation method defined by CompareDistancesWrapperT is used.
void testCompareDistances(CompareDistancesWrapperT)(
    S2Point x, S2Point a, S2Point b, int expected_sign, Precision expected_prec) {
  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(x)) x = x.normalize();
  if (!isUnitLength(a)) a = a.normalize();
  if (!isUnitLength(b)) b = b.normalize();

  int dbl_sign = CompareDistancesWrapperT.triage(x, a, b);
  int r_sign = CompareDistancesWrapperT.triage(
      Vector3_r.from(x), Vector3_r.from(a), Vector3_r.from(b));
  int exact_sign = s2pred.exactCompareDistances(
      s2pred.Vector3_xf.from(x), s2pred.Vector3_xf.from(a), s2pred.Vector3_xf.from(b));
  int actual_sign =
      exact_sign != 0 ? exact_sign : s2pred.symbolicCompareDistances(x, a, b);

  // Check that the signs are correct (if non-zero), and also that if dbl_sign
  // is non-zero then so is r_sign, etc.
  Assert.equal(expected_sign, actual_sign);
  if (exact_sign != 0) Assert.equal(exact_sign, actual_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);

  Precision actual_prec = dbl_sign ? Precision.DOUBLE
      : r_sign ? Precision.REAL
      : exact_sign ? Precision.EXACT
      : Precision.SYMBOLIC;
  Assert.equal(expected_prec, actual_prec);

  // Make sure that the top-level function returns the expected result.
  Assert.equal(expected_sign, s2pred.compareDistances(x, a, b));

  // Check that reversing the arguments negates the result.
  Assert.equal(-expected_sign, s2pred.compareDistances(x, b, a));
}

@("CompareDistances.Coverage")
unittest {
  // This test attempts to exercise all the code paths in all precisions.

  // Test TriageCompareSin2Distances.
  testCompareDistances!Sin2Distances(
      S2Point(1, 1, 1), S2Point(1, 1 - 1e-15, 1), S2Point(1, 1, 1 + 2e-15),
      -1, Precision.DOUBLE);
  testCompareDistances!Sin2Distances(
      S2Point(1, 1, 0), S2Point(1, 1 - 1e-15, 1e-21), S2Point(1, 1 - 1e-15, 0),
      1, Precision.DOUBLE);
  testCompareDistances!Sin2Distances(
      S2Point(2, 0, 0), S2Point(2, -1, 0), S2Point(2, 1, 1e-8),
      -1, Precision.REAL);
  testCompareDistances!Sin2Distances(
      S2Point(2, 0, 0), S2Point(2, -1, 0), S2Point(2, 1, 1e-100),
      -1, Precision.EXACT);
  testCompareDistances!Sin2Distances(
      S2Point(1, 0, 0), S2Point(1, -1, 0), S2Point(1, 1, 0),
      1, Precision.SYMBOLIC);
  testCompareDistances!Sin2Distances(
      S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(1, 0, 0),
      0, Precision.SYMBOLIC);

  // Test TriageCompareCosDistances.
  testCompareDistances!CosDistances(
      S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1, 1, 3e-15),
      1, Precision.DOUBLE);
  testCompareDistances!CosDistances(
      S2Point(1, 0, 0), S2Point(1, 1e-30, 0), S2Point(-1, 1e-40, 0),
      -1, Precision.DOUBLE);
  testCompareDistances!CosDistances(
      S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1, 1, 3e-18),
      1, Precision.REAL);
  testCompareDistances!CosDistances(
      S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1, 1, 1e-100),
      1, Precision.EXACT);
  testCompareDistances!CosDistances(
      S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(-1, 1, 0),
      -1, Precision.SYMBOLIC);
  testCompareDistances!CosDistances(
      S2Point(1, 1, 1), S2Point(1, -1, 0), S2Point(1, -1, 0),
      0, Precision.SYMBOLIC);

  // Test TriageCompareSin2Distances using distances greater than 90 degrees.
  testCompareDistances!MinusSin2Distances(
      S2Point(1, 1, 0), S2Point(-1, -1 + 1e-15, 0), S2Point(-1, -1, 0),
      -1, Precision.DOUBLE);
  testCompareDistances!MinusSin2Distances(
      S2Point(-1, -1, 0), S2Point(1, 1 - 1e-15, 0),
      S2Point(1, 1 - 1e-15, 1e-21), 1, Precision.DOUBLE);
  testCompareDistances!MinusSin2Distances(
      S2Point(-1, -1, 0), S2Point(2, 1, 0), S2Point(2, 1, 1e-8),
      1, Precision.REAL);
  testCompareDistances!MinusSin2Distances(
      S2Point(-1, -1, 0), S2Point(2, 1, 0), S2Point(2, 1, 1e-30),
      1, Precision.EXACT);
  testCompareDistances!MinusSin2Distances(
      S2Point(-1, -1, 0), S2Point(2, 1, 0), S2Point(1, 2, 0),
      -1, Precision.SYMBOLIC);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testCompareDistancesConsistency(CompareDistancesWrapperT)(
    in S2Point x, in S2Point a, in S2Point b) {
  int dbl_sign = CompareDistancesWrapperT.triage(x, a, b);
  int r_sign = CompareDistancesWrapperT.triage(
      Vector3_r.from(x), Vector3_r.from(a), Vector3_r.from(b));
  int exact_sign = s2pred.exactCompareDistances(
      s2pred.Vector3_xf.from(x), s2pred.Vector3_xf.from(a), s2pred.Vector3_xf.from(b));
  if (dbl_sign != 0) {
    Assert.equal(r_sign, dbl_sign);
  }
  if (r_sign != 0) {
    Assert.equal(exact_sign, r_sign);
  }
  if (exact_sign != 0) {
    Assert.equal(exact_sign, s2pred.compareDistances(x, a, b));
    return (r_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.REAL : Precision.DOUBLE;
  } else {
    // Unlike the other methods, SymbolicCompareDistances has the
    // precondition that the exact sign must be zero.
    int symbolic_sign = s2pred.symbolicCompareDistances(x, a, b);
    Assert.equal(symbolic_sign, s2pred.compareDistances(x, a, b));
    return Precision.SYMBOLIC;
  }
}

@("CompareDistances.Consistency")
unittest {
  // This test chooses random point pairs that are nearly equidistant from a
  // target point, and then checks that the answer given by a method at one
  // level of precision is consistent with the answer given at the next higher
  // level of precision.
  //
  // The way the .cc file is structured, we can only do comparisons using a
  // specific precision if we also choose the specific distance calculation
  // method.  The code below checks that the Cos, Sin2, and MinusSin2 methods
  // are consistent across their entire valid range of inputs, and also
  // simulates the logic in CompareDistance that chooses which method to use
  // in order to gather statistics about how often each precision is needed.
  // (These statistics are only useful for coverage purposes, not benchmarks,
  // since the input points are chosen to be pathological worst cases.)
  testCompareDistancesConsistency!CosDistances(
      S2Point(1, 0, 0), S2Point(0, -1, 0), S2Point(0, 1, 0));
  auto rnd = S2Testing.rnd;
  PrecisionStats sin2_stats, cos_stats, minus_sin2_stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point x = choosePoint();
    S2Point dir = choosePoint();
    S1Angle r = S1Angle.fromRadians(math.PI_2 * math.pow(1e-30, rnd.randDouble()));
    if (rnd.oneIn(2)) r = S1Angle.fromRadians(math.PI_2) - r;
    if (rnd.oneIn(2)) r = S1Angle.fromRadians(math.PI_2) + r;
    S2Point a = interpolateAtDistance(r, x, dir);
    S2Point b = interpolateAtDistance(r, x, -dir);
    Precision prec = testCompareDistancesConsistency!CosDistances(x, a, b);
    if (r.degrees() >= 45 && r.degrees() <= 135) cos_stats.tally(prec);
    // The Sin2 method is only valid if both distances are less than 90
    // degrees, and similarly for the MinusSin2 method.  (In the actual
    // implementation these methods are only used if both distances are less
    // than 45 degrees or greater than 135 degrees respectively.)
    if (r.radians() < math.PI_2 - 1e-14) {
      prec = testCompareDistancesConsistency!Sin2Distances(x, a, b);
      if (r.degrees() < 45) {
        // Don't skew the statistics by recording degenerate inputs.
        if (a == b) {
          Assert.equal(Precision.SYMBOLIC, prec);
        } else {
          sin2_stats.tally(prec);
        }
      }
    } else if (r.radians() > math.PI_2 + 1e-14) {
      prec = testCompareDistancesConsistency!MinusSin2Distances(x, a, b);
      if (r.degrees() > 135) minus_sin2_stats.tally(prec);
    }
  }
  writeln("\nsin2:  ", sin2_stats.toString(), "\ncos:   ", cos_stats.toString(),
      "\n-sin2: ", minus_sin2_stats.toString());
}

// Helper classes for testing the various distance calculation methods.
class Sin2Distance {
public:
  static int triage(T)(in Vector!(T, 3) x, in Vector!(T, 3) y, S1ChordAngle r) {
    return s2pred.triageCompareSin2Distance(x, y, cast(T) r.length2());
  }
}

class CosDistance {
public:
  static int triage(T)(in Vector!(T, 3) x, in Vector!(T, 3) y, S1ChordAngle r) {
    return s2pred.triageCompareCosDistance(x, y, cast(T) r.length2());
  }
}

// Verifies that CompareDistance(x, y, r) == expected_sign, and furthermore
// checks that the minimum required precision is "expected_prec" when the
// distance calculation method defined by CompareDistanceWrapper is used.
void testCompareDistance(CompareDistanceWrapperT)(
    S2Point x, S2Point y, S1ChordAngle r, int expected_sign, Precision expected_prec) {
  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(x)) x = x.normalize();
  if (!isUnitLength(y)) y = y.normalize();

  int dbl_sign = CompareDistanceWrapperT.triage(x, y, r);
  int r_sign = CompareDistanceWrapperT.triage(Vector3_r.from(x), Vector3_r.from(y), r);
  int exact_sign = s2pred.exactCompareDistance(
      s2pred.Vector3_xf.from(x), s2pred.Vector3_xf.from(y), ExactFloat(r.length2()));

  // Check that the signs are correct (if non-zero), and also that if dbl_sign
  // is non-zero then so is ld_sign, etc.
  Assert.equal(expected_sign, exact_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);

  Precision actual_prec = dbl_sign ? Precision.DOUBLE : r_sign ? Precision.REAL : Precision.EXACT;
  Assert.equal(expected_prec, actual_prec);

  // Make sure that the top-level function returns the expected result.
  Assert.equal(expected_sign, s2pred.compareDistance(x, y, r));

  // Mathematically, if d(X, Y) < r then d(-X, Y) > (Pi - r).  Unfortunately
  // there can be rounding errors when computing the supplementary distance,
  // so to ensure the two distances are exactly supplementary we need to do
  // the following.
  S1ChordAngle r_supp = S1ChordAngle.straight() - r;
  r = S1ChordAngle.straight() - r_supp;
  Assert.equal(-s2pred.compareDistance(x, y, r), s2pred.compareDistance(-x, y, r_supp));
}

@("CompareDistance.Coverage")
unittest {
  // Test TriageCompareSin2Distance.
  testCompareDistance!Sin2Distance(
      S2Point(1, 1, 1), S2Point(1, 1 - 1e-15, 1),
      S1ChordAngle.fromRadians(1e-15), -1, Precision.DOUBLE);
  testCompareDistance!Sin2Distance(
      S2Point(1, 0, 0), S2Point(1, 1, 0),
      S1ChordAngle.fromRadians(math.PI_4), -1, Precision.REAL);
  testCompareDistance!Sin2Distance(
      S2Point(1, 1e-40, 0), S2Point(1 + double.epsilon, 1e-40, 0),
      S1ChordAngle.fromRadians(0.9 * double.epsilon * 1e-40), 1, Precision.EXACT);
  testCompareDistance!Sin2Distance(
      S2Point(1, 1e-40, 0), S2Point(1 + double.epsilon, 1e-40, 0),
      S1ChordAngle.fromRadians(1.1 * double.epsilon * 1e-40), -1, Precision.EXACT);
  testCompareDistance!Sin2Distance(
      S2Point(1, 0, 0), S2Point(1 + double.epsilon, 0, 0),
      S1ChordAngle.zero(), 0, Precision.EXACT);

  // Test TriageCompareCosDistance.
  testCompareDistance!CosDistance(
      S2Point(1, 0, 0), S2Point(1, 1e-8, 0),
      S1ChordAngle.fromRadians(1e-7), -1, Precision.DOUBLE);
  testCompareDistance!CosDistance(
      S2Point(1, 0, 0), S2Point(-1, 1e-8, 0),
      S1ChordAngle.fromRadians(math.PI - 1e-7), 1, Precision.DOUBLE);
  testCompareDistance!CosDistance(
      S2Point(1, 1, 0), S2Point(1, -1 - 2 * double.epsilon, 0),
      S1ChordAngle.right(), 1, Precision.DOUBLE);
  testCompareDistance!CosDistance(
      S2Point(1, 1, 0), S2Point(1, -1 - double.epsilon, 0),
      S1ChordAngle.right(), 1, Precision.REAL);
  testCompareDistance!CosDistance(
      S2Point(1, 1, 0), S2Point(1, -1, 1e-30),
      S1ChordAngle.right(), 0, Precision.EXACT);
  // The angle between these two points is exactly 60 degrees.
  testCompareDistance!CosDistance(
      S2Point(1, 1, 0), S2Point(0, 1, 1),
      S1ChordAngle.fromLength2(1), 0, Precision.EXACT);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testCompareDistanceConsistency(CompareDistanceWrapperT)(
    in S2Point x, in S2Point y, S1ChordAngle r) {
  int dbl_sign = CompareDistanceWrapperT.triage(x, y, r);
  int r_sign = CompareDistanceWrapperT.triage(Vector3_r.from(x), Vector3_r.from(y), r);
  int exact_sign = s2pred.exactCompareDistance(
      s2pred.Vector3_xf.from(x), s2pred.Vector3_xf.from(y), ExactFloat(r.length2()));
  Assert.equal(exact_sign, s2pred.compareDistance(x, y, r));
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  return (r_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.REAL : Precision.DOUBLE;
}

@("CompareDistance.Consistency")
unittest {
  // This test chooses random inputs such that the distance between points X
  // and Y is very close to the threshold distance "r".  It then checks that
  // the answer given by a method at one level of precision is consistent with
  // the answer given at the next higher level of precision.  See also the
  // comments in the CompareDistances consistency test.
  auto rnd = S2Testing.rnd;
  PrecisionStats sin2_stats, cos_stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point x = choosePoint();
    S2Point dir = choosePoint();
    S1Angle r = S1Angle.fromRadians(math.PI_2 * math.pow(1e-30, rnd.randDouble()));
    if (rnd.oneIn(2)) r = S1Angle.fromRadians(math.PI_2) - r;
    if (rnd.oneIn(5)) r = S1Angle.fromRadians(math.PI_2) + r;
    S2Point y = interpolateAtDistance(r, x, dir);
    Precision prec = testCompareDistanceConsistency!CosDistance(x, y, S1ChordAngle(r));
    if (r.degrees() >= 45) cos_stats.tally(prec);
    if (r.radians() < math.PI_2 - 1e-14) {
      prec = testCompareDistanceConsistency!Sin2Distance(
          x, y, S1ChordAngle(r));
      if (r.degrees() < 45) sin2_stats.tally(prec);
    }
  }
  writeln(
      "\nsin2:  ", sin2_stats.toString(),
      "\ncos:   ", cos_stats.toString());
}

// Verifies that CompareEdgeDistance(x, a0, a1, r) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
void testCompareEdgeDistance(S2Point x, S2Point a0, S2Point a1, S1ChordAngle r,
                             int expected_sign, Precision expected_prec) {

  scope(failure) writeln(
      "testCompareEdgeDistance failure::",
      "\n  x=", x, ", a0=", a0, ", a1=", a1, ", r=", r,
      "\n  expected_sign=", expected_sign, ", expected_prec=", expected_prec);

  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(x)) x = x.normalize();
  if (!isUnitLength(a0)) a0 = a0.normalize();
  if (!isUnitLength(a1)) a1 = a1.normalize();

  int dbl_sign = s2pred.triageCompareEdgeDistance(x, a0, a1, r.length2());
  int r_sign = s2pred.triageCompareEdgeDistance(
      Vector3_r.from(x), Vector3_r.from(a0), Vector3_r.from(a1), r.length2());
  int exact_sign = s2pred.exactCompareEdgeDistance(x, a0, a1, r);

  // Check that the signs are correct (if non-zero), and also that if dbl_sign
  // is non-zero then so is ld_sign, etc.
  Assert.equal(expected_sign, exact_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);

  Precision actual_prec = dbl_sign ? Precision.DOUBLE : r_sign ? Precision.REAL : Precision.EXACT;
  Assert.equal(actual_prec, expected_prec);

  // Make sure that the top-level function returns the expected result.
  Assert.equal(s2pred.compareEdgeDistance(x, a0, a1, r), expected_sign);
}

@("CompareEdgeDistance.Coverage")
unittest {
  // Test TriageCompareLineSin2Distance.
  testCompareEdgeDistance(
      S2Point(1, 1e-10, 1e-15), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.fromRadians(1e-15 + double.epsilon), -1, Precision.DOUBLE);
  testCompareEdgeDistance(
      S2Point(1, 1, 1e-15), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.fromRadians(1e-15 + double.epsilon), -1, Precision.REAL);
  testCompareEdgeDistance(
      S2Point(1, 1, 1e-40), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.fromRadians(1e-40), -1, Precision.EXACT);
  testCompareEdgeDistance(
      S2Point(1, 1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.zero(), 0, Precision.EXACT);

  // Test TriageCompareLineCos2Distance.
  testCompareEdgeDistance(
      S2Point(1e-15, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.fromRadians(cast(double) math.PI_2 - 1e-15 - 5 * double.epsilon),
      1, Precision.DOUBLE);
  testCompareEdgeDistance(
      S2Point(1e-15, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.fromRadians(cast(double) math.PI_2 - 1e-15 - double.epsilon),
      1, Precision.REAL);
  testCompareEdgeDistance(
      S2Point(1e-40, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.right(), -1, Precision.EXACT);
  testCompareEdgeDistance(
      S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
      S1ChordAngle.right(), 0, Precision.EXACT);

  // Test cases where the closest point is an edge endpoint.
  testCompareEdgeDistance(
      S2Point(1e-15, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
      S1ChordAngle.right(), -1, Precision.DOUBLE);
  testCompareEdgeDistance(
      S2Point(1e-18, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
      S1ChordAngle.right(), -1, Precision.REAL);
  testCompareEdgeDistance(
      S2Point(1e-100, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
      S1ChordAngle.right(), -1, Precision.EXACT);
  testCompareEdgeDistance(
      S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
      S1ChordAngle.right(), 0, Precision.EXACT);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testCompareEdgeDistanceConsistency(
    in S2Point x, in S2Point a0, in S2Point a1, S1ChordAngle r) {
  int dbl_sign = s2pred.triageCompareEdgeDistance(x, a0, a1, r.length2());
  int r_sign = s2pred.triageCompareEdgeDistance(
      Vector3_r.from(x), Vector3_r.from(a0), Vector3_r.from(a1), r.length2());
  int exact_sign = s2pred.exactCompareEdgeDistance(x, a0, a1, r);
  Assert.equal(exact_sign, s2pred.compareEdgeDistance(x, a0, a1, r));
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  return (r_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.REAL : Precision.DOUBLE;
}

@("CompareEdgeDistance.Consistency")
unittest {
  // This test chooses random inputs such that the distance between "x" and
  // the line (a0, a1) is very close to the threshold distance "r".  It then
  // checks that the answer given by a method at one level of precision is
  // consistent with the answer given at the next higher level of precision.
  // See also the comments in the CompareDistances consistency test.
  auto rnd = S2Testing.rnd;
  PrecisionStats stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point a0 = choosePoint();
    S1Angle len = S1Angle.fromRadians(math.PI * math.pow(1e-20, rnd.randDouble()));
    S2Point a1 = interpolateAtDistance(len, a0, choosePoint());
    if (rnd.oneIn(2)) a1 = -a1;
    if (a0 == -a1) continue;  // Not allowed by API.
    S2Point n = robustCrossProd(a0, a1).normalize();
    double f = math.pow(1e-20, rnd.randDouble());
    S2Point a = ((1 - f) * a0 + f * a1).normalize();
    S1Angle r = S1Angle.fromRadians(math.PI_2 * math.pow(1e-20, rnd.randDouble()));
    if (rnd.oneIn(2)) r = S1Angle.fromRadians(math.PI_2) - r;
    S2Point x = interpolateAtDistance(r, a, n);
    if (rnd.oneIn(5)) {
      // Replace "x" with a random point that is closest to an edge endpoint.
      do {
        x = choosePoint();
      } while (s2pred.compareEdgeDirections(a0, x, a0, a1) > 0 &&
               s2pred.compareEdgeDirections(x, a1, a0, a1) > 0);
      r = algorithm.min(S1Angle(x, a0), S1Angle(x, a1));
    }
    Precision prec = testCompareEdgeDistanceConsistency(x, a0, a1, S1ChordAngle(r));
    stats.tally(prec);
  }
  writeln(stats.toString());
}

// Verifies that CompareEdgeDirections(a0, a1, b0, b1) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
void testCompareEdgeDirections(
    S2Point a0, S2Point a1, S2Point b0, S2Point b1, int expected_sign, Precision expected_prec) {
  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(a0)) a0 = a0.normalize();
  if (!isUnitLength(a1)) a1 = a1.normalize();
  if (!isUnitLength(b0)) b0 = b0.normalize();
  if (!isUnitLength(b1)) b1 = b1.normalize();

  int dbl_sign = s2pred.triageCompareEdgeDirections(a0, a1, b0, b1);
  int r_sign = s2pred.triageCompareEdgeDirections(
      Vector3_r.from(a0), Vector3_r.from(a1), Vector3_r.from(b0), Vector3_r.from(b1));
  int exact_sign = s2pred.exactCompareEdgeDirections(
      s2pred.Vector3_xf.from(a0), s2pred.Vector3_xf.from(a1), s2pred.Vector3_xf.from(b0),
      s2pred.Vector3_xf.from(b1));

  // Check that the signs are correct (if non-zero), and also that if dbl_sign
  // is non-zero then so is ld_sign, etc.
  Assert.equal(expected_sign, exact_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);

  Precision actual_prec = dbl_sign ? Precision.DOUBLE : r_sign ? Precision.REAL : Precision.EXACT;
  Assert.equal(expected_prec, actual_prec);

  // Make sure that the top-level function returns the expected result.
  Assert.equal(expected_sign, s2pred.compareEdgeDirections(a0, a1, b0, b1));

  // Check various identities involving swapping or negating arguments.
  Assert.equal(expected_sign, s2pred.compareEdgeDirections(b0, b1, a0, a1));
  Assert.equal(expected_sign, s2pred.compareEdgeDirections(-a0, -a1, b0, b1));
  Assert.equal(expected_sign, s2pred.compareEdgeDirections(a0, a1, -b0, -b1));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(a1, a0, b0, b1));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(a0, a1, b1, b0));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(-a0, a1, b0, b1));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(a0, -a1, b0, b1));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(a0, a1, -b0, b1));
  Assert.equal(-expected_sign, s2pred.compareEdgeDirections(a0, a1, b0, -b1));
}

@("CompareEdgeDirections.Coverage")
unittest {
  testCompareEdgeDirections(
      S2Point(1, 0, 0), S2Point(1, 1, 0), S2Point(1, -1, 0), S2Point(1, 0, 0),
      1, Precision.DOUBLE);
  testCompareEdgeDirections(
      S2Point(1, 0, 1.5e-15), S2Point(1, 1, 0), S2Point(0, -1, 0), S2Point(0, 0, 1),
      1, Precision.DOUBLE);
  testCompareEdgeDirections(
      S2Point(1, 0, 1e-18), S2Point(1, 1, 0), S2Point(0, -1, 0), S2Point(0, 0, 1),
      1, Precision.REAL);
  testCompareEdgeDirections(
      S2Point(1, 0, 1e-50), S2Point(1, 1, 0), S2Point(0, -1, 0), S2Point(0, 0, 1),
      1, Precision.EXACT);
  testCompareEdgeDirections(
      S2Point(1, 0, 0), S2Point(1, 1, 0), S2Point(0, -1, 0), S2Point(0, 0, 1),
      0, Precision.EXACT);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testCompareEdgeDirectionsConsistency(
    in S2Point a0, in S2Point a1, in S2Point b0, in S2Point b1) {
  int dbl_sign = s2pred.triageCompareEdgeDirections(a0, a1, b0, b1);
  int r_sign = s2pred.triageCompareEdgeDirections(
      Vector3_r.from(a0), Vector3_r.from(a1), Vector3_r.from(b0), Vector3_r.from(b1));
  int exact_sign = s2pred.exactCompareEdgeDirections(
      s2pred.Vector3_xf.from(a0), s2pred.Vector3_xf.from(a1), s2pred.Vector3_xf.from(b0),
      s2pred.Vector3_xf.from(b1));
  Assert.equal(exact_sign, s2pred.compareEdgeDirections(a0, a1, b0, b1));
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  return (r_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.REAL : Precision.DOUBLE;
}

@("CompareEdgeDirections.Consistency")
unittest {
  // This test chooses random pairs of edges that are nearly perpendicular,
  // then checks that the answer given by a method at one level of precision
  // is consistent with the answer given at the next higher level of
  // precision.  See also the comments in the CompareDistances test.
  auto rnd = S2Testing.rnd;
  PrecisionStats stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point a0 = choosePoint();
    S1Angle a_len = S1Angle.fromRadians(math.PI * math.pow(1e-20, rnd.randDouble()));
    S2Point a1 = interpolateAtDistance(a_len, a0, choosePoint());
    S2Point a_norm = robustCrossProd(a0, a1).normalize();
    S2Point b0 = choosePoint();
    S1Angle b_len = S1Angle.fromRadians(math.PI * math.pow(1e-20, rnd.randDouble()));
    S2Point b1 = interpolateAtDistance(b_len, b0, a_norm);
    if (a0 == -a1 || b0 == -b1) continue;  // Not allowed by API.
    Precision prec = testCompareEdgeDirectionsConsistency(a0, a1, b0, b1);
    // Don't skew the statistics by recording degenerate inputs.
    if (a0 == a1 || b0 == b1) {
      Assert.equal(Precision.EXACT, prec);
    } else {
      stats.tally(prec);
    }
  }
  writeln(stats.toString());
}

// Verifies that EdgeCircumcenterSign(x0, x1, a, b, c) == expected_sign, and
// furthermore checks that the minimum required precision is "expected_prec".
void testEdgeCircumcenterSign(
    S2Point x0, S2Point x1, S2Point a, S2Point b, S2Point c,
    int expected_sign, Precision expected_prec) {
  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(x0)) x0 = x0.normalize();
  if (!isUnitLength(x1)) x1 = x1.normalize();
  if (!isUnitLength(a)) a = a.normalize();
  if (!isUnitLength(b)) b = b.normalize();
  if (!isUnitLength(c)) c = c.normalize();

  int abc_sign = s2pred.sign(a, b, c);
  int dbl_sign = s2pred.triageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign);
  int r_sign = s2pred.triageEdgeCircumcenterSign(
      Vector3_r.from(x0), Vector3_r.from(x1), Vector3_r.from(a), Vector3_r.from(b),
      Vector3_r.from(c), abc_sign);
  int exact_sign = s2pred.exactEdgeCircumcenterSign(
      s2pred.Vector3_xf.from(x0), s2pred.Vector3_xf.from(x1), s2pred.Vector3_xf.from(a),
      s2pred.Vector3_xf.from(b), s2pred.Vector3_xf.from(c), abc_sign);
  int actual_sign = (exact_sign != 0 ? exact_sign :
                     s2pred.symbolicEdgeCircumcenterSign(x0, x1, a, b, c));

  // Check that the signs are correct (if non-zero), and also that if dbl_sign
  // is non-zero then so is ld_sign, etc.
  Assert.equal(expected_sign, actual_sign);
  if (exact_sign != 0) Assert.equal(exact_sign, actual_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);

  Precision actual_prec = (dbl_sign ? Precision.DOUBLE :
                           r_sign ? Precision.REAL :
                           exact_sign ? Precision.EXACT : Precision.SYMBOLIC);
  Assert.equal(expected_prec, actual_prec);

  // Make sure that the top-level function returns the expected result.
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, a, b, c));

  // Check various identities involving swapping or negating arguments.
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, a, c, b));
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, b, a, c));
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, b, c, a));
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, c, a, b));
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(x0, x1, c, b, a));
  Assert.equal(-expected_sign, s2pred.edgeCircumcenterSign(x1, x0, a, b, c));
  Assert.equal(expected_sign, s2pred.edgeCircumcenterSign(-x0, -x1, a, b, c));
  if (actual_sign == exact_sign) {
    // Negating the input points may not preserve the result when symbolic
    // perturbations are used, since -X is not an exact multiple of X.
    Assert.equal(-expected_sign, s2pred.edgeCircumcenterSign(x0, x1, -a, -b, -c));
  }
}

@("EdgeCircumcenterSign.Coverage")
unittest {
  testEdgeCircumcenterSign(
      S2Point(1, 0, 0), S2Point(1, 1, 0),
      S2Point(0, 0, 1), S2Point(1, 0, 1), S2Point(0, 1, 1),
      1, Precision.DOUBLE);
  testEdgeCircumcenterSign(
      S2Point(1, 0, 0), S2Point(1, 1, 0),
      S2Point(0, 0, -1), S2Point(1, 0, -1), S2Point(0, 1, -1),
      -1, Precision.DOUBLE);
  testEdgeCircumcenterSign(
      S2Point(1, -1, 0), S2Point(1, 1, 0),
      S2Point(1, -1e-5, 1), S2Point(1, 1e-5, -1), S2Point(1, 1 - 1e-5, 1e-5),
      -1, Precision.DOUBLE);
  testEdgeCircumcenterSign(
      S2Point(1, -1, 0), S2Point(1, 1, 0),
      S2Point(1, -1e-5, 1), S2Point(1, 1e-5, -1), S2Point(1, 1 - 1e-9, 1e-5),
      -1, Precision.REAL);
  testEdgeCircumcenterSign(
      S2Point(1, -1, 0), S2Point(1, 1, 0),
      S2Point(1, -1e-5, 1), S2Point(1, 1e-5, -1), S2Point(1, 1 - 1e-15, 1e-5),
      -1, Precision.EXACT);
  testEdgeCircumcenterSign(
      S2Point(1, -1, 0), S2Point(1, 1, 0),
      S2Point(1, -1e-5, 1), S2Point(1, 1e-5, -1), S2Point(1, 1, 1e-5),
      1, Precision.SYMBOLIC);

  // This test falls back to the second symbolic perturbation:
  testEdgeCircumcenterSign(
      S2Point(1, -1, 0), S2Point(1, 1, 0),
      S2Point(0, -1, 0), S2Point(0, 0, -1), S2Point(0, 0, 1),
      -1, Precision.SYMBOLIC);

  // This test falls back to the third symbolic perturbation:
  testEdgeCircumcenterSign(
      S2Point(0, -1, 1), S2Point(0, 1, 1),
      S2Point(0, 1, 0), S2Point(0, -1, 0), S2Point(1, 0, 0),
      -1, Precision.SYMBOLIC);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testEdgeCircumcenterSignConsistency(
    in S2Point x0, in S2Point x1,
    in S2Point a, in S2Point b, in S2Point c) {
  int abc_sign = s2pred.sign(a, b, c);
  int dbl_sign = s2pred.triageEdgeCircumcenterSign(x0, x1, a, b, c, abc_sign);
  int r_sign = s2pred.triageEdgeCircumcenterSign(
      Vector3_r.from(x0), Vector3_r.from(x1), Vector3_r.from(a), Vector3_r.from(b),
      Vector3_r.from(c), abc_sign);
  int exact_sign = s2pred.exactEdgeCircumcenterSign(
      s2pred.Vector3_xf.from(x0), s2pred.Vector3_xf.from(x1), s2pred.Vector3_xf.from(a),
      s2pred.Vector3_xf.from(b), s2pred.Vector3_xf.from(c), abc_sign);
  if (dbl_sign != 0) Assert.equal(r_sign, dbl_sign);
  if (r_sign != 0) Assert.equal(exact_sign, r_sign);
  if (exact_sign != 0) {
    Assert.equal(exact_sign, s2pred.edgeCircumcenterSign(x0, x1, a, b, c));
    return (r_sign == 0) ? Precision.EXACT : (dbl_sign == 0) ? Precision.REAL : Precision.DOUBLE;
  } else {
    // Unlike the other methods, SymbolicEdgeCircumcenterSign has the
    // precondition that the exact sign must be zero.
    int symbolic_sign = s2pred.symbolicEdgeCircumcenterSign(x0, x1, a, b, c);
    Assert.equal(symbolic_sign, s2pred.edgeCircumcenterSign(x0, x1, a, b, c));
    return Precision.SYMBOLIC;
  }
}

@("EdgeCircumcenterSign.Consistency")
unittest {
  // This test chooses random a random edge X, then chooses a random point Z
  // on the great circle through X, and finally choose three points A, B, C
  // that are nearly equidistant from X.  It then checks that the answer given
  // by a method at one level of precision is consistent with the answer given
  // at the next higher level of precision.
  auto rnd = S2Testing.rnd;
  PrecisionStats stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point x0 = choosePoint();
    S2Point x1 = choosePoint();
    if (x0 == -x1) continue;  // Not allowed by API.
    double c0 = (rnd.oneIn(2) ? -1 : 1) * math.pow(1e-20, rnd.randDouble());
    double c1 = (rnd.oneIn(2) ? -1 : 1) * math.pow(1e-20, rnd.randDouble());
    S2Point z = (c0 * x0 + c1 * x1).normalize();
    S1Angle r = S1Angle.fromRadians(math.PI * math.pow(1e-30, rnd.randDouble()));
    S2Point a = interpolateAtDistance(r, z, choosePoint());
    S2Point b = interpolateAtDistance(r, z, choosePoint());
    S2Point c = interpolateAtDistance(r, z, choosePoint());
    Precision prec = testEdgeCircumcenterSignConsistency(x0, x1, a, b, c);
    // Don't skew the statistics by recording degenerate inputs.
    if (x0 == x1) {
      // This precision would be SYMBOLIC if we handled this degeneracy.
      Assert.equal(Precision.EXACT, prec);
    } else if (a == b || b == c || c == a) {
      Assert.equal(Precision.SYMBOLIC, prec);
    } else {
      stats.tally(prec);
    }
  }
  writeln(stats.toString());
}

// Verifies that VoronoiSiteExclusion(a, b, x0, x1, r) == expected_result, and
// furthermore checks that the minimum required precision is "expected_prec".
void testVoronoiSiteExclusion(
    S2Point a, S2Point b, S2Point x0, S2Point x1, S1ChordAngle r,
    s2pred.Excluded expected_result, Precision expected_prec) {
  s2pred.Excluded UNCERTAIN = s2pred.Excluded.UNCERTAIN;

  // Don't normalize the arguments unless necessary (to allow testing points
  // that differ only in magnitude).
  if (!isUnitLength(a)) a = a.normalize();
  if (!isUnitLength(b)) b = b.normalize();
  if (!isUnitLength(x0)) x0 = x0.normalize();
  if (!isUnitLength(x1)) x1 = x1.normalize();

  // The internal methods (Triage, Exact, etc) require that site A is closer
  // to X0 and site B is closer to X1.  GetVoronoiSiteExclusion has special
  // code to handle the case where this is not true.  We need to duplicate
  // that code here.  Essentially, since the API requires site A to be closer
  // than site B to X0, then if site A is also closer to X1 then site B must
  // be excluded.
  if (s2pred.compareDistances(x1, a, b) < 0) {
    Assert.equal(expected_result, s2pred.Excluded.SECOND);
    // We don't know what precision was used by CompareDistances(), but we
    // arbitrarily require the test to specify it as DOUBLE.
    Assert.equal(expected_prec, Precision.DOUBLE);
  } else {
    s2pred.Excluded dbl_result = s2pred.triageVoronoiSiteExclusion(a, b, x0, x1, r.length2());
    s2pred.Excluded r_result = s2pred.triageVoronoiSiteExclusion(
        Vector3_r.from(a), Vector3_r.from(b), Vector3_r.from(x0), Vector3_r.from(x1),
        r.length2());
    s2pred.Excluded exact_result = s2pred.exactVoronoiSiteExclusion(
        s2pred.Vector3_xf.from(a), s2pred.Vector3_xf.from(b), s2pred.Vector3_xf.from(x0),
        s2pred.Vector3_xf.from(x1), ExactFloat(r.length2()));

    // Check that the results are correct (if not UNCERTAIN), and also that if
    // dbl_result is not UNCERTAIN then so is r_result, etc.
    Assert.equal(expected_result, exact_result);
    if (r_result != UNCERTAIN) Assert.equal(exact_result, r_result);
    if (dbl_result != UNCERTAIN) Assert.equal(r_result, dbl_result);

    Precision actual_prec = (dbl_result != UNCERTAIN ? Precision.DOUBLE :
                             r_result != UNCERTAIN ? Precision.REAL : Precision.EXACT);
    Assert.equal(expected_prec, actual_prec);
  }
  // Make sure that the top-level function returns the expected result.
  Assert.equal(expected_result, s2pred.getVoronoiSiteExclusion(a, b, x0, x1, r));

  // If site B is closer to X1, then the same site should be excluded (if any)
  // when we swap the sites and the edge direction.
  s2pred.Excluded swapped_result =
      expected_result == s2pred.Excluded.FIRST ? s2pred.Excluded.SECOND :
      expected_result == s2pred.Excluded.SECOND ? s2pred.Excluded.FIRST : expected_result;
  if (s2pred.compareDistances(x1, b, a) < 0) {
    Assert.equal(swapped_result, s2pred.getVoronoiSiteExclusion(b, a, x1, x0, r));
  }
}

@("VoronoiSiteExclusion.Coverage")
unittest {
  // Both sites are closest to edge endpoint X0.
  testVoronoiSiteExclusion(
      S2Point(1, -1e-5, 0), S2Point(1, -2e-5, 0),
      S2Point(1, 0, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-3),
      s2pred.Excluded.SECOND, Precision.DOUBLE);

  // Both sites are closest to edge endpoint X1.
  testVoronoiSiteExclusion(
      S2Point(1, 1, 1e-30), S2Point(1, 1, -1e-20),
      S2Point(1, 0, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-10),
      s2pred.Excluded.SECOND, Precision.DOUBLE);

  // Test cases where neither site is excluded.
  testVoronoiSiteExclusion(
      S2Point(1, -1e-10, 1e-5), S2Point(1, 1e-10, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-4),
      s2pred.Excluded.NEITHER, Precision.DOUBLE);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-10, 1e-5), S2Point(1, 1e-10, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-5),
      s2pred.Excluded.NEITHER, Precision.REAL);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-17, 1e-5), S2Point(1, 1e-17, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-4),
      s2pred.Excluded.NEITHER, Precision.REAL);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-20, 1e-5), S2Point(1, 1e-20, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1e-5),
      s2pred.Excluded.NEITHER, Precision.EXACT);

  // Test cases where the first site is excluded.  (Tests where the second
  // site is excluded are constructed by TestVoronoiSiteExclusion.)
  testVoronoiSiteExclusion(
      S2Point(1, -1e-6, 1.0049999999e-5), S2Point(1, 0, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1.005e-5),
      s2pred.Excluded.FIRST, Precision.DOUBLE);
  testVoronoiSiteExclusion(
      S2Point(1, -1.00105e-6, 1.0049999999e-5), S2Point(1, 0, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1.005e-5),
      s2pred.Excluded.FIRST, Precision.REAL);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-6, 1.005e-5), S2Point(1, 0, -1e-5),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1.005e-5),
      s2pred.Excluded.FIRST, Precision.REAL);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-31, 1.005e-30), S2Point(1, 0, -1e-30),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1.005e-30),
      s2pred.Excluded.FIRST, Precision.EXACT);
  testVoronoiSiteExclusion(
      S2Point(1, -1e-31, 1.005e-30), S2Point(1, 0, -1e-30),
      S2Point(1, -1, 0), S2Point(1, 1, 0), S1ChordAngle.fromRadians(1.005e-30),
      s2pred.Excluded.FIRST, Precision.EXACT);

  // These two sites are exactly 60 degrees away from the point (1, 1, 0),
  // which is the midpoint of edge X.  This case requires symbolic
  // perturbations to resolve correctly.  Site A is closer to every point in
  // its coverage interval except for (1, 1, 0), but site B is considered
  // closer to that point symbolically.
  testVoronoiSiteExclusion(
      S2Point(0, 1, 1), S2Point(1, 0, 1),
      S2Point(0, 1, 1), S2Point(1, 0, -1), S1ChordAngle.fromLength2(1),
      s2pred.Excluded.NEITHER, Precision.EXACT);

  // This test is similar except that site A is considered closer to the
  // equidistant point (-1, 1, 0), and therefore site B is excluded.
  testVoronoiSiteExclusion(
      S2Point(0, 1, 1), S2Point(-1, 0, 1),
      S2Point(0, 1, 1), S2Point(-1, 0, -1), S1ChordAngle.fromLength2(1),
      s2pred.Excluded.SECOND, Precision.EXACT);
}

// Checks that the result at one level of precision is consistent with the
// result at the next higher level of precision.  Returns the minimum
// precision that yielded a non-zero result.
Precision testVoronoiSiteExclusionConsistency(
    in S2Point a, in S2Point b, in S2Point x0, in S2Point x1, S1ChordAngle r) {
  s2pred.Excluded UNCERTAIN = s2pred.Excluded.UNCERTAIN;

  // The internal methods require this (see TestVoronoiSiteExclusion).
  if (s2pred.compareDistances(x1, a, b) < 0) return Precision.DOUBLE;

  s2pred.Excluded dbl_result = s2pred.triageVoronoiSiteExclusion(a, b, x0, x1, r.length2());
  s2pred.Excluded r_result = s2pred.triageVoronoiSiteExclusion(
      Vector3_r.from(a), Vector3_r.from(b), Vector3_r.from(x0), Vector3_r.from(x1), r.length2());
  s2pred.Excluded exact_result = s2pred.exactVoronoiSiteExclusion(
      s2pred.Vector3_xf.from(a), s2pred.Vector3_xf.from(b), s2pred.Vector3_xf.from(x0),
      s2pred.Vector3_xf.from(x1), ExactFloat(r.length2()));
  Assert.equal(exact_result, s2pred.getVoronoiSiteExclusion(a, b, x0, x1, r));

  Assert.notEqual(UNCERTAIN, exact_result);
  if (r_result == UNCERTAIN) {
    Assert.equal(UNCERTAIN, dbl_result);
    return Precision.EXACT;
  }
  Assert.equal(exact_result, r_result);
  if (dbl_result == UNCERTAIN) {
    return Precision.REAL;
  }
  Assert.equal(exact_result, dbl_result);
  return Precision.DOUBLE;
}

@("VoronoiSiteExclusion.Consistency")
unittest {
  // This test chooses random a random edge X, a random point P on that edge,
  // and a random threshold distance "r".  It then choose two sites A and B
  // whose distance to P is almost exactly "r".  This ensures that the
  // coverage intervals for A and B will (almost) share a common endpoint.  It
  // then checks that the answer given by a method at one level of precision
  // is consistent with the answer given at higher levels of precision.
  auto rnd = S2Testing.rnd;
  PrecisionStats stats;
  for (int iter = 0; iter < CONSISTENCY_ITERS; ++iter) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point x0 = choosePoint();
    S2Point x1 = choosePoint();
    if (x0 == -x1) continue;  // Not allowed by API.
    double f = math.pow(1e-20, rnd.randDouble());
    S2Point p = ((1 - f) * x0 + f * x1).normalize();
    S1Angle r1 = S1Angle.fromRadians(math.PI_2 * math.pow(1e-20, rnd.randDouble()));
    S2Point a = interpolateAtDistance(r1, p, choosePoint());
    S2Point b = interpolateAtDistance(r1, p, choosePoint());
    // Check that the other API requirements are met.
    S1ChordAngle r = S1ChordAngle(r1);
    if (s2pred.compareEdgeDistance(a, x0, x1, r) > 0) continue;
    if (s2pred.compareEdgeDistance(b, x0, x1, r) > 0) continue;
    if (s2pred.compareDistances(x0, a, b) > 0) algorithm.swap(a, b);
    if (a == b) continue;

    Precision prec = testVoronoiSiteExclusionConsistency(a, b, x0, x1, r);
    // Don't skew the statistics by recording degenerate inputs.
    if (x0 == x1) {
      Assert.equal(Precision.DOUBLE, prec);
    } else {
      stats.tally(prec);
    }
  }
  writeln(stats.toString());
}
