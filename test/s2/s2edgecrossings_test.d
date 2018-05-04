module s2.s2edgecrossings_test;

// Original author: ericv@google.com (Eric Veach)

import fluent.asserts;
import s2.s2point;
import s2.s1angle;
import s2.s2edgecrossings;
import s2.s2edgedistances;
import s2.s2predicates;
import s2.s2testing;
import s2.s2edgedistances;
import s2.util.math.exactfloat;
import s2.util.math.vector;
import std.stdio;
import algorithm = std.algorithm;
import math = std.math;
import traits = std.traits;

// CrossingSign, VertexCrossing, and EdgeOrVertexCrossing are tested in
// s2edge_crosser_test.cc.

// This class records statistics about the intersection methods used by
// GetIntersection().
class GetIntersectionStats {
private:
  int[IntersectionMethod] _tally;

public:

  this() {
    foreach (method; traits.EnumMembers!IntersectionMethod) {
      _tally[method] = 0;
    }
    intersectionMethodTally = &_tally;
  }

  ~this() {
    intersectionMethodTally = null;
  }

  void print() {
    int total = 0;
    int[IntersectionMethod] totals;
    foreach (method; traits.EnumMembers!IntersectionMethod) {
      total += _tally[method];
      totals[method] = total;
    }
    writefln("%10s %16s %16s  %6s", "Method", "Successes", "Attempts", "Rate");
    foreach (method; traits.EnumMembers!IntersectionMethod) {
      if (_tally[method] == 0) continue;
      writefln("%10s %9d %5.1f%% %9d %5.1f%%  %5.1f%%",
          method,
          _tally[method], 100.0 * _tally[method] / total,
          totals[method], 100.0 * totals[method] / total,
          100.0 * _tally[method] / totals[method]);
    }
    foreach (method; traits.EnumMembers!IntersectionMethod) {
      _tally[method] = 0;
    }
  }
}

// This returns the true intersection point of two line segments (a0,a1) and
// (b0,b1), with a relative error of at most DBL_EPSILON in each coordinate
// (i.e., one ulp, or twice the double precision rounding error).
private S2Point getFixedIntersectionExact(
    in S2Point a0, in S2Point a1, in S2Point b0, in S2Point b1) {
  S2Point x = getIntersectionExact(a0, a1, b0, b1);
  if (x.dotProd((a0 + a1) + (b0 + b1)) < 0) x = -x;
  return x;
}

// The approximate maximum error in GetDistance() for small distances.
immutable S1Angle GET_DISTANCE_ABS_ERROR = S1Angle.fromRadians(3 * double.epsilon);

@("S2EdgeUtil.IntersectionError")
unittest {
  // We repeatedly construct two edges that cross near a random point "p", and
  // measure the distance from the actual intersection point "x" to the
  // exact intersection point and also to the edges.

  GetIntersectionStats stats = new GetIntersectionStats();
  S1Angle max_point_dist, max_edge_dist;
  auto rnd = S2Testing.rnd;
  for (int iter = 0; iter < 5000; ++iter) {
    // We construct two edges AB and CD that intersect near "p".  The angle
    // between AB and CD (expressed as a slope) is chosen randomly between
    // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
    // Similarly, two edge lengths approximately between 1e-15 and 1 are
    // chosen.  The edge endpoints are chosen such that they are often very
    // close to the other edge (i.e., barely crossing).  Taken together this
    // ensures that we test both long and very short edges that intersect at
    // both large and very small angles.
    //
    // Sometimes the edges we generate will not actually cross, in which case
    // we simply try again.
    Vector3_d p, d1, d2;
    S2Testing.getRandomFrame(p, d1, d2);
    double slope = 1e-15 * math.pow(1e30, rnd.randDouble());
    d2 = (d1 + slope * d2).normalize();
    S2Point a, b, c, d;
    do {
      double ab_len = math.pow(1e-15, rnd.randDouble());
      double cd_len = math.pow(1e-15, rnd.randDouble());
      double a_fraction = math.pow(1e-5, rnd.randDouble());
      if (rnd.oneIn(2)) a_fraction = 1 - a_fraction;
      double c_fraction = math.pow(1e-5, rnd.randDouble());
      if (rnd.oneIn(2)) c_fraction = 1 - c_fraction;
      a = (p - a_fraction * ab_len * d1).normalize();
      b = (p + (1 - a_fraction) * ab_len * d1).normalize();
      c = (p - c_fraction * cd_len * d2).normalize();
      d = (p + (1 - c_fraction) * cd_len * d2).normalize();
    } while (crossingSign(a, b, c, d) <= 0);

    // Each constructed edge should be at most 1.5 * DBL_EPSILON away from the
    // original point P.
    Assert.notGreaterThan(
        getDistance(p, a, b),
        S1Angle.fromRadians(1.5 * double.epsilon) + GET_DISTANCE_ABS_ERROR);
    Assert.notGreaterThan(
        getDistance(p, c, d),
        S1Angle.fromRadians(1.5 * double.epsilon) + GET_DISTANCE_ABS_ERROR);

    // Verify that the expected intersection point is close to both edges and
    // also close to the original point P.  (It might not be very close to P
    // if the angle between the edges is very small.)
    S2Point expected = getFixedIntersectionExact(a, b, c, d);
    Assert.notGreaterThan(
        getDistance(expected, a, b),
        S1Angle.fromRadians(3 * double.epsilon) + GET_DISTANCE_ABS_ERROR);
    Assert.notGreaterThan(
        getDistance(expected, c, d),
        S1Angle.fromRadians(3 * double.epsilon) + GET_DISTANCE_ABS_ERROR);
    Assert.notGreaterThan(
        S1Angle(expected, p),
        S1Angle.fromRadians(3 * double.epsilon / slope) + INTERSECTION_ERROR);

    // Now we actually test the GetIntersection() method.
    S2Point actual = getIntersection(a, b, c, d);
    S1Angle dist_ab = getDistance(actual, a, b);
    S1Angle dist_cd = getDistance(actual, c, d);
    Assert.notGreaterThan(dist_ab, INTERSECTION_ERROR + GET_DISTANCE_ABS_ERROR);
    Assert.notGreaterThan(dist_cd, INTERSECTION_ERROR + GET_DISTANCE_ABS_ERROR);
    max_edge_dist = algorithm.max(max_edge_dist, algorithm.max(dist_ab, dist_cd));
    S1Angle point_dist = S1Angle(expected, actual);
    Assert.notGreaterThan(point_dist, INTERSECTION_ERROR);
    max_point_dist = algorithm.max(max_point_dist, point_dist);
  }
  stats.print();
  writeln("Max distance to either edge being intersected: ", max_edge_dist.radians());
  writeln("Maximum distance to expected intersection point: ", max_point_dist.radians());
}

/+
// Chooses a point in the XY plane that is separated from X by at least 1e-15
// (to avoid choosing too many duplicate points) and by at most Pi/2 - 1e-3
// (to avoid nearly-diametric edges, since the test below is not sophisticated
// enough to test such edges).
static S2Point ChooseSemicirclePoint(const S2Point& x, const S2Point& y) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  double sign = (2 * rnd->Uniform(2)) - 1;
  return (x + sign * 1e3 * pow(1e-18, rnd->RandDouble()) * y).Normalize();
}

TEST(S2EdgeUtil, GrazingIntersections) {
  // This test choose 5 points along a great circle (i.e., as collinear as
  // possible), and uses them to construct an edge AB and a triangle CDE such
  // that CD and CE both cross AB.  It then checks that the intersection
  // points returned by GetIntersection() have the correct relative ordering
  // along AB (to within kIntersectionError).
  GetIntersectionStats stats;
  for (int iter = 0; iter < 1000; ++iter) {
    Vector3_d x, y, z;
    S2Testing::GetRandomFrame(&x, &y, &z);
    S2Point a, b, c, d, e, ab;
    do {
      a = ChooseSemicirclePoint(x, y);
      b = ChooseSemicirclePoint(x, y);
      c = ChooseSemicirclePoint(x, y);
      d = ChooseSemicirclePoint(x, y);
      e = ChooseSemicirclePoint(x, y);
      ab = (a - b).CrossProd(a + b);
    } while (ab.Norm() < 50 * DBL_EPSILON ||
             S2::CrossingSign(a, b, c, d) <= 0 ||
             S2::CrossingSign(a, b, c, e) <= 0);
    S2Point xcd = S2::GetIntersection(a, b, c, d);
    S2Point xce = S2::GetIntersection(a, b, c, e);
    // Essentially this says that if CDE and CAB have the same orientation,
    // then CD and CE should intersect along AB in that order.
    ab = ab.Normalize();
    if (S1Angle(xcd, xce) > 2 * S2::kIntersectionError) {
      EXPECT_EQ(s2pred::Sign(c, d, e) == s2pred::Sign(c, a, b),
                s2pred::Sign(ab, xcd, xce) > 0);
    }
  }
  stats.Print();
}

TEST(S2EdgeUtil, ExactIntersectionUnderflow) {
  // Tests that a correct intersection is computed even when two edges are
  // exactly collinear and the normals of both edges underflow in double
  // precision when normalized (see S2PointFromExact function for details).
  S2Point a0(1, 0, 0), a1(1, 2e-300, 0);
  S2Point b0(1, 1e-300, 0), b1(1, 3e-300, 0);
  EXPECT_EQ(S2Point(1, 1e-300, 0), S2::GetIntersection(a0, a1, b0, b1));
}

TEST(S2EdgeUtil, GetIntersectionInvariants) {
  // Test that the result of GetIntersection does not change when the edges
  // are swapped and/or reversed.  The number of iterations is high because it
  // is difficult to generate test cases that show that CompareEdges() is
  // necessary and correct, for example.
  const int kIters = google::DEBUG_MODE ? 5000 : 50000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Point a, b, c, d;
    do {
      // GetIntersectionStable() sorts the two edges by length, so construct
      // edges (a,b) and (c,d) that cross and have exactly the same length.
      // This can be done by swapping the "x" and "y" coordinates.
      // [Swapping other coordinate pairs doesn't work because it changes the
      // order of addition in Norm2() == (x**2 + y**2) + z**2.]
      a = c = S2Testing::RandomPoint();
      b = d = S2Testing::RandomPoint();
      swap(c[0], c[1]);
      swap(d[0], d[1]);
    } while (S2::CrossingSign(a, b, c, d) <= 0);
    EXPECT_EQ((a - b).Norm2(), (c - d).Norm2());

    // Now verify that GetIntersection returns exactly the same result when
    // the edges are swapped and/or reversed.
    S2Point result = S2::GetIntersection(a, b, c, d);
    if (S2Testing::rnd.OneIn(2)) { swap(a, b); }
    if (S2Testing::rnd.OneIn(2)) { swap(c, d); }
    if (S2Testing::rnd.OneIn(2)) { swap(a, c); swap(b, d); }
    EXPECT_EQ(result, S2::GetIntersection(a, b, c, d));
  }
}
+/
