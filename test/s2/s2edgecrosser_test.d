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

module s2.s2edgecrosser_test;

import fluent.asserts;
import s2.s2edgecrosser;
import s2.s2edgecrossings;
import s2.s2edgedistances;
import s2.s2point;
import s2.s2pointutil;
import s2.s2testing;
import math = std.math;

// In non-debug builds, check that default-constructed and/or NaN S2Point
// arguments don't cause crashes, especially on the very first method call
// (since S2CopyingEdgeCrosser checks whether the first vertex of each edge is
// the same as the last vertex of the previous edged when deciding whether or
// not to call Restart).

void testCrossingSignInvalid(in S2Point point, int expected) {
  auto crosser = new S2EdgeCrosser(point, point);
  Assert.equal(expected, crosser.crossingSign(point, point));
  auto crosser2 = new S2CopyingEdgeCrosser(point, point);
  Assert.equal(expected, crosser2.crossingSign(point, point));
}

void testEdgeOrVertexCrossingInvalid(in S2Point point, bool expected) {
  auto crosser = new S2EdgeCrosser(point, point);
  Assert.equal(expected, crosser.edgeOrVertexCrossing(point, point));
  auto crosser2 = new S2CopyingEdgeCrosser(point, point);
  Assert.equal(expected, crosser2.edgeOrVertexCrossing(point, point));
}

@("S2EdgeUtil.InvalidDefaultPoints")
unittest {
  // Check that default-constructed S2Point arguments don't cause crashes.
  auto point = S2Point(0, 0, 0);
  testCrossingSignInvalid(point, 0);
  testEdgeOrVertexCrossingInvalid(point, false);
}

@("S2EdgeUtil.InvalidNanPoints")
unittest {
  // Check that NaN S2Point arguments don't cause crashes.
  auto point = S2Point(double.nan, double.nan, double.nan);
  testCrossingSignInvalid(point, -1);
  testEdgeOrVertexCrossingInvalid(point, false);
}

void testCrossing(
    in S2Point a, in S2Point b, in S2Point c, in S2Point d, int robust, bool edge_or_vertex) {
  // Modify the expected result if two vertices from different edges match.
  if (a == c || a == d || b == c || b == d) robust = 0;
  Assert.equal(robust, crossingSign(a, b, c, d));
  scope auto crosser = new S2EdgeCrosser(a, b, c);
  Assert.equal(robust, crosser.crossingSign(d));
  Assert.equal(robust, crosser.crossingSign(c));
  Assert.equal(robust, crosser.crossingSign(d, c));
  Assert.equal(robust, crosser.crossingSign(c, d));

  Assert.equal(edge_or_vertex, edgeOrVertexCrossing(a, b, c, d));
  crosser.restartAt(c);
  Assert.equal(edge_or_vertex, crosser.edgeOrVertexCrossing(d));
  Assert.equal(edge_or_vertex, crosser.edgeOrVertexCrossing(c));
  Assert.equal(edge_or_vertex, crosser.edgeOrVertexCrossing(d, c));
  Assert.equal(edge_or_vertex, crosser.edgeOrVertexCrossing(c, d));

  // Check that the crosser can be re-used.
  crosser.init(c, d);
  crosser.restartAt(a);
  Assert.equal(robust, crosser.crossingSign(b));
  Assert.equal(robust, crosser.crossingSign(a));

  // Now try all the same tests with CopyingEdgeCrosser.
  scope auto crosser2 = new S2CopyingEdgeCrosser(a, b, c);
  Assert.equal(robust, crosser2.crossingSign(d));
  Assert.equal(robust, crosser2.crossingSign(c));
  Assert.equal(robust, crosser2.crossingSign(d, c));
  Assert.equal(robust, crosser2.crossingSign(c, d));

  Assert.equal(edge_or_vertex, edgeOrVertexCrossing(a, b, c, d));
  crosser2.restartAt(c);
  Assert.equal(edge_or_vertex, crosser2.edgeOrVertexCrossing(d));
  Assert.equal(edge_or_vertex, crosser2.edgeOrVertexCrossing(c));
  Assert.equal(edge_or_vertex, crosser2.edgeOrVertexCrossing(d, c));
  Assert.equal(edge_or_vertex, crosser2.edgeOrVertexCrossing(c, d));

  // Check that the crosser can be re-used.
  crosser2.init(c, d);
  crosser2.restartAt(a);
  Assert.equal(robust, crosser2.crossingSign(b));
  Assert.equal(robust, crosser2.crossingSign(a));
}

void testCrossings(
    S2Point a, S2Point b, S2Point c, S2Point d, int robust, bool edge_or_vertex) {
  a = a.normalize();
  b = b.normalize();
  c = c.normalize();
  d = d.normalize();
  testCrossing(a, b, c, d, robust, edge_or_vertex);
  testCrossing(b, a, c, d, robust, edge_or_vertex);
  testCrossing(a, b, d, c, robust, edge_or_vertex);
  testCrossing(b, a, d, c, robust, edge_or_vertex);
  testCrossing(a, a, c, d, -1, false);
  testCrossing(a, b, c, c, -1, false);
  testCrossing(a, a, c, c, -1, false);
  testCrossing(a, b, a, b, 0, true);
  testCrossing(c, d, a, b, robust, edge_or_vertex != (robust == 0));
}

@("S2EdgeUtil.Crossings")
unittest {
  // The real tests of edge crossings are in s2{loop,polygon}_test,
  // but we do a few simple tests here.

  // Two regular edges that cross.
  testCrossings(
      S2Point(1, 2, 1), S2Point(1, -3, 0.5), S2Point(1, -0.5, -3), S2Point(0.1, 0.5, 3), 1, true);

  // Two regular edges that intersect antipodal points.
  testCrossings(
      S2Point(1, 2, 1), S2Point(1, -3, 0.5), S2Point(-1, 0.5, 3), S2Point(-0.1, -0.5, -3),
      -1, false);

  // Two edges on the same great circle that start at antipodal points.
  testCrossings(
      S2Point(0, 0, -1), S2Point(0, 1, 0), S2Point(0, 0, 1), S2Point(0, 1, 1), -1, false);

  // Two edges that cross where one vertex is S2::Origin().
  testCrossings(
      S2Point(1, 0, 0), origin(), S2Point(1, -0.1, 1), S2Point(1, 1, -0.1), 1, true);

  // Two edges that intersect antipodal points where one vertex is
  // S2::Origin().
  testCrossings(
      S2Point(1, 0, 0), origin(), S2Point(-1, 0.1, -1), S2Point(-1, -1, 0.1), -1, false);

  // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
  // and edge CD is further CCW around (2,3,4) than AB.
  testCrossings(
      S2Point(2, 3, 4), S2Point(-1, 2, 5), S2Point(7, -2, 3), S2Point(2, 3, 4), 0, false);

  // Two edges that barely cross each other near the middle of one edge.  The
  // edge AB is approximately in the x=y plane, while CD is approximately
  // perpendicular to it and ends exactly at the x=y plane.
  testCrossings(
      S2Point(1, 1, 1), S2Point(1, math.nextafter(1.0, 0.0), -1),
      S2Point(11, -12, -1), S2Point(10, 10, 1), 1, true);

  // In this version, the edges are separated by a distance of about 1e-15.
  testCrossings(
      S2Point(1, 1, 1), S2Point(1, math.nextafter(1.0, 2.0), -1),
      S2Point(1, -1, 0), S2Point(1, 1, 0), -1, false);

  // Two edges that barely cross each other near the end of both edges.  This
  // example cannot be handled using regular double-precision arithmetic due
  // to floating-point underflow.
  testCrossings(
      S2Point(0, 0, 1), S2Point(2, -double.min_normal, 1),
      S2Point(1, -1, 1), S2Point(double.min_normal, 0, 1), 1, true);

  // In this version, the edges are separated by a distance of about 1e-640.
  testCrossings(
      S2Point(0, 0, 1), S2Point(2, double.min_normal, 1),
      S2Point(1, -1, 1), S2Point(double.min_normal, 0, 1), -1, false);

  // Two edges that barely cross each other near the middle of one edge.
  // Computing the exact determinant of some of the triangles in this test
  // requires more than 2000 bits of precision.
  testCrossings(
      S2Point(1, -double.min_normal, -double.min_normal),
      S2Point(double.min_normal, 1, double.min_normal),
      S2Point(1, -1, double.min_normal), S2Point(1, 1, 0),
      1, true);

  // In this version, the edges are separated by a distance of about 1e-640.
  testCrossings(
      S2Point(1, double.min_normal, -double.min_normal),
      S2Point(-double.min_normal, 1, double.min_normal),
      S2Point(1, -1, double.min_normal), S2Point(1, 1, 0), -1, false);
}

@("S2EdgeUtil.CollinearEdgesThatDontTouch")
unittest {
  const int kIters = 500;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Point a = S2Testing.randomPoint();
    S2Point d = S2Testing.randomPoint();
    S2Point b = interpolate(0.05, a, d);
    S2Point c = interpolate(0.95, a, d);
    Assert.greaterThan(0, crossingSign(a, b, c, d));
    Assert.greaterThan(0, crossingSign(a, b, c, d));
    auto crosser = new S2EdgeCrosser(a, b, c);
    Assert.greaterThan(0, crosser.crossingSign(d));
    Assert.greaterThan(0, crosser.crossingSign(c));
  }
}

@("S2EdgeUtil.CoincidentZeroLengthEdgesThatDontTouch")
unittest {
  // It is important that the edge primitives can handle vertices that exactly
  // exactly proportional to each other, i.e. that are not identical but are
  // nevertheless exactly coincident when projected onto the unit sphere.
  // There are various ways that such points can arise.  For example,
  // Normalize() itself is not idempotent: there exist distinct points A,B
  // such that Normalize(A) == B  and Normalize(B) == A.  Another issue is
  // that sometimes calls to Normalize() are skipped when the result of a
  // calculation "should" be unit length mathematically (e.g., when computing
  // the cross product of two orthonormal vectors).
  //
  // This test checks pairs of edges AB and CD where A,B,C,D are exactly
  // coincident on the sphere and the norms of A,B,C,D are monotonically
  // increasing.  Such edge pairs should never intersect.  (This is not
  // obvious, since it depends on the particular symbolic perturbations used
  // by s2pred::Sign().  It would be better to replace this with a test that
  // says that the CCW results must be consistent with each other.)
  const int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    // Construct a point P where every component is zero or a power of 2.
    S2Point p;
    for (int i = 0; i < 3; ++i) {
      int binary_exp = S2Testing.rnd.skewed(11);
      p[i] = (binary_exp > 1022) ? 0 : math.pow(2.0, -binary_exp);
    }
    // If all components were zero, try again.  Note that normalization may
    // convert a non-zero point into a zero one due to underflow (!)
    p = p.normalize();
    if (p == S2Point(0, 0, 0)) { --iter; continue; }

    // Now every non-zero component should have exactly the same mantissa.
    // This implies that if we scale the point by an arbitrary factor, every
    // non-zero component will still have the same mantissa.  Scale the points
    // so that they are all distinct and are still very likely to satisfy
    // S2::IsUnitLength (which allows for a small amount of error in the norm).

    S2Point a = (1-3e-16) * p;
    S2Point b = (1-1e-16) * p;
    S2Point c = p;
    S2Point d = (1+2e-16) * p;
    if (!isUnitLength(a) || !isUnitLength(d)) {
      --iter;
      continue;
    }
    // Verify that the expected edges do not cross.
    Assert.greaterThan(0, crossingSign(a, b, c, d));
    auto crosser = new S2EdgeCrosser(a, b, c);
    Assert.greaterThan(0, crosser.crossingSign(d));
    Assert.greaterThan(0, crosser.crossingSign(c));
  }
}
