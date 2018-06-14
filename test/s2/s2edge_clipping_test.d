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

module s2.s2edge_clipping_test;

import fluent.asserts;
import s2.s2edge_clipping;
import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2coords;
import s2.s2point;
import s2.s2pointutil;
import s2.s2testing;

import std.algorithm : max;
import std.format : format;
import std.math : atan2, fabs, nextafter, pow;
import std.stdio;

enum double DOUBLE_ERR = 0.001;

void testFaceClipping(in S2Point a_raw, in S2Point b_raw) {
  S2Point a = a_raw.normalize();
  S2Point b = b_raw.normalize();
  // TODO(ericv): Remove the following line once S2::RobustCrossProd is
  // extended to use simulation of simplicity.
  if (a == -b) return;

  // First we test GetFaceSegments.
  FaceSegmentVector segments;
  getFaceSegments(a, b, segments);
  auto n = segments.length;
  Assert.notLessThan(n, 1);

  string msg;
  msg ~= format("\nA=%s\nB=%s", a_raw, b_raw);
  msg ~= format("\nN=%s\nSegments:\n", robustCrossProd(a, b));
  int i = 0;
  foreach (const FaceSegment s; segments) {
    msg ~= format("%d: face=%d, a=%s, b=%s\n", i++, s.face, s.a, s.b);
  }
  writeln(msg);

  auto biunit = R2Rect(R1Interval(-1, 1), R1Interval(-1, 1));
  const double kErrorRadians = FACE_CLIP_ERROR_RADIANS;

  // The first and last vertices should approximately equal A and B.
  Assert.notGreaterThan(a.angle(FaceUVtoXYZ(segments[0].face, segments[0].a)), kErrorRadians);
  Assert.notGreaterThan(b.angle(FaceUVtoXYZ(segments[n-1].face, segments[n-1].b)), kErrorRadians);

  S2Point norm = robustCrossProd(a, b).normalize();
  S2Point a_tangent = norm.crossProd(a);
  S2Point b_tangent = b.crossProd(norm);
  for (i = 0; i < n; ++i) {
    // Vertices may not protrude outside the biunit square.
    Assert.equal(biunit.contains(segments[i].a), true);
    Assert.equal(biunit.contains(segments[i].b), true);
    if (i == 0) continue;

    // The two representations of each interior vertex (on adjacent faces)
    // must correspond to exactly the same S2Point.
    Assert.notEqual(segments[i-1].face, segments[i].face);
    Assert.equal(
        FaceUVtoXYZ(segments[i-1].face, segments[i-1].b),
        FaceUVtoXYZ(segments[i].face, segments[i].a));

    // Interior vertices should be in the plane containing A and B, and should
    // be contained in the wedge of angles between A and B (i.e., the dot
    // products with a_tangent and b_tangent should be non-negative).
    S2Point p = FaceUVtoXYZ(segments[i].face, segments[i].a).normalize();
    Assert.notGreaterThan(fabs(p.dotProd(norm)), kErrorRadians);
    Assert.notLessThan(p.dotProd(a_tangent), -kErrorRadians);
    Assert.notLessThan(p.dotProd(b_tangent), -kErrorRadians);
  }

  // Now we test ClipToPaddedFace (sometimes with a padding of zero).  We do
  // this by defining an (x,y) coordinate system for the plane containing AB,
  // and converting points along the great circle AB to angles in the range
  // [-Pi, Pi].  We then accumulate the angle intervals spanned by each
  // clipped edge; the union over all 6 faces should approximately equal the
  // interval covered by the original edge.
  auto rnd = S2Testing.rnd;
  double padding = rnd.oneIn(10) ? 0.0 : 1e-10 * pow(1e-5, rnd.randDouble());
  S2Point x_axis = a, y_axis = a_tangent;
  auto expected_angles = S1Interval(0, a.angle(b));
  S1Interval max_angles = expected_angles.expanded(kErrorRadians);
  S1Interval actual_angles;
  for (int face = 0; face < 6; ++face) {
    R2Point a_uv, b_uv;
    if (clipToPaddedFace(a, b, face, padding, a_uv, b_uv)) {
      S2Point a_clip = FaceUVtoXYZ(face, a_uv).normalize();
      S2Point b_clip = FaceUVtoXYZ(face, b_uv).normalize();
      Assert.notGreaterThan(fabs(a_clip.dotProd(norm)), kErrorRadians);
      Assert.notGreaterThan(fabs(b_clip.dotProd(norm)), kErrorRadians);
      if (a_clip.angle(a) > kErrorRadians) {
        Assert.approximately(1 + padding, max(fabs(a_uv[0]), fabs(a_uv[1])), DOUBLE_ERR);
      }
      if (b_clip.angle(b) > kErrorRadians) {
        Assert.approximately(1 + padding, max(fabs(b_uv[0]), fabs(b_uv[1])), DOUBLE_ERR);
      }
      double a_angle = atan2(a_clip.dotProd(y_axis), a_clip.dotProd(x_axis));
      double b_angle = atan2(b_clip.dotProd(y_axis), b_clip.dotProd(x_axis));
      // Rounding errors may cause b_angle to be slightly less than a_angle.
      // We handle this by constructing the interval with FromPointPair(),
      // which is okay since the interval length is much less than M_PI.
      S1Interval face_angles = S1Interval.fromPointPair(a_angle, b_angle);
      Assert.equal(max_angles.contains(face_angles), true);
      actual_angles = actual_angles.unite(face_angles);
    }
  }
  Assert.equal(actual_angles.expanded(kErrorRadians).contains(expected_angles), true);
}

void testFaceClippingEdgePair(in S2Point a, in S2Point b) {
  testFaceClipping(a, b);
  testFaceClipping(b, a);
}

// This function is designed to choose line segment endpoints that are
// difficult to handle correctly.  Given two adjacent cube vertices P and Q,
// it returns either an edge midpoint, face midpoint, or corner vertex that is
// in the plane of PQ and that has been perturbed slightly.  It also sometimes
// returns a random point from anywhere on the sphere.
S2Point perturbedCornerOrMidpoint(in S2Point p, in S2Point q) {
  auto rnd = S2Testing.rnd;
  S2Point a = (rnd.uniform(3) - 1) * p + (rnd.uniform(3) - 1) * q;
  if (rnd.oneIn(10)) {
    // This perturbation often has no effect except on coordinates that are
    // zero, in which case the perturbed value is so small that operations on
    // it often result in underflow.
    a += pow(1e-300, rnd.randDouble()) * S2Testing.randomPoint();
  } else if (rnd.oneIn(2)) {
    // For coordinates near 1 (say > 0.5), this perturbation yields values
    // that are only a few representable values away from the initial value.
    a += 4 * double.epsilon * S2Testing.randomPoint();
  } else {
    // A perturbation whose magnitude is in the range [1e-25, 1e-10].
    a += 1e-10 * pow(1e-15, rnd.randDouble()) * S2Testing.randomPoint();
  }
  if (a.norm2() < double.min_normal) {
    // If a.Norm2() is denormalized, Normalize() loses too much precision.
    return perturbedCornerOrMidpoint(p, q);
  }
  return a;
}

@("S2EdgeUtil.FaceClipping")
unittest {
  // Start with a few simple cases.
  // An edge that is entirely contained within one cube face:
  testFaceClippingEdgePair(S2Point(1, -0.5, -0.5), S2Point(1, 0.5, 0.5));
  // An edge that crosses one cube edge:
  testFaceClippingEdgePair(S2Point(1, 0, 0), S2Point(0, 1, 0));
  // An edge that crosses two opposite edges of face 0:
  testFaceClippingEdgePair(S2Point(0.75, 0, -1), S2Point(0.75, 0, 1));
  // An edge that crosses two adjacent edges of face 2:
  testFaceClippingEdgePair(S2Point(1, 0, 0.75), S2Point(0, 1, 0.75));
  // An edges that crosses three cube edges (four faces):
  testFaceClippingEdgePair(S2Point(1, 0.9, 0.95), S2Point(-1, 0.95, 0.9));

  // Comprehensively test edges that are difficult to handle, especially those
  // that nearly follow one of the 12 cube edges.
  auto rnd = S2Testing.rnd;
  auto biunit = R2Rect(R1Interval(-1, 1), R1Interval(-1, 1));
  const int kIters = 1000;  // Test passes with 1e6 iterations
  for (int iter = 0; iter < kIters; ++iter) {
    writeln("Iteration ", iter);
    // Choose two adjacent cube corners P and Q.
    int face = rnd.uniform(6);
    int i = rnd.uniform(4);
    int j = (i + 1) & 3;
    S2Point p = FaceUVtoXYZ(face, biunit.getVertex(i));
    S2Point q = FaceUVtoXYZ(face, biunit.getVertex(j));

    // Now choose two points that are nearly in the plane of PQ, preferring
    // points that are near cube corners, face midpoints, or edge midpoints.
    S2Point a = perturbedCornerOrMidpoint(p, q);
    S2Point b = perturbedCornerOrMidpoint(p, q);
    testFaceClipping(a, b);
  }
}

// Choose a random point in the rectangle defined by points A and B, sometimes
// returning a point on the edge AB or the points A and B themselves.
R2Point chooseRectPoint(in R2Point a, in R2Point b) {
  auto rnd = S2Testing.rnd;
  if (rnd.oneIn(5)) {
    return rnd.oneIn(2) ? a : b;
  } else if (rnd.oneIn(3)) {
    return a + rnd.randDouble() * (b - a);
  } else {
    return R2Point(rnd.uniformDouble(a[0], b[0]), rnd.uniformDouble(a[1], b[1]));
  }
}

// Given a point X on the line AB (which is checked), return the fraction "t"
// such that x = (1-t)*a + t*b.  Return 0 if A = B.
double getFraction(in R2Point x, in R2Point a, in R2Point b) {
  // A bound for the error in edge clipping plus the error in the calculation
  // below (which is similar to IntersectsRect).
  const double kError = EDGE_CLIP_ERROR_UV_DIST + INTERSECTS_RECT_ERROR_UV_DIST;
  if (a == b) return 0.0;
  R2Point dir = (b - a).normalize();
  Assert.notGreaterThan(fabs((x - a).dotProd(dir.ortho())), kError);
  return (x - a).dotProd(dir);
}

// Given a point P representing a possibly clipped endpoint A of an edge AB,
// verify that "clip" contains P, and that if clipping occurred (i.e., P != A)
// then P is on the boundary of "clip".
void checkPointOnBoundary(in R2Point p, in R2Point a, in R2Rect clip) {
  Assert.equal(clip.contains(p), true);
  if (p != a) {
    Assert.equal(clip.contains(R2Point(nextafter(p[0], a[0]), nextafter(p[1], a[1]))), false);
  }
}

// Given an edge AB and a rectangle "clip", verify that IntersectsRect(),
// ClipEdge(), and ClipEdgeBound() produce consistent results.
void testClipEdge(in R2Point a, in R2Point b, in R2Rect clip) {
  // A bound for the error in edge clipping plus the error in the
  // IntersectsRect calculation below.
  const double kError = EDGE_CLIP_ERROR_UV_DIST + INTERSECTS_RECT_ERROR_UV_DIST;
  R2Point a_clipped, b_clipped;
  if (!clipEdge(a, b, clip, a_clipped, b_clipped)) {
    Assert.equal(intersectsRect(a, b, clip.expanded(-kError)), false);
  } else {
    Assert.equal(intersectsRect(a, b, clip.expanded(kError)), true);
    // Check that the clipped points lie on the edge AB, and that the points
    // have the expected order along the segment AB.
    Assert.notGreaterThan(getFraction(a_clipped, a, b), getFraction(b_clipped, a, b));
    // Check that the clipped portion of AB is as large as possible.
    checkPointOnBoundary(a_clipped, a, clip);
    checkPointOnBoundary(b_clipped, b, clip);
  }
  // Choose a random initial bound to pass to ClipEdgeBound.
  R2Rect initial_clip = R2Rect.fromPointPair(chooseRectPoint(a, b), chooseRectPoint(a, b));
  R2Rect bound = getClippedEdgeBound(a, b, initial_clip);
  if (bound.isEmpty()) return;  // Precondition of ClipEdgeBound not met
  R2Rect max_bound = bound.intersection(clip);
  if (!clipEdgeBound(a, b, clip, bound)) {
    Assert.equal(intersectsRect(a, b, max_bound.expanded(-kError)), false);
  } else {
    Assert.equal(intersectsRect(a, b, max_bound.expanded(kError)), true);
    // Check that the bound is as large as possible.
    int ai = (a[0] > b[0]), aj = (a[1] > b[1]);
    checkPointOnBoundary(bound.getVertex(ai, aj), a, max_bound);
    checkPointOnBoundary(bound.getVertex(1-ai, 1-aj), b, max_bound);
  }
}

// Given an interval "clip", randomly choose either a value in the interval, a
// value outside the interval, or one of the two interval endpoints, ensuring
// that all cases have reasonable probability for any interval "clip".
double chooseEndpoint(const R1Interval clip) {
  auto rnd = S2Testing.rnd;
  if (rnd.oneIn(5)) {
    return rnd.oneIn(2) ? clip.lo() : clip.hi();
  } else {
    switch (rnd.uniform(3)) {
      case 0:  return clip.lo() - rnd.randDouble();
      case 1:  return clip.hi() + rnd.randDouble();
      default: return clip.lo() + rnd.randDouble() * clip.getLength();
    }
  }
}

// Given a rectangle "clip", choose a point that may lie in the rectangle
// interior, along an extended edge, exactly at a vertex, or in one of the
// eight regions exterior to "clip" that are separated by its extended edges.
// Also sometimes return points that are exactly on one of the extended
// diagonals of "clip".  All cases are reasonably likely to occur for any
// given rectangle "clip".
R2Point chooseEndpoint(in R2Rect clip) {
  if (S2Testing.rnd.oneIn(10)) {
    // Return a point on one of the two extended diagonals.
    int diag = S2Testing.rnd.uniform(2);
    double t = S2Testing.rnd.uniformDouble(-1, 2);
    return (1 - t) * clip.getVertex(diag) + t * clip.getVertex(diag + 2);
  } else {
    return R2Point(chooseEndpoint(clip[0]), chooseEndpoint(clip[1]));
  }
}

// Given a rectangle "clip", test the S2EdgeUtil edge clipping methods using
// many edges that are randomly constructed to trigger special cases.
void testEdgeClipping(in R2Rect clip) {
  const int kIters = 1000;  // Test passes with 1e6 iterations
  for (int iter = 0; iter < kIters; ++iter) {
    writeln("Iteration ", iter);
    testClipEdge(chooseEndpoint(clip), chooseEndpoint(clip), clip);
  }
}

@("S2EdgeUtil.EdgeClipping")
unittest {
  auto rnd = S2Testing.rnd;
  // Test clipping against random rectangles.
  for (int i = 0; i < 5; ++i) {
    testEdgeClipping(R2Rect.fromPointPair(
        R2Point(rnd.uniformDouble(-1, 1), rnd.uniformDouble(-1, 1)),
        R2Point(rnd.uniformDouble(-1, 1), rnd.uniformDouble(-1, 1))));
  }
  // Also clip against one-dimensional, singleton, and empty rectangles.
  testEdgeClipping(R2Rect(R1Interval(-0.7, -0.7), R1Interval(0.3, 0.35)));
  testEdgeClipping(R2Rect(R1Interval(0.2, 0.5), R1Interval(0.3, 0.3)));
  testEdgeClipping(R2Rect(R1Interval(-0.7, 0.3), R1Interval(0, 0)));
  testEdgeClipping(R2Rect.fromPoint(R2Point(0.3, 0.8)));
  testEdgeClipping(R2Rect.empty());
}
