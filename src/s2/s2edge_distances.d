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

module s2.s2edge_distances;

// Defines a collection of functions for computing the distance to an edge,
// interpolating along an edge, projecting points onto edges, etc.

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2point;
import s2pointutil = s2.s2pointutil;
import s2.util.math.vector;
import algorithm = std.algorithm;
import math = std.math;

/////////////////////////////////////////////////////////////////////////////
///////////////            (point, edge) functions            ///////////////

// Returns the minimum distance from X to any point on the edge AB.  All
// arguments should be unit length.  The result is very accurate for small
// distances but may have some numerical error if the distance is large
// (approximately Pi/2 or greater).  The case A == B is handled correctly.
//
// If you want to compare a distance against a fixed threshold, e.g.
//    if (S2::GetDistance(x, a, b) < limit)
// then it is significantly faster to use UpdateMinDistance() below.
S1Angle getDistance(in S2Point x, in S2Point a, in S2Point b) {
  S1ChordAngle min_dist;
  alwaysUpdateMinDistance!(true)(x, a, b, min_dist);
  return min_dist.toS1Angle();
}

// Returns true if the distance from X to the edge AB is less than "limit".
// (Specify limit.Successor() for "less than or equal to".)  This method is
// significantly faster than GetDistance().  If you want to compare against a
// fixed S1Angle, you should convert it to an S1ChordAngle once and save the
// value, since this step is relatively expensive.
//
// See s2pred::CompareEdgeDistance() for an exact version of this predicate.
//bool IsDistanceLess(const S2Point& x, const S2Point& a, const S2Point& b,
//                    S1ChordAngle limit);

// If the distance from X to the edge AB is less than "min_dist", this
// method updates "min_dist" and returns true.  Otherwise it returns false.
// The case A == B is handled correctly.
//
// Use this method when you want to compute many distances and keep track of
// the minimum.  It is significantly faster than using GetDistance(),
// because (1) using S1ChordAngle is much faster than S1Angle, and (2) it
// can save a lot of work by not actually computing the distance when it is
// obviously larger than the current minimum.
bool updateMinDistance(in S2Point x, in S2Point a, in S2Point b, ref S1ChordAngle min_dist) {
  return alwaysUpdateMinDistance!(false)(x, a, b, min_dist);
}

// This function computes the distance from a point X to a line segment AB.
// If the distance is less than "min_dist" or "always_update" is true, it
// updates "min_dist" and returns true.  Otherwise it returns false.
//
// The "Always" in the function name refers to the template argument, i.e.
// AlwaysUpdateMinDistance<true> always updates the given distance, while
// AlwaysUpdateMinDistance<false> does not.  This optimization increases the
// speed of GetDistance() by about 10% without creating code duplication.
bool alwaysUpdateMinDistance(bool alwaysUpdate)(
    in S2Point x, in S2Point a, in S2Point b, ref S1ChordAngle min_dist)
in {
  assert(s2pointutil.isUnitLength(x) && s2pointutil.isUnitLength(a)
      && s2pointutil.isUnitLength(b));
} body {
  double xa2 = (x-a).norm2();
  double xb2 = (x-b).norm2();
  if (alwaysUpdateMinInteriorDistance!(alwaysUpdate)(x, a, b, xa2, xb2, min_dist)) {
    return true;  // Minimum distance is attained along the edge interior.
  }
  // Otherwise the minimum distance is to one of the endpoints.
  double dist2 = algorithm.min(xa2, xb2);
  if (!alwaysUpdate && dist2 >= min_dist.length2()) {
    return false;
  }
  min_dist = S1ChordAngle.fromLength2(dist2);
  return true;
}

// If the distance from X to the edge AB is greater than "max_dist", this
// method updates "max_dist" and returns true.  Otherwise it returns false.
// The case A == B is handled correctly.
bool updateMaxDistance(in S2Point x, in S2Point a, in S2Point b, out S1ChordAngle max_dist) {
  auto dist = algorithm.max(S1ChordAngle(x, a), S1ChordAngle(x, b));
  if (dist > S1ChordAngle.right()) {
    alwaysUpdateMinDistance!true(-x, a, b, dist);
    dist = S1ChordAngle.straight() - dist;
  }
  if (max_dist < dist) {
    max_dist = dist;
    return true;
  }

  return false;
}

// Returns the maximum error in the result of UpdateMinDistance (and
// associated functions such as UpdateMinInteriorDistance, IsDistanceLess,
// etc), assuming that all input points are normalized to within the bounds
// guaranteed by S2Point::Normalize().  The error can be added or subtracted
// from an S1ChordAngle "x" using x.PlusError(error).
//
// Note that accuracy goes down as the distance approaches 0 degrees or 180
// degrees (for different reasons).  Near 0 degrees the error is acceptable
// for all practical purposes (about 1.2e-15 radians ~= 8 nanometers).  For
// exactly antipodal points the maximum error is quite high (0.5 meters), but
// this error drops rapidly as the points move away from antipodality
// (approximately 1 millimeter for points that are 50 meters from antipodal,
// and 1 micrometer for points that are 50km from antipodal).
//
// TODO(ericv): Currently the error bound does not hold for edges whose
// endpoints are antipodal to within about 1e-15 radians (less than 1 micron).
// This could be fixed by extending S2::RobustCrossProd to use higher
// precision when necessary.
//double GetUpdateMinDistanceMaxError(S1ChordAngle dist);

// Returns true if the minimum distance from X to the edge AB is attained at
// an interior point of AB (i.e., not an endpoint), and that distance is less
// than "limit".  (Specify limit.Successor() for "less than or equal to".)
//bool IsInteriorDistanceLess(const S2Point& x,
//                            const S2Point& a, const S2Point& b,
//                            S1ChordAngle limit);

// If the minimum distance from X to AB is attained at an interior point of AB
// (i.e., not an endpoint), and that distance is less than "min_dist", then
// this method updates "min_dist" and returns true.  Otherwise returns false.
//bool UpdateMinInteriorDistance(const S2Point& x,
//                               const S2Point& a, const S2Point& b,
//                               S1ChordAngle* min_dist);

// Returns the point along the edge AB that is closest to the point X.
// The fractional distance of this point along the edge AB can be obtained
// using GetDistanceFraction() above.  Requires that all vectors have
// unit length.
//S2Point Project(const S2Point& x, const S2Point& a, const S2Point& b);

// A slightly more efficient version of Project() where the cross product of
// the two endpoints has been precomputed.  The cross product does not need to
// be normalized, but should be computed using S2::RobustCrossProd() for the
// most accurate results.  Requires that x, a, and b have unit length.
//S2Point Project(const S2Point& x, const S2Point& a, const S2Point& b,
//                const Vector3_d& a_cross_b);


/////////////////////////////////////////////////////////////////////////////
///////////////         (point along edge) functions          ///////////////


// Given a point X and an edge AB, returns the distance ratio AX / (AX + BX).
// If X happens to be on the line segment AB, this is the fraction "t" such
// that X == Interpolate(t, A, B).  Requires that A and B are distinct.
//double GetDistanceFraction(const S2Point& x,
//                           const S2Point& a, const S2Point& b);

// Returns the point X along the line segment AB whose distance from A is the
// given fraction "t" of the distance AB.  Does NOT require that "t" be
// between 0 and 1.  Note that all distances are measured on the surface of
// the sphere, so this is more complicated than just computing (1-t)*a + t*b
// and normalizing the result.
S2Point interpolate(double t, in S2Point a, in S2Point b) {
  if (t == 0) return a;
  if (t == 1) return b;
  auto ab = S1Angle(a, b);
  return interpolateAtDistance(t * ab, a, b);
}

// Like Interpolate(), except that the parameter "ax" represents the desired
// distance from A to the result X rather than a fraction between 0 and 1.
S2Point interpolateAtDistance(S1Angle ax_angle, in S2Point a, in S2Point b)
in {
  assert(s2pointutil.isUnitLength(a));
  assert(s2pointutil.isUnitLength(b));
} body {
  double ax = ax_angle.radians();

  // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
  // result is always perpendicular to A, even if A=B or A=-B, but it is not
  // necessarily unit length.  (We effectively normalize it below.)
  Vector3_d normal = s2pointutil.robustCrossProd(a, b);
  Vector3_d tangent = normal.crossProd(a);
  assert(tangent != S2Point(0, 0, 0));

  // Now compute the appropriate linear combination of A and "tangent".  With
  // infinite precision the result would always be unit length, but we
  // normalize it anyway to ensure that the error is within acceptable bounds.
  // (Otherwise errors can build up when the result of one interpolation is
  // fed into another interpolation.)
  return (math.cos(ax) * a + (math.sin(ax) / tangent.norm()) * tangent).normalize();
}


/////////////////////////////////////////////////////////////////////////////
///////////////            (edge, edge) functions             ///////////////


// Like UpdateMinDistance(), but computes the minimum distance between the
// given pair of edges.  (If the two edges cross, the distance is zero.)
// The cases a0 == a1 and b0 == b1 are handled correctly.
// bool UpdateEdgePairMinDistance(const S2Point& a0, const S2Point& a1,
//                                const S2Point& b0, const S2Point& b1,
//                                S1ChordAngle* min_dist);

// As above, but for maximum distances. If one edge crosses the antipodal
// reflection of the other, the distance is Pi.
// bool UpdateEdgePairMaxDistance(const S2Point& a0, const S2Point& a1,
//                                const S2Point& b0, const S2Point& b1,
//                                S1ChordAngle* max_dist);

// Returns the pair of points (a, b) that achieves the minimum distance
// between edges a0a1 and b0b1, where "a" is a point on a0a1 and "b" is a
// point on b0b1.  If the two edges intersect, "a" and "b" are both equal to
// the intersection point.  Handles a0 == a1 and b0 == b1 correctly.
// std::pair<S2Point, S2Point> GetEdgePairClosestPoints(
//     const S2Point& a0, const S2Point& a1,
//     const S2Point& b0, const S2Point& b1);

// Returns true if every point on edge B=b0b1 is no further than "tolerance"
// from some point on edge A=a0a1.  Equivalently, returns true if the directed
// Hausdorff distance from B to A is no more than "tolerance".
// Requires that tolerance is less than 90 degrees.
// bool IsEdgeBNearEdgeA(const S2Point& a0, const S2Point& a1,
//                       const S2Point& b0, const S2Point& b1,
//                       S1Angle tolerance);


//////////////////   Implementation details follow   ////////////////////


// inline bool IsDistanceLess(const S2Point& x, const S2Point& a,
//                            const S2Point& b, S1ChordAngle limit) {
//   return UpdateMinDistance(x, a, b, &limit);
// }

// inline bool IsInteriorDistanceLess(const S2Point& x, const S2Point& a,
//                                    const S2Point& b, S1ChordAngle limit) {
//   return UpdateMinInteriorDistance(x, a, b, &limit);
// }

// If the minimum distance from X to AB is attained at an interior point of AB
// (i.e., not an endpoint), and that distance is less than "min_dist" or
// "always_update" is true, then update "min_dist" and return true.  Otherwise
// return false.
//
// The "Always" in the function name refers to the template argument, i.e.
// AlwaysUpdateMinInteriorDistance<true> always updates the given distance,
// while AlwaysUpdateMinInteriorDistance<false> does not.  This optimization
// increases the speed of GetDistance() by about 10% without creating code
// duplication.
bool alwaysUpdateMinInteriorDistance(bool alwaysUpdate)(
    in S2Point x, in S2Point a, in S2Point b,
    in double xa2, in double xb2, ref S1ChordAngle min_dist)
in {
  assert(s2pointutil.isUnitLength(x) && s2pointutil.isUnitLength(a)
      && s2pointutil.isUnitLength(b));
  assert(xa2 == (x-a).norm2());
  assert(xb2 == (x-b).norm2());
} body {

  // The closest point on AB could either be one of the two vertices (the
  // "vertex case") or in the interior (the "interior case").  Let C = A x B.
  // If X is in the spherical wedge extending from A to B around the axis
  // through C, then we are in the interior case.  Otherwise we are in the
  // vertex case.
  //
  // Check whether we might be in the interior case.  For this to be true, XAB
  // and XBA must both be acute angles.  Checking this condition exactly is
  // expensive, so instead we consider the planar triangle ABX (which passes
  // through the sphere's interior).  The planar angles XAB and XBA are always
  // less than the corresponding spherical angles, so if we are in the
  // interior case then both of these angles must be acute.
  //
  // We check this by computing the squared edge lengths of the planar
  // triangle ABX, and testing acuteness using the law of cosines:
  //
  //             max(XA^2, XB^2) < min(XA^2, XB^2) + AB^2
  //
  if (algorithm.max(xa2, xb2) >= algorithm.min(xa2, xb2) + (a-b).norm2()) {
    return false;
  }
  // The minimum distance might be to a point on the edge interior.  Let R
  // be closest point to X that lies on the great circle through AB.  Rather
  // than computing the geodesic distance along the surface of the sphere,
  // instead we compute the "chord length" through the sphere's interior.
  // If the squared chord length exceeds min_dist.length2() then we can
  // return "false" immediately.
  //
  // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
  // is the point X projected onto the plane through the great circle AB.
  // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
  // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
  // is faster and the corresponding distance on the Earth's surface is
  // accurate to within 1% for distances up to about 1800km.
  S2Point c = s2pointutil.robustCrossProd(a, b);
  double c2 = c.norm2();
  double x_dot_c = x.dotProd(c);
  double x_dot_c2 = x_dot_c * x_dot_c;
  if (!alwaysUpdate && x_dot_c2 >= c2 * min_dist.length2()) {
    // The closest point on the great circle AB is too far away.
    return false;
  }
  // Otherwise we do the exact, more expensive test for the interior case.
  // This test is very likely to succeed because of the conservative planar
  // test we did initially.
  S2Point cx = c.crossProd(x);
  if (a.dotProd(cx) >= 0 || b.dotProd(cx) <= 0) {
    return false;
  }
  // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
  // This calculation has good accuracy for all chord lengths since it
  // is based on both the dot product and cross product (rather than
  // deriving one from the other).  However, note that the chord length
  // representation itself loses accuracy as the angle approaches Pi.
  double qr = 1 - math.sqrt(cx.norm2() / c2);
  double dist2 = (x_dot_c2 / c2) + (qr * qr);
  if (!alwaysUpdate && dist2 >= min_dist.length2()) {
    return false;
  }
  min_dist = S1ChordAngle.fromLength2(dist2);
  return true;
}
