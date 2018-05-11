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

module s2.s2edgecrossings;

// Defines functions related to determining whether two geodesic edges cross
// and for computing intersection points.
//
// The predicates CrossingSign(), VertexCrossing(), and EdgeOrVertexCrossing()
// are robust, meaning that they produce correct, consistent results even in
// pathological cases.  See s2predicates.h for additional robust predicates.
//
// See also S2EdgeCrosser (which efficiently tests an edge against a sequence
// of other edges) and S2CrossingEdgeQuery (which uses an index to speed up
// the process).

import s2.s1angle;
import s2.s1chordangle;
import s2.s1interval;
import s2.s2latlng;
import s2.s2latlng;
import s2.s2point;
import s2.s2pointutil;
import s2.s2predicates;
import s2.s2edgecrosser;
import s2.util.math.exactfloat;
import s2.util.math.vector;
import std.stdio;
import algorithm = std.algorithm;
import math = std.math;
import traits = std.traits;

package int[IntersectionMethod]* intersectionMethodTally;

/**
 * This function determines whether the edge AB intersects the edge CD.
 * Returns +1 if AB crosses CD at a point that is interior to both edges.
 * Returns  0 if any two vertices from different edges are the same.
 * Returns -1 otherwise.
 *
 * Note that if an edge is degenerate (A == B or C == D), the return value
 * is 0 if two vertices from different edges are the same and -1 otherwise.
 *
 * Properties of CrossingSign:
 *
 *  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
 *  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
 *  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
 *  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
 *
 * This function implements an exact, consistent perturbation model such
 * that no three points are ever considered to be collinear.  This means
 * that even if you have 4 points A, B, C, D that lie exactly in a line
 * (say, around the equator), C and D will be treated as being slightly to
 * one side or the other of AB.  This is done in a way such that the
 * results are always consistent (see s2pred::Sign).
 *
 * Note that if you want to check an edge against a collection of other edges,
 * it is much more efficient to use an S2EdgeCrosser (see s2edge_crosser.h).
 */
int crossingSign(in S2Point a, in S2Point b, in S2Point c, in S2Point d) {
  S2EdgeCrosser crosser = new S2EdgeCrosser(a, b, c);
  return crosser.crossingSign(d);
}

/**
 * Given two edges AB and CD where at least two vertices are identical
 * (i.e. CrossingSign(a,b,c,d) == 0), this function defines whether the
 * two edges "cross" in a such a way that point-in-polygon containment tests
 * can be implemented by counting the number of edge crossings.  The basic
 * rule is that a "crossing" occurs if AB is encountered after CD during a
 * CCW sweep around the shared vertex starting from a fixed reference point.
 *
 * Note that according to this rule, if AB crosses CD then in general CD
 * does not cross AB.  However, this leads to the correct result when
 * counting polygon edge crossings.  For example, suppose that A,B,C are
 * three consecutive vertices of a CCW polygon.  If we now consider the edge
 * crossings of a segment BP as P sweeps around B, the crossing number
 * changes parity exactly when BP crosses BA or BC.
 *
 * Useful properties of VertexCrossing (VC):
 *
 *  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
 *  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
 *  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
 *  (3) If exactly one of a,b equals one of c,d, then exactly one of
 *      VC(a,b,c,d) and VC(c,d,a,b) is true
 *
 * It is an error to call this method with 4 distinct vertices.
 */
bool vertexCrossing(in S2Point a, in S2Point b, in S2Point c, in S2Point d) {
  // If A == B or C == D there is no intersection.  We need to check this
  // case first in case 3 or more input points are identical.
  if (a == b || c == d) {
    return false;
  }

  // If any other pair of vertices is equal, there is a crossing if and only
  // if OrderedCCW() indicates that the edge AB is further CCW around the
  // shared vertex O (either A or B) than the edge CD, starting from an
  // arbitrary fixed reference point.
  //
  // Optimization: if AB=CD or AB=DC, we can avoid most of the calculations.
  if (a == c) {
    return (b == d) || orderedCCW(ortho(a), d, b, a);
  }
  if (b == d) {
    return orderedCCW(ortho(b), c, a, b);
  }

  if (a == d) {
    return (b == c) || orderedCCW(ortho(a), c, b, a);
  }
  if (b == c) {
    return orderedCCW(ortho(b), d, a, b);
  }

  return false;
}


/**
 * A convenience function that calls CrossingSign() to handle cases
 * where all four vertices are distinct, and VertexCrossing() to handle
 * cases where two or more vertices are the same.  This defines a crossing
 * function such that point-in-polygon containment tests can be implemented
 * by simply counting edge crossings.
 */
bool edgeOrVertexCrossing(in S2Point a, in S2Point b, in S2Point c, in S2Point d) {
  int crossing = crossingSign(a, b, c, d);
  if (crossing < 0) {
    return false;
  }
  if (crossing > 0) {
    return true;
  }
  return vertexCrossing(a, b, c, d);
}

/**
 * Returns whether (a0,a1) is less than (b0,b1) with respect to a total
 * ordering on edges that is invariant under edge reversals.
 */
private bool compareEdges(T)(
    in Vector!(T, 3) a0, in Vector!(T, 3) a1, in Vector!(T, 3) b0, in Vector!(T, 3) b1)
if (traits.isFloatingPoint!T) {
  const(Vector3!(T))* pa0 = &a0;
  const(Vector3!(T))* pa1 = &a1;
  const(Vector3!(T))* pb0 = &b0;
  const(Vector3!(T))* pb1 = &b1;
  if (*pa0 >= *pa1) algorithm.swap(pa0, pa1);
  if (*pb0 >= *pb1) algorithm.swap(pb0, pb1);
  return *pa0 < *pb0 || (*pa0 == *pb0 && *pb0 < *pb1);
}

/**
 * If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
 * to within an error of at most INTERSECTION_ERROR by this function, then set
 * "result" to the intersection point and return true.
 *
 * The intersection point is not guaranteed to have the correct sign
 * (i.e., it may be either "result" or "-result").
 */
private bool getIntersectionStable(T)(
    in Vector!(T, 3) a0, in Vector!(T, 3) a1,
    in Vector!(T, 3) b0, in Vector!(T, 3) b1,
    out Vector!(T, 3) result)
if (traits.isFloatingPoint!T) {
  // Sort the two edges so that (a0,a1) is longer, breaking ties in a
  // deterministic way that does not depend on the ordering of the endpoints.
  // This is desirable for two reasons:
  //  - So that the result doesn't change when edges are swapped or reversed.
  //  - It reduces error, since the first edge is used to compute the edge
  //    normal (where a longer edge means less error), and the second edge
  //    is used for interpolation (where a shorter edge means less error).
  T a_len2 = (a1 - a0).norm2();
  T b_len2 = (b1 - b0).norm2();
  if (a_len2 < b_len2 || (a_len2 == b_len2 && compareEdges(a0, a1, b0, b1))) {
    return getIntersectionStableSorted(b0, b1, a0, a1, result);
  } else {
    return getIntersectionStableSorted(a0, a1, b0, b1, result);
  }
}

/**
 * Given a point X and a vector "a_norm" (not necessarily unit length),
 * compute x.DotProd(a_norm) and return a bound on the error in the result.
 * The remaining parameters allow this dot product to be computed more
 * accurately and efficiently.  They include the length of "a_norm"
 * ("a_norm_len") and the edge endpoints "a0" and "a1".
 */
private T getProjection(T)(in Vector!(T, 3) x, in Vector!(T, 3) a_norm, in T a_norm_len,
    in Vector!(T, 3) a0, in Vector!(T, 3) a1, out T error)
if (traits.isFloatingPoint!T) {
  // The error in the dot product is proportional to the lengths of the input
  // vectors, so rather than using "x" itself (a unit-length vector) we use
  // the vectors from "x" to the closer of the two edge endpoints.  This
  // typically reduces the error by a huge factor.
  Vector3!T x0 = x - a0;
  Vector3!T x1 = x - a1;
  T x0_dist2 = x0.norm2();
  T x1_dist2 = x1.norm2();

  // If both distances are the same, we need to be careful to choose one
  // endpoint deterministically so that the result does not change if the
  // order of the endpoints is reversed.
  T dist, result;
  if (x0_dist2 < x1_dist2 || (x0_dist2 == x1_dist2 && x0 < x1)) {
    dist = math.sqrt(x0_dist2);
    result = x0.dotProd(a_norm);
  } else {
    dist = math.sqrt(x1_dist2);
    result = x1.dotProd(a_norm);
  }
  // This calculation bounds the error from all sources: the computation of
  // the normal, the subtraction of one endpoint, and the dot product itself.
  // (DBL_ERR appears because the input points are assumed to be normalized in
  // double precision rather than in the given type T.)
  //
  // For reference, the bounds that went into this calculation are:
  // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * DBL_ERR) * T_ERR
  // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * T_ERR
  // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * T_ERR
  T T_ERR = roundingEpsilon!T();
  error = (((3.5 + 2 * math.sqrt(3.0)) * a_norm_len + 32 * math.sqrt(3.0) * DBL_ERR)
            * dist + 1.5 * math.abs(result)) * T_ERR;
  return result;
}

/**
 * Helper function for GetIntersectionStable().  It expects that the edges
 * (a0,a1) and (b0,b1) have been sorted so that the first edge is longer.
 */
private bool getIntersectionStableSorted(T)(
    in Vector!(T, 3) a0, in Vector!(T, 3) a1,
    in Vector!(T, 3) b0, in Vector!(T, 3) b1,
    out Vector!(T, 3) result)
in {
  assert((a1 - a0).norm2() >= (b1 - b0).norm2());
} body {

  // Compute the normal of the plane through (a0, a1) in a stable way.
  Vector3!T a_norm = (a0 - a1).crossProd(a0 + a1);
  T a_norm_len = a_norm.norm();
  T b_len = (b1 - b0).norm();

  // Compute the projection (i.e., signed distance) of b0 and b1 onto the
  // plane through (a0, a1).  Distances are scaled by the length of a_norm.
  T b0_error, b1_error;
  T b0_dist = getProjection(b0, a_norm, a_norm_len, a0, a1, b0_error);
  T b1_dist = getProjection(b1, a_norm, a_norm_len, a0, a1, b1_error);

  // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
  // |b0_dist - b1_dist|.  Note that b0_dist and b1_dist generally have
  // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
  // code below finds the intersection point by interpolating along the edge
  // (b0, b1) to a fractional distance of b0_dist / (b0_dist - b1_dist).
  //
  // It can be shown that the maximum error in the interpolation fraction is
  //
  //     (b0_dist * b1_error - b1_dist * b0_error) /
  //        (dist_sum * (dist_sum - error_sum))
  //
  // We save ourselves some work by scaling the result and the error bound by
  // "dist_sum", since the result is normalized to be unit length anyway.
  T dist_sum = math.abs(b0_dist - b1_dist);
  T error_sum = b0_error + b1_error;
  if (dist_sum <= error_sum) {
    return false;  // Error is unbounded in this case.
  }
  Vector3!T x = b0_dist * b1 - b1_dist * b0;
  T T_ERR = roundingEpsilon!T();
  T error = b_len * math.abs(b0_dist * b1_error - b1_dist * b0_error) /
      (dist_sum - error_sum) + 2 * T_ERR * dist_sum;

  // Finally we normalize the result, compute the corresponding error, and
  // check whether the total error is acceptable.
  T x_len = x.norm();
  const T kMaxError = INTERSECTION_ERROR.radians();
  if (error > (kMaxError - T_ERR) * x_len) {
    return false;
  }
  result = (1 / x_len) * x;
  return true;
}

static bool getIntersectionStableR(
    in S2Point a0, in S2Point a1, in S2Point b0, in S2Point b1, out S2Point result) {
  Vector3_r result_r;
  if (getIntersectionStable(
          Vector3_r.from(a0), Vector3_r.from(a1),
          Vector3_r.from(b0), Vector3_r.from(b1),
          result_r)) {
    result = S2Point.from(result_r);
    return true;
  }
  return false;
}

package enum IntersectionMethod {
  SIMPLE,
  SIMPLE_R,
  STABLE,
  STABLE_R,
  EXACT,
  NUM_METHODS
}

/**
 * Given three points "a", "x", "b", returns true if these three points occur
 * in the given order along the edge (a,b) to within the given tolerance.
 * More precisely, either "x" must be within "tolerance" of "a" or "b", or
 * when "x" is projected onto the great circle through "a" and "b" it must lie
 * along the edge (a,b) (i.e., the shortest path from "a" to "b").
 */
private bool approximatelyOrdered(
    in S2Point a, in S2Point x, in S2Point b, double tolerance) {
  if ((x - a).norm2() <= tolerance * tolerance) {
    return true;
  }
  if ((x - b).norm2() <= tolerance * tolerance) {
    return true;
  }
  return orderedCCW(a, x, b, robustCrossProd(a, b).normalize());
}

/**
 * Given two edges AB and CD such that CrossingSign(A, B, C, D) > 0, returns
 * their intersection point.  Useful properties of GetIntersection (GI):
 *
 *  (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
 *  (2) GI(c,d,a,b) == GI(a,b,c,d)
 *
 * The returned intersection point X is guaranteed to be very close to the
 * true intersection point of AB and CD, even if the edges intersect at a
 * very small angle.  See "INTERSECTION_ERROR" below for details.
 */
S2Point getIntersection(in S2Point a0, in S2Point a1, in S2Point b0, in S2Point b1)
in {
  assert(crossingSign(a0, a1, b0, b1) > 0);
} out (result) {
  // Make sure that the intersection point lies on both edges.
  assert(approximatelyOrdered(a0, result, a1, INTERSECTION_ERROR.radians()));
  assert(approximatelyOrdered(b0, result, b1, INTERSECTION_ERROR.radians()));
} body {

  // It is difficult to compute the intersection point of two edges accurately
  // when the angle between the edges is very small.  Previously we handled
  // this by only guaranteeing that the returned intersection point is within
  // INTERSECTION_ERROR of each edge.  However, this means that when the edges
  // cross at a very small angle, the computed result may be very far from the
  // true intersection point.
  //
  // Instead this function now guarantees that the result is always within
  // INTERSECTION_ERROR of the true intersection.  This requires using more
  // sophisticated techniques and in some cases extended precision.
  //
  // Three different techniques are implemented, but only two are used:
  //
  //  - GetIntersectionSimple() computes the intersection point using
  //    numerically stable cross products in "long double" precision.
  //
  //  - GetIntersectionStable() computes the intersection point using
  //    projection and interpolation, taking care to minimize cancellation
  //    error.  This method exists in "double" and "long double" versions.
  //
  //  - GetIntersectionExact() computes the intersection point using exact
  //    arithmetic and converts the final result back to an S2Point.
  //
  // We don't actually use the first method (GetIntersectionSimple) because it
  // turns out that GetIntersectionStable() is twice as fast and also much
  // more accurate (even in double precision).  The "long double" version
  // (only available on Intel platforms) uses 80-bit precision and is about
  // twice as slow.  The exact arithmetic version is about 100x slower.
  //
  // So our strategy is to first call GetIntersectionStable() in double
  // precision; if that doesn't work and this platform supports "long double",
  // then we try again in "long double"; if that doesn't work then we fall
  // back to exact arithmetic.

  S2Point result;
  IntersectionMethod method;
  if (getIntersectionStable(a0, a1, b0, b1, result)) {
    method = IntersectionMethod.STABLE;
  } else if (getIntersectionStableR(a0, a1, b0, b1, result)) {
    method = IntersectionMethod.STABLE_R;
  } else {
    result = getIntersectionExact(a0, a1, b0, b1);
    method = IntersectionMethod.EXACT;
    method = IntersectionMethod.STABLE_R;
  }

  if (intersectionMethodTally) {
    ++(*intersectionMethodTally)[method];
  }

  // Make sure the intersection point is on the correct side of the sphere.
  // Since all vertices are unit length, and edges are less than 180 degrees,
  // (a0 + a1) and (b0 + b1) both have positive dot product with the
  // intersection point.  We use the sum of all vertices to make sure that the
  // result is unchanged when the edges are swapped or reversed.
  if (result.dotProd((a0 + a1) + (b0 + b1)) < 0) result = -result;

  return result;
}

/**
 * INTERSECTION_ERROR is an upper bound on the distance from the intersection
 * point returned by GetIntersection() to the true intersection point.
 */
immutable S1Angle INTERSECTION_ERROR = S1Angle.fromRadians(8 * DBL_ERR);

/**
 * This value can be used as the S2Builder snap_radius() to ensure that edges
 * that have been displaced by up to INTERSECTION_ERROR are merged back
 * together again.  For example this can happen when geometry is intersected
 * with a set of tiles and then unioned.  It is equal to twice the
 * intersection error because input edges might have been displaced in
 * opposite directions.
 */
immutable S1Angle INTERSECTION_MERGE_RADIUS = 2 * INTERSECTION_ERROR;


// Compute the intersection point of (a0, a1) and (b0, b1) using exact
// arithmetic.  Note that the result is not exact because it is rounded to
// double precision.  Also, the intersection point is not guaranteed to have
// the correct sign (i.e., the return value may need to be negated).
package S2Point getIntersectionExact(
    in S2Point a0, in S2Point a1, in S2Point b0, in S2Point b1)
out(x) {
  assert(isUnitLength(x));
} body {
  // Since we are using exact arithmetic, we don't need to worry about
  // numerical stability.
  Vector3_xf a0_xf = Vector3_xf.from(a0);
  Vector3_xf a1_xf = Vector3_xf.from(a1);
  Vector3_xf b0_xf = Vector3_xf.from(b0);
  Vector3_xf b1_xf = Vector3_xf.from(b1);
  Vector3_xf a_norm_xf = a0_xf.crossProd(a1_xf);
  Vector3_xf b_norm_xf = b0_xf.crossProd(b1_xf);
  Vector3_xf x_xf = a_norm_xf.crossProd(b_norm_xf);

  // The final Normalize() call is done in double precision, which creates a
  // directional error of up to 2 * DBL_ERR.  (ToDouble() and Normalize() each
  // contribute up to DBL_ERR of directional error.)
  S2Point x = S2PointFromExact(x_xf);

  if (x == S2Point(0, 0, 0)) {
    // The two edges are exactly collinear, but we still consider them to be
    // "crossing" because of simulation of simplicity.  Out of the four
    // endpoints, exactly two lie in the interior of the other edge.  Of
    // those two we return the one that is lexicographically smallest.
    x = S2Point(10, 10, 10);  // Greater than any valid S2Point
    S2Point a_norm = S2PointFromExact(a_norm_xf);
    S2Point b_norm = S2PointFromExact(b_norm_xf);
    if (a_norm == S2Point(0, 0, 0) || b_norm == S2Point(0, 0, 0)) {
      // TODO(ericv): To support antipodal edges properly, we would need to
      // add an s2pred::CrossProd() function that computes the cross product
      // using simulation of simplicity and rounds the result to the nearest
      // floating-point representation.
      writeln("Exactly antipodal edges not supported by GetIntersection");
    }
    if (orderedCCW(b0, a0, b1, b_norm) && a0 < x) x = a0;
    if (orderedCCW(b0, a1, b1, b_norm) && a1 < x) x = a1;
    if (orderedCCW(a0, b0, a1, a_norm) && b0 < x) x = b0;
    if (orderedCCW(a0, b1, a1, a_norm) && b1 < x) x = b1;
  }
  return x;
}

private S2Point S2PointFromExact(in Vector3_xf xf) {
  // If all components of "x" have absolute value less than about 1e-154,
  // then x.Norm2() is zero in double precision due to underflow.  Therefore
  // we need to scale "x" by an appropriate power of 2 before the conversion.
  S2Point x = S2Point(xf[0].toDouble(), xf[1].toDouble(), xf[2].toDouble());
  if (x.norm2() > 0) return x.normalize();

  // Scale so that the largest component magnitude is in the range [0.5, 1).
  int exp = ExactFloat.MIN_EXP - 1;
  for (int i = 0; i < 3; ++i) {
    if (xf[i].isNormal()) exp = algorithm.max(exp, xf[i].exp());
  }
  if (exp < ExactFloat.MIN_EXP) {
    return S2Point(0, 0, 0);
  }
  return S2Point(ldexp(xf[0], -exp).toDouble(),
                 ldexp(xf[1], -exp).toDouble(),
                 ldexp(xf[2], -exp).toDouble()).normalize();
}

