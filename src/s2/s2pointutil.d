module s2.s2pointutil;

// Original author: ericv@google.com (Eric Veach)
//
// Defines additional operations for points on the unit sphere (in addition to
// the standard vector operations defined in "util/math/vector.h").

import s2.s1angle;
import s2.s2point;
import s2.util.math.matrix3x3;
import s2.util.math.vector;

import math = std.math;

// Uncomment the following line for testing purposes only.
// version = S2_TEST_DEGENERACIES;

// Return a unique "origin" on the sphere for operations that need a fixed
// reference point.  In particular, this is the "point at infinity" used for
// point-in-polygon testing (by counting the number of edge crossings).
S2Point origin() {
  version(S2_TEST_DEGENERACIES) {
    // This value makes polygon operations much slower, because it greatly
    // increases the number of degenerate cases that need to be handled using
    // s2pred::ExpensiveSign().
    return S2Point(0, 0, 1);
  } else {
    // The origin should not be a point that is commonly used in edge tests in
    // order to avoid triggering code to handle degenerate cases.  (This rules
    // out the north and south poles.)  It should also not be on the boundary of
    // any low-level S2Cell for the same reason.
    //
    // The point chosen here is about 66km from the north pole towards the East
    // Siberian Sea.  See the unittest for more details.  It is written out
    // explicitly using floating-point literals because the optimizer doesn't
    // seem willing to evaluate Normalize() at compile time.
    return S2Point(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195);
  }
}

// Return true if the given point is approximately unit length
// (this is mainly useful for assertions).
bool isUnitLength(in S2Point p) {
  // Normalize() is guaranteed to return a vector whose L2-norm differs from 1
  // by less than 2 * DBL_EPSILON.  Thus the squared L2-norm differs by less
  // than 4 * DBL_EPSILON.  The actual calculated Norm2() can have up to 1.5 *
  // DBL_EPSILON of additional error.  The total error of 5.5 * DBL_EPSILON
  // can then be rounded down since the result must be a representable
  // double-precision value.
  return math.fabs(p.norm2() - 1) <= 5 * double.epsilon;  // About 1.11e-15
}


// Return true if two points are within the given distance of each other
// (this is mainly useful for testing).
bool approxEquals(in S2Point a, in S2Point b, S1Angle max_error = S1Angle.fromRadians(1e-15)) {
  return S1Angle(a, b) <= max_error;
}

// Return a unit-length vector that is orthogonal to "a".  Satisfies
// Ortho(-a) = -Ortho(a) for all a.
//
// Note that Vector3_d also defines an "Ortho" method, but this one is
// preferred for use in S2 code because it explicitly tries to avoid result
// result coordinates that are zero.  (This is a performance optimization that
// reduces the amount of time spent in functions which handle degeneracies.)
S2Point ortho(in S2Point a) {
  version (S2_TEST_DEGENERACIES) {
    // Vector3.ortho() always returns a point on the X-Y, Y-Z, or X-Z planes.
    // This leads to many more degenerate cases in polygon operations.
    return a.ortho();
  } else {
    int k = a.largestAbsComponent() - 1;
    if (k < 0) {
      k = 2;
    }
    S2Point temp = S2Point(0.012, 0.0053, 0.00457);
    temp[k] = 1;
    return a.crossProd(temp).normalize();
  }
}

// Return a vector "c" that is orthogonal to the given unit-length vectors
// "a" and "b".  This function is similar to a.CrossProd(b) except that it
// does a better job of ensuring orthogonality when "a" is nearly parallel
// to "b", and it returns a non-zero result even when a == b or a == -b.
//
// It satisfies the following properties (RCP == RobustCrossProd):
//
//   (1) RCP(a,b) != 0 for all a, b
//   (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b
//   (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
//   (4) RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
//
// The result is not guaranteed to be unit length.
S2Point robustCrossProd(in S2Point a, in S2Point b)
in {
  assert(isUnitLength(a));
  assert(isUnitLength(b));
} body {
  // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
  // approaches zero.  This leads to situations where a.CrossProd(b) is not
  // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
  // but we also want b.RobustCrossProd(a) == -a.RobustCrossProd(b).
  //
  // The easiest fix is to just compute the cross product of (b+a) and (b-a).
  // Mathematically, this cross product is exactly twice the cross product of
  // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
  // are always perpendicular (since "a" and "b" are unit length).  This
  // yields a result that is nearly orthogonal to both "a" and "b" even if
  // these two values differ only in the lowest bit of one component.
  Vector3_d x = (b + a).crossProd(b - a);
  if (x != S2Point(0, 0, 0)) {
    return x;
  }

  // The only result that makes sense mathematically is to return zero, but
  // we find it more convenient to return an arbitrary orthogonal vector.
  return ortho(a);
}

// Rotate the given point about the given axis by the given angle.  "p" and
// "axis" must be unit length; "angle" has no restrictions (e.g., it can be
// positive, negative, greater than 360 degrees, etc).
S2Point rotate(in S2Point p, in S2Point axis, in S1Angle angle)
in {
  assert(isUnitLength(p));
  assert(isUnitLength(axis));
} body {
  // Let M be the plane through P that is perpendicular to "axis", and let
  // "center" be the point where M intersects "axis".  We construct a
  // right-handed orthogonal frame (dx, dy, center) such that "dx" is the
  // vector from "center" to P, and "dy" has the same length as "dx".  The
  // result can then be expressed as (cos(angle)*dx + sin(angle)*dy + center).
  S2Point center = p.dotProd(axis) * axis;
  S2Point dx = p - center;
  S2Point dy = axis.crossProd(p);
  // Mathematically the result is unit length, but normalization is necessary
  // to ensure that numerical errors don't accumulate.
  return (angle.cos() * dx + angle.sin() * dy + center).normalize();
}

// Extend the given point "z" on the unit sphere into a right-handed
// coordinate frame of unit-length column vectors m = (x,y,z).  Note that the
// vectors (x,y) are an orthonormal frame for the tangent space at "z", while
// "z" itself is an orthonormal frame for the normal space at "z".
Matrix3x3_d getFrame(in S2Point z) {
  Matrix3x3_d m;
  getFrame(z, m);
  return m;
}

void getFrame(in S2Point z, out Matrix3x3_d m)
in {
  assert(isUnitLength(z));
} body {
  m.setCol(2, z);
  m.setCol(1, ortho(z));
  m.setCol(0, m.col(1).crossProd(z));  // Already unit-length.
}

// Given an orthonormal basis "m" of column vectors and a point "p", return
// the coordinates of "p" with respect to the basis "m".  The resulting
// point "q" satisfies the identity (m * q == p).
S2Point toFrame(in Matrix3x3_d m, in S2Point p) {
  // The inverse of an orthonormal matrix is its transpose.
  return m.transpose() * p;
}

// Given an orthonormal basis "m" of column vectors and a point "q" with
// respect to that basis, return the equivalent point "p" with respect to
// the standard axis-aligned basis.  The result satisfies (p == m * q).
S2Point fromFrame(in Matrix3x3_d m, in S2Point q) {
  return m * q;
}

// Return true if the points A, B, C are strictly counterclockwise.  Return
// false if the points are clockwise or collinear (i.e. if they are all
// contained on some great circle).
//
// Due to numerical errors, situations may arise that are mathematically
// impossible, e.g. ABC may be considered strictly CCW while BCA is not.
// However, the implementation guarantees the following:
//
//   If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
deprecated("Use s2pred::Sign instead.")
bool simpleCCW(in S2Point a, in S2Point b, in S2Point c) {
  // We compute the signed volume of the parallelepiped ABC.  The usual
  // formula for this is (AxB).C, but we compute it here using (CxA).B
  // in order to ensure that ABC and CBA are not both CCW.  This follows
  // from the following identities (which are true numerically, not just
  // mathematically):
  //
  //     (1) x.CrossProd(y) == -(y.CrossProd(x))
  //     (2) (-x).DotProd(y) == -(x.DotProd(y))

  return c.crossProd(a).dotProd(b) > 0;
}
