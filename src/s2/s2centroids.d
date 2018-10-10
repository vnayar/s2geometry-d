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
//
// There are several notions of the "centroid" of a triangle.  First, there
// is the planar centroid, which is simply the centroid of the ordinary
// (non-spherical) triangle defined by the three vertices.  Second, there is
// the surface centroid, which is defined as the intersection of the three
// medians of the spherical triangle.  It is possible to show that this
// point is simply the planar centroid projected to the surface of the
// sphere.  Finally, there is the true centroid (mass centroid), which is
// defined as the surface integral over the spherical triangle of (x,y,z)
// divided by the triangle area.  This is the point that the triangle would
// rotate around if it was spinning in empty space.
//
// The best centroid for most purposes is the true centroid.  Unlike the
// planar and surface centroids, the true centroid behaves linearly as
// regions are added or subtracted.  That is, if you split a triangle into
// pieces and compute the average of their centroids (weighted by triangle
// area), the result equals the centroid of the original triangle.  This is
// not true of the other centroids.
//
// Also note that the surface centroid may be nowhere near the intuitive
// "center" of a spherical triangle.  For example, consider the triangle
// with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
// The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
// within a distance of 2*eps of the vertex B.  Note that the median from A
// (the segment connecting A to the midpoint of BC) passes through S, since
// this is the shortest path connecting the two endpoints.  On the other
// hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
// the surface is a much more reasonable interpretation of the "center" of
// this triangle.

module s2.s2centroids;

import s2.s2point;
import s2.s2pointutil : isUnitLength;

import std.math;

// Return the centroid of the planar triangle ABC.  This can be normalized
// to unit length to obtain the "surface centroid" of the corresponding
// spherical triangle, i.e. the intersection of the three medians.  However,
// note that for large spherical triangles the surface centroid may be
// nowhere near the intuitive "center" (see example above).
S2Point planarCentroid(in S2Point a, in S2Point b, in S2Point c) {
  return (1.0 / 3.0) * (a + b + c);
}

// Returns the true centroid of the spherical triangle ABC multiplied by the
// signed area of spherical triangle ABC.  The reasons for multiplying by
// the signed area are (1) this is the quantity that needs to be summed to
// compute the centroid of a union or difference of triangles, and (2) it's
// actually easier to calculate this way.  All points must have unit length.
S2Point trueCentroid(in S2Point a, in S2Point b, in S2Point c)
in {
  assert(isUnitLength(a));
  assert(isUnitLength(b));
  assert(isUnitLength(c));
} body {

  // I couldn't find any references for computing the true centroid of a
  // spherical triangle...  I have a truly marvellous demonstration of this
  // formula which this margin is too narrow to contain :)

  // Use Angle() in order to get accurate results for small triangles.
  double angle_a = b.angle(c);
  double angle_b = c.angle(a);
  double angle_c = a.angle(b);
  double ra = (angle_a == 0) ? 1 : (angle_a / sin(angle_a));
  double rb = (angle_b == 0) ? 1 : (angle_b / sin(angle_b));
  double rc = (angle_c == 0) ? 1 : (angle_c / sin(angle_c));

  // Now compute a point M such that:
  //
  //  [Ax Ay Az] [Mx]                       [ra]
  //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
  //  [Cx Cy Cz] [Mz]                       [rc]
  //
  // To improve the numerical stability we subtract the first row (A) from the
  // other two rows; this reduces the cancellation error when A, B, and C are
  // very close together.  Then we solve it using Cramer's rule.
  //
  // TODO(ericv): This code still isn't as numerically stable as it could be.
  // The biggest potential improvement is to compute B-A and C-A more
  // accurately so that (B-A)x(C-A) is always inside triangle ABC.
  auto x = S2Point(a.x(), b.x() - a.x(), c.x() - a.x());
  auto y = S2Point(a.y(), b.y() - a.y(), c.y() - a.y());
  auto z = S2Point(a.z(), b.z() - a.z(), c.z() - a.z());
  auto r = S2Point(ra, rb - ra, rc - ra);
  return 0.5 * S2Point(
      y.crossProd(z).dotProd(r), z.crossProd(x).dotProd(r), x.crossProd(y).dotProd(r));
}
