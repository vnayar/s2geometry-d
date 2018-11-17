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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)
//
// Defines a collection of functions for:
//
//   (1) Robustly clipping geodesic edges to the faces of the S2 biunit cube
//       (see s2coords.h), and
//
//   (2) Robustly clipping 2D edges against 2D rectangles.
//
// These functions can be used to efficiently find the set of S2CellIds that
// are intersected by a geodesic edge (e.g., see S2CrossingEdgeQuery).

module s2.s2edge_clipping;

import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s2point;

import s2.s2pointutil : isUnitLength, robustCrossProd;
import s2.s2coords
    : GetFace, GetUVWFace, FaceUVtoXYZ, FaceXYZtoUVW, ValidFaceXYZtoUV, XYZtoFaceUV;

import std.exception : enforce;
import std.algorithm : min, max;
import std.math : fabs, ldexp, signbit, SQRT1_2, SQRT2;

// FaceSegment represents an edge AB clipped to an S2 cube face.  It is
// represented by a face index and a pair of (u,v) coordinates.
struct FaceSegment {
  int face;
  R2Point a, b;
}

alias FaceSegmentVector = FaceSegment[];

// S2PointUVW is used to document that a given S2Point is expressed in the
// (u,v,w) coordinates of some cube face.
alias S2PointUVW = S2Point;

// Subdivides the given edge AB at every point where it crosses the boundary
// between two S2 cube faces and returns the corresponding FaceSegments.  The
// segments are returned in order from A toward B.  The input points must be
// unit length.
//
// This method guarantees that the returned segments form a continuous path
// from A to B, and that all vertices are within kFaceClipErrorUVDist of the
// line AB.  All vertices lie within the [-1,1]x[-1,1] cube face rectangles.
// The results are consistent with s2pred::Sign(), i.e. the edge is
// well-defined even its endpoints are antipodal.  [TODO(ericv): Extend the
// implementation of S2::RobustCrossProd so that this statement is true.]
void getFaceSegments(in S2Point a, in S2Point b, out FaceSegmentVector segments)
in {
  assert(isUnitLength(a));
  assert(isUnitLength(b));
} body {
  // Fast path: both endpoints are on the same face.
  FaceSegment segment;
  int a_face = XYZtoFaceUV(a, segment.a);
  int b_face = XYZtoFaceUV(b, segment.b);
  if (a_face == b_face) {
    segment.face = a_face;
    segments ~= segment;
    return;
  }

  // Starting at A, we follow AB from face to face until we reach the face
  // containing B.  The following code is designed to ensure that we always
  // reach B, even in the presence of numerical errors.
  //
  // First we compute the normal to the plane containing A and B.  This normal
  // becomes the ultimate definition of the line AB; it is used to resolve all
  // questions regarding where exactly the line goes.  Unfortunately due to
  // numerical errors, the line may not quite intersect the faces containing
  // the original endpoints.  We handle this by moving A and/or B slightly if
  // necessary so that they are on faces intersected by the line AB.
  S2Point ab = robustCrossProd(a, b);
  a_face = moveOriginToValidFace(a_face, a, ab, segment.a);
  b_face = moveOriginToValidFace(b_face, b, -ab, segment.b);

  // Now we simply follow AB from face to face until we reach B.
  segment.face = a_face;
  R2Point b_saved = segment.b;
  for (int face = a_face; face != b_face; ) {
    // Complete the current segment by finding the point where AB exits the
    // current face.
    S2PointUVW n = FaceXYZtoUVW(face, ab);
    int exit_axis = getExitAxis(n);
    segment.b = getExitPoint(n, exit_axis);
    segments ~= segment;

    // Compute the next face intersected by AB, and translate the exit point
    // of the current segment into the (u,v) coordinates of the next face.
    // This becomes the first point of the next segment.
    S2Point exit_xyz = FaceUVtoXYZ(face, segment.b);
    face = getNextFace(face, segment.b, exit_axis, n, b_face);
    S2PointUVW exit_uvw = FaceXYZtoUVW(face, exit_xyz);
    segment.face = face;
    segment.a = R2Point(exit_uvw[0], exit_uvw[1]);
  }
  // Finish the last segment.
  segment.b = b_saved;
  segments ~= segment;
}

// This helper function does two things.  First, it clips the line segment AB
// to find the clipped destination B' on a given face.  (The face is specified
// implicitly by expressing *all arguments* in the (u,v,w) coordinates of that
// face.)  Second, it partially computes whether the segment AB intersects
// this face at all.  The actual condition is fairly complicated, but it turns
// out that it can be expressed as a "score" that can be computed
// independently when clipping the two endpoints A and B.  This function
// returns the score for the given endpoint, which is an integer ranging from
// 0 to 3.  If the sum of the two scores is 3 or more, then AB does not
// intersect this face.  See the calling function for the meaning of the
// various parameters.
private int clipDestination(
    in S2PointUVW a, in S2PointUVW b, in S2PointUVW scaled_n,
    in S2PointUVW a_tangent, in S2PointUVW b_tangent, double scale_uv,
    out R2Point uv)
in {
  assert(intersectsFace(scaled_n));
} body {

  // Optimization: if B is within the safe region of the face, use it.
  const double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
  if (b[2] > 0) {
    uv = R2Point(b[0] / b[2], b[1] / b[2]);
    if (max(fabs(uv[0]), fabs(uv[1])) <= kMaxSafeUVCoord)
      return 0;
  }
  // Otherwise find the point B' where the line AB exits the face.
  uv = scale_uv * getExitPoint(scaled_n, getExitAxis(scaled_n));
  auto p = S2PointUVW(uv[0], uv[1], 1.0);

  // Determine if the exit point B' is contained within the segment.  We do this
  // by computing the dot products with two inward-facing tangent vectors at A
  // and B.  If either dot product is negative, we say that B' is on the "wrong
  // side" of that point.  As the point B' moves around the great circle AB past
  // the segment endpoint B, it is initially on the wrong side of B only; as it
  // moves further it is on the wrong side of both endpoints; and then it is on
  // the wrong side of A only.  If the exit point B' is on the wrong side of
  // either endpoint, we can't use it; instead the segment is clipped at the
  // original endpoint B.
  //
  // We reject the segment if the sum of the scores of the two endpoints is 3
  // or more.  Here is what that rule encodes:
  //  - If B' is on the wrong side of A, then the other clipped endpoint A'
  //    must be in the interior of AB (otherwise AB' would go the wrong way
  //    around the circle).  There is a similar rule for A'.
  //  - If B' is on the wrong side of either endpoint (and therefore we must
  //    use the original endpoint B instead), then it must be possible to
  //    project B onto this face (i.e., its w-coordinate must be positive).
  //    This rule is only necessary to handle certain zero-length edges (A=B).
  int score = 0;
  if ((p - a).dotProd(a_tangent) < 0) {
    score = 2;  // B' is on wrong side of A.
  } else if ((p - b).dotProd(b_tangent) < 0) {
    score = 1;  // B' is on wrong side of B.
  }
  if (score > 0) {  // B' is not in the interior of AB.
    if (b[2] <= 0) {
      score = 3;    // B cannot be projected onto this face.
    } else {
      uv = R2Point(b[0] / b[2], b[1] / b[2]);
    }
  }
  return score;
}

// Given an edge AB and a face, returns the (u,v) coordinates for the portion
// of AB that intersects that face.  This method guarantees that the clipped
// vertices lie within the [-1,1]x[-1,1] cube face rectangle and are within
// kFaceClipErrorUVDist of the line AB, but the results may differ from
// those produced by GetFaceSegments.  Returns false if AB does not
// intersect the given face.
bool clipToFace(in S2Point a, in S2Point b, int face, out R2Point a_uv, out R2Point b_uv) {
  return clipToPaddedFace(a, b, face, 0.0, a_uv, b_uv);
}

// Like ClipToFace, but rather than clipping to the square [-1,1]x[-1,1]
// in (u,v) space, this method clips to [-R,R]x[-R,R] where R=(1+padding).
bool clipToPaddedFace(
    in S2Point a_xyz, in S2Point b_xyz, int face, double padding,
    out R2Point a_uv, out R2Point b_uv)
in {
  assert(padding >= 0);
} body {
  // Fast path: both endpoints are on the given face.
  if (GetFace(a_xyz) == face && GetFace(b_xyz) == face) {
    ValidFaceXYZtoUV(face, a_xyz, a_uv);
    ValidFaceXYZtoUV(face, b_xyz, b_uv);
    return true;
  }
  // Convert everything into the (u,v,w) coordinates of the given face.  Note
  // that the cross product *must* be computed in the original (x,y,z)
  // coordinate system because RobustCrossProd (unlike the mathematical cross
  // product) can produce different results in different coordinate systems
  // when one argument is a linear multiple of the other, due to the use of
  // symbolic perturbations.
  S2PointUVW n = FaceXYZtoUVW(face, robustCrossProd(a_xyz, b_xyz));
  S2PointUVW a = FaceXYZtoUVW(face, a_xyz);
  S2PointUVW b = FaceXYZtoUVW(face, b_xyz);

  // Padding is handled by scaling the u- and v-components of the normal.
  // Letting R=1+padding, this means that when we compute the dot product of
  // the normal with a cube face vertex (such as (-1,-1,1)), we will actually
  // compute the dot product with the scaled vertex (-R,-R,1).  This allows
  // methods such as IntersectsFace(), GetExitAxis(), etc, to handle padding
  // with no further modifications.
  const double scale_uv = 1 + padding;
  auto scaled_n = S2PointUVW(scale_uv * n[0], scale_uv * n[1], n[2]);
  if (!intersectsFace(scaled_n)) return false;

  // TODO(ericv): This is a temporary hack until I rewrite S2::RobustCrossProd;
  // it avoids loss of precision in Normalize() when the vector is so small
  // that it underflows.
  if (max(fabs(n[0]), max(fabs(n[1]), fabs(n[2]))) < ldexp(1.0, -511)) {
    n *= ldexp(1.0, 563);
  }  // END OF HACK
  n = n.normalize();
  S2PointUVW a_tangent = n.crossProd(a);
  S2PointUVW b_tangent = b.crossProd(n);
  // As described above, if the sum of the scores from clipping the two
  // endpoints is 3 or more, then the segment does not intersect this face.
  int a_score = clipDestination(b, a, -scaled_n, b_tangent, a_tangent, scale_uv, a_uv);
  int b_score = clipDestination(a, b, scaled_n, a_tangent, b_tangent, scale_uv, b_uv);
  return a_score + b_score < 3;
}

// The maximum error in the vertices returned by GetFaceSegments and
// ClipToFace (compared to an exact calculation):
//
//  - kFaceClipErrorRadians is the maximum angle between a returned vertex
//    and the nearest point on the exact edge AB.  It is equal to the
//    maximum directional error in S2::RobustCrossProd, plus the error when
//    projecting points onto a cube face.
//
//  - kFaceClipErrorDist is the same angle expressed as a maximum distance
//    in (u,v)-space.  In other words, a returned vertex is at most this far
//    from the exact edge AB projected into (u,v)-space.

//  - kFaceClipErrorUVCoord is the same angle expressed as the maximum error
//    in an individual u- or v-coordinate.  In other words, for each
//    returned vertex there is a point on the exact edge AB whose u- and
//    v-coordinates differ from the vertex by at most this amount.

enum double FACE_CLIP_ERROR_RADIANS = 3.0 * double.epsilon;
enum double FACE_CLIP_ERROR_UV_DIST = 9.0 * double.epsilon;
enum double FACE_CLIP_ERROR_UV_COORD = 9.0 * SQRT1_2 * double.epsilon;

// Returns true if the edge AB intersects the given (closed) rectangle to
// within the error bound below.
bool intersectsRect(in R2Point a, in R2Point b, in R2Rect rect) {
  // First check whether the bound of AB intersects "rect".
  R2Rect bound = R2Rect.fromPointPair(a, b);
  if (!rect.intersects(bound)) return false;

  // Otherwise AB intersects "rect" if and only if all four vertices of "rect"
  // do not lie on the same side of the extended line AB.  We test this by
  // finding the two vertices of "rect" with minimum and maximum projections
  // onto the normal of AB, and computing their dot products with the edge
  // normal.
  R2Point n = (b - a).ortho();
  int i = (n[0] >= 0) ? 1 : 0;
  int j = (n[1] >= 0) ? 1 : 0;
  double max = n.dotProd(rect.getVertex(i, j) - a);
  double min = n.dotProd(rect.getVertex(1-i, 1-j) - a);
  return (max >= 0) && (min <= 0);
}

private bool updateEndpoint(ref R1Interval bound, int end, double value) {
  if (end == 0) {
    if (bound.hi() < value) return false;
    if (bound.lo() < value) bound.setLo(value);
  } else {
    if (bound.lo() > value) return false;
    if (bound.hi() > value) bound.setHi(value);
  }
  return true;
}

// The maximum error in IntersectRect.  If some point of AB is inside the
// rectangle by at least this distance, the result is guaranteed to be true;
// if all points of AB are outside the rectangle by at least this distance,
// the result is guaranteed to be false.  This bound assumes that "rect" is
// a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it
// (e.g., by 1e-10 or less).
enum double INTERSECTS_RECT_ERROR_UV_DIST = 3 * SQRT2 * double.epsilon;

// Given an edge AB, returns the portion of AB that is contained by the given
// rectangle "clip".  Returns false if there is no intersection.
bool clipEdge(
    in R2Point a, in R2Point b, in R2Rect clip, out R2Point a_clipped, out R2Point b_clipped) {
  // Compute the bounding rectangle of AB, clip it, and then extract the new
  // endpoints from the clipped bound.
  R2Rect bound = R2Rect.fromPointPair(a, b);
  if (clipEdgeBound(a, b, clip, bound)) {
    int ai = (a[0] > b[0]), aj = (a[1] > b[1]);
    a_clipped = bound.getVertex(ai, aj);
    b_clipped = bound.getVertex(1 - ai, 1 - aj);
    return true;
  }
  return false;
}

// Given an edge AB and a rectangle "clip", returns the bounding rectangle of
// the portion of AB intersected by "clip".  The resulting bound may be
// empty.  This is a convenience function built on top of ClipEdgeBound.
R2Rect getClippedEdgeBound(in R2Point a, in R2Point b, in R2Rect clip) {
  R2Rect bound = R2Rect.fromPointPair(a, b);
  if (clipEdgeBound(a, b, clip, bound)) return bound;
  return R2Rect.empty();
}

// This function can be used to clip an edge AB to sequence of rectangles
// efficiently.  It represents the clipped edges by their bounding boxes
// rather than as a pair of endpoints.  Specifically, let A'B' be some
// portion of an edge AB, and let "bound" be a tight bound of A'B'.  This
// function updates "bound" (in place) to be a tight bound of A'B'
// intersected with a given rectangle "clip".  If A'B' does not intersect
// "clip", returns false and does not necessarily update "bound".
//
// REQUIRES: "bound" is a tight bounding rectangle for some portion of AB.
// (This condition is automatically satisfied if you start with the bounding
// box of AB and clip to a sequence of rectangles, stopping when the method
// returns false.)
bool clipEdgeBound(in R2Point a, in R2Point b, in R2Rect clip, ref R2Rect bound) {
  // "diag" indicates which diagonal of the bounding box is spanned by AB: it
  // is 0 if AB has positive slope, and 1 if AB has negative slope.  This is
  // used to determine which interval endpoints need to be updated each time
  // the edge is clipped.
  int diag = (a[0] > b[0]) != (a[1] > b[1]);
  return clipBoundAxis(a[0], b[0], bound[0], a[1], b[1], bound[1], diag, clip[0])
      && clipBoundAxis(a[1], b[1], bound[1], a[0], b[0], bound[0], diag, clip[1]);
}

// Given a line segment from (a0,a1) to (b0,b1) and a bounding interval for
// each axis, clip the segment further if necessary so that "bound0" does not
// extend outside the given interval "clip".  "diag" is a a precomputed helper
// variable that indicates which diagonal of the bounding box is spanned by AB:
// it is 0 if AB has positive slope, and 1 if AB has negative slope.
private bool clipBoundAxis(
    double a0, double b0, ref R1Interval bound0, double a1, double b1, ref R1Interval bound1,
    int diag, in R1Interval clip0) {
  if (bound0.lo() < clip0.lo()) {
    if (bound0.hi() < clip0.lo()) return false;
    bound0[0] = clip0.lo();
    if (!updateEndpoint(bound1, diag, interpolateDouble(clip0.lo(), a0, b0, a1, b1)))
      return false;
  }
  if (bound0.hi() > clip0.hi()) {
    if (bound0.lo() > clip0.hi()) return false;
    bound0[1] = clip0.hi();
    if (!updateEndpoint(bound1, 1 - diag, interpolateDouble(clip0.hi(), a0, b0, a1, b1)))
      return false;
  }
  return true;
}


// The maximum error in the vertices generated by ClipEdge and the bounds
// generated by ClipEdgeBound (compared to an exact calculation):
//
//  - kEdgeClipErrorUVCoord is the maximum error in a u- or v-coordinate
//    compared to the exact result, assuming that the points A and B are in
//    the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
//
//  - kEdgeClipErrorUVDist is the maximum distance from a clipped point to
//    the corresponding exact result.  It is equal to the error in a single
//    coordinate because at most one coordinate is subject to error.

enum double EDGE_CLIP_ERROR_UV_COORD = 2.25 * double.epsilon;
enum double EDGE_CLIP_ERROR_UV_DIST = 2.25 * double.epsilon;

// Given a value x that is some linear combination of a and b, returns the
// value x1 that is the same linear combination of a1 and b1.  This function
// makes the following guarantees:
//  - If x == a, then x1 = a1 (exactly).
//  - If x == b, then x1 = b1 (exactly).
//  - If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
// REQUIRES: a != b
double interpolateDouble(double x, double a, double b, double a1, double b1)
in {
  assert(a != b);
} body {
  // To get results that are accurate near both A and B, we interpolate
  // starting from the closer of the two points.
  if (fabs(a - x) <= fabs(b - x)) {
    return a1 + (b1 - a1) * (x - a) / (b - a);
  } else {
    return b1 + (a1 - b1) * (x - b) / (a - b);
  }
}


// Given a line segment AB whose origin A has been projected onto a given cube
// face, determine whether it is necessary to project A onto a different face
// instead.  This can happen because the normal of the line AB is not computed
// exactly, so that the line AB (defined as the set of points perpendicular to
// the normal) may not intersect the cube face containing A.  Even if it does
// intersect the face, the "exit point" of the line from that face may be on
// the wrong side of A (i.e., in the direction away from B).  If this happens,
// we reproject A onto the adjacent face where the line AB approaches A most
// closely.  This moves the origin by a small amount, but never more than the
// error tolerances documented in the header file.
private int moveOriginToValidFace(int face, in S2Point a, in S2Point ab, ref R2Point a_uv) {
  // Fast path: if the origin is sufficiently far inside the face, it is
  // always safe to use it.
  const double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
  if (max(fabs(a_uv[0]), fabs(a_uv[1])) <= kMaxSafeUVCoord) {
    return face;
  }
  // Otherwise check whether the normal AB even intersects this face.
  S2PointUVW n = FaceXYZtoUVW(face, ab);
  if (intersectsFace(n)) {
    // Check whether the point where the line AB exits this face is on the
    // wrong side of A (by more than the acceptable error tolerance).
    S2Point exit = FaceUVtoXYZ(face, getExitPoint(n, getExitAxis(n)));
    S2Point a_tangent = ab.normalize().crossProd(a);
    if ((exit - a).dotProd(a_tangent) >= -FACE_CLIP_ERROR_RADIANS) {
      return face;  // We can use the given face.
    }
  }
  // Otherwise we reproject A to the nearest adjacent face.  (If line AB does
  // not pass through a given face, it must pass through all adjacent faces.)
  if (fabs(a_uv[0]) >= fabs(a_uv[1])) {
    face = GetUVWFace(face, 0 /*U axis*/, a_uv[0] > 0);
  } else {
    face = GetUVWFace(face, 1 /*V axis*/, a_uv[1] > 0);
  }
  enforce(intersectsFace(FaceXYZtoUVW(face, ab)));
  ValidFaceXYZtoUV(face, a, a_uv);
  a_uv[0] = max(-1.0, min(1.0, a_uv[0]));
  a_uv[1] = max(-1.0, min(1.0, a_uv[1]));
  return face;
}

// Given cube face F and a directed line L (represented by its CCW normal N in
// the (u,v,w) coordinates of F), compute the axis of the cube face edge where
// L exits the face: return 0 if L exits through the u=-1 or u=+1 edge, and 1
// if L exits through the v=-1 or v=+1 edge.  Either result is acceptable if L
// exits exactly through a corner vertex of the cube face.
private int getExitAxis(in S2PointUVW n)
in {
  assert(intersectsFace(n));
} body {
  if (intersectsOppositeEdges(n)) {
    // The line passes through through opposite edges of the face.
    // It exits through the v=+1 or v=-1 edge if the u-component of N has a
    // larger absolute magnitude than the v-component.
    return (fabs(n[0]) >= fabs(n[1])) ? 1 : 0;
  } else {
    // The line passes through through two adjacent edges of the face.
    // It exits the v=+1 or v=-1 edge if an even number of the components of N
    // are negative.  We test this using signbit() rather than multiplication
    // to avoid the possibility of underflow.
    enforce(n[0] != 0 && n[1] != 0  && n[2] != 0);
    return ((signbit(n[0]) ^ signbit(n[1]) ^ signbit(n[2])) == 0) ? 1 : 0;
  }
}

// Given a cube face F, a directed line L (represented by its CCW normal N in
// the (u,v,w) coordinates of F), and result of GetExitAxis(N), return the
// (u,v) coordinates of the point where L exits the cube face.
private R2Point getExitPoint(in S2PointUVW n, int axis) {
  if (axis == 0) {
    double u = (n[1] > 0) ? 1.0 : -1.0;
    return R2Point(u, (-u * n[0] - n[2]) / n[1]);
  } else {
    double v = (n[0] < 0) ? 1.0 : -1.0;
    return R2Point((-v * n[1] - n[2]) / n[0], v);
  }
}

// Return the next face that should be visited by GetFaceSegments, given that
// we have just visited "face" and we are following the line AB (represented
// by its normal N in the (u,v,w) coordinates of that face).  The other
// arguments include the point where AB exits "face", the corresponding
// exit axis, and the "target face" containing the destination point B.
private int getNextFace(int face, in R2Point exit, int axis, in S2PointUVW n, int target_face) {
  // We return the face that is adjacent to the exit point along the given
  // axis.  If line AB exits *exactly* through a corner of the face, there are
  // two possible next faces.  If one is the "target face" containing B, then
  // we guarantee that we advance to that face directly.
  //
  // The three conditions below check that (1) AB exits approximately through
  // a corner, (2) the adjacent face along the non-exit axis is the target
  // face, and (3) AB exits *exactly* through the corner.  (The SumEquals()
  // code checks whether the dot product of (u,v,1) and "n" is exactly zero.)
  if (fabs(exit[1 - axis]) == 1
      && GetUVWFace(face, 1 - axis, exit[1 - axis] > 0) == target_face
      && sumEquals(exit[0] * n[0], exit[1] * n[1], -n[2])) {
    return target_face;
  }
  // Otherwise return the face that is adjacent to the exit point in the
  // direction of the exit axis.
  return GetUVWFace(face, axis, exit[axis] > 0);
}

// The three functions below all compare a sum (u + v) to a third value w.
// They are implemented in such a way that they produce an exact result even
// though all calculations are done with ordinary floating-point operations.
// Here are the principles on which these functions are based:
//
// A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
//
// B. If u + v < w in exact arithmetic, then at least one of the following
//    expressions is true in floating-point:
//       u + v < w
//       u < w - v
//       v < w - u
//
//    Proof: By rearranging terms and substituting ">" for "<", we can assume
//    that all values are non-negative.  Now clearly "w" is not the smallest
//    value, so assume WLOG that "u" is the smallest.  We want to show that
//    u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
//    exact since the result is smaller in magnitude than either input value,
//    so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
//    (even in floating point), so the result also holds.

// Return true if u + v == w exactly.
private bool sumEquals(double u, double v, double w) {
  return (u + v == w) && (u == w - v) && (v == w - u);
}

// Return true if a given directed line L intersects the cube face F.  The
// line L is defined by its normal N in the (u,v,w) coordinates of F.
private bool intersectsFace(in S2PointUVW n) {
  // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
  // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
  // and (-1,1,1) do not all have the same sign.  This is true exactly when
  // |Nu| + |Nv| >= |Nw|.  The code below evaluates this expression exactly
  // (see comments above).
  double u = fabs(n[0]), v = fabs(n[1]), w = fabs(n[2]);
  // We only need to consider the cases where u or v is the smallest value,
  // since if w is the smallest then both expressions below will have a
  // positive LHS and a negative RHS.
  return (v >= w - u) && (u >= w - v);
}

// Given a directed line L intersecting a cube face F, return true if L
// intersects two opposite edges of F (including the case where L passes
// exactly through a corner vertex of F).  The line L is defined by its
// normal N in the (u,v,w) coordinates of F.
private bool intersectsOppositeEdges(in S2PointUVW n) {
  // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
  // and only exactly two of the corner vertices lie on each side of L.  This
  // is true exactly when ||Nu| - |Nv|| >= |Nw|.  The code below evaluates this
  // expression exactly (see comments above).
  double u = fabs(n[0]), v = fabs(n[1]), w = fabs(n[2]);
  // If w is the smallest, the following line returns an exact result.
  if (fabs(u - v) != w) return fabs(u - v) >= w;
  // Otherwise u - v = w exactly, or w is not the smallest value.  In either
  // case the following line returns the correct result.
  return (u >= v) ? (u - w >= v) : (v - w >= u);
}
