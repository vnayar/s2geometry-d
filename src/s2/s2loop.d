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

module s2.s2loop;

import s2.logger;
import s2.mutable_s2shape_index;
import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2centroids : trueCentroid;
import s2.s2closest_edge_query;
import s2.s2crossing_edge_query;
import s2.s2debug;
import s2.s2edge_crosser;
import s2.s2edge_distances : getDistance;
import s2.s2error;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2measures : signedArea, turnAngle;
import s2.s2padded_cell;
import s2.s2point;
import s2.s2pointutil : approxEquals, fromFrame, getFrame, origin, robustCrossProd;
import s2.s2predicates : orderedCCW;
import s2.s2region;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2wedge_relations : wedgeContains, wedgeIntersects;
import s2.shapeutil.visit_crossing_edge_pairs : findSelfIntersection;
import s2.util.math.matrix3x3;
import s2.util.math.vector;

import std.algorithm : copy, min, max, reverse, swap;
import std.array : appender;
import std.conv : to;
import std.exception;
import std.math;
import std.range : empty, back, isInputRange, popBack;
import core.atomic;


// Build the S2ShapeIndex only when it is first needed.  This can save
// significant amounts of memory and time when geometry is constructed but
// never queried, for example when loops are passed directly to S2Polygon,
// or when geometry is being converted from one format to another.
enum bool LAZY_INDEXING = true;

enum double M_PI = cast(double) PI;
enum double M_PI_2 = cast(double) PI_2;

/**
 * An S2Loop represents a simple spherical polygon.  It consists of a single
 * chain of vertices where the first vertex is implicitly connected to the
 * last. All loops are defined to have a CCW orientation, i.e. the interior of
 * the loop is on the left side of the edges.  This implies that a clockwise
 * loop enclosing a small area is interpreted to be a CCW loop enclosing a
 * very large area.
 *
 * Loops are not allowed to have any duplicate vertices (whether adjacent or
 * not).  Non-adjacent edges are not allowed to intersect, and furthermore edges
 * of length 180 degrees are not allowed (i.e., adjacent vertices cannot be
 * antipodal).  Loops must have at least 3 vertices (except for the "empty" and
 * "full" loops discussed below).  Although these restrictions are not enforced
 * in optimized code, you may get unexpected results if they are violated.
 *
 * There are two special loops: the "empty" loop contains no points, while the
 * "full" loop contains all points.  These loops do not have any edges, but to
 * preserve the invariant that every loop can be represented as a vertex
 * chain, they are defined as having exactly one vertex each (see kEmpty and
 * kFull).
 *
 * Point containment of loops is defined such that if the sphere is subdivided
 * into faces (loops), every point is contained by exactly one face.  This
 * implies that loops do not necessarily contain their vertices.
 *
 * Note: The reason that duplicate vertices and intersecting edges are not
 * allowed is that they make it harder to define and implement loop
 * relationships, e.g. whether one loop contains another.  If your data does
 * not satisfy these restrictions, you can use S2Builder to normalize it.
 *
 * TODO: Convert logic to use a ForwardRange rather than making a copy of its vertices.
 */
class S2Loop : S2Region {
public:
  // Default constructor.  The loop must be initialized by calling Init() or
  // Decode() before it is used.
  this() {
    _index = new MutableS2ShapeIndex();
  }

  this(in S2Point[] vertices) {
    this(vertices, S2Debug.ALLOW);
  }

  // Convenience constructor that calls Init() with the given vertices.
  this(in S2Point[] vertices, S2Debug s2DebugOverride) {
    this();
    _s2DebugOverride = s2DebugOverride;
    initialize(vertices);
  }

  // Initialize a loop with given vertices.  The last vertex is implicitly
  // connected to the first.  All points should be unit length.  Loops must
  // have at least 3 vertices (except for the "empty" and "full" loops, see
  // kEmpty and kFull).  This method may be called multiple times.
  void initialize(RangeT)(RangeT vertexRange)
  if (isInputRange!RangeT /*&& is(ElementType!RangeT == Vector!(double, 3))*/) {
    clearIndex();
    _vertices.length = 0;
    _vertices ~= vertexRange;
    //copy(vertexRange, appender(_vertices));
    initOriginAndBound();
  }

  // A special vertex chain of length 1 that creates an empty loop (i.e., a
  // loop with no edges that contains no points).  Example usage:
  //
  //    S2Loop emptyLoop = new S2Loop(S2Loop.empty());
  //
  // The loop may be safely encoded lossily (e.g. by snapping it to an S2Cell
  // center) as long as its position does not move by 90 degrees or more.
  static S2Point[] empty() {
    return [S2Loop.emptyVertex()];
  }

  // A special vertex chain of length 1 that creates a full loop (i.e., a loop
  // with no edges that contains all points).  See kEmpty() for details.
  static S2Point[] full() {
    return [S2Loop.fullVertex()];
  }

  // Construct a loop corresponding to the given cell.
  //
  // Note that the loop and cell *do not* contain exactly the same set of
  // points, because S2Loop and S2Cell have slightly different definitions of
  // point containment.  For example, an S2Cell vertex is contained by all
  // four neighboring S2Cells, but it is contained by exactly one of four
  // S2Loops constructed from those cells.  As another example, the S2Cell
  // coverings of "cell" and "S2Loop(cell)" will be different, because the
  // loop contains points on its boundary that actually belong to other cells
  // (i.e., the covering will include a layer of neighboring cells).
  this(in S2Cell cell) {
    this();
    _depth = 0;
    _vertices.length = 4;
    _s2DebugOverride = S2Debug.ALLOW;
    _unindexedContainsCalls = 0;
    foreach (i; 0 .. 4) {
      _vertices[i] = cell.getVertex(i);
    }
    // We recompute the bounding rectangle ourselves, since S2Cell uses a
    // different method and we need all the bounds to be consistent.
    initOriginAndBound();
  }

  // Allows overriding the automatic validity checks controlled by the
  // --s2debug flag.  If this flag is true, then loops are automatically
  // checked for validity as they are initialized.  The main reason to disable
  // this flag is if you intend to call IsValid() explicitly, like this:
  //
  //   S2Loop loop;
  //   loop.set_s2debug_override(S2Debug::DISABLE);
  //   loop.Init(...);
  //   if (!loop.IsValid()) { ... }
  //
  // Without the call to set_s2debug_override(), invalid data would cause a
  // fatal error in Init() whenever the --s2debug flag is enabled.
  //
  // This setting is preserved across calls to Init() and Decode().
  @property
  void s2DebugOverride(S2Debug s2DebugOverride) {
    _s2DebugOverride = s2DebugOverride;
  }

  @property
  S2Debug s2DebugOverride() const {
    return _s2DebugOverride;
  }

  // Returns true if this is a valid loop.  Note that validity is checked
  // automatically during initialization when --s2debug is enabled (true by
  // default in debug binaries).
  bool isValid() {
    S2Error error;
    if (findValidationError(error)) {
      if (_s2DebugOverride == S2Debug.ALLOW) logger.logError(error);
      return false;
    }
    return true;
  }

  // Returns true if this is *not* a valid loop and sets "error"
  // appropriately.  Otherwise returns false and leaves "error" unchanged.
  bool findValidationError(out S2Error error) {
    return (findValidationErrorNoIndex(error) || findSelfIntersection(_index, error));
  }

  // Like FindValidationError(), but skips any checks that would require
  // building the S2ShapeIndex (i.e., self-intersection tests).  This is used
  // by the S2Polygon implementation, which uses its own index to check for
  // loop self-intersections.
  bool findValidationErrorNoIndex(S2Error error) const
  in {
    // subregion_bound_ must be at least as large as bound_.  (This is an
    // internal consistency check rather than a test of client data.)
    assert(_subregionBound.contains(_bound));
  } body {
    import s2.s2pointutil : isUnitLength;
    // All vertices must be unit length.  (Unfortunately this check happens too
    // late in debug mode, because S2Loop construction calls s2pred::Sign which
    // expects vertices to be unit length.  But it is still a useful check in
    // optimized builds.)
    for (int i = 0; i < numVertices(); ++i) {
      if (!isUnitLength(vertex(i))) {
        error.initialize(S2Error.Code.NOT_UNIT_LENGTH, "Vertex %d is not unit length", i);
        return true;
      }
    }
    // Loops must have at least 3 vertices (except for "empty" and "full").
    if (numVertices() < 3) {
      if (isEmptyOrFull()) {
        return false;  // Skip remaining tests.
      }
      error.initialize(S2Error.Code.LOOP_NOT_ENOUGH_VERTICES,
          "Non-empty, non-full loops must have at least 3 vertices");
      return true;
    }
    // Loops are not allowed to have any duplicate vertices or edge crossings.
    // We split this check into two parts.  First we check that no edge is
    // degenerate (identical endpoints).  Then we check that there are no
    // intersections between non-adjacent edges (including at vertices).  The
    // second part needs the S2ShapeIndex, so it does not fall within the scope
    // of this method.
    for (int i = 0; i < numVertices(); ++i) {
      if (vertex(i) == vertex(i+1)) {
        error.initialize(
            S2Error.Code.DUPLICATE_VERTICES, "Edge %d is degenerate (duplicate vertex)", i);
        return true;
      }
      if (vertex(i) == -vertex(i + 1)) {
        error.initialize(S2Error.Code.ANTIPODAL_VERTICES,
            "Vertices %d and %d are antipodal", i, (i + 1) % numVertices());
        return true;
      }
    }
    return false;
  }

  int numVertices() const {
    return cast(int) _vertices.length;
  }

  // For convenience, we make two entire copies of the vertex list available:
  // vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == num_vertices().
  //
  // REQUIRES: 0 <= i < 2 * num_vertices()
  S2Point vertex(int i) const
  in {
    assert(i >= 0);
    assert(i < 2 * numVertices());
  } body {
    int j = i - numVertices();
    return _vertices[j < 0 ? i : j];
  }

  const(S2Point[]) vertices() const {
    return _vertices;
  }

  // Like vertex(), but this method returns vertices in reverse order if the
  // loop represents a polygon hole.  For example, arguments 0, 1, 2 are
  // mapped to vertices n-1, n-2, n-3, where n == num_vertices().  This
  // ensures that the interior of the polygon is always to the left of the
  // vertex chain.
  //
  // REQUIRES: 0 <= i < 2 * num_vertices()
  S2Point orientedVertex(int i) const
  in {
    assert(i >= 0);
    assert(i < 2 * numVertices());
  } body {
    int j = i - numVertices();
    if (j < 0) j = i;
    if (isHole()) j = numVertices() - 1 - j;
    return _vertices[j];
  }

  // Return true if this is the special "empty" loop that contains no points.
  bool isEmpty() const {
    return isEmptyOrFull() && !containsOrigin();
  }

  // Return true if this is the special "full" loop that contains all points.
  bool isFull() const {
    return isEmptyOrFull() && containsOrigin();
  }

  // Return true if this loop is either "empty" or "full".
  bool isEmptyOrFull() const {
    return numVertices() == 1;
  }

  // The depth of a loop is defined as its nesting level within its containing
  // polygon.  "Outer shell" loops have depth 0, holes within those loops have
  // depth 1, shells within those holes have depth 2, etc.  This field is only
  // used by the S2Polygon implementation.
  int depth() const {
    return _depth;
  }

  void setDepth(int depth) {
    _depth = depth;
  }

  // Return true if this loop represents a hole in its containing polygon.
  bool isHole() const {
    return (_depth & 1) != 0;
  }

  // The sign of a loop is -1 if the loop represents a hole in its containing
  // polygon, and +1 otherwise.
  int sign() const {
    return isHole() ? -1 : 1;
  }

  // Return true if the loop area is at most 2*Pi.  Degenerate loops are
  // handled consistently with s2pred::Sign(), i.e., if a loop can be
  // expressed as the union of degenerate or nearly-degenerate CCW triangles,
  // then it will always be considered normalized.
  bool isNormalized() const {
    // Optimization: if the longitude span is less than 180 degrees, then the
    // loop covers less than half the sphere and is therefore normalized.
    if (_bound.lng().getLength() < M_PI) return true;

    // We allow some error so that hemispheres are always considered normalized.
    // TODO(ericv): This is no longer required by the S2Polygon implementation,
    // so alternatively we could create the invariant that a loop is normalized
    // if and only if its complement is not normalized.
    return getTurningAngle() >= -getTurningAngleMaxError();
  }

  // Invert the loop if necessary so that the area enclosed by the loop is at
  // most 2*Pi.
  void normalize()
  out {
    assert(isNormalized());
  } body {
    if (!isNormalized()) invert();
  }

  // Reverse the order of the loop vertices, effectively complementing the
  // region represented by the loop.  For example, the loop ABCD (with edges
  // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
  // Notice that the last edge is the same in both cases except that its
  // direction has been reversed.
  void invert() {
    clearIndex();
    if (isEmptyOrFull()) {
      _vertices[0] = isFull() ? emptyVertex() : fullVertex();
    } else {
      reverse(_vertices);
    }
    // origin_inside_ must be set correctly before building the S2ShapeIndex.
    _originInside ^= true;
    if (_bound.lat().lo() > -M_PI_2 && _bound.lat().hi() < M_PI_2) {
      // The complement of this loop contains both poles.
      _subregionBound = _bound = S2LatLngRect.full();
    } else {
      initBound();
    }
    initIndex();
  }

  // Return the area of the loop interior, i.e. the region on the left side of
  // the loop.  The return value is between 0 and 4*Pi.  (Note that the return
  // value is not affected by whether this loop is a "hole" or a "shell".)
  double getArea() const {
    // It is suprisingly difficult to compute the area of a loop robustly.  The
    // main issues are (1) whether degenerate loops are considered to be CCW or
    // not (i.e., whether their area is close to 0 or 4*Pi), and (2) computing
    // the areas of small loops with good relative accuracy.
    //
    // With respect to degeneracies, we would like GetArea() to be consistent
    // with S2Loop::Contains(S2Point) in that loops that contain many points
    // should have large areas, and loops that contain few points should have
    // small areas.  For example, if a degenerate triangle is considered CCW
    // according to s2pred::Sign(), then it will contain very few points and
    // its area should be approximately zero.  On the other hand if it is
    // considered clockwise, then it will contain virtually all points and so
    // its area should be approximately 4*Pi.

    // More precisely, let U be the set of S2Points for which S2::IsUnitLength()
    // is true, let P(U) be the projection of those points onto the mathematical
    // unit sphere, and let V(P(U)) be the Voronoi diagram of the projected
    // points.  Then for every loop x, we would like GetArea() to approximately
    // equal the sum of the areas of the Voronoi regions of the points p for
    // which x.Contains(p) is true.
    //
    // The second issue is that we want to compute the area of small loops
    // accurately.  This requires having good relative precision rather than
    // good absolute precision.  For example, if the area of a loop is 1e-12 and
    // the error is 1e-15, then the area only has 3 digits of accuracy.  (For
    // reference, 1e-12 is about 40 square meters on the surface of the earth.)
    // We would like to have good relative accuracy even for small loops.
    //
    // To achieve these goals, we combine two different methods of computing the
    // area.  This first method is based on the Gauss-Bonnet theorem, which says
    // that the area enclosed by the loop equals 2*Pi minus the total geodesic
    // curvature of the loop (i.e., the sum of the "turning angles" at all the
    // loop vertices).  The big advantage of this method is that as long as we
    // use s2pred::Sign() to compute the turning angle at each vertex, then
    // degeneracies are always handled correctly.  In other words, if a
    // degenerate loop is CCW according to the symbolic perturbations used by
    // s2pred::Sign(), then its turning angle will be approximately 2*Pi.
    //
    // The disadvantage of the Gauss-Bonnet method is that its absolute error is
    // about 2e-15 times the number of vertices (see GetTurningAngleMaxError).
    // So, it cannot compute the area of small loops accurately.
    //
    // The second method is based on splitting the loop into triangles and
    // summing the area of each triangle.  To avoid the difficulty and expense
    // of decomposing the loop into a union of non-overlapping triangles,
    // instead we compute a signed sum over triangles that may overlap (see the
    // comments for S2Loop::GetSurfaceIntegral).  The advantage of this method
    // is that the area of each triangle can be computed with much better
    // relative accuracy (using l'Huilier's theorem).  The disadvantage is that
    // the result is a signed area: CCW loops may yield a small positive value,
    // while CW loops may yield a small negative value (which is converted to a
    // positive area by adding 4*Pi).  This means that small errors in computing
    // the signed area may translate into a very large error in the result (if
    // the sign of the sum is incorrect).
    //
    // So, our strategy is to combine these two methods as follows.  First we
    // compute the area using the "signed sum over triangles" approach (since it
    // is generally more accurate).  We also estimate the maximum error in this
    // result.  If the signed area is too close to zero (i.e., zero is within
    // the error bounds), then we double-check the sign of the result using the
    // Gauss-Bonnet method.  (In fact we just call IsNormalized(), which is
    // based on this method.)  If the two methods disagree, we return either 0
    // or 4*Pi based on the result of IsNormalized().  Otherwise we return the
    // area that we computed originally.

    if (isEmptyOrFull()) {
      return containsOrigin() ? (4 * M_PI) : 0;
    }
    double area = getSurfaceIntegral(&signedArea);

    // TODO(ericv): This error estimate is very approximate.  There are two
    // issues: (1) SignedArea needs some improvements to ensure that its error
    // is actually never higher than GirardArea, and (2) although the number of
    // triangles in the sum is typically N-2, in theory it could be as high as
    // 2*N for pathological inputs.  But in other respects this error bound is
    // very conservative since it assumes that the maximum error is achieved on
    // every triangle.
    double max_error = getTurningAngleMaxError();

    // The signed area should be between approximately -4*Pi and 4*Pi.
    enforce(fabs(area) <= 4 * M_PI + max_error);
    if (area < 0) {
      // We have computed the negative of the area of the loop exterior.
      area += 4 * M_PI;
    }
    area = max(0.0, min(4 * M_PI, area));

    // If the area is close enough to zero or 4*Pi so that the loop orientation
    // is ambiguous, then we compute the loop orientation explicitly.
    if (area < max_error && !isNormalized()) {
      return 4 * M_PI;
    } else if (area > (4 * M_PI - max_error) && isNormalized()) {
      return 0.0;
    } else {
      return area;
    }

  }

  // Return the true centroid of the loop multiplied by the area of the loop
  // (see s2centroids.h for details on centroids).  The result is not unit
  // length, so you may want to normalize it.  Also note that in general, the
  // centroid may not be contained by the loop.
  //
  // We prescale by the loop area for two reasons: (1) it is cheaper to
  // compute this way, and (2) it makes it easier to compute the centroid of
  // more complicated shapes (by splitting them into disjoint regions and
  // adding their centroids).
  //
  // Note that the return value is not affected by whether this loop is a
  // "hole" or a "shell".
  S2Point getCentroid() const {
    // GetSurfaceIntegral() returns either the integral of position over loop
    // interior, or the negative of the integral of position over the loop
    // exterior.  But these two values are the same (!), because the integral of
    // position over the entire sphere is (0, 0, 0).
    return getSurfaceIntegral(&trueCentroid);
  }

  // Return the sum of the turning angles at each vertex.  The return value is
  // positive if the loop is counter-clockwise, negative if the loop is
  // clockwise, and zero if the loop is a great circle.  Degenerate and
  // nearly-degenerate loops are handled consistently with s2pred::Sign().
  // So for example, if a loop has zero area (i.e., it is a very small CCW
  // loop) then the turning angle will always be negative.
  //
  // This quantity is also called the "geodesic curvature" of the loop.
  double getTurningAngle() const {
    // For empty and full loops, we return the limit value as the loop area
    // approaches 0 or 4*Pi respectively.
    if (isEmptyOrFull()) {
      return containsOrigin() ? (-2 * M_PI) : (2 * M_PI);
    }
    // Don't crash even if the loop is not well-defined.
    if (numVertices() < 3) return 0;

    // To ensure that we get the same result when the vertex order is rotated,
    // and that the result is negated when the vertex order is reversed, we need
    // to add up the individual turn angles in a consistent order.  (In general,
    // adding up a set of numbers in a different order can change the sum due to
    // rounding errors.)
    //
    // Furthermore, if we just accumulate an ordinary sum then the worst-case
    // error is quadratic in the number of vertices.  (This can happen with
    // spiral shapes, where the partial sum of the turning angles can be linear
    // in the number of vertices.)  To avoid this we use the Kahan summation
    // algorithm (http://en.wikipedia.org/wiki/Kahan_summation_algorithm).

    int n = numVertices();
    int dir, i = getCanonicalFirstVertex(dir);
    double sum = turnAngle(vertex((i + n - dir) % n), vertex(i), vertex((i + dir) % n));
    double compensation = 0;  // Kahan summation algorithm
    while (--n > 0) {
      i += dir;
      double angle = turnAngle(vertex(i - dir), vertex(i), vertex(i + dir));
      double old_sum = sum;
      angle += compensation;
      sum += angle;
      compensation = (old_sum - sum) + angle;
    }
    return dir * (sum + compensation);
  }

  // Return the maximum error in GetTurningAngle().  The return value is not
  // constant; it depends on the loop.
  double getTurningAngleMaxError() const {
    // The maximum error can be bounded as follows:
    //   2.24 * DBL_EPSILON    for RobustCrossProd(b, a)
    //   2.24 * DBL_EPSILON    for RobustCrossProd(c, b)
    //   3.25 * DBL_EPSILON    for Angle()
    //   2.00 * DBL_EPSILON    for each addition in the Kahan summation
    //   ------------------
    //   9.73 * DBL_EPSILON
    const double kMaxErrorPerVertex = 9.73 * double.epsilon;
    return kMaxErrorPerVertex * numVertices();
  }

  // Return the distance from the given point to the loop interior.  If the
  // loop is empty, return S1Angle::Infinity().  "x" should be unit length.
  S1Angle getDistance(in S2Point x) {
    // Note that S2Loop::Contains(S2Point) is slightly more efficient than the
    // generic version used by S2ClosestEdgeQuery.
    if (contains(x)) return S1Angle.zero();
    return getDistanceToBoundary(x);
  }

  // Return the distance from the given point to the loop boundary.  If the
  // loop is empty or full, return S1Angle::Infinity() (since the loop has no
  // boundary).  "x" should be unit length.
  S1Angle getDistanceToBoundary(in S2Point x) {
    auto options = new S2ClosestEdgeQuery.Options();
    options.setIncludeInteriors(false);
    auto t = new S2ClosestEdgeQuery.PointTarget(x);
    return new S2ClosestEdgeQuery(_index, options).getDistance(t).toS1Angle();
  }

  // If the given point is contained by the loop, return it.  Otherwise return
  // the closest point on the loop boundary.  If the loop is empty, return the
  // input argument.  Note that the result may or may not be contained by the
  // loop.  "x" should be unit length.
  S2Point project(in S2Point x) {
    if (contains(x)) return x;
    return projectToBoundary(x);
  }

  // Return the closest point on the loop boundary to the given point.  If the
  // loop is empty or full, return the input argument (since the loop has no
  // boundary).  "x" should be unit length.
  S2Point projectToBoundary(in S2Point x) {
    auto options = new S2ClosestEdgeQuery.Options();
    options.setIncludeInteriors(false);
    auto q = new S2ClosestEdgeQuery(_index, options);
    auto target = new S2ClosestEdgeQuery.PointTarget(x);
    S2ClosestEdgeQuery.Result edge = q.findClosestEdge(target);
    return q.project(x, edge);
  }

  // Return true if the region contained by this loop is a superset of the
  // region contained by the given other loop.
  bool contains(S2Loop b) {
    // For this loop A to contains the given loop B, all of the following must
    // be true:
    //
    //  (1) There are no edge crossings between A and B except at vertices.
    //
    //  (2) At every vertex that is shared between A and B, the local edge
    //      ordering implies that A contains B.
    //
    //  (3) If there are no shared vertices, then A must contain a vertex of B
    //      and B must not contain a vertex of A.  (An arbitrary vertex may be
    //      chosen in each case.)
    //
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.
    if (!_subregionBound.contains(b._bound)) return false;

    // Special cases to handle either loop being empty or full.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return isFull() || b.isEmpty();
    }

    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    auto relation = new ContainsRelation();
    if (hasCrossingRelation(this, b, relation)) return false;

    // There are no crossings, and if there are any shared vertices then A
    // contains B locally at each shared vertex.
    if (relation.foundSharedVertex()) return true;

    // Since there are no edge intersections or shared vertices, we just need to
    // test condition (3) above.  We can skip this test if we discovered that A
    // contains at least one point of B while checking for edge crossings.
    if (!contains(b.vertex(0))) return false;

    // We still need to check whether (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if ((b._subregionBound.contains(_bound) || b._bound.unite(_bound).isFull())
        && b.contains(vertex(0))) {
      return false;
    }
    return true;
  }

  // Return true if the region contained by this loop intersects the region
  // contained by the given other loop.
  bool intersects(S2Loop b) {
    // a->Intersects(b) if and only if !a->Complement()->Contains(b).
    // This code is similar to Contains(), but is optimized for the case
    // where both loops enclose less than half of the sphere.
    if (!_bound.intersects(b._bound)) return false;

    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    auto relation = new IntersectsRelation();
    if (hasCrossingRelation(this, b, relation)) return true;
    if (relation.foundSharedVertex()) return false;

    // Since there are no edge intersections or shared vertices, the loops
    // intersect only if A contains B, B contains A, or the two loops contain
    // each other's boundaries.  These checks are usually cheap because of the
    // bounding box preconditions.  Note that neither loop is empty (because of
    // the bounding box check above), so it is safe to access vertex(0).

    // Check whether A contains B, or A and B contain each other's boundaries.
    // (Note that A contains all the vertices of B in either case.)
    if (_subregionBound.contains(b._bound) || _bound.unite(b._bound).isFull()) {
      if (contains(b.vertex(0))) return true;
    }
    // Check whether B contains A.
    if (b._subregionBound.contains(_bound)) {
      if (b.contains(vertex(0))) return true;
    }
    return false;
  }

  // Return true if two loops have the same vertices in the same linear order
  // (i.e., cyclic rotations are not allowed).
  bool equals(in S2Loop b) const {
    if (numVertices() != b.numVertices()) return false;
    for (int i = 0; i < numVertices(); ++i) {
      if (vertex(i) != b.vertex(i)) return false;
    }
    return true;
  }

  // Return true if two loops have the same boundary.  This is true if and
  // only if the loops have the same vertices in the same cyclic order (i.e.,
  // the vertices may be cyclically rotated).  The empty and full loops are
  // considered to have different boundaries.
  bool boundaryEquals(in S2Loop b) const {
    if (numVertices() != b.numVertices()) return false;

    // Special case to handle empty or full loops.  Since they have the same
    // number of vertices, if one loop is empty/full then so is the other.
    if (isEmptyOrFull()) return isEmpty() == b.isEmpty();

    for (int offset = 0; offset < numVertices(); ++offset) {
      if (vertex(offset) == b.vertex(0)) {
        // There is at most one starting offset since loop vertices are unique.
        for (int i = 0; i < numVertices(); ++i) {
          if (vertex(i + offset) != b.vertex(i)) return false;
        }
        return true;
      }
    }
    return false;
  }

  // Return true if two loops have the same boundary except for vertex
  // perturbations.  More precisely, the vertices in the two loops must be in
  // the same cyclic order, and corresponding vertex pairs must be separated
  // by no more than "max_error".
  bool boundaryApproxEquals(in S2Loop b, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    if (numVertices() != b.numVertices()) return false;

    // Special case to handle empty or full loops.  Since they have the same
    // number of vertices, if one loop is empty/full then so is the other.
    if (isEmptyOrFull()) return isEmpty() == b.isEmpty();

    for (int offset = 0; offset < numVertices(); ++offset) {
      if (approxEquals(vertex(offset), b.vertex(0), max_error)) {
        bool success = true;
        for (int i = 0; i < numVertices(); ++i) {
          if (!approxEquals(vertex(i + offset), b.vertex(i), max_error)) {
            success = false;
            break;
          }
        }
        if (success) return true;
        // Otherwise continue looping.  There may be more than one candidate
        // starting offset since vertices are only matched approximately.
      }
    }
    return false;
  }

  // Return true if the two loop boundaries are within "max_error" of each
  // other along their entire lengths.  The two loops may have different
  // numbers of vertices.  More precisely, this method returns true if the two
  // loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
  // distance(a(t), b(t)) <= max_error for all t.  You can think of this as
  // testing whether it is possible to drive two cars all the way around the
  // two loops such that no car ever goes backward and the cars are always
  // within "max_error" of each other.
  bool boundaryNear(in S2Loop b, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    // Special case to handle empty or full loops.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return (isEmpty() && b.isEmpty()) || (isFull() && b.isFull());
    }

    for (int a_offset = 0; a_offset < numVertices(); ++a_offset) {
      if (matchBoundaries(this, b, a_offset, max_error)) return true;
    }
    return false;
  }

  // This method computes the oriented surface integral of some quantity f(x)
  // over the loop interior, given a function f_tri(A,B,C) that returns the
  // corresponding integral over the spherical triangle ABC.  Here "oriented
  // surface integral" means:
  //
  // (1) f_tri(A,B,C) must be the integral of f if ABC is counterclockwise,
  //     and the integral of -f if ABC is clockwise.
  //
  // (2) The result of this function is *either* the integral of f over the
  //     loop interior, or the integral of (-f) over the loop exterior.
  //
  // Note that there are at least two common situations where it easy to work
  // around property (2) above:
  //
  //  - If the integral of f over the entire sphere is zero, then it doesn't
  //    matter which case is returned because they are always equal.
  //
  //  - If f is non-negative, then it is easy to detect when the integral over
  //    the loop exterior has been returned, and the integral over the loop
  //    interior can be obtained by adding the integral of f over the entire
  //    unit sphere (a constant) to the result.
  //
  // Also requires that the default constructor for T must initialize the
  // value to zero.  (This is true for built-in types such as "double".)
  T getSurfaceIntegral(T)(
      T function(in Vector!(double, 3), in Vector!(double, 3), in Vector!(double, 3)) fTri) const {
    // We sum "f_tri" over a collection T of oriented triangles, possibly
    // overlapping.  Let the sign of a triangle be +1 if it is CCW and -1
    // otherwise, and let the sign of a point "x" be the sum of the signs of the
    // triangles containing "x".  Then the collection of triangles T is chosen
    // such that either:
    //
    //  (1) Each point in the loop interior has sign +1, and sign 0 otherwise; or
    //  (2) Each point in the loop exterior has sign -1, and sign 0 otherwise.
    //
    // The triangles basically consist of a "fan" from vertex 0 to every loop
    // edge that does not include vertex 0.  These triangles will always satisfy
    // either (1) or (2).  However, what makes this a bit tricky is that
    // spherical edges become numerically unstable as their length approaches
    // 180 degrees.  Of course there is not much we can do if the loop itself
    // contains such edges, but we would like to make sure that all the triangle
    // edges under our control (i.e., the non-loop edges) are stable.  For
    // example, consider a loop around the equator consisting of four equally
    // spaced points.  This is a well-defined loop, but we cannot just split it
    // into two triangles by connecting vertex 0 to vertex 2.
    //
    // We handle this type of situation by moving the origin of the triangle fan
    // whenever we are about to create an unstable edge.  We choose a new
    // location for the origin such that all relevant edges are stable.  We also
    // create extra triangles with the appropriate orientation so that the sum
    // of the triangle signs is still correct at every point.

    // The maximum length of an edge for it to be considered numerically stable.
    // The exact value is fairly arbitrary since it depends on the stability of
    // the "f_tri" function.  The value below is quite conservative but could be
    // reduced further if desired.
    const auto kMaxLength = S1ChordAngle.fromRadians(M_PI - 1e-5);

    // The default constructor for T must initialize the value to zero.
    // (This is true for built-in types such as "double".)
    T sum = to!T(0);
    S2Point origin = vertex(0);
    for (int i = 1; i + 1 < numVertices(); ++i) {
      // Let V_i be vertex(i), let O be the current origin, and let length(A,B)
      // be the length of edge (A,B).  At the start of each loop iteration, the
      // "leading edge" of the triangle fan is (O,V_i), and we want to extend
      // the triangle fan so that the leading edge is (O,V_i+1).
      //
      // Invariants:
      //  1. length(O,V_i) < kMaxLength for all (i > 1).
      //  2. Either O == V_0, or O is approximately perpendicular to V_0.
      //  3. "sum" is the oriented integral of f over the area defined by
      //     (O, V_0, V_1, ..., V_i).
      enforce(i == 1 || S1ChordAngle(origin, vertex(i)) < kMaxLength);
      enforce(origin == vertex(0) || fabs(origin.dotProd(vertex(0))) < 1e-15);

      if (S1ChordAngle(vertex(i + 1), origin) > kMaxLength) {
        // We are about to create an unstable edge, so choose a new origin O'
        // for the triangle fan.
        S2Point old_origin = origin;
        if (origin == vertex(0)) {
          // The following point is well-separated from V_i and V_0 (and
          // therefore V_i+1 as well).
          origin = robustCrossProd(vertex(0), vertex(i)).normalize();
        } else if (S1ChordAngle(vertex(i), vertex(0)) < kMaxLength) {
          // All edges of the triangle (O, V_0, V_i) are stable, so we can
          // revert to using V_0 as the origin.
          origin = vertex(0);
        } else {
          // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
          // perpendicular.  Therefore V_0.CrossProd(O) is approximately
          // perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
          // this point O' as the new origin.
          origin = vertex(0).crossProd(old_origin);

          // Advance the edge (V_0,O) to (V_0,O').
          sum += fTri(vertex(0), old_origin, origin);
        }
        // Advance the edge (O,V_i) to (O',V_i).
        sum += fTri(old_origin, vertex(i), origin);
      }
      // Advance the edge (O,V_i) to (O,V_i+1).
      sum += fTri(origin, vertex(i), vertex(i+1));
    }
    // If the origin is not V_0, we need to sum one more triangle.
    if (origin != vertex(0)) {
      // Advance the edge (O,V_n-1) to (O,V_0).
      sum += fTri(origin, vertex(numVertices() - 1), vertex(0));
    }
    return sum;
  }

  // Constructs a regular polygon with the given number of vertices, all
  // located on a circle of the specified radius around "center".  The radius
  // is the actual distance from "center" to each vertex.
  static S2Loop makeRegularLoop(in S2Point center, S1Angle radius, int num_vertices) {
    Matrix3x3_d m;
    getFrame(center, m);  // TODO(ericv): Return by value
    return makeRegularLoop(m, radius, num_vertices);
  }

  // Like the function above, but this version constructs a loop centered
  // around the z-axis of the given coordinate frame, with the first vertex in
  // the direction of the positive x-axis.  (This allows the loop to be
  // rotated for testing purposes.)
  static S2Loop makeRegularLoop(in Matrix3x3_d frame, S1Angle radius, int num_vertices) {
    // We construct the loop in the given frame coordinates, with the center at
    // (0, 0, 1).  For a loop of radius "r", the loop vertices have the form
    // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r).  The distance on the
    // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
    double z = cos(radius.radians());
    double r = sin(radius.radians());
    double radian_step = 2 * M_PI / num_vertices;
    S2Point[] vertices;
    for (int i = 0; i < num_vertices; ++i) {
      double angle = i * radian_step;
      auto p = S2Point(r * cos(angle), r * sin(angle), z);
      vertices ~= fromFrame(frame, p).normalize();
    }
    return new S2Loop(vertices);
  }

  // Returns the total number of bytes used by the loop.
  size_t spaceUsed() const {
    size_t size = this.classinfo.m_init.sizeof;
    size += numVertices() * S2Point.sizeof;
    // index_ itself is already included in sizeof(*this).
    size += _index.spaceUsed() - _index.sizeof;
    return size;
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  override
  S2Loop clone() const {
    return new S2Loop(this);
  }

  // GetRectBound() returns essentially tight results, while GetCapBound()
  // might have a lot of extra padding.  Both bounds are conservative in that
  // if the loop contains a point P, then the bound contains P also.
  override
  S2Cap getCapBound() const {
    return _bound.getCapBound();
  }

  override
  S2LatLngRect getRectBound() {
    return _bound;
  }

  override
  void getCellUnionBound(out S2CellId[] cell_ids) const {
    return getCapBound().getCellUnionBound(cell_ids);
  }

  override
  bool contains(in S2Cell target) {
    auto it = new MutableS2ShapeIndex.Iterator(_index);
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If "target" is disjoint from all index cells, it is not contained.
    // Similarly, if "target" is subdivided into one or more index cells then it
    // is not contained, since index cells are subdivided only if they (nearly)
    // intersect a sufficient number of edges.  (But note that if "target" itself
    // is an index cell then it may be contained, since it could be a cell with
    // no edges in the loop interior.)
    if (relation != S2ShapeIndex.CellRelation.INDEXED) return false;

    // Otherwise check if any edges intersect "target".
    if (boundaryApproxIntersects(it, target)) return false;

    // Otherwise check if the loop contains the center of "target".
    return contains(it, target.getCenter());
  }

  override
  bool mayIntersect(in S2Cell target) {
    auto it = new MutableS2ShapeIndex.Iterator(_index);
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If "target" does not overlap any index cell, there is no intersection.
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) return false;

    // If "target" is subdivided into one or more index cells, there is an
    // intersection to within the S2ShapeIndex error bound (see Contains).
    if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) return true;

    // If "target" is an index cell, there is an intersection because index cells
    // are created only if they have at least one edge or they are entirely
    // contained by the loop.
    if (it.id() == target.id()) return true;

    // Otherwise check if any edges intersect "target".
    if (boundaryApproxIntersects(it, target)) return true;

    // Otherwise check if the loop contains the center of "target".
    return contains(it, target.getCenter());
  }

  // The point 'p' does not need to be normalized.
  override
  bool contains(in S2Point p) {
    // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
    // worthwhile only when it might allow us to delay building the index.
    if (!_index.isFresh() && !_bound.contains(p)) return false;

    // For small loops it is faster to just check all the crossings.  We also
    // use this method during loop initialization because InitOriginAndBound()
    // calls Contains() before InitIndex().  Otherwise, we keep track of the
    // number of calls to Contains() and only build the index when enough calls
    // have been made so that we think it is worth the effort.  Note that the
    // code below is structured so that if many calls are made in parallel only
    // one thread builds the index, while the rest continue using brute force
    // until the index is actually available.
    //
    // The constants below were tuned using the benchmarks.  It turns out that
    // building the index costs roughly 50x as much as Contains().  (The ratio
    // increases slowly from 46x with 64 points to 61x with 256k points.)  The
    // textbook approach to this problem would be to wait until the cumulative
    // time we would have saved with an index approximately equals the cost of
    // building the index, and then build it.  (This gives the optimal
    // competitive ratio of 2; look up "competitive algorithms" for details.)
    // We set the limit somewhat lower than this (20 rather than 50) because
    // building the index may be forced anyway by other API calls, and so we
    // want to err on the side of building it too early.

    enum int kMaxBruteForceVertices = 32;
    enum int kMaxUnindexedContainsCalls = 20;  // See notes above.
    if (_index.numShapeIds() == 0  // InitIndex() not called yet
        || numVertices() <= kMaxBruteForceVertices
        || (!_index.isFresh() && atomicOp!"+="(_unindexedContainsCalls, 1) != kMaxUnindexedContainsCalls)) {
      return bruteForceContains(p);
    }
    // Otherwise we look up the S2ShapeIndex cell containing this point.  Note
    // the index is built automatically the first time an iterator is created.
    auto it = new MutableS2ShapeIndex.Iterator(_index);
    if (!it.locate(p)) return false;
    return contains(it, p);
  }

  // Appends a serialized representation of the S2Loop to "encoder".
  //
  // Generally clients should not use S2Loop::Encode().  Instead they should
  // encode an S2Polygon, which unlike this method supports (lossless)
  // compression.
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  // void Encode(Encoder* const encoder) const;

  // Decodes a loop encoded with Encode() or the private method
  // EncodeCompressed() (used by the S2Polygon encoder).  Returns true on
  // success.
  //
  // This method may be called with loops that have already been initialized.
  // bool Decode(Decoder* const decoder);

  // Provides the same functionality as Decode, except that decoded regions
  // are allowed to point directly into the Decoder's memory buffer rather
  // than copying the data.  This can be much faster, but the decoded loop is
  // only valid within the scope (lifetime) of the Decoder's memory buffer.
  // bool DecodeWithinScope(Decoder* const decoder);

  ////////////////////////////////////////////////////////////////////////
  // Methods intended primarily for use by the S2Polygon implementation:

  // Given two loops of a polygon, return true if A contains B.  This version
  // of Contains() is cheap because it does not test for edge intersections.
  // The loops must meet all the S2Polygon requirements; for example this
  // implies that their boundaries may not cross or have any shared edges
  // (although they may have shared vertices).
  bool containsNested(S2Loop b) {
    if (!_subregionBound.contains(b._bound)) return false;

    // Special cases to handle either loop being empty or full.  Also bail out
    // when B has no vertices to avoid heap overflow on the vertex(1) call
    // below.  (This method is called during polygon initialization before the
    // client has an opportunity to call IsValid().)
    if (isEmptyOrFull() || b.numVertices() < 2) {
      return isFull() || b.isEmpty();
    }

    // We are given that A and B do not share any edges, and that either one
    // loop contains the other or they do not intersect.
    int m = findVertex(b.vertex(1));
    if (m < 0) {
      // Since b->vertex(1) is not shared, we can check whether A contains it.
      return contains(b.vertex(1));
    }
    // Check whether the edge order around b->vertex(1) is compatible with
    // A containing B.
    return wedgeContains(vertex(m-1), vertex(m), vertex(m+1), b.vertex(0), b.vertex(2));
  }

  // Return +1 if A contains the boundary of B, -1 if A excludes the boundary
  // of B, and 0 if the boundaries of A and B cross.  Shared edges are handled
  // as follows: If XY is a shared edge, define Reversed(XY) to be true if XY
  // appears in opposite directions in A and B.  Then A contains XY if and
  // only if Reversed(XY) == B->is_hole().  (Intuitively, this checks whether
  // A contains a vanishingly small region extending from the boundary of B
  // toward the interior of the polygon to which loop B belongs.)
  //
  // This method is used for testing containment and intersection of
  // multi-loop polygons.  Note that this method is not symmetric, since the
  // result depends on the direction of loop A but not on the direction of
  // loop B (in the absence of shared edges).
  //
  // REQUIRES: neither loop is empty.
  // REQUIRES: if b->is_full(), then !b->is_hole().
  int compareBoundary(S2Loop b)
  in {
    assert(!isEmpty() && !b.isEmpty());
    assert(!b.isFull() || !b.isHole());
  } body {
    // The bounds must intersect for containment or crossing.
    if (!_bound.intersects(b._bound)) return -1;

    // Full loops are handled as though the loop surrounded the entire sphere.
    if (isFull()) return 1;
    if (b.isFull()) return -1;

    // Check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    auto relation = new CompareBoundaryRelation(b.isHole());
    if (hasCrossingRelation(this, b, relation)) return 0;
    if (relation.foundSharedVertex()) {
      return relation.containsEdge() ? 1 : -1;
    }

    // There are no edge intersections or shared vertices, so we can check
    // whether A contains an arbitrary vertex of B.
    return contains(b.vertex(0)) ? 1 : -1;
  }

  // Given two loops whose boundaries do not cross (see CompareBoundary),
  // return true if A contains the boundary of B.  If "reverse_b" is true, the
  // boundary of B is reversed first (which only affects the result when there
  // are shared edges).  This method is cheaper than CompareBoundary() because
  // it does not test for edge intersections.
  //
  // REQUIRES: neither loop is empty.
  // REQUIRES: if b->is_full(), then reverse_b == false.
  package bool containsNonCrossingBoundary(in S2Loop b, bool reverse_b)
  in {
    assert(!isEmpty() && !b.isEmpty());
    assert(!b.isFull() || !reverse_b);
  } body {
    // The bounds must intersect for containment.
    if (!_bound.intersects(b._bound)) return false;

    // Full loops are handled as though the loop surrounded the entire sphere.
    if (isFull()) return true;
    if (b.isFull()) return false;

    int m = findVertex(b.vertex(0));
    if (m < 0) {
      // Since vertex b0 is not shared, we can check whether A contains it.
      return contains(b.vertex(0));
    }
    // Otherwise check whether the edge (b0, b1) is contained by A.
    return wedgeContainsSemiwedge(vertex(m-1), vertex(m), vertex(m+1), b.vertex(1), reverse_b);
  }

  // Wrapper class for indexing a loop (see S2ShapeIndex).  Once this object
  // is inserted into an S2ShapeIndex it is owned by that index, and will be
  // automatically deleted when no longer needed by the index.  Note that this
  // class does not take ownership of the loop itself (see OwningShape below).
  // You can also subtype this class to store additional data (see S2Shape for
  // details).
  static class Shape : S2Shape {
  public:
    // Must call Init().
    this() {
      _loop = null;
    }

    // Initialize the shape.  Does not take ownership of "loop".
    this(S2Loop loop) {
      initialize(loop);
    }

    void initialize(S2Loop loop) {
      _loop = loop;
    }

    const(S2Loop) loop() const {
      return _loop;
    }

    // S2Shape interface:
    final override
    int numEdges() const {
      return _loop.isEmptyOrFull() ? 0 : _loop.numVertices();
    }

    final override
    Edge edge(int e) const {
      return S2Shape.Edge(_loop.vertex(e), _loop.vertex(e + 1));
    }

    final override
    int dimension() const {
      return 2;
    }

    final override
    S2Shape.ReferencePoint getReferencePoint() const {
      return S2Shape.ReferencePoint(origin(), _loop.containsOrigin());
    }

    final override
    int numChains() const {
      return _loop.isEmptyOrFull() ? 0 : 1;
    }

    final override
    Chain chain(int i) const
    in {
      assert(i == 0);
    } body {
      return Chain(0, numEdges());
    }

    final override
    Edge chainEdge(int i, int j) const
    in {
      assert(i == 0);
    } body {
      return Edge(_loop.vertex(j), _loop.vertex(j + 1));
    }

    final override
    ChainPosition chainPosition(int e) const {
      return ChainPosition(0, e);
    }

   private:
    S2Loop _loop;
  }

  // Like Shape, except that the S2Loop is automatically deleted when this
  // object is deleted by the S2ShapeIndex.  This is useful when an S2Loop
  // is constructed solely for the purpose of indexing it.
  //
  // Not needed in GC languages.
  // class OwningShape : public Shape { ... }

  int opCmp(in S2Loop other) const {
    import std.algorithm : cmp;
    return cmp(_vertices, other._vertices);
  }

package:
  // Return true if this loop contains S2::Origin().
  bool containsOrigin() const {
    return _originInside;
  }

private:
  // Internal copy constructor used only by Clone() that makes a deep copy of
  // its argument.
  this(in S2Loop src) {
    this();
    _depth = src._depth;
    _vertices = src._vertices.dup;
    _s2DebugOverride = src._s2DebugOverride;
    _originInside = src._originInside;
    _unindexedContainsCalls = 0;
    _bound = src._bound.clone();
    _subregionBound = src._subregionBound.clone();
    initIndex();
  }

  // Any single-vertex loop is interpreted as being either the empty loop or the
  // full loop, depending on whether the vertex is in the northern or southern
  // hemisphere respectively.

  // The single vertex in the "empty loop" vertex chain.
  static S2Point emptyVertex() {
    return S2Point(0, 0, 1);
  }

  // The single vertex in the "full loop" vertex chain.
  static S2Point fullVertex() {
    return S2Point(0, 0, -1);
  }

  void initOriginAndBound() {
    import s2.s2predicates : orderedCCW;
    import s2.s2pointutil : ortho;
    if (numVertices() < 3) {
      // Check for the special "empty" and "full" loops (which have one vertex).
      if (!isEmptyOrFull()) {
        _originInside = false;
        return;  // Bail out without trying to access non-existent vertices.
      }
      // If the vertex is in the southern hemisphere then the loop is full,
      // otherwise it is empty.
      _originInside = (vertex(0).z() < 0);
    } else {
      // Point containment testing is done by counting edge crossings starting
      // at a fixed point on the sphere (S2::Origin()).  Historically this was
      // important, but it is now no longer necessary, and it may be worthwhile
      // experimenting with using a loop vertex as the reference point.  In any
      // case, we need to know whether the reference point (S2::Origin) is
      // inside or outside the loop before we can construct the S2ShapeIndex.
      // We do this by first guessing that it is outside, and then seeing
      // whether we get the correct containment result for vertex 1.  If the
      // result is incorrect, the origin must be inside the loop.
      //
      // A loop with consecutive vertices A,B,C contains vertex B if and only if
      // the fixed vector R = S2::Ortho(B) is contained by the wedge ABC.  The
      // wedge is closed at A and open at C, i.e. the point B is inside the loop
      // if A=R but not if C=R.  This convention is required for compatibility
      // with S2::VertexCrossing.  (Note that we can't use S2::Origin()
      // as the fixed vector because of the possibility that B == S2::Origin().)
      //
      // TODO(ericv): Investigate using vertex(0) as the reference point.

      _originInside = false;  // Initialize before calling Contains().
      bool v1_inside = orderedCCW(ortho(vertex(1)), vertex(0), vertex(2), vertex(1));
      // Note that Contains(S2Point) only does a bounds check once InitIndex()
      // has been called, so it doesn't matter that bound_ is undefined here.
      if (v1_inside != contains(vertex(1))) {
        _originInside = true;
      }
    }
    // We *must* call InitBound() before InitIndex(), because InitBound() calls
    // Contains(S2Point), and Contains(S2Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    //
    // TODO(ericv): When fewer S2Loop methods depend on internal bounds checks,
    // consider computing the bound on demand as well.
    initBound();
    initIndex();
  }

  void initBound() {
    // Check for the special "empty" and "full" loops.
    if (isEmptyOrFull()) {
      if (isEmpty()) {
        _subregionBound = _bound = S2LatLngRect.empty();
      } else {
        _subregionBound = _bound = S2LatLngRect.full();
      }
      return;
    }

    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices.  First, the maximal latitude may be
    // attained along the interior of an edge.  Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe).  Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.

    auto bounder = new S2LatLngRectBounder();
    for (int i = 0; i <= numVertices(); ++i) {
      bounder.addPoint(vertex(i));
    }
    S2LatLngRect b = bounder.getBound();
    if (contains(S2Point(0, 0, 1))) {
      b = new S2LatLngRect(R1Interval(b.lat().lo(), M_PI_2), S1Interval.full());
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.lng().is_full() due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.lng().is_full().
    if (b.lng().isFull() && contains(S2Point(0, 0, -1))) {
      b.mutableLat().setLo(-M_PI_2);
    }
    _bound = b;
    _subregionBound = S2LatLngRectBounder.expandForSubregions(_bound);
  }

  void initIndex() {
    import s2.s2debug : flagsS2Debug;
    _index.add(new Shape(this));
    if (!LAZY_INDEXING) {
      _index.forceBuild();
    }
    if (flagsS2Debug && _s2DebugOverride == S2Debug.ALLOW) {
      // Note that FLAGS_s2debug is false in optimized builds (by default).
      enforce(isValid());
    }
  }

  // A version of Contains(S2Point) that does not use the S2ShapeIndex.
  // Used by the S2Polygon implementation.
  package bool bruteForceContains(in S2Point p) const {
    // Empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if (numVertices() < 3) return _originInside;

    S2Point origin = origin();
    auto crosser = new S2CopyingEdgeCrosser(origin, p, vertex(0));
    bool inside = _originInside;
    for (int i = 1; i <= numVertices(); ++i) {
      inside ^= crosser.edgeOrVertexCrossing(vertex(i));
    }
    return inside;
  }

  // Internal implementation of the Decode and DecodeWithinScope methods above.
  // If within_scope is true, memory is allocated for vertices_ and data
  // is copied from the decoder using std::copy. If it is false, vertices_
  // will point to the memory area inside the decoder, and the field
  // owns_vertices_ is set to false.
  //bool DecodeInternal(Decoder* const decoder, bool within_scope);

  // Converts the loop vertices to the S2XYZFaceSiTi format and store the result
  // in the given array, which must be large enough to store all the vertices.
  // void getXYZFaceSiTiVertices(S2XYZFaceSiTi vertices) const {
  //   for (int i = 0; i < num_vertices(); ++i) {
  //     vertices[i].xyz = vertex(i);
  //     vertices[i].cell_level = S2::XYZtoFaceSiTi(vertices[i].xyz,
  //         &vertices[i].face, &vertices[i].si, &vertices[i].ti);
  //   }
  // }

  // Encode the loop's vertices using S2EncodePointsCompressed.  Uses
  // approximately 8 bytes for the first vertex, going down to less than 4 bytes
  // per vertex on Google's geographic repository, plus 24 bytes per vertex that
  // does not correspond to the center of a cell at level 'snap_level'. The loop
  // vertices must first be converted to the S2XYZFaceSiTi format with
  // GetXYZFaceSiTiVertices.
  //
  // REQUIRES: the loop is initialized and valid.
  // void EncodeCompressed(Encoder* encoder, const S2XYZFaceSiTi* vertices,
  //                       int snap_level) const;

  // Decode a loop encoded with EncodeCompressed. The parameters must be the
  // same as the one used when EncodeCompressed was called.
  // bool DecodeCompressed(Decoder* decoder, int snap_level);

  // Returns a bitset of properties used by EncodeCompressed
  // to efficiently encode boolean values.  Properties are
  // origin_inside and whether the bound was encoded.
  // std::bitset<2> GetCompressedEncodingProperties() const;

  // Given an iterator that is already positioned at the S2ShapeIndexCell
  // containing "p", returns Contains(p).
  bool contains(in MutableS2ShapeIndex.Iterator it, in S2Point p) const {
    // Test containment by drawing a line segment from the cell center to the
    // given point and counting edge crossings.
    const(S2ClippedShape) a_clipped = it.cell().clipped(0);
    bool inside = a_clipped.containsCenter();
    int a_num_edges = a_clipped.numEdges();
    if (a_num_edges > 0) {
      S2Point center = it.center();
      auto crosser = new S2CopyingEdgeCrosser(center, p);
      int ai_prev = -2;
      for (int i = 0; i < a_num_edges; ++i) {
        int ai = a_clipped.edge(i);
        if (ai != ai_prev + 1) crosser.restartAt(vertex(ai));
        ai_prev = ai;
        inside ^= crosser.edgeOrVertexCrossing(vertex(ai+1));
      }
    }
    return inside;
  }

  // Return true if the loop boundary intersects "target".  It may also
  // return true when the loop boundary does not intersect "target" but
  // some edge comes within the worst-case error tolerance.
  //
  // REQUIRES: it.id().contains(target.id())
  // [This condition is true whenever it.Locate(target) returns INDEXED.]
  bool boundaryApproxIntersects(in MutableS2ShapeIndex.Iterator it, in S2Cell target) const
  in {
    assert(it.id().contains(target.id()));
  } body {
    import s2.s2edge_clipping :
        clipToPaddedFace, intersectsRect, FACE_CLIP_ERROR_UV_COORD, INTERSECTS_RECT_ERROR_UV_DIST;

    const S2ClippedShape a_clipped = it.cell().clipped(0);
    int a_num_edges = a_clipped.numEdges();

    // If there are no edges, there is no intersection.
    if (a_num_edges == 0) return false;

    // We can save some work if "target" is the index cell itself.
    if (it.id() == target.id()) return true;

    // Otherwise check whether any of the edges intersect "target".
    static const double kMaxError = FACE_CLIP_ERROR_UV_COORD + INTERSECTS_RECT_ERROR_UV_DIST;
    R2Rect bound = target.getBoundUV().expanded(kMaxError);
    for (int i = 0; i < a_num_edges; ++i) {
      int ai = a_clipped.edge(i);
      R2Point v0, v1;
      if (clipToPaddedFace(vertex(ai), vertex(ai+1), target.face(), kMaxError, v0, v1)
          && intersectsRect(v0, v1, bound)) {
        return true;
      }
    }
    return false;
  }

  // Return an index "first" and a direction "dir" (either +1 or -1) such that
  // the vertex sequence (first, first+dir, ..., first+(n-1)*dir) does not
  // change when the loop vertex order is rotated or inverted.  This allows
  // the loop vertices to be traversed in a canonical order.  The return
  // values are chosen such that (first, ..., first+n*dir) are in the range
  // [0, 2*n-1] as expected by the vertex() method.
  package int getCanonicalFirstVertex(out int dir) const {
    int first = 0;
    int n = numVertices();
    for (int i = 1; i < n; ++i) {
      if (vertex(i) < vertex(first)) first = i;
    }
    if (vertex(first + 1) < vertex(first + n - 1)) {
      dir = 1;
      // 0 <= first <= n-1, so (first+n*dir) <= 2*n-1.
    } else {
      dir = -1;
      first += n;
      // n <= first <= 2*n-1, so (first+n*dir) >= 0.
    }
    return first;
  }

  // Return the index of a vertex at point "p", or -1 if not found.
  // The return value is in the range 1..num_vertices_ if found.
  int findVertex(in S2Point p) {
    if (numVertices() < 10) {
      // Exhaustive search.  Return value must be in the range [1..N].
      for (int i = 1; i <= numVertices(); ++i) {
        if (vertex(i) == p) return i;
      }
      return -1;
    }
    auto it = new MutableS2ShapeIndex.Iterator(_index);
    if (!it.locate(p)) return -1;

    const(S2ClippedShape) a_clipped = it.cell().clipped(0);
    for (int i = a_clipped.numEdges() - 1; i >= 0; --i) {
      int ai = a_clipped.edge(i);
      // Return value must be in the range [1..N].
      if (vertex(ai) == p) return (ai == 0) ? numVertices() : ai;
      if (vertex(ai+1) == p) return ai+1;
    }
    return -1;
  }

  // This method checks all edges of loop A for intersection against all edges
  // of loop B.  If there is any shared vertex, the wedges centered at this
  // vertex are sent to "relation".
  //
  // If the two loop boundaries cross, this method is guaranteed to return
  // true.  It also returns true in certain cases if the loop relationship is
  // equivalent to crossing.  For example, if the relation is Contains() and a
  // point P is found such that B contains P but A does not contain P, this
  // method will return true to indicate that the result is the same as though
  // a pair of crossing edges were found (since Contains() returns false in
  // both cases).
  //
  // See Contains(), Intersects() and CompareBoundary() for the three uses of
  // this function.
  static bool hasCrossingRelation(S2Loop a, S2Loop b, LoopRelation relation) {
    // We look for S2CellId ranges where the indexes of A and B overlap, and
    // then test those edges for crossings.
    auto ai = new RangeIterator(a._index);
    auto bi = new RangeIterator(b._index);
    auto ab = new LoopCrosser(a, b, relation, false);  // Tests edges of A against B
    auto ba = new LoopCrosser(b, a, relation, true);   // Tests edges of B against A
    int cnt = 0;
    while (!ai.done() || !bi.done()) {
      if (ai.rangeMax() < bi.rangeMin()) {
        // The A and B cells don't overlap, and A precedes B.
        ai.seekTo(bi);
      } else if (bi.rangeMax() < ai.rangeMin()) {
        // The A and B cells don't overlap, and B precedes A.
        bi.seekTo(ai);
      } else {
        // One cell contains the other.  Determine which cell is larger.
        long ab_relation = ai.id().lsb() - bi.id().lsb();
        if (ab_relation > 0) {
          // A's index cell is larger.
          if (ab.hasCrossingRelation(ai, bi)) return true;
        } else if (ab_relation < 0) {
          // B's index cell is larger.
          if (ba.hasCrossingRelation(bi, ai)) return true;
        } else {
          // The A and B cells are the same.  Since the two cells have the same
          // center point P, check whether P satisfies the crossing targets.
          if (ai.containsCenter() == ab.aCrossingTarget()
              && bi.containsCenter() == ab.bCrossingTarget()) {
            return true;
          }
          // Otherwise test all the edge crossings directly.
          if (ai.numEdges() > 0 && bi.numEdges() > 0
              && ab.cellCrossesCell(ai.clipped(), bi.clipped())) {
            return true;
          }
          ai.next();
          bi.next();
        }
      }
    }
    return false;
  }

  // When the loop is modified (Invert(), or Init() called again) then the
  // indexing structures need to be cleared since they become invalid.
  void clearIndex() {
    atomicStore!(MemoryOrder.raw)(_unindexedContainsCalls, 0);
    _index.clear();
  }

  // The nesting depth, if this field belongs to an S2Polygon.  We define it
  // here to optimize field packing.
  int _depth = 0;

  // We store the vertices in an array rather than a vector because we don't
  // need any STL methods, and computing the number of vertices using size()
  // would be relatively expensive (due to division by sizeof(S2Point) == 24).
  // When DecodeWithinScope is used to initialize the loop, we do not
  // take ownership of the memory for vertices_, and the owns_vertices_ field
  // is used to prevent deallocation and overwriting.
  S2Point[] _vertices;

  S2Debug _s2DebugOverride = S2Debug.ALLOW;
  bool _originInside = false;  // Does the loop contain S2::Origin()?

  // In general we build the index the first time it is needed, but we make an
  // exception for Contains(S2Point) because this method has a simple brute
  // force implementation that is also relatively cheap.  For this one method
  // we keep track of the number of calls made and only build the index once
  // enough calls have been made that we think an index would be worthwhile.
  shared int _unindexedContainsCalls;

  // "bound_" is a conservative bound on all points contained by this loop:
  // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
  S2LatLngRect _bound;

  // Since "bound_" is not exact, it is possible that a loop A contains
  // another loop B whose bounds are slightly larger.  "subregion_bound_"
  // has been expanded sufficiently to account for this error, i.e.
  // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
  S2LatLngRect _subregionBound;

  // Spatial index for this loop.
  MutableS2ShapeIndex _index;
}

// Loop relation for Contains().
class ContainsRelation : LoopRelation {
public:
  this() {
    _foundSharedVertex = false;
  }

  bool foundSharedVertex() const {
    return _foundSharedVertex;
  }

  // If A.Contains(P) == false && B.Contains(P) == true, it is equivalent to
  // having an edge crossing (i.e., Contains returns false).
  override
  int aCrossingTarget() const {
    return false;
  }

  override
  int bCrossingTarget() const {
    return true;
  }

  override
  bool wedgesCross(
      in S2Point a0, in S2Point ab1, in S2Point a2,
      in S2Point b0, in S2Point b2) {
    _foundSharedVertex = true;
    return !wedgeContains(a0, ab1, a2, b0, b2);
  }

private:
  bool _foundSharedVertex;
}

// Loop relation for Intersects().
class IntersectsRelation : LoopRelation {
 public:
  this() {
    _foundSharedVertex = false;
  }

  bool foundSharedVertex() const {
    return _foundSharedVertex;
  }

  // If A.Contains(P) == true && B.Contains(P) == true, it is equivalent to
  // having an edge crossing (i.e., Intersects returns true).
  final override
  int aCrossingTarget() const {
    return true;
  }

  final override
  int bCrossingTarget() const {
    return true;
  }

  override
  bool wedgesCross(
      in S2Point a0, in S2Point ab1, in S2Point a2,
      in S2Point b0, in S2Point b2) {
    _foundSharedVertex = true;
    return wedgeIntersects(a0, ab1, a2, b0, b2);
  }

private:
  bool _foundSharedVertex;
}

// Returns true if the wedge (a0, ab1, a2) contains the "semiwedge" defined as
// any non-empty open set of rays immediately CCW from the edge (ab1, b2).  If
// "reverse_b" is true, then substitute "clockwise" for "CCW"; this simulates
// what would happen if the direction of loop B was reversed.
private bool wedgeContainsSemiwedge(
    in S2Point a0, in S2Point ab1, in S2Point a2, in S2Point b2, bool reverse_b) {
  if (b2 == a0 || b2 == a2) {
    // We have a shared or reversed edge.
    return (b2 == a0) == reverse_b;
  } else {
    return orderedCCW(a0, a2, b2, ab1);
  }
}

/// Loop relation for CompareBoundary().
class CompareBoundaryRelation : LoopRelation {
public:
  this(bool reverse_b) {
    _reverseB = reverse_b;
    _foundSharedVertex = false;
    _containsEdge = false;
    _excludesEdge = false;
  }

  bool foundSharedVertex() const {
    return _foundSharedVertex;
  }

  bool containsEdge() const {
    return _containsEdge;
  }

  // The CompareBoundary relation does not have a useful early-exit condition,
  // so we return -1 for both crossing targets.
  //
  // Aside: A possible early exit condition could be based on the following.
  //   If A contains a point of both B and ~B, then A intersects Boundary(B).
  //   If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
  //   So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
  //   the return value is 0, i.e., Boundary(A) intersects Boundary(B).
  // Unfortunately it isn't worth detecting this situation because by the
  // time we have seen a point in all four intersection regions, we are also
  // guaranteed to have seen at least one pair of crossing edges.
  override
  int aCrossingTarget() const {
    return -1;
  }

  override
  int bCrossingTarget() const {
    return -1;
  }

  override
  bool wedgesCross(
      in S2Point a0, in S2Point ab1, in S2Point a2,
      in S2Point b0, in S2Point b2) {
    // Because we don't care about the interior of B, only its boundary, it is
    // sufficient to check whether A contains the semiwedge (ab1, b2).
    _foundSharedVertex = true;
    if (wedgeContainsSemiwedge(a0, ab1, a2, b2, _reverseB)) {
      _containsEdge = true;
    } else {
      _excludesEdge = true;
    }
    return _containsEdge & _excludesEdge;
  }

protected:
  const bool _reverseB;      // True if loop B should be reversed.
  bool _foundSharedVertex;   // True if any wedge was processed.
  bool _containsEdge;        // True if any edge of B is contained by A.
  bool _excludesEdge;        // True if any edge of B is excluded by A.
}

// LoopRelation is an abstract class that defines a relationship between two
// loops (Contains, Intersects, or CompareBoundary).
abstract class LoopRelation {
public:
  this() {}

  // Optionally, a_target() and b_target() can specify an early-exit condition
  // for the loop relation.  If any point P is found such that
  //
  //   A.Contains(P) == a_crossing_target() &&
  //   B.Contains(P) == b_crossing_target()
  //
  // then the loop relation is assumed to be the same as if a pair of crossing
  // edges were found.  For example, the Contains() relation has
  //
  //   a_crossing_target() == 0
  //   b_crossing_target() == 1
  //
  // because if A.Contains(P) == 0 (false) and B.Contains(P) == 1 (true) for
  // any point P, then it is equivalent to finding an edge crossing (i.e.,
  // since Contains() returns false in both cases).
  //
  // Loop relations that do not have an early-exit condition of this form
  // should return -1 for both crossing targets.
  int aCrossingTarget() const;
  int bCrossingTarget() const;

  // Given a vertex "ab1" that is shared between the two loops, return true if
  // the two associated wedges (a0, ab1, b2) and (b0, ab1, b2) are equivalent
  // to an edge crossing.  The loop relation is also allowed to maintain its
  // own internal state, and can return true if it observes any sequence of
  // wedges that are equivalent to an edge crossing.
  bool wedgesCross(
      in S2Point a0, in S2Point ab1,
      in S2Point a2, in S2Point b0,
      in S2Point b2);
}

// RangeIterator is a wrapper over MutableS2ShapeIndex::Iterator with extra
// methods that are useful for merging the contents of two or more
// S2ShapeIndexes.
class RangeIterator {
public:
  // Construct a new RangeIterator positioned at the first cell of the index.
  this(MutableS2ShapeIndex index) {
    _it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
    refresh();
  }

  // The current S2CellId and cell contents.
  S2CellId id() const {
    return _it.id();
  }

  const(S2ShapeIndexCell) cell() const {
    return _it.cell();
  }

  // The min and max leaf cell ids covered by the current cell.  If Done() is
  // true, these methods return a value larger than any valid cell id.
  S2CellId rangeMin() const {
    return _rangeMin;
  }

  S2CellId rangeMax() const {
    return _rangeMax;
  }

  // Various other convenience methods for the current cell.
  const(S2ClippedShape) clipped() const {
    return cell().clipped(0);
  }

  int numEdges() const {
    return clipped().numEdges();
  }

  bool containsCenter() const {
    return clipped().containsCenter();
  }

  void next() {
    _it.next();
    refresh();
  }

  bool done() {
    return _it.done();
  }

  // Position the iterator at the first cell that overlaps or follows
  // "target", i.e. such that range_max() >= target.range_min().
  void seekTo(in RangeIterator target) {
    _it.seek(target.rangeMin());
    // If the current cell does not overlap "target", it is possible that the
    // previous cell is the one we are looking for.  This can only happen when
    // the previous cell contains "target" but has a smaller S2CellId.
    if (_it.done() || _it.id().rangeMin() > target.rangeMax()) {
      if (_it.prev() && _it.id().rangeMax() < target.id()) _it.next();
    }
    refresh();
  }

  // Position the iterator at the first cell that follows "target", i.e. the
  // first cell such that range_min() > target.range_max().
  void seekBeyond(in RangeIterator target) {
    _it.seek(target.rangeMax().next());
    if (!_it.done() && _it.id().rangeMin() <= target.rangeMax()) {
      _it.next();
    }
    refresh();
  }

private:
  // Updates internal state after the iterator has been repositioned.
  void refresh() {
    _rangeMin = id().rangeMin();
    _rangeMax = id().rangeMax();
  }

  MutableS2ShapeIndex.Iterator _it;
  S2CellId _rangeMin, _rangeMax;
  S2ClippedShape _clipped;
}

// LoopCrosser is a helper class for determining whether two loops cross.
// It is instantiated twice for each pair of loops to be tested, once for the
// pair (A,B) and once for the pair (B,A), in order to be able to process
// edges in either loop nesting order.
class LoopCrosser {
public:
  // If "swapped" is true, the loops A and B have been swapped.  This affects
  // how arguments are passed to the given loop relation, since for example
  // A.Contains(B) is not the same as B.Contains(A).
  this(S2Loop a, S2Loop b, LoopRelation relation, bool swapped) {
    _a = a;
    _b = b;
    _relation = relation;
    _swapped = swapped;
    _aCrossingTarget = relation.aCrossingTarget();
    _bCrossingTarget = relation.bCrossingTarget();
    _bQuery = new S2CrossingEdgeQuery(b._index);
    if (swapped) swap(_aCrossingTarget, _bCrossingTarget);
    _crosser = new S2CopyingEdgeCrosser();
  }

  // Return the crossing targets for the loop relation, taking into account
  // whether the loops have been swapped.
  int aCrossingTarget() const {
    return _aCrossingTarget;
  }

  int bCrossingTarget() const {
    return _bCrossingTarget;
  }

  // Given two iterators positioned such that ai->id().Contains(bi->id()),
  // return true if there is a crossing relationship anywhere within ai->id().
  // Specifically, this method returns true if there is an edge crossing, a
  // wedge crossing, or a point P that matches both "crossing targets".
  // Advances both iterators past ai->id().
  bool hasCrossingRelation(RangeIterator ai, RangeIterator bi)
  in {
    assert(ai.id().contains(bi.id()));
  } body {
    if (ai.numEdges() == 0) {
      if (ai.containsCenter() == _aCrossingTarget) {
        // All points within ai->id() satisfy the crossing target for A, so it's
        // worth iterating through the cells of B to see whether any cell
        // centers also satisfy the crossing target for B.
        do {
          if (bi.containsCenter() == _bCrossingTarget) return true;
          bi.next();
        } while (bi.id() <= ai.rangeMax());
      } else {
        // The crossing target for A is not satisfied, so we skip over the cells
        // of B using binary search.
        bi.seekBeyond(ai);
      }
    } else {
      // The current cell of A has at least one edge, so check for crossings.
      if (hasCrossing(ai, bi)) return true;
    }
    ai.next();
    return false;
  }

  // Given two index cells, return true if there are any edge crossings or
  // wedge crossings within those cells.
  bool cellCrossesCell(in S2ClippedShape a_clipped, in S2ClippedShape b_clipped) {
    // Test all edges of "a_clipped" against all edges of "b_clipped".
    int a_num_edges = a_clipped.numEdges();
    for (int i = 0; i < a_num_edges; ++i) {
      startEdge(a_clipped.edge(i));
      if (edgeCrossesCell(b_clipped)) return true;
    }
    return false;
  }

private:
  // Given two iterators positioned such that ai->id().Contains(bi->id()),
  // return true if there is an edge crossing or wedge crosssing anywhere
  // within ai->id().  Advances "bi" (only) past ai->id().
  bool hasCrossing(RangeIterator ai, RangeIterator bi)
  in {
    assert(ai.id().contains(bi.id()));
  } body {
    // If ai->id() intersects many edges of B, then it is faster to use
    // S2CrossingEdgeQuery to narrow down the candidates.  But if it intersects
    // only a few edges, it is faster to check all the crossings directly.
    // We handle this by advancing "bi" and keeping track of how many edges we
    // would need to test.

    enum int kEdgeQueryMinEdges = 20;  // Tuned using benchmarks.
    int total_edges = 0;
    _bCells.length = 0;
    do {
      if (bi.numEdges() > 0) {
        total_edges += bi.numEdges();
        if (total_edges >= kEdgeQueryMinEdges) {
          // There are too many edges to test them directly, so use
          // S2CrossingEdgeQuery.
          if (cellCrossesAnySubcell(ai.clipped(), ai.id())) return true;
          bi.seekBeyond(ai);
          return false;
        }
        _bCells ~= bi.cell();
      }
      bi.next();
    } while (bi.id() <= ai.rangeMax());

    // Test all the edge crossings directly.
    for (int c = 0; c < _bCells.length; ++c) {
      if (cellCrossesCell(ai.clipped(), _bCells[c].clipped(0))) {
        return true;
      }
    }
    return false;
  }

  // Given an index cell of A, return true if there are any edge or wedge
  // crossings with any index cell of B contained within "b_id".
  bool cellCrossesAnySubcell(in S2ClippedShape a_clipped, S2CellId b_id) {
    // Test all edges of "a_clipped" against all edges of B.  The relevant B
    // edges are guaranteed to be children of "b_id", which lets us find the
    // correct index cells more efficiently.
    auto b_root = new S2PaddedCell(b_id, 0);
    int a_num_edges = a_clipped.numEdges();
    for (int i = 0; i < a_num_edges; ++i) {
      int aj = a_clipped.edge(i);
      // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
      // of B that might contain crossing edges.
      _bQuery.getCells(_a.vertex(aj), _a.vertex(aj+1), b_root, _bCells);
      if (_bCells.empty()) continue;
      startEdge(aj);
      for (int c = 0; c < _bCells.length; ++c) {
        if (edgeCrossesCell(_bCells[c].clipped(0))) return true;
      }
    }
    return false;
  }

  // Prepare to check the given edge of loop A for crossings.
  void startEdge(int aj) {
    // Start testing the given edge of A for crossings.
    _crosser.initialize(_a.vertex(aj), _a.vertex(aj+1));
    _aj = aj;
    _bjPrev = -2;
  }

  // Check the current edge of loop A for crossings with all edges of the
  // given index cell of loop B.
  bool edgeCrossesCell(in S2ClippedShape b_clipped) {
    // Test the current edge of A against all edges of "b_clipped".
    int b_num_edges = b_clipped.numEdges();
    for (int j = 0; j < b_num_edges; ++j) {
      int bj = b_clipped.edge(j);
      if (bj != _bjPrev + 1) _crosser.restartAt(_b.vertex(bj));
      _bjPrev = bj;
      int crossing = _crosser.crossingSign(_b.vertex(bj + 1));
      if (crossing < 0) continue;
      if (crossing > 0) return true;
      // We only need to check each shared vertex once, so we only
      // consider the case where a_vertex(aj_+1) == b_.vertex(bj+1).
      if (_a.vertex(_aj+1) == _b.vertex(bj+1)) {
        if (_swapped
            ? _relation.wedgesCross(
                _b.vertex(bj), _b.vertex(bj+1), _b.vertex(bj+2),
                _a.vertex(_aj), _a.vertex(_aj+2))
            : _relation.wedgesCross(
                _a.vertex(_aj), _a.vertex(_aj+1), _a.vertex(_aj+2),
                _b.vertex(bj), _b.vertex(bj+2))) {
          return true;
        }
      }
    }
    return false;
  }

  const(S2Loop) _a;
  const(S2Loop) _b;
  LoopRelation _relation;
  const(bool) _swapped;
  int _aCrossingTarget, _bCrossingTarget;

  // State maintained by StartEdge() and EdgeCrossesCell().
  S2CopyingEdgeCrosser _crosser;
  int _aj, _bjPrev;

  // Temporary data declared here to avoid repeated memory allocations.
  S2CrossingEdgeQuery _bQuery;
  const(S2ShapeIndexCell)[] _bCells;
}

private bool matchBoundaries(in S2Loop a, in S2Loop b, int a_offset, S1Angle max_error) {
  // The state consists of a pair (i,j).  A state transition consists of
  // incrementing either "i" or "j".  "i" can be incremented only if
  // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
  // applies to "j".  The function returns true iff we can proceed all the way
  // around both loops in this way.
  //
  // Note that when "i" and "j" can both be incremented, sometimes only one
  // choice leads to a solution.  We handle this using a stack and
  // backtracking.  We also keep track of which states have already been
  // explored to avoid duplicating work.

  struct Boundary {
    int first;
    int second;
  }
  Boundary[] pending;
  bool[Boundary] done;
  pending ~= Boundary(0, 0);
  while (!pending.empty()) {
    int i = pending.back().first;
    int j = pending.back().second;
    pending.popBack();
    if (i == a.numVertices() && j == b.numVertices()) {
      return true;
    }
    done[Boundary(i, j)] = true;

    // If (i == na && offset == na-1) where na == a->num_vertices(), then
    // then (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
    // So we reduce the range if necessary.
    int io = i + a_offset;
    if (io >= a.numVertices()) {
      io -= a.numVertices();
    }

    if (i < a.numVertices() && Boundary(i + 1, j) !in done
        && getDistance(a.vertex(io + 1), b.vertex(j), b.vertex(j + 1)) <= max_error) {
      pending ~= Boundary(i + 1, j);
    }
    if (j < b.numVertices() && Boundary(i, j + 1) !in done
        && getDistance(b.vertex(j + 1), a.vertex(io), a.vertex(io + 1)) <= max_error) {
      pending ~= Boundary(i, j + 1);
    }
  }
  return false;
}
