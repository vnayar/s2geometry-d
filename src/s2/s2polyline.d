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

module s2.s2polyline;

import s2.logger;
import s2.s1angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2debug;
import s2.s2edge_crosser;
import s2.s2edge_distances : getDistance, interpolateAtDistance, isEdgeBNearEdgeA;
import s2.s2error;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2point;
import s2.s2pointutil : approxEquals, getFrame, isUnitLength, toFrame;
import s2.s2predicates : orderedCCW, sign;
import s2.s2region;
import s2.s2shape;
import s2.util.math.matrix3x3;

import std.math;
import std.algorithm : copy, max, min;
import std.exception : enforce;
import std.range : back, empty, popBack;
import std.typecons : Rebindable;

/**
 * An S2Polyline represents a sequence of zero or more vertices connected by
 * straight edges (geodesics).  Edges of length 0 and 180 degrees are not
 * allowed, i.e. adjacent vertices should not be identical or antipodal.
 */
class S2Polyline : S2Region {
public:
  /// Creates an empty S2Polyline that should be initialized by calling Init() or Decode().
  this() {
    _s2debugOverride = S2Debug.ALLOW;
  }

  /// Convenience constructors that call Init() with the given vertices.
  this(S2Point[] vertices) {
    this(vertices, S2Debug.ALLOW);
  }

  this(S2LatLng[] vertices) {
    this(vertices, S2Debug.ALLOW);
  }

  /**
   * Convenience constructors to disable the automatic validity checking
   * controlled by the --s2debug flag.  Example:
   *
   *   S2Polyline* line = new S2Polyline(vertices, S2Debug::DISABLE);
   *
   * This is equivalent to:
   *
   *   S2Polyline* line = new S2Polyline;
   *   line->set_s2debug_override(S2Debug::DISABLE);
   *   line->Init(vertices);
   *
   * The main reason to use this constructors is if you intend to call
   * IsValid() explicitly.  See set_s2debug_override() for details.
   */
  this(S2Point[] vertices, S2Debug s2debugOverride) {
    _s2debugOverride = s2debugOverride;
    initialize(vertices);
  }

  this(S2LatLng[] vertices, S2Debug s2debugOverride) {
    _s2debugOverride = s2debugOverride;
    initialize(vertices);
  }

  /**
   * Initialize a polyline that connects the given vertices. Empty polylines are
   * allowed.  Adjacent vertices should not be identical or antipodal.  All
   * vertices should be unit length.
   */
  void initialize(S2Point[] vertices) {
    _vertices.length = vertices.length;
    copy(vertices, _vertices);
    if (flagsS2Debug && _s2debugOverride == S2Debug.ALLOW) {
      enforce(isValid());
    }
  }

  /**
   * Convenience initialization function that accepts latitude-longitude
   * coordinates rather than S2Points.
   */
  void initialize(S2LatLng[] vertices) {
    _vertices.length = vertices.length;
    for (int i = 0; i < _vertices.length; ++i) {
      _vertices[i] = vertices[i].toS2Point();
    }
    if (flagsS2Debug && _s2debugOverride == S2Debug.ALLOW) {
      enforce(isValid());
    }
  }

  /**
   * Allows overriding the automatic validity checks controlled by the
   * --s2debug flag.  If this flag is true, then polylines are automatically
   * checked for validity as they are initialized.  The main reason to disable
   * this flag is if you intend to call IsValid() explicitly, like this:
   *
   *   S2Polyline line;
   *   line.set_s2debug_override(S2Debug::DISABLE);
   *   line.Init(...);
   *   if (!line.IsValid()) { ... }
   *
   * Without the call to set_s2debug_override(), invalid data would cause a
   * fatal error in Init() whenever the --s2debug flag is enabled.
   */
  void setS2debugOverride(S2Debug s2debugOverride) {
    _s2debugOverride = s2debugOverride;
  }

  S2Debug s2debugOverride() const {
    return _s2debugOverride;
  }

  /// Return true if the given vertices form a valid polyline.
  bool isValid() const {
    S2Error error;
    if (findValidationError(error)) {
      if (flagsS2Debug) {
        logger.logError(error);
      }
      return false;
    }
    return true;
  }

  /**
   * Returns true if this is *not* a valid polyline and sets "error"
   * appropriately.  Otherwise returns false and leaves "error" unchanged.
   *
   * REQUIRES: error != nullptr
   */
  bool findValidationError(out S2Error error) const {
    // All vertices must be unit length.
    for (int i = 0; i < numVertices(); ++i) {
      if (!isUnitLength(vertex(i))) {
        error.initialize(
            S2Error.Code.NOT_UNIT_LENGTH, "Vertex %d is not unit length", i);
        return true;
      }
    }
    // Adjacent vertices must not be identical or antipodal.
    for (int i = 1; i < numVertices(); ++i) {
      if (vertex(i - 1) == vertex(i)) {
        error.initialize(
            S2Error.Code.DUPLICATE_VERTICES, "Vertices %d and %d are identical", i - 1, i);
        return true;
      }
      if (vertex(i - 1) == -vertex(i)) {
        error.initialize(
            S2Error.Code.ANTIPODAL_VERTICES, "Vertices %d and %d are antipodal", i - 1, i);
        return true;
      }
    }
    return false;
  }

  int numVertices() const {
    return cast(int) _vertices.length;
  }

  ref const(S2Point) vertex(int k) const
  in {
    assert(k >= 0);
    assert(k < _vertices.length);
  } do {
    return _vertices[k];
  }

  const(S2Point[]) vertices() const {
    return _vertices;
  }

  /// Return the length of the polyline.
  S1Angle getLength() const {
    S1Angle length;
    for (int i = 1; i < numVertices(); ++i) {
      length += S1Angle(vertex(i-1), vertex(i));
    }
    return length;
  }

  /**
   * Return the true centroid of the polyline multiplied by the length of the
   * polyline (see s2centroids.h for details on centroids).  The result is not
   * unit length, so you may want to normalize it.
   *
   * Prescaling by the polyline length makes it easy to compute the centroid
   * of several polylines (by simply adding up their centroids).
   */
  S2Point getCentroid() const {
    S2Point centroid;
    for (int i = 1; i < numVertices(); ++i) {
      // The centroid (multiplied by length) is a vector toward the midpoint
      // of the edge, whose length is twice the sin of half the angle between
      // the two vertices.  Defining theta to be this angle, we have:
      S2Point vsum = vertex(i-1) + vertex(i);    // Length == 2*cos(theta)
      S2Point vdiff = vertex(i-1) - vertex(i);   // Length == 2*sin(theta)
      double cos2 = vsum.norm2();
      double sin2 = vdiff.norm2();
      enforce(cos2 > 0.0);  // Otherwise edge is undefined, and result is NaN.
      centroid += sqrt(sin2 / cos2) * vsum;  // Length == 2*sin(theta)
    }
    return centroid;
  }

  /**
   * Return the point whose distance from vertex 0 along the polyline is the
   * given fraction of the polyline's total length.  Fractions less than zero
   * or greater than one are clamped.  The return value is unit length.  This
   * cost of this function is currently linear in the number of vertices.
   * The polyline must not be empty.
   */
  S2Point interpolate(double fraction) const {
    int next_vertex;
    return getSuffix(fraction, next_vertex);
  }

  /**
   * Like Interpolate(), but also return the index of the next polyline
   * vertex after the interpolated point P.  This allows the caller to easily
   * construct a given suffix of the polyline by concatenating P with the
   * polyline vertices starting at "next_vertex".  Note that P is guaranteed
   * to be different than vertex(*next_vertex), so this will never result in
   * a duplicate vertex.
   *
   * The polyline must not be empty.  Note that if "fraction" >= 1.0, then
   * "next_vertex" will be set to num_vertices() (indicating that no vertices
   * from the polyline need to be appended).  The value of "next_vertex" is
   * always between 1 and num_vertices().
   *
   * This method can also be used to construct a prefix of the polyline, by
   * taking the polyline vertices up to "next_vertex - 1" and appending the
   * returned point P if it is different from the last vertex (since in this
   * case there is no guarantee of distinctness).
   */
  S2Point getSuffix(double fraction, out int next_vertex) const
  in {
    assert(numVertices() > 0);
  } do {
    // We intentionally let the (fraction >= 1) case fall through, since
    // we need to handle it in the loop below in any case because of
    // possible roundoff errors.
    if (fraction <= 0) {
      next_vertex = 1;
      return vertex(0);
    }
    S1Angle length_sum;
    for (int i = 1; i < numVertices(); ++i) {
      length_sum += S1Angle(vertex(i-1), vertex(i));
    }
    S1Angle target = fraction * length_sum;
    for (int i = 1; i < numVertices(); ++i) {
      auto length = S1Angle(vertex(i-1), vertex(i));
      if (target < length) {
        // This interpolates with respect to arc length rather than
        // straight-line distance, and produces a unit-length result.
        S2Point result = interpolateAtDistance(target, vertex(i-1), vertex(i));
        // It is possible that (result == vertex(i)) due to rounding errors.
        next_vertex = (result == vertex(i)) ? (i + 1) : i;
        return result;
      }
      target -= length;
    }
    next_vertex = numVertices();
    return vertex(numVertices() - 1);
  }

  /**
   * The inverse operation of GetSuffix/Interpolate.  Given a point on the
   * polyline, returns the ratio of the distance to the point from the
   * beginning of the polyline over the length of the polyline.  The return
   * value is always betwen 0 and 1 inclusive.  See GetSuffix() for the
   * meaning of "next_vertex".
   *
   * The polyline should not be empty.  If it has fewer than 2 vertices, the
   * return value is zero.
   */
  double unInterpolate(in S2Point point, int next_vertex) const
  in {
    assert(numVertices() > 0);
  } do {
    if (numVertices() < 2) {
      return 0;
    }
    S1Angle length_sum;
    for (int i = 1; i < next_vertex; ++i) {
      length_sum += S1Angle(vertex(i-1), vertex(i));
    }
    S1Angle length_to_point = length_sum + S1Angle(vertex(next_vertex-1), point);
    for (int i = next_vertex; i < numVertices(); ++i) {
      length_sum += S1Angle(vertex(i-1), vertex(i));
    }
    // The ratio can be greater than 1.0 due to rounding errors or because the
    // point is not exactly on the polyline.
    return min(1.0, (length_to_point / length_sum).radians());
  }

  /**
   * Given a point, returns a point on the polyline that is closest to the given
   * point.  See GetSuffix() for the meaning of "next_vertex", which is chosen
   * here w.r.t. the projected point as opposed to the interpolated point in
   * GetSuffix().
   *
   * The polyline must be non-empty.
   */
  S2Point project(in S2Point point, out int next_vertex) const
  in {
    assert(numVertices() > 0);
  } do {
    import s2.s2edge_distances : project;

    if (numVertices() == 1) {
      // If there is only one vertex, it is always closest to any given point.
      next_vertex = 1;
      return vertex(0);
    }

    // Initial value larger than any possible distance on the unit sphere.
    S1Angle min_distance = S1Angle.fromRadians(10.0);
    int min_index = -1;

    // Find the line segment in the polyline that is closest to the point given.
    for (int i = 1; i < numVertices(); ++i) {
      S1Angle distance_to_segment = getDistance(point, vertex(i-1), vertex(i));
      if (distance_to_segment < min_distance) {
        min_distance = distance_to_segment;
        min_index = i;
      }
    }
    enforce(min_index != -1);

    // Compute the point on the segment found that is closest to the point given.
    S2Point closest_point = project(point, vertex(min_index - 1), vertex(min_index));

    next_vertex = min_index + (closest_point == vertex(min_index) ? 1 : 0);
    return closest_point;
  }

  /**
   * Returns true if the point given is on the right hand side of the polyline,
   * using a naive definition of "right-hand-sideness" where the point is on
   * the RHS of the polyline iff the point is on the RHS of the line segment in
   * the polyline which it is closest to.
   *
   * The polyline must have at least 2 vertices.
   */
  bool isOnRight(in S2Point point) const
  in {
    assert(numVertices() >= 2);
  } do {
    int next_vertex;
    S2Point closest_point = project(point, next_vertex);

    enforce(next_vertex >= 1);
    enforce(next_vertex <= numVertices());

    // If the closest point C is an interior vertex of the polyline, let B and D
    // be the previous and next vertices.  The given point P is on the right of
    // the polyline (locally) if B, P, D are ordered CCW around vertex C.
    if (closest_point == vertex(next_vertex - 1) && next_vertex > 1
        && next_vertex < numVertices()) {
      if (point == vertex(next_vertex-1))
        return false;  // Polyline vertices are not on the RHS.
      return orderedCCW(vertex(next_vertex-2), point, vertex(next_vertex), vertex(next_vertex-1));
    }

    // Otherwise, the closest point C is incident to exactly one polyline edge.
    // We test the point P against that edge.
    if (next_vertex == numVertices())
      --next_vertex;

    return sign(point, vertex(next_vertex), vertex(next_vertex - 1)) > 0;
  }

  /**
   * Return true if this polyline intersects the given polyline. If the
   * polylines share a vertex they are considered to be intersecting. When a
   * polyline endpoint is the only intersection with the other polyline, the
   * function may return true or false arbitrarily.
   *
   * The running time is quadratic in the number of vertices.  (To intersect
   * polylines more efficiently, or compute the actual intersection geometry,
   * use S2BooleanOperation.)
   */
  bool intersects(S2Polyline line) {
    if (numVertices() <= 0 || line.numVertices() <= 0) {
      return false;
    }

    if (!getRectBound().intersects(line.getRectBound())) {
      return false;
    }

    // TODO(ericv): Use S2ShapeIndex here.
    for (int i = 1; i < numVertices(); ++i) {
      auto crosser = new S2EdgeCrosser(vertex(i - 1), vertex(i), line.vertex(0));
      for (int j = 1; j < line.numVertices(); ++j) {
        if (crosser.crossingSign(line.vertex(j)) >= 0) {
          return true;
        }
      }
    }
    return false;
  }

  /// Reverse the order of the polyline vertices.
  void reverse() {
    import std.algorithm : reverse;
    reverse(_vertices);
  }

  /**
   * Return a subsequence of vertex indices such that the polyline connecting
   * these vertices is never further than "tolerance" from the original
   * polyline.  Provided the first and last vertices are distinct, they are
   * always preserved; if they are not, the subsequence may contain only a
   * single index.
   *
   * Some useful properties of the algorithm:
   *
   *  - It runs in linear time.
   *
   *  - The output is always a valid polyline.  In particular, adjacent
   *    output vertices are never identical or antipodal.
   *
   *  - The method is not optimal, but it tends to produce 2-3% fewer
   *    vertices than the Douglas-Peucker algorithm with the same tolerance.
   *
   *  - The output is *parametrically* equivalent to the original polyline to
   *    within the given tolerance.  For example, if a polyline backtracks on
   *    itself and then proceeds onwards, the backtracking will be preserved
   *    (to within the given tolerance).  This is different than the
   *    Douglas-Peucker algorithm, which only guarantees geometric equivalence.
   *
   * See also S2PolylineSimplifier, which uses the same algorithm but is more
   * efficient and supports more features, and also S2Builder, which can
   * simplify polylines and polygons, supports snapping (e.g. to E7 lat/lng
   * coordinates or S2CellId centers), and can split polylines at intersection
   * points.
   */
  void subsampleVertices(S1Angle tolerance, out int[] indices) const {
    if (numVertices() == 0) return;

    indices ~= 0;
    S1Angle clamped_tolerance = max(tolerance, S1Angle.fromRadians(0.0));
    for (int index = 0; index + 1 < numVertices(); ) {
      int next_index = findEndVertex(this, clamped_tolerance, index);
      // Don't create duplicate adjacent vertices.
      if (vertex(next_index) != vertex(index)) {
        indices ~= next_index;
      }
      index = next_index;
    }
  }

  /// Return true if two polylines are exactly the same.
  override
  bool opEquals(in Object o) const {
    const(S2Polyline) b = cast(S2Polyline) o;
    if (b is null) return false;
    if (numVertices() != b.numVertices()) return false;
    for (int offset = 0; offset < numVertices(); ++offset) {
      if (vertex(offset) != b.vertex(offset)) return false;
    }
    return true;
  }

  /**
   * Return true if two polylines have the same number of vertices, and
   * corresponding vertex pairs are separated by no more than "max_error".
   * (For testing purposes.)
   */
  bool approxEquals(in S2Polyline b, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    import s2.s2pointutil : approxEquals;

    if (numVertices() != b.numVertices()) return false;
    for (int offset = 0; offset < numVertices(); ++offset) {
      if (!approxEquals(vertex(offset), b.vertex(offset), max_error)) {
        return false;
      }
    }
    return true;
  }

  // Return true if "covered" is within "max_error" of a contiguous subpath of
  // this polyline over its entire length.  Specifically, this method returns
  // true if this polyline has parameterization a:[0,1] -> S^2, "covered" has
  // parameterization b:[0,1] -> S^2, and there is a non-decreasing function
  // f:[0,1] -> [0,1] such that distance(a(f(t)), b(t)) <= max_error for all t.
  //
  // You can think of this as testing whether it is possible to drive a car
  // along "covered" and a car along some subpath of this polyline such that no
  // car ever goes backward, and the cars are always within "max_error" of each
  // other.
  //
  // This function is well-defined for empty polylines:
  //    anything.covers(empty) = true
  //    empty.covers(nonempty) = false
  bool nearlyCovers(in S2Polyline covered, S1Angle max_error) const {
    import s2.s2edge_distances : project;

    // NOTE: This algorithm is described assuming that adjacent vertices in a
    // polyline are never at the same point.  That is, the ith and i+1th vertices
    // of a polyline are never at the same point in space.  The implementation
    // does not make this assumption.

    // DEFINITIONS:
    //   - edge "i" of a polyline is the edge from the ith to i+1th vertex.
    //   - covered_j is a polyline consisting of edges 0 through j of "covered."
    //   - this_i is a polyline consisting of edges 0 through i of this polyline.
    //
    // A search state is represented as an (int, int, bool) tuple, (i, j,
    // i_in_progress).  Using the "drive a car" analogy from the header comment, a
    // search state signifies that you can drive one car along "covered" from its
    // first vertex through a point on its jth edge, and another car along this
    // polyline from some point on or before its ith edge to a to a point on its
    // ith edge, such that no car ever goes backward, and the cars are always
    // within "max_error" of each other.  If i_in_progress is true, it means that
    // you can definitely drive along "covered" through the jth vertex (beginning
    // of the jth edge). Otherwise, you can definitely drive along "covered"
    // through the point on the jth edge of "covered" closest to the ith vertex of
    // this polyline.
    //
    // The algorithm begins by finding all edges of this polyline that are within
    // "max_error" of the first vertex of "covered," and adding search states
    // representing all of these possible starting states to the stack of
    // "pending" states.
    //
    // The algorithm proceeds by popping the next pending state,
    // (i,j,i_in_progress), off of the stack.  First it checks to see if that
    // state represents finding a valid covering of "covered" and returns true if
    // so.  Next, if the state represents reaching the end of this polyline
    // without finding a successful covering, the algorithm moves on to the next
    // state in the stack.  Otherwise, if state (i+1,j,false) is valid, it is
    // added to the stack of pending states.  Same for state (i,j+1,true).
    //
    // We need the stack because when "i" and "j" can both be incremented,
    // sometimes only one choice leads to a solution.  We use a set to keep track
    // of visited states to avoid duplicating work.  With the set, the worst-case
    // number of states examined is O(n+m) where n = this->num_vertices() and m =
    // covered.num_vertices().  Without it, the amount of work could be as high as
    // O((n*m)^2).  Using set, the running time is O((n*m) log (n*m)).
    //
    // TODO(user): Benchmark this, and see if the set is worth it.

    if (covered.numVertices() == 0) return true;
    if (numVertices() == 0) return false;

    SearchState[] pending;
    bool[SearchState] done;

    // Find all possible starting states.
    for (int i = 0, next_i = nextDistinctVertex(this, 0), next_next_i;
         next_i < numVertices(); i = next_i, next_i = next_next_i) {
      next_next_i = nextDistinctVertex(this, next_i);
      S2Point closest_point = project(covered.vertex(0), vertex(i), vertex(next_i));

      // In order to avoid duplicate starting states, we exclude the end vertex
      // of each edge *except* for the last non-degenerate edge.
      if ((next_next_i == numVertices() || closest_point != vertex(next_i))
          && S1Angle(closest_point, covered.vertex(0)) <= max_error) {
        pending ~= SearchState(i, 0, true);
      }
    }

    while (!pending.empty()) {
      const SearchState state = pending.back();
      pending.popBack();

      if (state in done) continue;
      done[state] = true;

      const int next_i = nextDistinctVertex(this, state.i);
      const int next_j = nextDistinctVertex(covered, state.j);
      if (next_j == covered.numVertices()) {
        return true;
      } else if (next_i == numVertices()) {
        continue;
      }

      S2Point i_begin, j_begin;
      if (state.iInProgress) {
        j_begin = covered.vertex(state.j);
        i_begin = project(j_begin, vertex(state.i), vertex(next_i));
      } else {
        i_begin = vertex(state.i);
        j_begin = project(i_begin, covered.vertex(state.j), covered.vertex(next_j));
      }

      if (isEdgeBNearEdgeA(j_begin, covered.vertex(next_j), i_begin, vertex(next_i), max_error)) {
        pending ~= SearchState(next_i, state.j, false);
      }
      if (isEdgeBNearEdgeA(i_begin, vertex(next_i), j_begin, covered.vertex(next_j), max_error)) {
        pending ~= SearchState(state.i, next_j, true);
      }
    }
    return false;
  }

  // Returns the total number of bytes used by the polyline.
  size_t spaceUsed() const {
    return this.classinfo.m_init.length + numVertices() * S2Point.sizeof;
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  override
  S2Polyline clone() {
    return new S2Polyline(this);
  }

  override
  void getCellUnionBound(out S2CellId[] cell_ids) {
    return getCapBound().getCellUnionBound(cell_ids);
  }

  override
  S2Cap getCapBound() {
    return getRectBound().getCapBound();
  }

  override
  S2LatLngRect getRectBound() {
    auto bounder = new S2LatLngRectBounder();
    for (int i = 0; i < numVertices(); ++i) {
      bounder.addPoint(vertex(i));
    }
    return bounder.getBound();
  }

  override
  bool contains(in S2Cell cell) const {
    return false;
  }

  override
  bool mayIntersect(in S2Cell cell) const {
    if (numVertices() == 0) return false;

    // We only need to check whether the cell contains vertex 0 for correctness,
    // but these tests are cheap compared to edge crossings so we might as well
    // check all the vertices.
    for (int i = 0; i < numVertices(); ++i) {
      if (cell.contains(vertex(i))) return true;
    }
    S2Point[4] cell_vertices;
    for (int i = 0; i < 4; ++i) {
      cell_vertices[i] = cell.getVertex(i);
    }
    for (int j = 0; j < 4; ++j) {
      auto crosser = new S2EdgeCrosser(cell_vertices[j], cell_vertices[(j + 1) & 3], vertex(0));
      for (int i = 1; i < numVertices(); ++i) {
        if (crosser.crossingSign(vertex(i)) >= 0) {
          // There is a proper crossing, or two vertices were the same.
          return true;
        }
      }
    }
    return false;
  }

  // Always return false, because "containment" is not numerically
  // well-defined except at the polyline vertices.
  override
  bool contains(in S2Point p) const {
    return false;
  }

  // TODO: Implement when Encoder is implemented.
  // Appends a serialized representation of the S2Polyline to "encoder".
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  //void Encode(Encoder* const encoder) const;

  // TODO: Implement when Decoder is implemented.
  // Decodes an S2Polyline encoded with Encode().  Returns true on success.
  //bool Decode(Decoder* const decoder);

  /**
   * Wrapper class for indexing a polyline (see S2ShapeIndex).  Once this
   * object is inserted into an S2ShapeIndex it is owned by that index, and
   * will be automatically deleted when no longer needed by the index.  Note
   * that this class does not take ownership of the polyline itself (see
   * OwningShape below).  You can also subtype this class to store additional
   * data (see S2Shape for details).
   */
  static class Shape : S2Shape {
  public:
    // Must call Init().
    this() {
      _polyline = null;
    }

    /**
     * Initialization.  Does not take ownership of "polyline".
     *
     * Note that a polyline with one vertex is defined to have no edges.  Use
     * S2LaxPolylineShape or S2LaxClosedPolylineShape if you want to define a
     * polyline consisting of a single degenerate edge.
     */
    this(in S2Polyline polyline) {
      init(polyline);
    }

    void init(in S2Polyline polyline) {
      if (polyline.numVertices() == 1) {
        logger.logWarn("S2Polyline::Shape with one vertex has no edges");
      }
      _polyline = polyline;
    }

    const(S2Polyline) polyline() const {
      return _polyline;
    }

    // S2Shape interface:

    override
    int numEdges() const {
      return max(0, _polyline.numVertices() - 1);
    }

    override
    Edge edge(int e) const {
      return Edge(_polyline.vertex(e), _polyline.vertex(e + 1));
    }

    override
    int dimension() const {
      return 1;
    }

    override
    S2Shape.ReferencePoint getReferencePoint() const {
      return S2Shape.ReferencePoint(false);
    }

    override
    int numChains() const {
      return min(1, numEdges());  // Avoid virtual call.
    }

    override
    Chain chain(int i) const
    in {
      assert(i == 0);
    } do {
      return Chain(0, numEdges());  // Avoid virtual call.
    }

    override
    Edge chainEdge(int i, int j) const
    in {
      assert(i == 0);
    } do {
      return Edge(_polyline.vertex(j), _polyline.vertex(j + 1));
    }

    override
    ChainPosition chainPosition(int e) const {
      return ChainPosition(0, e);
    }

   private:
    Rebindable!(const(S2Polyline)) _polyline;
  }

  // Like Shape, except that the S2Polyline is automatically deleted when this
  // object is deleted by the S2ShapeIndex.  This is useful when an S2Polyline
  // is constructed solely for the purpose of indexing it.
  //
  // class OwningShape : Shape -- Not needed in D with GC.

 private:
  // Internal copy constructor used only by Clone() that makes a deep copy of
  // its argument.
  this(in S2Polyline o) {
    _s2debugOverride = o._s2debugOverride;
    _vertices = o._vertices.dup;
  }

  // Allows overriding the automatic validity checking controlled by the
  // --s2debug flag.
  S2Debug _s2debugOverride;

  // We store the vertices in an array rather than a vector because we don't
  // need any STL methods, and computing the number of vertices using size()
  // would be relatively expensive (due to division by sizeof(S2Point) == 24).
  S2Point[] _vertices;

  // Given a polyline, a tolerance distance, and a start index, this function
  // returns the maximal end index such that the line segment between these two
  // vertices passes within "tolerance" of all interior vertices, in order.
  static int findEndVertex(in S2Polyline polyline, S1Angle tolerance, int index)
  in {
    assert(tolerance.radians() >= 0);
    assert((index + 1) < polyline.numVertices());
  } do {

    // The basic idea is to keep track of the "pie wedge" of angles from the
    // starting vertex such that a ray from the starting vertex at that angle
    // will pass through the discs of radius "tolerance" centered around all
    // vertices processed so far.

    // First we define a "coordinate frame" for the tangent and normal spaces
    // at the starting vertex.  Essentially this means picking three
    // orthonormal vectors X,Y,Z such that X and Y span the tangent plane at
    // the starting vertex, and Z is "up".  We use the coordinate frame to
    // define a mapping from 3D direction vectors to a one-dimensional "ray
    // angle" in the range (-Pi, Pi].  The angle of a direction vector is
    // computed by transforming it into the X,Y,Z basis, and then calculating
    // atan2(y,x).  This mapping allows us to represent a wedge of angles as a
    // 1D interval.  Since the interval wraps around, we represent it as an
    // S1Interval, i.e. an interval on the unit circle.
    Matrix3x3_d frame;
    const(S2Point) origin = polyline.vertex(index);
    getFrame(origin, frame);

    // As we go along, we keep track of the current wedge of angles and the
    // distance to the last vertex (which must be non-decreasing).
    S1Interval current_wedge = S1Interval.full();
    double last_distance = 0;

    for (++index; index < polyline.numVertices(); ++index) {
      const(S2Point) candidate = polyline.vertex(index);
      double distance = origin.angle(candidate);

      // We don't allow simplification to create edges longer than 90 degrees,
      // to avoid numeric instability as lengths approach 180 degrees.  (We do
      // need to allow for original edges longer than 90 degrees, though.)
      if (distance > M_PI_2 && last_distance > 0) break;

      // Vertices must be in increasing order along the ray, except for the
      // initial disc around the origin.
      if (distance < last_distance && last_distance > tolerance.radians()) break;
      last_distance = distance;

      // Points that are within the tolerance distance of the origin do not
      // constrain the ray direction, so we can ignore them.
      if (distance <= tolerance.radians()) continue;

      // If the current wedge of angles does not contain the angle to this
      // vertex, then stop right now.  Note that the wedge of possible ray
      // angles is not necessarily empty yet, but we can't continue unless we
      // are willing to backtrack to the last vertex that was contained within
      // the wedge (since we don't create new vertices).  This would be more
      // complicated and also make the worst-case running time more than linear.
      S2Point direction = toFrame(frame, candidate);
      double center = atan2(direction.y(), direction.x());
      if (!current_wedge.contains(center)) break;

      // To determine how this vertex constrains the possible ray angles,
      // consider the triangle ABC where A is the origin, B is the candidate
      // vertex, and C is one of the two tangent points between A and the
      // spherical cap of radius "tolerance" centered at B.  Then from the
      // spherical law of sines, sin(a)/sin(A) = sin(c)/sin(C), where "a" and
      // "c" are the lengths of the edges opposite A and C.  In our case C is a
      // 90 degree angle, therefore A = asin(sin(a) / sin(c)).  Angle A is the
      // half-angle of the allowable wedge.

      double half_angle = asin(sin(tolerance.radians()) / sin(distance));
      S1Interval target = S1Interval.fromPoint(center).expanded(half_angle);
      current_wedge = current_wedge.intersection(target);
      enforce(!current_wedge.isEmpty());
    }
    // We break out of the loop when we reach a vertex index that can't be
    // included in the line segment, so back up by one vertex.
    return index - 1;
  }

  // Return the first i > "index" such that the ith vertex of "pline" is not at
  // the same point as the "index"th vertex.  Returns pline.num_vertices() if
  // there is no such value of i.
  static int nextDistinctVertex(in S2Polyline pline, int index) {
    S2Point initial = pline.vertex(index);
    do {
      ++index;
    } while (index < pline.numVertices() && pline.vertex(index) == initial);
    return index;
  }

  // This struct represents a search state in the NearlyCovers algorithm
  // below.  See the description of the algorithm for details.
  struct SearchState {
    int i;
    int j;
    bool iInProgress;

    // This operator is needed for storing SearchStates in a set.  The ordering
    // chosen has no special meaning.
    int opCmp(SearchState b) const {
      if (i != b.i) return i - b.i;
      if (j != b.j) return j - b.j;
      return iInProgress - b.iInProgress;
    }
  }
}
