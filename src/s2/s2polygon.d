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

module s2.s2polygon;

import s2.builder.util.s2polygon_layer;
import s2.builder.util.s2polyline_layer;
import s2.builder.util.s2polyline_vector_layer;
import s2.builder.util.snap_functions;
import s2.logger;
import s2.mutable_s2shape_index;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s2boolean_operation;
import s2.s2builder;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2contains_point_query;
import s2.s2debug;
import s2.s2edge_crossings : INTERSECTION_ERROR, INTERSECTION_MERGE_RADIUS;
import s2.s2error;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2loop;
import s2.s2metrics;
import s2.s2point;
import s2.s2polyline;
import s2.s2region;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2shape_index_region : makeS2ShapeIndexRegion;
import s2.util.container.btree_map;

import core.atomic;
import std.algorithm;
import std.container.rbtree;
import std.conv;
import std.exception;
import std.math;
import std.range;
import std.typecons;

/**
 * Build the S2ShapeIndex only when it is first needed.  This can save
 * significant amounts of memory and time when geometry is constructed but
 * never queried, for example when converting from one format to another.
 */
enum bool S2POLYGON_LAZY_INDEXING = true;

/**
 * An S2Polygon is an S2Region object that represents a polygon.  A polygon is
 * defined by zero or more loops; recall that the interior of a loop is
 * defined to be its left-hand side (see S2Loop).  There are two different
 * conventions for creating an S2Polygon:
 *
 *   - InitNested() expects the input loops to be nested hierarchically.  The
 *     polygon interior then consists of the set of points contained by an odd
 *     number of loops.  So for example, a circular region with a hole in it
 *     would be defined as two CCW loops, with one loop containing the other.
 *     The loops can be provided in any order.
 *
 *     When the orientation of the input loops is unknown, the nesting
 *     requirement is typically met by calling S2Loop::Normalize() on each
 *     loop (which inverts the loop if necessary so that it encloses at most
 *     half the sphere).  But in fact any set of loops can be used as long as
 *     (1) there is no pair of loops that cross, and (2) there is no pair of
 *     loops whose union is the entire sphere.
 *
 *   - InitOriented() expects the input loops to be oriented such that the
 *     polygon interior is on the left-hand side of every loop.  So for
 *     example, a circular region with a hole in it would be defined using a
 *     CCW outer loop and a CW inner loop.  The loop orientations must all be
 *     consistent; for example, it is not valid to have one CCW loop nested
 *     inside another CCW loop, because the region between the two loops is on
 *     the left-hand side of one loop and the right-hand side of the other.
 *
 * Most clients will not call these methods directly; instead they should use
 * S2Builder, which has better support for dealing with imperfect data.
 *
 * When the polygon is initialized, the given loops are automatically
 * converted into a canonical form consisting of "shells" and "holes".  Shells
 * and holes are both oriented CCW, and are nested hierarchically.  The loops
 * are reordered to correspond to a preorder traversal of the nesting
 * hierarchy; InitOriented may also invert some loops. The set of input S2Loop
 * pointers is always preserved; the caller can use this to determine how the
 * loops were reordered if desired.
 *
 * Polygons may represent any region of the sphere with a polygonal boundary,
 * including the entire sphere (known as the "full" polygon).  The full
 * polygon consists of a single full loop (see S2Loop), whereas the empty
 * polygon has no loops at all.
 *
 * Polygons have the following restrictions:
 *
 *  - Loops may not cross, i.e. the boundary of a loop may not intersect
 *    both the interior and exterior of any other loop.
 *
 *  - Loops may not share edges, i.e. if a loop contains an edge AB, then
 *    no other loop may contain AB or BA.
 *
 *  - Loops may share vertices, however no vertex may appear twice in a
 *    single loop (see S2Loop).
 *
 *  - No loop may be empty.  The full loop may appear only in the full polygon.
 */
class S2Polygon : S2Region {
public:
  /**
   * The default constructor creates an empty polygon.  It can be made
   * non-empty by calling Init(), Decode(), etc.
   */
  this() {
    _bound = new S2LatLngRect();
    _subregionBound = new S2LatLngRect();
    _index = new MutableS2ShapeIndex();
    _s2debugOverride = S2Debug.ALLOW;
    _errorInconsistentLoopOrientations = false;
    _numVertices = 0;
    _unindexedContainsCalls = 0;
  }

  /**
   * Convenience constructor that calls InitNested() with the given loops.
   *
   * When called with override == S2Debug::ALLOW, the automatic validity
   * checking is controlled by --s2debug (which is true by default in
   * non-optimized builds).  When this flag is enabled, a fatal error is
   * generated whenever an invalid polygon is constructed.
   *
   * With override == S2Debug::DISABLE, the automatic validity checking
   * is disabled.  The main reason to do this is if you intend to call
   * IsValid() explicitly.  (See set_s2debug_override() for details.)
   * Example:
   *
   *   S2Polygon* polygon = new S2Polygon(loops, S2Debug::DISABLE);
   *
   * This is equivalent to:
   *
   *   S2Polygon* polygon = new S2Polygon;
   *   polygon->set_s2debug_override(S2Debug::DISABLE);
   *   polygon->InitNested(loops);
   */
  this(S2Loop[] loops, S2Debug s2debugOverride = S2Debug.ALLOW) {
    _bound = new S2LatLngRect();
    _subregionBound = new S2LatLngRect();
    _index = new MutableS2ShapeIndex();
    _s2debugOverride = s2debugOverride;
    initializeNested(loops);
  }

  /**
   * Convenience constructor that creates a polygon with a single loop
   * corresponding to the given cell.
   */
  this(in S2Cell cell) {
    _bound = new S2LatLngRect();
    _subregionBound = new S2LatLngRect();
    _index = new MutableS2ShapeIndex();
    _s2debugOverride = S2Debug.ALLOW;
    initialize(new S2Loop(cell));
  }

  /**
   * Convenience constructor that calls Init(S2Loop*).  Note that this method
   * automatically converts the special empty loop (see S2Loop) into an empty
   * polygon, unlike the vector-of-loops constructor which does not allow
   * empty loops at all.
   */
  this(S2Loop loop, S2Debug s2debugOverride = S2Debug.ALLOW) {
    _bound = new S2LatLngRect();
    _subregionBound = new S2LatLngRect();
    _index = new MutableS2ShapeIndex();
    _s2debugOverride = s2debugOverride;
    initialize(loop);
  }

  /**
   * Create a polygon from a set of hierarchically nested loops.  The polygon
   * interior consists of the points contained by an odd number of loops.
   * (Recall that a loop contains the set of points on its left-hand side.)
   *
   * This method figures out the loop nesting hierarchy and assigns every
   * loop a depth.  Shells have even depths, and holes have odd depths.  Note
   * that the loops are reordered so the hierarchy can be traversed more
   * easily (see GetParent(), GetLastDescendant(), and S2Loop::depth()).
   *
   * This method may be called more than once, in which case any existing
   * loops are deleted before being replaced by the input loops.
   */
  void initializeNested(S2Loop[] loops) {
    clearLoops();
    _loops = loops;

    if (numLoops() == 1) {
      initializeOneLoop();
      return;
    }
    LoopMap loop_map;
    foreach (int i; 0 .. numLoops()) {
      insertLoop(loop(i), null, loop_map);
    }
    _loops.length = 0;
    initializeLoops(loop_map);

    // Compute num_vertices_, bound_, subregion_bound_.
    initializeLoopProperties();
  }

  /**
   * Like InitNested(), but expects loops to be oriented such that the polygon
   * interior is on the left-hand side of all loops.  This implies that shells
   * and holes should have opposite orientations in the input to this method.
   * (During initialization, loops representing holes will automatically be
   * inverted.)
   */
  void initializeOriented(S2Loop[] loops) {
    // Here is the algorithm:
    //
    // 1. Remember which of the given loops contain S2::Origin().
    //
    // 2. Invert loops as necessary to ensure that they are nestable (i.e., no
    //    loop contains the complement of any other loop).  This may result in a
    //    set of loops corresponding to the complement of the given polygon, but
    //    we will fix that problem later.
    //
    //    We make the loops nestable by first normalizing all the loops (i.e.,
    //    inverting any loops whose turning angle is negative).  This handles
    //    all loops except those whose turning angle is very close to zero
    //    (within the maximum error tolerance).  Any such loops are inverted if
    //    and only if they contain S2::Origin().  (In theory this step is only
    //    necessary if there are at least two such loops.)  The resulting set of
    //    loops is guaranteed to be nestable.
    //
    // 3. Build the polygon.  This yields either the desired polygon or its
    //    complement.
    //
    // 4. If there is at least one loop, we find a loop L that is adjacent to
    //    S2::Origin() (where "adjacent" means that there exists a path
    //    connecting S2::Origin() to some vertex of L such that the path does
    //    not cross any loop).  There may be a single such adjacent loop, or
    //    there may be several (in which case they should all have the same
    //    contains_origin() value).  We choose L to be the loop containing the
    //    origin whose depth is greatest, or loop(0) (a top-level shell) if no
    //    such loop exists.
    //
    // 5. If (L originally contained origin) != (polygon contains origin), we
    //    invert the polygon.  This is done by inverting a top-level shell whose
    //    turning angle is minimal and then fixing the nesting hierarchy.  Note
    //    that because we normalized all the loops initially, this step is only
    //    necessary if the polygon requires at least one non-normalized loop to
    //    represent it.

    bool[S2Loop] contained_origin;
    for (int i = 0; i < loops.length; ++i) {
      S2Loop loop = loops[i];
      if (loop.containsOrigin()) {
        contained_origin[loop] = true;
      }
      double angle = loop.getTurningAngle();
      if (fabs(angle) > loop.getTurningAngleMaxError()) {
        // Normalize the loop.
        if (angle < 0) loop.invert();
      } else {
        // Ensure that the loop does not contain the origin.
        if (loop.containsOrigin()) loop.invert();
      }
    }
    initializeNested(loops);
    if (numLoops() > 0) {
      S2Loop origin_loop = loop(0);
      bool polygon_contains_origin = false;
      for (int i = 0; i < numLoops(); ++i) {
        if (loop(i).containsOrigin()) {
          polygon_contains_origin ^= true;
          origin_loop = loop(i);
        }
      }
      if (contained_origin.get(origin_loop, false) != polygon_contains_origin) {
        invert();
      }
    }
    // Verify that the original loops had consistent shell/hole orientations.
    // Each original loop L should have been inverted if and only if it now
    // represents a hole.
    for (int i = 0; i < _loops.length; ++i) {
      if ((contained_origin.get(loop(i), false) != loop(i).containsOrigin())
          != loop(i).isHole()) {
        // There is no point in saving the loop index, because the error is a
        // property of the entire set of loops.  In general there is no way to
        // determine which ones are incorrect.
        _errorInconsistentLoopOrientations = true;
        if (flagsS2Debug && _s2debugOverride == S2Debug.ALLOW) {
          // The FLAGS_s2debug validity checking usually happens in InitIndex(),
          // but this error is detected too late for that.
          enforce(isValid());  // Always fails.
        }
      }
    }
  }

  /**
   * Initialize a polygon from a single loop.  Note that this method
   * automatically converts the special empty loop (see S2Loop) into an empty
   * polygon, unlike the vector-of-loops InitNested() method which does not
   * allow empty loops at all.
   */
  void initialize(S2Loop loop) {
    // We don't allow empty loops in the other Init() methods because deleting
    // them changes the number of loops, which is awkward to handle.
    clearLoops();
    if (loop.isEmpty()) {
      initializeLoopProperties();
    } else {
      _loops ~= loop;
      initializeOneLoop();
    }
  }

  /**
   * Releases ownership of and returns the loops of this polygon, and resets
   * the polygon to be empty.
   */
  S2Loop[] release() {
    // Reset the polygon to be empty.
    S2Loop[] loops;
    loops.swap(_loops);
    clearLoops();
    _numVertices = 0;
    _bound = S2LatLngRect.empty();
    _subregionBound = S2LatLngRect.empty();
    return loops;
  }

  // Makes a deep copy of the given source polygon.  The destination polygon
  // will be cleared if necessary.
  void copy(in S2Polygon src) {
    clearLoops();
    for (int i = 0; i < src.numLoops(); ++i) {
      _loops ~= src.loop(i).clone();
    }
    _s2debugOverride = src._s2debugOverride;
    // Don't copy error_inconsistent_loop_orientations_, since this is not a
    // property of the polygon but only of the way the polygon was constructed.
    _numVertices = src.numVertices();
    atomicStore!(MemoryOrder.raw)(_unindexedContainsCalls, 0);
    _bound = new S2LatLngRect(src._bound);
    _subregionBound = new S2LatLngRect(src._subregionBound);
    initializeIndex();  // TODO(ericv): Copy the index efficiently.
  }

  /// Destroys the polygon and frees its loops.
  // ~this() {
  //   clearLoops();
  // }

  /**
   * Allows overriding the automatic validity checks controlled by
   * --s2debug (which is true by default in non-optimized builds).
   * When this flag is enabled, a fatal error is generated whenever
   * an invalid polygon is constructed.  The main reason to disable
   * this flag is if you intend to call IsValid() explicitly, like this:
   *
   *   S2Polygon polygon;
   *   polygon.set_s2debug_override(S2Debug::DISABLE);
   *   polygon.Init(...);
   *   if (!polygon.IsValid()) { ... }
   *
   * This setting is preserved across calls to Init() and Decode().
   */
  void setS2debugOverride(S2Debug s2debugOverride) {
    _s2debugOverride = s2debugOverride;
  }
  S2Debug s2debugOverride() const {
      return _s2debugOverride;
  }

  /**
   * Returns true if this is a valid polygon (including checking whether all
   * the loops are themselves valid).  Note that validity is checked
   * automatically during initialization when --s2debug is enabled (true by
   * default in debug binaries).
   */
  bool isValid() {
    S2Error error;
    if (findValidationError(error) && flagsS2Debug) {
      logger.logError(error);
      return false;
    }
    return true;
  }

  /**
   * Returns true if this is *not* a valid polygon and sets "error"
   * appropriately.  Otherwise returns false and leaves "error" unchanged.
   *
   * Note that in error messages, loops that represent holes have their edges
   * numbered in reverse order, starting from the last vertex of the loop.
   *
   * REQUIRES: error != nullptr
   */
  bool findValidationError(out S2Error error) {
    import s2.shapeutil.visit_crossing_edge_pairs : findSelfIntersection;

    for (int i = 0; i < numLoops(); ++i) {
      // Check for loop errors that don't require building an S2ShapeIndex.
      if (loop(i).findValidationErrorNoIndex(error)) {
        error.initialize(error.code(), "Loop %d: %s", i, error.text());
        return true;
      }
      // Check that no loop is empty, and that the full loop only appears in the
      // full polygon.
      if (loop(i).isEmpty()) {
        error.initialize(S2Error.Code.POLYGON_EMPTY_LOOP, "Loop %d: empty loops are not allowed", i);
        return true;
      }
      if (loop(i).isFull() && numLoops() > 1) {
        error.initialize(
            S2Error.Code.POLYGON_EXCESS_FULL_LOOP, "Loop %d: full loop appears in non-full polygon", i);
        return true;
      }
    }

    // Check for loop self-intersections and loop pairs that cross
    // (including duplicate edges and vertices).
    if (findSelfIntersection(_index, error)) return true;

    // Check whether InitOriented detected inconsistent loop orientations.
    if (_errorInconsistentLoopOrientations) {
      error.initialize(S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
          "Inconsistent loop orientations detected");
      return true;
    }

    // Finally, verify the loop nesting hierarchy.
    return findLoopNestingError(error);
  }

  /// Return true if this is the empty polygon (consisting of no loops).
  bool isEmpty() const {
    return _loops.empty();
  }

  /**
   * Return true if this is the full polygon (consisting of a single loop that
   * encompasses the entire sphere).
   */
  bool isFull() const {
    return numLoops() == 1 && loop(0).isFull();
  }

  /// Return the number of loops in this polygon.
  int numLoops() const {
    return cast(int)(_loops.length);
  }

  /// Total number of vertices in all loops.
  int numVertices() const {
    return _numVertices;
  }

  // Return the loop at the given index.  Note that during initialization, the
  // given loops are reordered according to a preorder traversal of the loop
  // nesting hierarchy.  This implies that every loop is immediately followed
  // by its descendants.  This hierarchy can be traversed using the methods
  // GetParent(), GetLastDescendant(), and S2Loop::depth().
  inout(S2Loop) loop(int k) inout {
    return _loops[k];
  }

  /// Return the index of the parent of loop k, or -1 if it has no parent.
  int getParent(int k) const {
    int depth = loop(k).depth();
    if (depth == 0) return -1;  // Optimization.
    while (--k >= 0 && loop(k).depth() >= depth) continue;
    return k;
  }

  /**
   * Return the index of the last loop that is contained within loop k.
   * Returns num_loops() - 1 if k < 0.  Note that loops are indexed according
   * to a preorder traversal of the nesting hierarchy, so the immediate
   * children of loop k can be found by iterating over loops
   * (k+1)..GetLastDescendant(k) and selecting those whose depth is equal to
   * (loop(k)->depth() + 1).
   */
  int getLastDescendant(int k) const {
    if (k < 0) return numLoops() - 1;
    int depth = loop(k).depth();
    while (++k < numLoops() && loop(k).depth() > depth) continue;
    return k - 1;
  }

  /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops.  The return value is between 0 and 4*Pi.
   */
  double getArea() const {
    double area = 0;
    for (int i = 0; i < numLoops(); ++i) {
      area += loop(i).sign() * loop(i).getArea();
    }
    return area;
  }

  /**
   * Return the true centroid of the polygon multiplied by the area of the
   * polygon (see s2centroids.h for details on centroids).  The result is not
   * unit length, so you may want to normalize it.  Also note that in general,
   * the centroid may not be contained by the polygon.
   *
   * We prescale by the polygon area for two reasons: (1) it is cheaper to
   * compute this way, and (2) it makes it easier to compute the centroid of
   * more complicated shapes (by splitting them into disjoint regions and
   * adding their centroids).
   */
  S2Point getCentroid() const {
    S2Point centroid;
    for (int i = 0; i < numLoops(); ++i) {
      centroid += loop(i).sign() * loop(i).getCentroid();
    }
    return centroid;
  }

  /**
   * If all of the polygon's vertices happen to be the centers of S2Cells at
   * some level, then return that level, otherwise return -1.  See also
   * InitToSnapped() and s2builderutil::S2CellIdSnapFunction.
   * Returns -1 if the polygon has no vertices.
   */
  int getSnapLevel() const {
    import s2.s2coords : XYZtoFaceSiTi;

    int snap_level = -1;
    foreach (child; _loops) {
      for (int j = 0; j < child.numVertices(); ++j) {
        int face;
        uint si, ti;
        int level = XYZtoFaceSiTi(child.vertex(j), face, si, ti);
        if (level < 0) return level;  // Vertex is not a cell center.
        if (level != snap_level) {
          if (snap_level < 0) {
            snap_level = level;  // First vertex.
          } else {
            return -1;  // Vertices at more than one cell level.
          }
        }
      }
    }
    return snap_level;
  }

  /**
   * Return the distance from the given point to the polygon interior.  If the
   * polygon is empty, return S1Angle::Infinity().  "x" should be unit length.
   */
  S1Angle getDistance(in S2Point x) {
    // Note that S2Polygon::Contains(S2Point) is slightly more efficient than
    // the generic version used by S2ClosestEdgeQuery.
    if (contains(x)) return S1Angle.zero();
    return getDistanceToBoundary(x);
  }

  /**
   * Return the distance from the given point to the polygon boundary.  If the
   * polygon is empty or full, return S1Angle::Infinity() (since the polygon
   * has no boundary).  "x" should be unit length.
   */
  S1Angle getDistanceToBoundary(in S2Point x) {
    import s2.s2closest_edge_query;

    auto options = new S2ClosestEdgeQuery.Options();
    options.setIncludeInteriors(false);
    auto t = new S2ClosestEdgeQuery.PointTarget(x);
    return new S2ClosestEdgeQuery(_index, options).getDistance(t).toS1Angle();
  }

  static struct OverlapFractions {
    double first;
    double second;
  }

  /**
   * Return the overlap fractions between two polygons, i.e. the ratios of the
   * area of intersection to the area of each polygon.
   */
  static OverlapFractions getOverlapFractions(S2Polygon a, S2Polygon b) {
    auto intersection = new S2Polygon();
    intersection.initializeToIntersection(a, b);
    double intersection_area = intersection.getArea();
    double a_area = a.getArea();
    double b_area = b.getArea();
    return OverlapFractions(
        intersection_area >= a_area ? 1.0 : intersection_area / a_area,
        intersection_area >= b_area ? 1.0 : intersection_area / b_area);
  }

  /**
   * If the given point is contained by the polygon, return it.  Otherwise
   * return the closest point on the polygon boundary.  If the polygon is
   * empty, return the input argument.  Note that the result may or may not be
   * contained by the polygon.  "x" should be unit length.
   */
  S2Point project(in S2Point x) {
    if (contains(x)) return x;
    return projectToBoundary(x);
  }

  /**
   * Return the closest point on the polygon boundary to the given point.  If
   * the polygon is empty or full, return the input argument (since the
   * polygon has no boundary).  "x" should be unit length.
   */
  S2Point projectToBoundary(in S2Point x) {
    import s2.s2closest_edge_query;

    auto options = new S2ClosestEdgeQuery.Options();
    options.setIncludeInteriors(false);
    auto q = new S2ClosestEdgeQuery(_index, options);
    auto target = new S2ClosestEdgeQuery.PointTarget(x);
    S2ClosestEdgeQuery.Result edge = q.findClosestEdge(target);
    return q.project(x, edge);
  }

  /**
   * Return true if this polygon contains the given other polygon, i.e.
   * if polygon A contains all points contained by polygon B.
   */
  bool contains(S2Polygon b) {
    // It's worth checking bounding rectangles, since they are precomputed.
    // Note that the first bound has been expanded to account for possible
    // numerical errors in the second bound.
    if (!_subregionBound.contains(b._bound)) {
      // It is possible that A contains B even though Bound(A) does not contain
      // Bound(B).  This can only happen when polygon B has at least two outer
      // shells and the union of the two bounds spans all longitudes.  For
      // example, suppose that B consists of two shells with a longitude gap
      // between them, while A consists of one shell that surrounds both shells
      // of B but goes the other way around the sphere (so that it does not
      // intersect the longitude gap).
      //
      // For convenience we just check whether B has at least two loops rather
      // than two outer shells.
      if (b.numLoops() == 1 || !_bound.lng().unite(b._bound.lng()).isFull()) {
        return false;
      }
    }

    // The following case is not handled by S2BooleanOperation because it only
    // determines whether the boundary of the result is empty (which does not
    // distinguish between the full and empty polygons).
    if (isEmpty() && b.isFull()) return false;

    return S2BooleanOperation.contains(_index, b._index);
  }

  /**
   * Returns true if this polgyon (A) approximately contains the given other
   * polygon (B). This is true if it is possible to move the vertices of B
   * no further than "tolerance" such that A contains the modified B.
   *
   * For example, the empty polygon will contain any polygon whose maximum
   * width is no more than "tolerance".
   */
  bool approxContains(S2Polygon b, S1Angle tolerance) {
    auto difference = new S2Polygon();
    difference.initializeToApproxDifference(b, this, tolerance);
    return difference.isEmpty();
  }

  /**
   * Return true if this polygon intersects the given other polygon, i.e.
   * if there is a point that is contained by both polygons.
   */
  bool intersects(S2Polygon b) {
    // It's worth checking bounding rectangles, since they are precomputed.
    if (!_bound.intersects(b._bound)) return false;

    // The following case is not handled by S2BooleanOperation because it only
    // determines whether the boundary of the result is empty (which does not
    // distinguish between the full and empty polygons).
    if (isFull() && b.isFull()) return true;

    return S2BooleanOperation.intersects(_index, b._index);
  }

  /**
   * Returns true if this polgyon (A) and the given polygon (B) are
   * approximately disjoint.  This is true if it is possible to ensure that A
   * and B do not intersect by moving their vertices no further than
   * "tolerance".  This implies that in borderline cases where A and B overlap
   * slightly, this method returns true (A and B are approximately disjoint).
   *
   * For example, any polygon is approximately disjoint from a polygon whose
   * maximum width is no more than "tolerance".
   */
  bool approxDisjoint(S2Polygon b, S1Angle tolerance) {
    auto intersection = new S2Polygon();
    intersection.initializeToApproxIntersection(b, this, tolerance);
    return intersection.isEmpty();
  }

  /**
   * Initialize this polygon to the intersection, union, difference (A - B),
   * or symmetric difference (XOR) of the given two polygons.
   *
   * "snap_function" allows you to specify a minimum spacing between output
   * vertices, and/or that the vertices should be snapped to a discrete set of
   * points (e.g. S2CellId centers or E7 lat/lng coordinates).  Any snap
   * function can be used, including the IdentitySnapFunction with a
   * snap_radius of zero (which preserves the input vertices exactly).
   *
   * The boundary of the output polygon before snapping is guaranteed to be
   * accurate to within S2::kIntersectionError of the exact result.
   * Snapping can move the boundary by an additional distance that depends on
   * the snap function.  Finally, any degenerate portions of the output
   * polygon are automatically removed (i.e., regions that do not contain any
   * points) since S2Polygon does not allow such regions.
   *
   * See S2Builder and s2builderutil for more details on snap functions.  For
   * example, you can snap to E7 coordinates by setting "snap_function" to
   * s2builderutil::IntLatLngSnapFunction(7).
   *
   * The default snap function is the IdentitySnapFunction with a snap radius
   * of S2::kIntersectionMergeRadius (equal to about 1.8e-15 radians
   * or 11 nanometers on the Earth's surface).  This means that vertices may
   * be positioned arbitrarily, but vertices that are extremely close together
   * can be merged together.  The reason for a non-zero default snap radius is
   * that it helps to eliminate narrow cracks and slivers when T-vertices are
   * present.  For example, adjacent S2Cells at different levels do not share
   * exactly the same boundary, so there can be a narrow crack between them.
   * If a polygon is intersected with those cells and the pieces are unioned
   * together, the result would have a narrow crack unless the snap radius is
   * set to a non-zero value.
   *
   * Note that if you want to encode the vertices in a lower-precision
   * representation (such as S2CellIds or E7), it is much better to use a
   * suitable SnapFunction rather than rounding the vertices yourself, because
   * this will create self-intersections unless you ensure that the vertices
   * and edges are sufficiently well-separated first.  In particular you need
   * to use a snap function whose min_edge_vertex_separation() is at least
   * twice the maximum distance that a vertex can move when rounded.
   */
  void initializeToIntersection(S2Polygon a, S2Polygon b) {
    initializeToApproxIntersection(a, b, INTERSECTION_MERGE_RADIUS);
  }

  void initializeToIntersection(
      S2Polygon a, S2Polygon b, S2Builder.SnapFunction snap_function) {
    if (!a._bound.intersects(b._bound)) return;
    initializeToOperation(S2BooleanOperation.OpType.INTERSECTION, snap_function, a, b);

    // If the boundary is empty then there are two possible results: the empty
    // polygon or the full polygon.  Note that the (approximate) intersection of
    // two non-full polygons may be full, because one or both polygons may have
    // tiny cracks or holes that are eliminated by snapping.  Similarly, the
    // (approximate) intersection of two polygons that contain a common point
    // may be empty, since the point might be contained by tiny loops that are
    // snapped away.
    //
    // So instead we fall back to heuristics.  First we compute the minimum and
    // maximum intersection area based on the areas of the two input polygons.
    // If only one of {0, 4*Pi} is possible then we return that result.  If
    // neither is possible (before snapping) then we return the one that is
    // closest to being possible.  (It never true that both are possible.)
    if (numLoops() == 0) {
      // We know that both polygons are non-empty due to the initial bounds
      // check.  By far the most common case is that the intersection is empty,
      // so we want to make that case fast.  The intersection area satisfies:
      //
      //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
      //
      // where A, B can refer to a polygon or its area.  Note that if either A
      // or B is at most 2*Pi, the result must be "empty".  We can use the
      // bounding rectangle areas as upper bounds on the polygon areas.
      if (a._bound.area() <= 2 * M_PI || b._bound.area() <= 2 * M_PI) return;
      double a_area = a.getArea();
      double b_area = b.getArea();
      double min_area = max(0.0, a_area + b_area - 4 * M_PI);
      double max_area = min(a_area, b_area);
      if (min_area > 4 * M_PI - max_area) {
        invert();
      }
    }
  }

  void initializeToUnion(S2Polygon a, S2Polygon b) {
      initializeToApproxUnion(a, b, INTERSECTION_MERGE_RADIUS);
  }

  void initializeToUnion(S2Polygon a, S2Polygon b, S2Builder.SnapFunction snap_function) {
    initializeToOperation(S2BooleanOperation.OpType.UNION, snap_function, a, b);
    if (numLoops() == 0) {
      // See comments in InitToApproxIntersection().  In this case, the union
      // area satisfies:
      //
      //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
      //
      // where A, B can refer to a polygon or its area.  The most common case is
      // that neither input polygon is empty, but the union is empty due to
      // snapping.
      if (a._bound.area() + b._bound.area() <= 2 * M_PI) return;
      double a_area = a.getArea();
      double b_area = b.getArea();
      double min_area = max(a_area, b_area);
      double max_area = min(4 * M_PI, a_area + b_area);
      if (min_area > 4 * M_PI - max_area) {
        invert();
      }
    }
  }

  void initializeToDifference(S2Polygon a, S2Polygon b) {
      initializeToApproxDifference(a, b, INTERSECTION_MERGE_RADIUS);
  }

  void initializeToDifference(
      S2Polygon a, S2Polygon b, S2Builder.SnapFunction snap_function) {
    initializeToOperation(S2BooleanOperation.OpType.DIFFERENCE, snap_function, a, b);
    if (numLoops() == 0) {
      // See comments in InitToApproxIntersection().  In this case, the
      // difference area satisfies:
      //
      //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
      //
      // where A, B can refer to a polygon or its area.  By far the most common
      // case is that result is empty.
      if (a._bound.area() <= 2 * M_PI || b._bound.area() >= 2 * M_PI) return;
      double a_area = a.getArea();
      double b_area = b.getArea();
      double min_area = max(0.0, a_area - b_area);
      double max_area = min(a_area, 4 * M_PI - b_area);
      if (min_area > 4 * M_PI - max_area) {
        invert();
      }
    }
  }

  void initializeToSymmetricDifference(S2Polygon a, S2Polygon b) {
    initializeToApproxSymmetricDifference(a, b, INTERSECTION_MERGE_RADIUS);
  }

  void initializeToSymmetricDifference(
      S2Polygon a, S2Polygon b, S2Builder.SnapFunction snap_function) {
    initializeToOperation(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snap_function, a, b);
    if (numLoops() == 0) {
      // See comments in InitToApproxIntersection().  In this case, the
      // difference area satisfies:
      //
      //   |A - B| <= SymmetricDifference(A, B) <= 4*Pi - |4*Pi - (A + B)|
      //
      // where A, B can refer to a polygon or its area.  By far the most common
      // case is that result is empty.
      if (a._bound.area() + b._bound.area() <= 2 * M_PI) return;
      double a_area = a.getArea();
      double b_area = b.getArea();
      double min_area = fabs(a_area - b_area);
      double max_area = 4 * M_PI - fabs(4 * M_PI - (a_area + b_area));
      // If both input polygons have area 2*Pi, the result could be either empty
      // or full.  We explicitly want to choose "empty" in this case since it is
      // much more likely that the user is computing the difference between two
      // nearly identical polygons.  Hence the bias below.
      enum double kBiasTowardsEmpty = 1e-14;
      if (min_area - kBiasTowardsEmpty > 4 * M_PI - max_area) {
        invert();
      }
    }
  }

  // Convenience functions that use the IdentitySnapFunction with the given
  // snap radius.  TODO(ericv): Consider deprecating these and require the
  // snap function to be specified explcitly?
  void initializeToApproxIntersection(S2Polygon a, S2Polygon b, S1Angle snap_radius) {
    initializeToIntersection(a, b, new IdentitySnapFunction(snap_radius));
  }
  void initializeToApproxUnion(S2Polygon a, S2Polygon b, S1Angle snap_radius) {
    initializeToUnion(a, b, new IdentitySnapFunction(snap_radius));
  }
  void initializeToApproxDifference(S2Polygon a, S2Polygon b, S1Angle snap_radius) {
    initializeToDifference(a, b, new IdentitySnapFunction(snap_radius));
  }
  void initializeToApproxSymmetricDifference(S2Polygon a, S2Polygon b, S1Angle snap_radius) {
    initializeToSymmetricDifference(a, b, new IdentitySnapFunction(snap_radius));
  }

  /**
   * Snaps the vertices of the given polygon using the given SnapFunction
   * (e.g., s2builderutil::IntLatLngSnapFunction(6) snaps to E6 coordinates).
   * This can change the polygon topology (merging loops, for example), but
   * the resulting polygon is guaranteed to be valid, and no vertex will move
   * by more than snap_function.snap_radius().  See S2Builder for other
   * guarantees (e.g., minimum edge-vertex separation).
   *
   * Note that this method is a thin wrapper over S2Builder, so if you are
   * starting with data that is not in S2Polygon format (e.g., integer E7
   * coordinates) then it is faster to just use S2Builder directly.
   */
  void initializeToSnapped(in S2Polygon a, in S2Builder.SnapFunction snap_function) {
    auto options = new S2Builder.Options(snap_function);
    options.setSimplifyEdgeChains(true);
    auto builder = new S2Builder(options);
    initializeFromBuilder(a, builder);
  }

  /**
   * Convenience function that snaps the vertices to S2CellId centers at the
   * given level (default level 30, which has S2CellId centers spaced about 1
   * centimeter apart).  Polygons can be efficiently encoded by Encode() after
   * they have been snapped.
   */
  void initializeToSnapped(in S2Polygon a, int snap_level = S2CellId.MAX_LEVEL) {
    auto builder = new S2Builder(new S2Builder.Options(new S2CellIdSnapFunction(snap_level)));
    initializeFromBuilder(a, builder);
  }

  /**
   * Snaps the input polygon according to the given "snap_function" and
   * reduces the number of vertices if possible, while ensuring that no vertex
   * moves further than snap_function.snap_radius().
   *
   * Simplification works by replacing nearly straight chains of short edges
   * with longer edges, in a way that preserves the topology of the input
   * polygon up to the creation of degeneracies.  This means that loops or
   * portions of loops may become degenerate, in which case they are removed.
   * For example, if there is a very small island in the original polygon, it
   * may disappear completely.  (Even if there are dense islands, they could
   * all be removed rather than being replaced by a larger simplified island
   * if more area is covered by water than land.)
   */
  void initializeToSimplified(in S2Polygon a, S2Builder.SnapFunction snap_function) {
    auto options = new S2Builder.Options(snap_function);
    options.setSimplifyEdgeChains(true);
    auto builder = new S2Builder(options);
    initializeFromBuilder(a, builder);
  }

  /**
   * Like InitToSimplified, except that any vertices or edges on the boundary
   * of the given S2Cell are preserved if possible.  This method requires that
   * the polygon has already been clipped so that it does not extend outside
   * the cell by more than "boundary_tolerance".  In other words, it operates
   * on polygons that have already been intersected with a cell.
   *
   * Typically this method is used in geometry-processing pipelines that
   * intersect polygons with a collection of S2Cells and then process those
   * cells in parallel, where each cell generates some geometry that needs to
   * be simplified.  In contrast, if you just need to simplify the *input*
   * geometry then it is easier and faster to do the simplification before
   * computing the intersection with any S2Cells.
   *
   * "boundary_tolerance" specifies how close a vertex must be to the cell
   * boundary to be kept.  The default tolerance is large enough to handle any
   * reasonable way of interpolating points along the cell boundary, such as
   * S2::GetIntersection(), S2::Interpolate(), or direct (u,v)
   * interpolation using S2::FaceUVtoXYZ().  However, if the vertices have
   * been snapped to a lower-precision representation (e.g., S2CellId centers
   * or E7 coordinates) then you will need to set this tolerance explicitly.
   * For example, if the vertices were snapped to E7 coordinates then
   * "boundary_tolerance" should be set to
   *
   *   s2builderutil::IntLatLngSnapFunction::MinSnapRadiusForExponent(7)
   *
   * Degenerate portions of loops are always removed, so if a vertex on the
   * cell boundary belongs only to degenerate regions then it will not be
   * kept.  For example, if the input polygon is a narrow strip of width less
   * than "snap_radius" along one side of the cell, then the entire loop may
   * become degenerate and be removed.
   *
   * REQUIRES: all vertices of "a" are within "boundary_tolerance" of "cell".
   */
  void initilizeToSimplifiedInCell(
      in S2Polygon a, in S2Cell cell, S1Angle snap_radius,
      S1Angle boundary_tolerance = S1Angle.fromRadians(1e-15)) {
    // The polygon to be simplified consists of "boundary edges" that follow the
    // cell boundary and "interior edges" that do not.  We want to simplify the
    // interior edges while leaving the boundary edges unchanged.  It's not
    // sufficient to call S2Builder::ForceVertex() on all boundary vertices.
    // For example, suppose the polygon includes a triangle ABC where all three
    // vertices are on the cell boundary and B is a cell corner.  Then if
    // interior edge AC snaps to vertex B, this loop would become degenerate and
    // be removed.  Similarly, we don't want boundary edges to snap to interior
    // vertices, since this also would cause portions of the polygon along the
    // boundary to be removed.
    //
    // Instead we use a two-pass algorithm.  In the first pass, we simplify
    // *only* the interior edges, using ForceVertex() to ensure that any edge
    // endpoints on the cell boundary do not move.  In the second pass, we add
    // the boundary edges (which are guaranteed to still form loops with the
    // interior edges) and build the output polygon.
    //
    // Note that in theory, simplifying the interior edges could create an
    // intersection with one of the boundary edges, since if two interior edges
    // intersect very near the boundary then the intersection point could be
    // slightly outside the cell (by at most S2::kIntersectionError).
    // This is the *only* way that a self-intersection can be created, and it is
    // expected to be extremely rare.  Nevertheless we use a small snap radius
    // in the second pass in order to eliminate any such self-intersections.
    //
    // We also want to preserve the cyclic vertex order of loops, so that the
    // original polygon can be reconstructed when no simplification is possible
    // (i.e., idempotency).  In order to do this, we group the input edges into
    // a sequence of polylines.  Each polyline contains only one type of edge
    // (interior or boundary).  We use S2Builder to simplify the interior
    // polylines, while the boundary polylines are passed through unchanged.
    // Each interior polyline is in its own S2Builder layer in order to keep the
    // edges in sequence.  This lets us ensure that in the second pass, the
    // edges are added in their original order so that S2PolygonLayer can
    // reconstruct the original loops.

    // We want an upper bound on how much "u" or "v" can change when a point on
    // the boundary of the S2Cell is moved away by up to "boundary_tolerance".
    // Inverting this, instead we could compute a lower bound on how far a point
    // can move away from an S2Cell edge when "u" or "v" is changed by a given
    // amount.  The latter quantity is simply (S2::kMinWidth.deriv() / 2)
    // under the S2_LINEAR_PROJECTION model, where we divide by 2 because we
    // want the bound in terms of (u = 2 * s - 1) rather than "s" itself.
    // Consulting s2metrics.cc, this value is sqrt(2/3)/2 = sqrt(1/6).
    // Going back to the original problem, this gives:
    double boundary_tolerance_uv = sqrt(6.0) * boundary_tolerance.radians();

    // The first pass yields a collection of simplified polylines that preserve
    // the original cyclic vertex order.
    auto polylines = simplifyEdgesInCell(a, cell, boundary_tolerance_uv, snap_radius);

    // The second pass eliminates any intersections between interior edges and
    // boundary edges, and then assembles the edges into a polygon.
    auto options = new S2Builder.Options(new IdentitySnapFunction(INTERSECTION_ERROR));
    options.setIdempotent(false);  // Force snapping up to the given radius
    auto builder = new S2Builder(options);
    builder.startLayer(new S2PolygonLayer(this));
    foreach (polyline; polylines) {
      builder.addPolyline(polyline);
    }
    S2Error error;
    if (!builder.build(error)) {
      logger.logFatal("Could not build polygon: ", error);
      return;
    }
    // If there are no loops, check whether the result should be the full
    // polygon rather than the empty one.  (See InitToApproxIntersection.)
    if (numLoops() == 0) {
      if (a._bound.area() > 2 * PI && a.getArea() > 2 * PI) invert();
    }
  }

  /// Initialize this polygon to the complement of the given polygon.
  void initializeToComplement(in S2Polygon a) {
    copy(a);
    invert();
  }

  /// Invert the polygon (replace it by its complement).
  void invert() {
    // Inverting any one loop will invert the polygon.  The best loop to invert
    // is the one whose area is largest, since this yields the smallest area
    // after inversion.  The loop with the largest area is always at depth 0.
    // The descendents of this loop all have their depth reduced by 1, while the
    // former siblings of this loop all have their depth increased by 1.

    // The empty and full polygons are handled specially.
    if (isEmpty()) {
      _loops ~= new S2Loop(S2Loop.full());
    } else if (isFull()) {
      clearLoops();
    } else {
      // Find the loop whose area is largest (i.e., whose turning angle is
      // smallest), minimizing calls to GetTurningAngle().  In particular, for
      // polygons with a single shell at level 0 there is not need to call
      // GetTurningAngle() at all.  (This method is relatively expensive.)
      int best = 0;
      const double kNone = 10.0;  // Flag that means "not computed yet"
      double best_angle = kNone;
      for (int i = 1; i < numLoops(); ++i) {
        if (loop(i).depth() == 0) {
          // We defer computing the turning angle of loop 0 until we discover
          // that the polygon has another top-level shell.
          if (best_angle == kNone) best_angle = loop(best).getTurningAngle();
          double angle = loop(i).getTurningAngle();
          // We break ties deterministically in order to avoid having the output
          // depend on the input order of the loops.
          if (angle < best_angle ||
              (angle == best_angle && compareLoops(loop(i), loop(best)) < 0)) {
            best = i;
            best_angle = angle;
          }
        }
      }
      // Build the new loops vector, starting with the inverted loop.
      loop(best).invert();
      S2Loop[] new_loops;
      new_loops.reserve(numLoops());
      // Add the former siblings of this loop as descendants.
      int last_best = getLastDescendant(best);
      new_loops ~= _loops[best]; // TODO: Remove from _loops?
      for (int i = 0; i < numLoops(); ++i) {
        if (i < best || i > last_best) {
          loop(i).setDepth(loop(i).depth() + 1);
          new_loops ~= _loops[i]; // TODO: Remove from _loops?
        }
      }
      // Add the former children of this loop as siblings.
      for (int i = 0; i < numLoops(); ++i) {
        if (i > best && i <= last_best) {
          loop(i).setDepth(loop(i).depth() - 1);
          new_loops ~= _loops[i]; // TODO: Remove from _loops?
        }
      }
      _loops.swap(new_loops);
      enforce(new_loops.length == numLoops());
    }
    clearIndex();
    initializeLoopProperties();
  }

  /**
   * Return true if this polygon contains the given polyline.  This method
   * returns an exact result, according to the following model:
   *
   *  - All edges are geodesics (of course).
   *
   *  - Vertices are ignored for the purposes of defining containment.
   *    (This is because polygons often do not contain their vertices, in
   *    order to that when a set of polygons tiles the sphere then every point
   *    is contained by exactly one polygon.)
   *
   *  - Points that lie exactly on geodesic edges are resolved using symbolic
   *    perturbations (i.e., they are considered to be infinitesmally offset
   *    from the edge).
   *
   *  - If the polygon and polyline share an edge, it is handled as follows.
   *    First, the polygon edges are oriented so that the interior is always
   *    on the left.  Then the shared polyline edge is contained if and only
   *    if it is in the same direction as the corresponding polygon edge.
   *    (This model ensures that when a polyline is intersected with a polygon
   *    and its complement, the edge only appears in one of the two results.)
   *
   * TODO(ericv): Update the implementation to correspond to the model above.
   */
  bool contains(in S2Polyline b) {
    return approxContains(b, INTERSECTION_MERGE_RADIUS);
  }

  /**
   * Returns true if this polgyon approximately contains the given polyline
   * This is true if it is possible to move the polyline vertices no further
   * than "tolerance" such that the polyline is now contained.
   */
  bool approxContains(in S2Polyline b, S1Angle tolerance) {
    auto difference = approxSubtractFromPolyline(b, tolerance);
    return difference.empty();
  }

  /**
   * Return true if this polygon intersects the given polyline.  This method
   * returns an exact result; see Contains(S2Polyline) for details.
   */
  bool intersects(in S2Polyline b) {
    return !approxDisjoint(b, INTERSECTION_MERGE_RADIUS);
  }

  /**
   * Returns true if this polgyon is approximately disjoint from the given
   * polyline.  This is true if it is possible to avoid intersection by moving
   * their vertices no further than "tolerance".
   *
   * This implies that in borderline cases where there is a small overlap,
   * this method returns true (i.e., they are approximately disjoint).
   */
  bool approxDisjoint(const S2Polyline b, S1Angle tolerance) {
    auto intersection = approxIntersectWithPolyline(b, tolerance);
    return intersection.empty();
  }

  /**
   * Intersect this polygon with the polyline "in" and return the resulting
   * zero or more polylines.  The polylines are returned in the order they
   * would be encountered by traversing "in" from beginning to end.
   * Note that the output may include polylines with only one vertex,
   * but there will not be any zero-vertex polylines.
   *
   * This is equivalent to calling ApproxIntersectWithPolyline() with the
   * "snap_radius" set to S2::kIntersectionMergeRadius.
   */
  S2Polyline[] intersectWithPolyline(in S2Polyline a) {
      return approxIntersectWithPolyline(a, INTERSECTION_MERGE_RADIUS);
  }

  /**
   * Similar to IntersectWithPolyline(), except that vertices will be
   * dropped as necessary to ensure that all adjacent vertices in the
   * sequence obtained by concatenating the output polylines will be
   * farther than "snap_radius" apart.  Note that this can change
   * the number of output polylines and/or yield single-vertex polylines.
   */
  S2Polyline[] approxIntersectWithPolyline(in S2Polyline a, S1Angle snap_radius) {
    return intersectWithPolyline(a, new IdentitySnapFunction(snap_radius));
  }

  // TODO(ericv): Update documentation.
  S2Polyline[] intersectWithPolyline(
      in S2Polyline a, S2Builder.SnapFunction snap_function) {
    return operationWithPolyline(S2BooleanOperation.OpType.INTERSECTION, snap_function, a);
  }

  /**
   * Same as IntersectWithPolyline, but subtracts this polygon from
   * the given polyline.
   */
  S2Polyline[] subtractFromPolyline(in S2Polyline a) {
    return approxSubtractFromPolyline(a, INTERSECTION_MERGE_RADIUS);
  }

  /**
   * Same as ApproxIntersectWithPolyline, but subtracts this polygon
   * from the given polyline.
   */
  S2Polyline[] approxSubtractFromPolyline(in S2Polyline a, S1Angle snap_radius) {
    return subtractFromPolyline(a, new IdentitySnapFunction(snap_radius));
  }

  S2Polyline[] subtractFromPolyline(in S2Polyline a, S2Builder.SnapFunction snap_function) {
    return operationWithPolyline(S2BooleanOperation.OpType.DIFFERENCE, snap_function, a);
  }

  /// Return a polygon which is the union of the given polygons.
  static S2Polygon destructiveUnion(S2Polygon[] polygons) {
    return destructiveApproxUnion(polygons, INTERSECTION_MERGE_RADIUS);
  }

  static S2Polygon destructiveApproxUnion(S2Polygon[] polygons, S1Angle snap_radius) {
    // Effectively create a priority queue of polygons in order of number of
    // vertices.  Repeatedly union the two smallest polygons and add the result
    // to the queue until we have a single polygon to return.
    alias QueueType = BTreeMap!(int, S2Polygon);
    auto queue = new QueueType();  // Map from # of vertices to polygon.
    foreach (polygon; polygons)
      queue.insert(polygon.numVertices(), polygon);

    while (queue.length > 1) {
      // Pop two simplest polygons from queue.
      QueueType.Iterator smallest_it = queue.begin();
      int a_size = smallest_it.getValue().key;
      auto a_polygon = smallest_it.getValue().value;
      queue.remove(smallest_it.getValue().key);
      smallest_it = queue.begin();
      int b_size = smallest_it.getValue().key;
      auto b_polygon = smallest_it.getValue().value;
      queue.remove(smallest_it.getValue().key);

      // Union and add result back to queue.
      auto union_polygon = new S2Polygon();
      union_polygon.initializeToApproxUnion(a_polygon, b_polygon, snap_radius);
      queue.insert(a_size + b_size, union_polygon);
      // We assume that the number of vertices in the union polygon is the
      // sum of the number of vertices in the original polygons, which is not
      // always true, but will almost always be a decent approximation, and
      // faster than recomputing.
    }

    if (queue.empty())
      return new S2Polygon();
    else
      return queue.begin().getValue().value;
  }

  /**
   * Initialize this polygon to the outline of the given cell union.
   * In principle this polygon should exactly contain the cell union and
   * this polygon's inverse should not intersect the cell union, but rounding
   * issues may cause this not to be the case.
   */
  void initializeToCellUnionBorder(in S2CellUnion cells) {
    // We use S2Builder to compute the union.  Due to rounding errors, we can't
    // compute an exact union - when a small cell is adjacent to a larger cell,
    // the shared edges can fail to line up exactly.  Two cell edges cannot come
    // closer then kMinWidth, so if we have S2Builder snap edges within half
    // that distance, then we should always merge shared edges without merging
    // different edges.
    double snap_radius = 0.5 * MIN_WIDTH.getValue(S2CellId.MAX_LEVEL);
    auto builder = new S2Builder(new S2Builder.Options(
            new IdentitySnapFunction(S1Angle.fromRadians(snap_radius))));
    builder.startLayer(new S2PolygonLayer(this));
    foreach (S2CellId id; cells.cellIds) {
      builder.addLoop(new S2Loop(new S2Cell(id)));
    }
    S2Error error;
    if (!builder.build(error)) {
      logger.logFatal("InitToCellUnionBorder failed: ", error);
    }

    // If there are no loops, check whether the result should be the full
    // polygon rather than the empty one.  There are only two ways that this can
    // happen: either the cell union is empty, or it consists of all six faces.
    if (numLoops() == 0) {
      if (cells.empty()) return;
      enforce(6uL << (2 * S2CellId.MAX_LEVEL) == cells.leafCellsCovered());
      invert();
    }
  }

  /**
   * Return true if every loop of this polygon shares at most one vertex with
   * its parent loop.  Every polygon has a unique normalized form.  A polygon
   * can be normalized by passing it through S2Builder (with no snapping) in
   * order to reconstruct the polygon from its edges.
   *
   * Generally there is no reason to convert polygons to normalized form.  It
   * is mainly useful for testing in order to compare whether two polygons
   * have exactly the same interior, even when they have a different loop
   * structure.  For example, a diamond nested within a square (touching at
   * four points) could be represented as a square with a diamond-shaped hole,
   * or as four triangles.  Methods such as BoundaryApproxEquals() will report
   * these polygons as being different (because they have different
   * boundaries) even though they contain the same points.  However if they
   * are both converted to normalized form (the "four triangles" version) then
   * they can be compared more easily.
   *
   * Also see ApproxEquals(), which can determine whether two polygons contain
   * approximately the same set of points without any need for normalization.
   */
  bool isNormalized() const {
    import std.algorithm : count;

    // TODO(ericv): The condition tested here is insufficient.  The correct
    // condition is that each *connected component* of child loops can share at
    // most one vertex with their parent loop.  Example: suppose loop A has
    // children B, C, D, and the following pairs are connected: AB, BC, CD, DA.
    // Then the polygon is not normalized.
    auto vertices = new RedBlackTree!S2Point();
    Rebindable!(const(S2Loop)) last_parent = null;
    for (int i = 0; i < numLoops(); ++i) {
      const(S2Loop) child = loop(i);
      if (child.depth() == 0) continue;
      Rebindable!(const(S2Loop)) parent = loop(getParent(i));
      if (parent != last_parent) {
        vertices.clear();
        for (int j = 0; j < parent.numVertices(); ++j) {
          vertices.insert(parent.vertex(j));
        }
        last_parent = parent;
      }
      int cnt = 0;
      for (int j = 0; j < child.numVertices(); ++j) {
        if (count(vertices.equalRange(child.vertex(j))) > 0) ++cnt;
      }
      if (cnt > 1) return false;
    }
    return true;
  }

  /**
   * Return true if two polygons have exactly the same loops.  The loops must
   * appear in the same order, and corresponding loops must have the same
   * linear vertex ordering (i.e., cyclic rotations are not allowed).
   */
  override
  bool opEquals(in Object o) const {
    S2Polygon b = cast(S2Polygon) o;
    if (b is null) return false;
    if (numLoops() != b.numLoops()) return false;
    for (int i = 0; i < numLoops(); ++i) {
      const S2Loop a_loop = loop(i);
      const S2Loop b_loop = b.loop(i);
      if ((b_loop.depth() != a_loop.depth()) || !b_loop.equals(a_loop)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if two polygons are approximately equal to within the given
   * tolerance.  This is true if it is possible to move the vertices of the
   * two polygons so that they contain the same set of points.
   *
   * Note that according to this model, small regions less than "tolerance" in
   * width do not need to be considered, since these regions can be collapsed
   * into degenerate loops (which contain no points) by moving their vertices.
   *
   * This model is not as strict as using the Hausdorff distance would be, and
   * it is also not as strict as BoundaryNear (defined below).  However, it is
   * a good choice for comparing polygons that have been snapped, simplified,
   * unioned, etc, since these operations use a model similar to this one
   * (i.e., degenerate loops or portions of loops are automatically removed).
   */
  bool approxEquals(S2Polygon b, S1Angle tolerance) {
    // TODO(ericv): This can be implemented more cheaply with S2Builder, by
    // simply adding all the edges from one polygon, adding the reversed edge
    // from the other polygon, and turning on the options to split edges and
    // discard sibling pairs.  Then the polygons are approximately equal if the
    // output graph has no edges.
    auto symmetric_difference = new S2Polygon();
    symmetric_difference.initializeToApproxSymmetricDifference(b, this, tolerance);
    return symmetric_difference.isEmpty();
  }

  /**
   * Returns true if two polygons have the same boundary.  More precisely,
   * this method requires that both polygons have loops with the same cyclic
   * vertex order and the same nesting hierarchy.  (This implies that vertices
   * may be cyclically rotated between corresponding loops, and the loop
   * ordering may be different between the two polygons as long as the nesting
   * hierarchy is the same.)
   */
  bool boundaryEquals(in S2Polygon b) const {
    if (numLoops() != b.numLoops()) return false;

    for (int i = 0; i < numLoops(); ++i) {
      const S2Loop a_loop = loop(i);
      bool success = false;
      for (int j = 0; j < numLoops(); ++j) {
        const S2Loop b_loop = b.loop(j);
        if (b_loop.depth() == a_loop.depth() && b_loop.boundaryEquals(a_loop)) {
          success = true;
          break;
        }
      }
      if (!success) return false;
    }
    return true;
  }

  /**
   * Return true if two polygons have the same boundary except for vertex
   * perturbations.  Both polygons must have loops with the same cyclic vertex
   * order and the same nesting hierarchy, but the vertex locations are
   * allowed to differ by up to "max_error".
   */
  bool boundaryApproxEquals(in S2Polygon b, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    if (numLoops() != b.numLoops()) return false;

    // For now, we assume that there is at most one candidate match for each
    // loop.  (So far this method is just used for testing.)

    for (int i = 0; i < numLoops(); ++i) {
      const S2Loop a_loop = loop(i);
      bool success = false;
      for (int j = 0; j < numLoops(); ++j) {
        const S2Loop b_loop = b.loop(j);
        if (b_loop.depth() == a_loop.depth() && b_loop.boundaryApproxEquals(a_loop, max_error)) {
          success = true;
          break;
        }
      }
      if (!success) return false;
    }
    return true;
  }

  /**
   * Return true if two polygons have boundaries that are within "max_error"
   * of each other along their entire lengths.  More precisely, there must be
   * a bijection between the two sets of loops such that for each pair of
   * loops, "a_loop->BoundaryNear(b_loop)" is true.
   */
  bool boundaryNear(in S2Polygon b, S1Angle max_error = S1Angle.fromRadians(1e-15)) const {
    if (numLoops() != b.numLoops()) return false;

    // For now, we assume that there is at most one candidate match for each
    // loop.  (So far this method is just used for testing.)

    for (int i = 0; i < numLoops(); ++i) {
      const S2Loop a_loop = loop(i);
      bool success = false;
      for (int j = 0; j < numLoops(); ++j) {
        const S2Loop b_loop = b.loop(j);
        if (b_loop.depth() == a_loop.depth() && b_loop.boundaryNear(a_loop, max_error)) {
          success = true;
          break;
        }
      }
      if (!success) return false;
    }
    return true;
  }

  // TODO: Fix this calculation, it is currently inaccurate.
  /// Returns the total number of bytes used by the polygon.
  size_t spaceUsed() const {
    size_t size = this.sizeof;
    for (int i = 0; i < numLoops(); ++i) {
      size += loop(i).spaceUsed();
    }
    size += _index.spaceUsed() - _index.sizeof;
    return size;
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  /**
   * GetRectBound() returns essentially tight results, while GetCapBound()
   * might have a lot of extra padding.  Both bounds are conservative in that
   * if the loop contains a point P, then the bound contains P also.
   */
  override
  S2Polygon clone() const {
    S2Polygon result = new S2Polygon();
    result.copy(this);
    return result;
  }

  /// Cap surrounding rect bound.
  override
  S2Cap getCapBound() const {
    return _bound.getCapBound();
  }

  override
  S2LatLngRect getRectBound() const {
    return _bound.clone();
  }

  override
  void getCellUnionBound(out S2CellId[] cell_ids) {
    return makeS2ShapeIndexRegion(_index).getCellUnionBound(cell_ids);
  }

  override
  bool contains(in S2Cell target) {
    return makeS2ShapeIndexRegion(_index).contains(target);
  }

  override
  bool mayIntersect(in S2Cell target) {
    return makeS2ShapeIndexRegion(_index).mayIntersect(target);
  }

  /// The point 'p' does not need to be normalized.
  override
  bool contains(in S2Point p) {
    // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
    // worthwhile only when it might allow us to delay building the index.
    if (!_index.isFresh() && !_bound.contains(p)) return false;

    // For small polygons it is faster to just check all the crossings.
    // Otherwise we keep track of the number of calls to Contains() and only
    // build the index once enough calls have been made so that we think it is
    // worth the effort.  See S2Loop::Contains(S2Point) for detailed comments.
    enum int kMaxBruteForceVertices = 32;
    enum int kMaxUnindexedContainsCalls = 20;
    if (numVertices() <= kMaxBruteForceVertices
        || !_index.isFresh() && atomicOp!"+="(_unindexedContainsCalls, 1) != kMaxUnindexedContainsCalls) {
      bool inside = false;
      for (int i = 0; i < numLoops(); ++i) {
        // Use brute force to avoid building the loop's S2ShapeIndex.
        inside ^= loop(i).bruteForceContains(p);
      }
      return inside;
    }
    // Otherwise we look up the S2ShapeIndex cell containing this point.
    return makeS2ContainsPointQuery(_index).contains(p);
  }

  // TODO: Implement when encode/decode support is added.
  // Appends a serialized representation of the S2Polygon to "encoder".
  //
  // The encoding uses about 4 bytes per vertex for typical polygons in
  // Google's geographic repository, assuming that most vertices have been
  // snapped to the centers of S2Cells at some fixed level (typically using
  // InitToSnapped). The remaining vertices are stored using 24 bytes.
  // Decoding a polygon encoded this way always returns the original polygon,
  // without any loss of precision.
  //
  // The snap level is chosen to be the one that has the most vertices snapped
  // to S2Cells at that level.  If most vertices need 24 bytes, then all
  // vertices are encoded this way (this method automatically chooses the
  // encoding that has the best chance of giving the smaller output size).
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  // void Encode(Encoder* const encoder) const;

  // Decodes a polygon encoded with Encode().  Returns true on success.
  // bool Decode(Decoder* const decoder);

  // Decodes a polygon by pointing the S2Loop vertices directly into the
  // decoder's memory buffer (which needs to persist for the lifetime of the
  // decoded S2Polygon).  It is much faster than Decode(), but requires that
  // all the polygon vertices were encoded exactly using 24 bytes per vertex.
  // This essentially requires that the polygon was not snapped beforehand to
  // a given S2Cell level; otherwise this method falls back to Decode().
  //
  // Returns true on success.
  // bool DecodeWithinScope(Decoder* const decoder);

  /**
   * Wrapper class for indexing a polygon (see S2ShapeIndex).  Once this
   * object is inserted into an S2ShapeIndex it is owned by that index, and
   * will be automatically deleted when no longer needed by the index.  Note
   * that this class does not take ownership of the polygon itself (see
   * OwningShape below).  You can also subtype this class to store additional
   * data (see S2Shape for details).
   *
   * Note that unlike S2Polygon, the edges of S2Polygon::Shape are directed
   * such that the polygon interior is always on the left.
   */
  static class Shape : S2Shape {
  public:
    this() {
      _polygon = null;
      _cumulativeEdges = null;
    }

    /**
     * Initialization.  Does not take ownership of "polygon".  May be called
     * more than once.
     * TODO(ericv/jrosenstock): Make "polygon" a const reference.
     */
    this(in S2Polygon polygon) {
      initialize(polygon);
    }

    void initialize(in S2Polygon polygon) {
      _polygon = polygon;
      _cumulativeEdges = null;
      _numEdges = 0;
      if (!polygon.isFull()) {
        const int kMaxLinearSearchLoops = 12;  // From benchmarks.
        int num_loops = polygon.numLoops();
        if (num_loops > kMaxLinearSearchLoops) {
          _cumulativeEdges = new int[num_loops];
        }
        for (int i = 0; i < num_loops; ++i) {
          if (_cumulativeEdges) _cumulativeEdges[i] = _numEdges;
          _numEdges += polygon.loop(i).numVertices();
        }
      }
    }

    const(S2Polygon) polygon() const {
      return _polygon;
    }

    /// S2Shape interface:
    final override
    int numEdges() const {
      return _numEdges;
    }

    final override
    Edge edge(int e) const
    in {
      assert(e < numEdges());
    } body {
      import std.algorithm : findSplitBefore;

      const S2Polygon p = polygon();
      int i;
      if (_cumulativeEdges) {
        // "upper_bound" finds the loop just beyond the one we want.
        auto r = findSplitAfter!"a > b"(_cumulativeEdges[0 .. p.numLoops()], [e])[0];
        int start = r.back();
        i = cast(int) r.length - 1;
        e -= start;
      } else {
        // When the number of loops is small, linear search is faster.  Most often
        // there is exactly one loop and the code below executes zero times.
        for (i = 0; e >= p.loop(i).numVertices(); ++i) {
          e -= p.loop(i).numVertices();
        }
      }
      return Edge(p.loop(i).orientedVertex(e), p.loop(i).orientedVertex(e + 1));
    }

    final override
    int dimension() const {
      return 2;
    }

    final override
    ReferencePoint getReferencePoint() const {
      import s2.s2pointutil : origin;
      const S2Polygon p = polygon();
      bool contains_origin = false;
      for (int i = 0; i < p.numLoops(); ++i) {
        contains_origin ^= p.loop(i).containsOrigin();
      }
      return ReferencePoint(origin(), contains_origin);
    }

    final override
    int numChains() const {
      return _polygon.isFull() ? 0 : _polygon.numLoops();
    }

    final override
    Chain chain(int i) const
    in {
      assert(i < numChains());
    } body {
      if (_cumulativeEdges) {
        return Chain(_cumulativeEdges[i], _polygon.loop(i).numVertices());
      } else {
        int e = 0;
        for (int j = 0; j < i; ++j) e += _polygon.loop(j).numVertices();
        return Chain(e, _polygon.loop(i).numVertices());
      }
    }

    final override
    Edge chainEdge(int i, int j) const
    in {
      assert(i < numChains());
      assert(j < _polygon.loop(i).numVertices());
    } body {
      return Edge(polygon().loop(i).orientedVertex(j), polygon().loop(i).orientedVertex(j + 1));
    }

    final override
    ChainPosition chainPosition(int e) const
    in {
      // TODO(ericv): Make inline to remove code duplication with GetEdge.
      assert(e < numEdges());
    } body {
      const S2Polygon p = polygon();
      int i;
      if (_cumulativeEdges) {
        // "upper_bound" finds the loop just beyond the one we want.
        auto r = findSplitAfter!"a < b"(_cumulativeEdges[0 .. p.numLoops()], [e])[0];
        int start = r.back();
        i = cast(int) r.length - 1;
        e -= start;
      } else {
        // When the number of loops is small, linear search is faster.  Most often
        // there is exactly one loop and the code below executes zero times.
        for (i = 0; e >= p.loop(i).numVertices(); ++i) {
          e -= p.loop(i).numVertices();
        }
      }
      return ChainPosition(i, e);
    }

  private:
    /**
     * The total number of edges in the polygon.  This is the same as
     * polygon_->num_vertices() except in one case (polygon_->is_full()).  On
     * the other hand this field doesn't take up any extra space due to field
     * packing with S2Shape::id_.
     *
     * TODO(ericv): Consider using this field instead as an atomic<int> hint to
     * speed up edge location when there are a large number of loops.  Also
     * consider changing S2Polygon::num_vertices to num_edges instead.
     */
    int _numEdges;

    Rebindable!(const(S2Polygon)) _polygon;

    // An array where element "i" is the total number of edges in loops 0..i-1.
    // This field is only used for polygons that have a large number of loops.
    int[] _cumulativeEdges;
  }

  // Like Shape, except that the S2Polygon is automatically deleted when this
  // object is deleted by the S2ShapeIndex.  This is useful when an S2Polygon
  // is constructed solely for the purpose of indexing it.
  // TODO: Not needed in garbage-collected languages.
  // class OwningShape : public Shape {
  //  public:
  //   OwningShape() {}  // Must call Init().
  //   explicit OwningShape(std::unique_ptr<const S2Polygon> polygon)
  //       : Shape(polygon.release()) {
  //   }
  //   void Init(std::unique_ptr<const S2Polygon> polygon) {
  //     Shape::Init(polygon.release());
  //   }
  //   ~OwningShape() override { delete polygon(); }
  // };

  /**
   * Returns the built-in S2ShapeIndex associated with every S2Polygon.  This
   * can be used in conjunction with the various S2ShapeIndex query classes
   * (S2ClosestEdgeQuery, S2BooleanOperation, etc) to do things beyond what is
   * possible with S2Polygon built-in convenience methods.
   *
   * For example, to measure the distance from one S2Polygon to another, you
   * can write:
   *   S2ClosestEdgeQuery query(&polygon1.index());
   *   S2ClosestEdgeQuery::ShapeIndexTarget target(&polygon2.index());
   *   S1ChordAngle distance = query.GetDistance(&target);
   *
   * The index contains a single S2Polygon::Shape object.
   */
  const(MutableS2ShapeIndex) index() const {
    return _index;
  }

private:

  /// Given that loops_ contains a single loop, initialize all other fields.
  void initializeOneLoop()
  in {
    assert(numLoops() == 1);
  } body {
    S2Loop loop = _loops[0];
    loop.setDepth(0);
    _errorInconsistentLoopOrientations = false;
    _numVertices = loop.numVertices();
    _bound = loop.getRectBound();
    _subregionBound = S2LatLngRectBounder.expandForSubregions(_bound);
    initializeIndex();
  }

  /// Compute num_vertices_, bound_, subregion_bound_.
  void initializeLoopProperties() {
    _numVertices = 0;
    _bound = S2LatLngRect.empty();
    for (int i = 0; i < numLoops(); ++i) {
      if (loop(i).depth() == 0) {
        _bound = _bound.unite(loop(i).getRectBound());
      }
      _numVertices += loop(i).numVertices();
    }
    _subregionBound = S2LatLngRectBounder.expandForSubregions(_bound);
    initializeIndex();
  }

  /// Deletes the contents of the loops_ vector and resets the polygon state.
  void clearLoops() {
    clearIndex();
    _loops.length = 0;
    _errorInconsistentLoopOrientations = false;
  }

  /// Return true if there is an error in the loop nesting hierarchy.
  bool findLoopNestingError(ref S2Error error) {
    // First check that the loop depths make sense.
    for (int last_depth = -1, i = 0; i < numLoops(); ++i) {
      int depth = loop(i).depth();
      if (depth < 0 || depth > last_depth + 1) {
        error.initialize(S2Error.Code.POLYGON_INVALID_LOOP_DEPTH,
            "Loop %d: invalid loop depth (%d)", i, depth);
        return true;
      }
      last_depth = depth;
    }
    // Then check that they correspond to the actual loop nesting.  This test
    // is quadratic in the number of loops but the cost per iteration is small.
    for (int i = 0; i < numLoops(); ++i) {
      int last = getLastDescendant(i);
      for (int j = 0; j < numLoops(); ++j) {
        if (i == j) continue;
        bool nested = (j >= i + 1) && (j <= last);
        const bool reverse_b = false;
        if (loop(i).containsNonCrossingBoundary(loop(j), reverse_b) != nested) {
          error.initialize(S2Error.Code.POLYGON_INVALID_LOOP_NESTING,
              "Invalid nesting: loop %d should %scontain loop %d",
              i, nested ? "" : "not ", j);
          return true;
        }
      }
    }
    return false;
  }

  /**
   * A map from each loop to its immediate children with respect to nesting.
   * This map is built during initialization of multi-loop polygons to
   * determine which are shells and which are holes, and then discarded.
   */
  alias LoopMap = S2Loop[][S2Loop];

  void insertLoop(S2Loop new_loop, S2Loop parent, ref LoopMap loop_map) {
    S2Loop[]* children;
    for (bool done = false; !done; ) {
      children = &loop_map.require(parent, new S2Loop[0]);
      done = true;
      foreach (S2Loop child; *children) {
        if (child.containsNested(new_loop)) {
          parent = child;
          done = false;
          break;
        }
      }
    }

    // Some of the children of the parent loop may now be children of
    // the new loop.
    S2Loop[]* new_children = &loop_map.require(new_loop, new S2Loop[0]);
    for (int i = 0; i < children.length;) {
      S2Loop child = (*children)[i];
      if (new_loop.containsNested(child)) {
        *new_children ~= child;
        *children = (*children).remove(i);
      } else {
        ++i;
      }
    }
    *children ~= new_loop;
  }

  void initializeLoops(LoopMap loop_map) {
    S2Loop[] loop_stack = [null];
    int depth = -1;
    while (!loop_stack.empty()) {
      S2Loop loop = loop_stack.back();
      loop_stack.popBack();
      if (loop !is null) {
        depth = loop.depth();
        _loops ~= loop;
      }
      if (loop !in loop_map) {
        continue;
      }
      S2Loop[] children = loop_map[loop];
      for (int i = cast(int) children.length - 1; i >= 0; --i) {
        S2Loop child = children[i];
        enforce(child !is null);
        child.setDepth(depth + 1);
        loop_stack ~= child;
      }
    }
  }

  /**
   * Add the polygon's loops to the S2ShapeIndex.  (The actual work of
   * building the index only happens when the index is first used.)
   */
  void initializeIndex()
  in {
    assert(_index.numShapeIds() == 0);
  } body {
    _index.add(new Shape(this));
    if (!S2POLYGON_LAZY_INDEXING) {
      _index.forceBuild();
    }
    if (flagsS2Debug && _s2debugOverride == S2Debug.ALLOW) {
      // Note that FLAGS_s2debug is false in optimized builds (by default).
      enforce(isValid());
    }
  }

  /**
   * When the loop is modified (Invert(), or Init() called again) then the
   * indexing structures need to be cleared since they become invalid.
   */
  void clearIndex() {
    atomicStore!(MemoryOrder.raw)(_unindexedContainsCalls, 0);
    _index.clear();
  }

  /// Initializes the polygon to the result of the given boolean operation.
  bool initializeToOperation(S2BooleanOperation.OpType op_type,
      S2Builder.SnapFunction snap_function, S2Polygon a, S2Polygon b) {
    S2BooleanOperation.Options options;
    options.setSnapFunction(snap_function);
    auto op = new S2BooleanOperation(op_type, new S2PolygonLayer(this), options);
    S2Error error;
    if (!op.build(a._index, b._index, error)) {
      logger.logFatal(op_type.to!string ~ " operation failed: ", error);
      return false;
    }
    return true;
  }

  /**
   * Initializes the polygon from input polygon "a" using the given S2Builder.
   * If the result has an empty boundary (no loops), also decides whether the
   * result should be the full polygon rather than the empty one based on the
   * area of the input polygon.  (See comments in InitToApproxIntersection.)
   */
  void initializeFromBuilder(in S2Polygon a, S2Builder builder) {
    builder.startLayer(new S2PolygonLayer(this));
    builder.addPolygon(a);
    S2Error error;
    if (!builder.build(error)) {
      logger.logFatal("Could not build polygon: ", error);
    }
    // If there are no loops, check whether the result should be the full
    // polygon rather than the empty one.  (See InitToApproxIntersection.)
    if (numLoops() == 0) {
      if (a._bound.area() > 2 * M_PI && a.getArea() > 2 * M_PI) invert();
    }
  }

  S2Polyline[] operationWithPolyline(
      S2BooleanOperation.OpType op_type, S2Builder.SnapFunction snap_function, in S2Polyline a) {
    S2BooleanOperation.Options options;
    options.setSnapFunction(snap_function);
    S2Polyline[] result;
    S2PolylineVectorLayer.Options layer_options;
    layer_options.setPolylineType(S2PolylineVectorLayer.Options.PolylineType.WALK);
    auto op = new S2BooleanOperation(
        op_type, new S2PolylineVectorLayer(&result, layer_options), options);
    auto a_index = new MutableS2ShapeIndex();
    a_index.add(new S2Polyline.Shape(a));
    S2Error error;
    if (!op.build(a_index, _index, error)) {
      logger.logFatal("Polyline ", op_type.to!string, " operation failed: ", error);
    }
    return result;
  }

  // TODO: Add when decode/encode are ready.
  // Encode the polygon's S2Points directly as three doubles using
  // (40 + 43 * num_loops + 24 * num_vertices) bytes.
  // void EncodeLossless(Encoder* encoder) const;

  // Decode a polygon encoded with EncodeLossless().  Used by the Decode and
  // DecodeWithinScope methods above.  The within_scope parameter specifies
  // whether to call DecodeWithinScope on the loops.
  // bool DecodeLossless(Decoder* const decoder, bool within_scope);

  // Encode the polygon's vertices using about 4 bytes / vertex plus 24 bytes /
  // unsnapped vertex. All the loop vertices must be converted first to the
  // S2XYZFaceSiTi format using S2Loop::GetXYZFaceSiTiVertices, and concatenated
  // in the all_vertices array.
  //
  // REQUIRES: snap_level >= 0.
  // void EncodeCompressed(Encoder* encoder, const S2XYZFaceSiTi* all_vertices,
  //                       int snap_level) const;

  // Decode a polygon encoded with EncodeCompressed().
  // bool DecodeCompressed(Decoder* decoder);

  /// See comments in InitToSimplifiedInCell.
  static S2Polyline[] simplifyEdgesInCell(
      in S2Polygon a, in S2Cell cell, double tolerance_uv, S1Angle snap_radius) {
    auto options = new S2Builder.Options(new IdentitySnapFunction(snap_radius));
    options.setSimplifyEdgeChains(true);
    auto builder = new S2Builder(options);
    // The output consists of a sequence of polylines.  Polylines consisting of
    // interior edges are simplified using S2Builder, while polylines consisting
    // of boundary edges are returned unchanged.
    S2Polyline[] polylines;
    for (int i = 0; i < a.numLoops(); ++i) {
      const S2Loop a_loop = a.loop(i);
      S2Point v0 = a_loop.orientedVertex(0);
      ubyte mask0 = getCellEdgeIncidenceMask(cell, v0, tolerance_uv);
      bool in_interior = false;  // Was the last edge an interior edge?
      for (int j = 1; j <= a_loop.numVertices(); ++j) {
        const S2Point v1 = a_loop.orientedVertex(j);
        ubyte mask1 = getCellEdgeIncidenceMask(cell, v1, tolerance_uv);
        if ((mask0 & mask1) != 0) {
          // This is an edge along the cell boundary.  Such edges do not get
          // simplified; we add them directly to the output.  (We create a
          // separate polyline for each edge to keep things simple.)  We call
          // ForceVertex on all boundary vertices to ensure that they don't
          // move, and so that nearby interior edges are snapped to them.
          enforce(!in_interior);
          builder.forceVertex(v1);
          polylines ~= new S2Polyline([v0, v1]);
        } else {
          // This is an interior edge.  If this is the first edge of an interior
          // chain, then start a new S2Builder layer.  Also ensure that any
          // polyline vertices on the boundary do not move, so that they will
          // still connect with any boundary edge(s) there.
          if (!in_interior) {
            S2Polyline polyline = new S2Polyline();
            builder.startLayer(new S2PolylineLayer(polyline));
            polylines ~= polyline;
            in_interior = true;
          }
          builder.addEdge(v0, v1);
          if (mask1 != 0) {
            builder.forceVertex(v1);
            in_interior = false;  // Terminate this polyline.
          }
        }
        v0 = v1;
        mask0 = mask1;
      }
    }
    S2Error error;
    if (!builder.build(error)) {
      logger.logFatal("InitToSimplifiedInCell failed: ", error);
    }
    return polylines;
  }

  // Internal implementation of intersect/subtract polyline functions above.
  // S2Polyline[] InternalClipPolyline(
  //     bool invert, const S2Polyline& a, S1Angle snap_radius) const;

  /**
   * Defines a total ordering on S2Loops that does not depend on the cyclic
   * order of loop vertices.  This function is used to choose which loop to
   * invert in the case where several loops have exactly the same area.
   */
  static int compareLoops(in S2Loop a, in S2Loop b) {
    if (a.numVertices() != b.numVertices()) {
      return a.numVertices() - b.numVertices();
    }
    int a_dir;
    int ai = a.getCanonicalFirstVertex(a_dir);
    int b_dir;
    int bi = b.getCanonicalFirstVertex(b_dir);
    if (a_dir != b_dir) return a_dir - b_dir;
    for (int n = a.numVertices(); --n >= 0; ai += a_dir, bi += b_dir) {
      if (a.vertex(ai) < b.vertex(bi)) return -1;
      if (a.vertex(ai) > b.vertex(bi)) return 1;
    }
    return 0;
  }

  S2Loop[] _loops;

  /**
   * Allows overriding the automatic validity checking controlled by the
   * --s2debug flag.
   */
  S2Debug _s2debugOverride;

  /**
   * True if InitOriented() was called and the given loops had inconsistent
   * orientations (i.e., it is not possible to construct a polygon such that
   * the interior is on the left-hand side of all loops).  We need to remember
   * this error so that it can be returned later by FindValidationError(),
   * since it is not possible to detect this error once the polygon has been
   * initialized.  This field is not preserved by Encode/Decode.
   */
  ubyte _errorInconsistentLoopOrientations;

  /// Cache for num_vertices().
  int _numVertices;

  /**
   * In general we build the index the first time it is needed, but we make an
   * exception for Contains(S2Point) because this method has a simple brute
   * force implementation that is also relatively cheap.  For this one method
   * we keep track of the number of calls made and only build the index once
   * enough calls have been made that we think an index would be worthwhile.
   */
  shared int _unindexedContainsCalls;

  /**
   * "bound_" is a conservative bound on all points contained by this polygon:
   * if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
   */
  S2LatLngRect _bound;

  /**
   * Since "bound_" is not exact, it is possible that a polygon A contains
   * another polygon B whose bounds are slightly larger.  "subregion_bound_"
   * has been expanded sufficiently to account for this error, i.e.
   * if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
   */
  S2LatLngRect _subregionBound;

  /// Spatial index containing this polygon.
  MutableS2ShapeIndex _index;

}

// Given a point "p" inside an S2Cell or on its boundary, return a mask
// indicating which of the S2Cell edges the point lies on.  All boundary
// comparisons are to within a maximum "u" or "v" error of "tolerance_uv".
// Bit "i" in the result is set if and only "p" is incident to the edge
// corresponding to S2Cell::edge(i).
private ubyte getCellEdgeIncidenceMask(in S2Cell cell, in S2Point p, double tolerance_uv) {
  import s2.s2coords : FaceXYZtoUV;
  ubyte mask = 0;
  R2Point uv;
  if (FaceXYZtoUV(cell.face(), p, uv)) {
    R2Rect bound = cell.getBoundUV();
    if (flagsS2Debug) debug enforce(bound.expanded(tolerance_uv).contains(uv));
    if (fabs(uv[1] - bound[1][0]) <= tolerance_uv) mask |= 1;
    if (fabs(uv[0] - bound[0][1]) <= tolerance_uv) mask |= 2;
    if (fabs(uv[1] - bound[1][1]) <= tolerance_uv) mask |= 4;
    if (fabs(uv[0] - bound[0][0]) <= tolerance_uv) mask |= 8;
  }
  return mask;
}
