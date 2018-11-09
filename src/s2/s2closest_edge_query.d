// Copyright 2013 Google Inc. All Rights Reserved.
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

module s2.s2closest_edge_query;

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2closest_edge_query_base;
import s2.s2edge_distances : getUpdateMinDistanceMaxError, project;
import s2.s2min_distance_targets;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;


/**
 * S2ClosestEdgeQuery is a helper class for finding the closest edge(s) to a
 * given point, edge, S2Cell, or geometry collection.  For example, given a
 * set of polylines, the following code efficiently finds the closest 5 edges
 * to a query point:
 *
 * void Test(const vector<S2Polyline*>& polylines, const S2Point& point) {
 *   MutableS2ShapeIndex index;
 *   for (S2Polyline* polyline : polylines) {
 *     index.Add(new S2Polyline::Shape(polyline));
 *   }
 *   S2ClosestEdgeQuery query(&index);
 *   query.mutable_options()->set_max_edges(5);
 *   S2ClosestEdgeQuery::PointTarget target(point);
 *   for (const auto& result : query.FindClosestEdges(&target)) {
 *     // The Result struct contains the following fields:
 *     //   "distance" is the distance to the edge.
 *     //   "shape_id" identifies the S2Shape containing the edge.
 *     //   "edge_id" identifies the edge with the given shape.
 *     // The following convenience methods may also be useful:
 *     //   query.GetEdge(result) returns the endpoints of the edge.
 *     //   query.Project(point, result) computes the closest point on the
 *     //       result edge to the given target point.
 *     int polyline_index = result.shape_id;
 *     int edge_index = result.edge_id;
 *     S1ChordAngle distance = result.distance;  // Use ToAngle() for S1Angle.
 *     S2Shape::Edge edge = query.GetEdge(result);
 *     S2Point closest_point = query.Project(point, result);
 *   }
 * }
 *
 * You can find either the k closest edges, or all edges within a given
 * radius, or both (i.e., the k closest edges up to a given maximum radius).
 * E.g. to find all the edges within 5 kilometers, call
 *
 *   query.mutable_options()->set_max_distance(
 *       S2Earth::ToAngle(util::units::Kilometers(5)));
 *
 * By default *all* edges are returned, so you should always specify either
 * max_edges() or max_distance() or both.  There is also a FindClosestEdge()
 * convenience method that returns only the closest edge.
 *
 * Note that by default, distances are measured to the boundary and interior
 * of polygons.  For example, if a point is inside a polygon then its distance
 * is zero.  To change this behavior, call set_include_interiors(false).
 *
 * If you only need to test whether the distance is above or below a given
 * threshold (e.g., 10 km), you can use the IsDistanceLess() method.  This is
 * much faster than actually calculating the distance with FindClosestEdge(),
 * since the implementation can stop as soon as it can prove that the minimum
 * distance is either above or below the threshold.
 *
 * To find the closest edges to a query edge rather than a point, use:
 *
 *   S2ClosestEdgeQuery::EdgeTarget target(v0, v1);
 *   query.FindClosestEdges(&target);
 *
 * Similarly you can find the closest edges to an S2Cell by using an
 * S2ClosestEdgeQuery::CellTarget, and you can find the closest edges to an
 * arbitrary collection of points, polylines, and polygons by using an
 * S2ClosestEdgeQuery::ShapeIndexTarget.
 *
 * The implementation is designed to be fast for both simple and complex
 * geometric objects.
 */
class S2ClosestEdgeQuery {
public:
  // See S2ClosestEdgeQueryBase for full documentation.

  // S2MinDistance is a thin wrapper around S1ChordAngle that implements the
  // Distance concept required by S2ClosestPointQueryBase.
  alias Distance = S2MinDistance;
  alias Base = S2ClosestEdgeQueryBase!Distance;

  // Each "Result" object represents a closest edge.  It has the following
  // fields:
  //
  //   S1ChordAngle distance;  // The distance from the target to this edge.
  //   int32 shape_id;         // Identifies an indexed shape.
  //   int32 edge_id;          // Identifies an edge within the shape.
  alias Result = Base.Result;

  /**
   * Options that control the set of edges returned.  Note that by default
   * *all* edges are returned, so you will always want to set either the
   * max_edges() option or the max_distance() option (or both).
   */
  final class Options : Base.Options {
  public:
    // See S2ClosestEdgeQueryBase::Options for the full set of options.

    this() { }

    this(Options options) {
      super(options);
    }

    /**
     * Specifies that only edges whose distance to the target is less than
     * "max_distance" should be returned.
     *
     * Note that edges whose distance is exactly equal to "max_distance" are
     * not returned.  Normally this doesn't matter, because distances are not
     * computed exactly in the first place, but if such edges are needed then
     * see set_inclusive_max_distance() below.
     *
     * DEFAULT: Distance::Infinity()
     */
    void setMaxDistance(S1ChordAngle max_distance) {
      super.setMaxDistance(Distance(max_distance));
    }

    /**
     * Like set_max_distance(), except that edges whose distance is exactly
     * equal to "max_distance" are also returned.  Equivalent to calling
     * set_max_distance(max_distance.Successor()).
     */
    void setInclusiveMaxDistance(S1ChordAngle max_distance) {
      setMaxDistance(max_distance.successor());
    }

    // Like set_inclusive_max_distance(), except that "max_distance" is also
    // increased by the maximum error in the distance calculation.  This
    // ensures that all edges whose true distance is less than or equal to
    // "max_distance" will be returned (along with some edges whose true
    // distance is slightly greater).
    //
    // Algorithms that need to do exact distance comparisons can use this
    // option to find a set of candidate edges that can then be filtered
    // further (e.g., using s2pred::CompareDistance).
    void setConservativeMaxDistance(S1ChordAngle max_distance) {
      setMaxDistance(
          Distance(
              max_distance.plusError(
                  getUpdateMinDistanceMaxError(max_distance))
              .successor()));

    }

    // Versions of set_max_distance that take an S1Angle argument.  (Note that
    // these functions require a conversion, and that the S1ChordAngle versions
    // are preferred.)
    void setMaxDistance(S1Angle max_distance) {
      super.setMaxDistance(Distance(S1ChordAngle(max_distance)));
    }

    void setInclusiveMaxDistance(S1Angle max_distance) {
      setInclusiveMaxDistance(S1ChordAngle(max_distance));
    }

    void setConservativeMaxDistance(S1Angle max_distance) {
      setConservativeMaxDistance(S1ChordAngle(max_distance));
    }

    // See S2ClosestEdgeQueryBase::Options for documentation.
    void setMaxError(S1Angle max_error) {
      super.setMaxError(S1ChordAngle(max_error));
    }

  }

  // "Target" represents the geometry to which the distance is measured.
  // There are subtypes for measuring the distance to a point, an edge, an
  // S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
  alias Target = S2MinDistanceTarget;

  /// Target subtype that computes the closest distance to a point.
  final static class PointTarget : S2MinDistancePointTarget {
  public:
    this(in S2Point point) {
      super(point);
    }

    override
    int maxBruteForceIndexSize() const {
      // Using BM_FindClosest (which finds the single closest edge), the
      // break-even points are approximately 80, 100, and 250 edges for point
      // cloud, fractal, and regular loop geometry respectively.
      return 120;
    }
  }

  // Target subtype that computes the closest distance to an edge.
  final static class EdgeTarget : S2MinDistanceEdgeTarget {
  public:
    this(in S2Point a, in S2Point b) {
      super(a, b);
    }

    override
    int maxBruteForceIndexSize() const {
      // Using BM_FindClosestToEdge (which finds the single closest edge), the
      // break-even points are approximately 40, 50, and 100 edges for point
      // cloud, fractal, and regular loop geometry respectively.
      return 60;
    }
  }

  // Target subtype that computes the closest distance to an S2Cell
  // (including the interior of the cell).
  final static class CellTarget : S2MinDistanceCellTarget {
  public:
    this(in S2Cell cell) {
      super(cell);
    }

    override
    int maxBruteForceIndexSize() const {
      // Using BM_FindClosestToCell (which finds the single closest edge), the
      // break-even points are approximately 20, 25, and 40 edges for point cloud,
      // fractal, and regular loop geometry respectively.
      return 30;
    }
  }

  // Target subtype that computes the closest distance to an S2ShapeIndex
  // (an arbitrary collection of points, polylines, and/or polygons).
  //
  // By default, distances are measured to the boundary and interior of
  // polygons in the S2ShapeIndex rather than to polygon boundaries only.
  // If you wish to change this behavior, you may call
  //
  //   target.set_include_interiors(false);
  //
  // (see S2MinDistanceShapeIndexTarget for details).
  final static class ShapeIndexTarget : S2MinDistanceShapeIndexTarget {
  public:
    this(S2ShapeIndex index) {
      super(index);
    }

    override
    int maxBruteForceIndexSize() const {
      // For BM_FindClosestToSameSizeAbuttingIndex (which uses two nearby indexes
      // with similar edge counts), the break-even points are approximately 20,
      // 30, and 40 edges for point cloud, fractal, and regular loop geometry
      // respectively.
      return 25;
    }
  }

  // Convenience constructor that calls Init().  Options may be specified here
  // or changed at any time using the mutable_options() accessor method.
  this(S2ShapeIndex index, Options options = null) {
    if (options is null)
      options = new Options();
    initialize(index, options);
  }

  // Default constructor; requires Init() to be called.
  this() { }

  // Initializes the query.  Options may be specified here or changed at any
  // time using the mutable_options() accessor method.
  //
  // REQUIRES: "index" must persist for the lifetime of this object.
  // REQUIRES: ReInit() must be called if "index" is modified.
  void initialize(S2ShapeIndex index, Options options) {
    _options = options is null ? new Options() : options;
    _base.initialize(index);
  }

  // Reinitializes the query.  This method must be called whenever the
  // underlying S2ShapeIndex is modified.
  void reInititialize() {
    _base.reInitialize();
  }

  // Returns a reference to the underlying S2ShapeIndex.
  const(S2ShapeIndex) index() const {
    return _base.index();
  }

  // Returns the query options.  Options can be modifed between queries.
  const(Options) options() const {
    return _options;
  }

  Options mutableOptions() {
    return _options;
  }

  // Returns the closest edges to the given target that satisfy the given
  // options.  This method may be called multiple times.
  //
  // Note that if options().include_interiors() is true, the result vector may
  // include some entries with edge_id == -1.  This indicates that the target
  // intersects the indexed polygon with the given shape_id.
  Result[] findClosestEdges(Target target) {
    return _base.findClosestEdges(target, _options);
  }

  // This version can be more efficient when this method is called many times,
  // since it does not require allocating a new vector on each call.
  void findClosestEdges(Target target, out Result[] results) {
    _base.findClosestEdges(target, _options, results);
  }

  //////////////////////// Convenience Methods ////////////////////////

  // Returns the closest edge to the target.  If no edge satisfies the search
  // criteria, then the Result object will have distance == Infinity() and
  // shape_id == edge_id == -1.
  //
  // Note that if options.include_interiors() is true, edge_id == -1 is also
  // used to indicate that the target intersects an indexed polygon (but in
  // that case distance == Zero() and shape_id >= 0).
  Result findClosestEdge(Target target) {
    import std.conv;
    static assert(
        __traits(classInstanceSize, Options) <= 48,
        "Consider not copying Options (" ~ to!string(__traits(classInstanceSize, Options))
            ~ " bytes) here");
    auto tmp_options = new Options(_options);
    tmp_options.setMaxEdges(1);
    return _base.findClosestEdge(target, tmp_options);
  }

  // Returns the minimum distance to the target.  If the index or target is
  // empty, returns S1ChordAngle::Infinity().
  //
  // Use IsDistanceLess() if you only want to compare the distance against a
  // threshold value, since it is often much faster.
  S1ChordAngle getDistance(Target target) {
    return findClosestEdge(target).distance;
  }

  // Returns true if the distance to "target" is less than "limit".
  //
  // This method is usually much faster than GetDistance(), since it is much
  // less work to determine whether the minimum distance is above or below a
  // threshold than it is to calculate the actual minimum distance.
  bool isDistanceLess(Target target, S1ChordAngle limit) {
    import std.conv;
    static assert(
        __traits(classInstanceSize, Options) <= 48,
        "Consider not copying Options (" ~ to!string(__traits(classInstanceSize, Options))
            ~ " bytes) here");
    Options tmp_options = new Options(_options);
    tmp_options.setMaxEdges(1);
    tmp_options.setMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight().toS1Angle());
    return _base.findClosestEdge(target, tmp_options).shapeId >= 0;
  }

  // Like IsDistanceLess(), but also returns true if the distance to "target"
  // is exactly equal to "limit".
  bool isDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    import std.conv;
    static assert(
        __traits(classInstanceSize, Options) <= 48,
        "Consider not copying Options (" ~ to!string(__traits(classInstanceSize, Options))
            ~ " bytes) here");
    Options tmp_options = new Options(_options);
    tmp_options.setMaxEdges(1);
    tmp_options.setInclusiveMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight().toS1Angle());
    return _base.findClosestEdge(target, tmp_options).shapeId >= 0;
  }

  // Like IsDistanceLessOrEqual(), except that "limit" is increased by the
  // maximum error in the distance calculation.  This ensures that this
  // function returns true whenever the true, exact distance is less than
  // or equal to "limit".
  //
  // For example, suppose that we want to test whether two geometries might
  // intersect each other after they are snapped together using S2Builder
  // (using the IdentitySnapFunction with a given "snap_radius").  Since
  // S2Builder uses exact distance predicates (s2predicates.h), we need to
  // measure the distance between the two geometries conservatively.  If the
  // distance is definitely greater than "snap_radius", then the geometries
  // are guaranteed to not intersect after snapping.
  bool isConservativeDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    import std.conv;
    static assert(
        __traits(classInstanceSize, Options) <= 48,
        "Consider not copying Options (" ~ to!string(__traits(classInstanceSize, Options))
            ~ " bytes) here");
    Options tmp_options = _options;
    tmp_options.setMaxEdges(1);
    tmp_options.setConservativeMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight().toS1Angle());
    return _base.findClosestEdge(target, tmp_options).shapeId >= 0;
  }

  // Returns the endpoints of the given result edge.
  //
  // CAVEAT: If options().include_interiors() is true, then clients must not
  // pass this method any Result objects that correspond to shape interiors,
  // i.e. those where result.edge_id < 0.
  //
  // REQUIRES: result.edge_id >= 0
  S2Shape.Edge getEdge(in Result result) const {
    return index().shape(result.shapeId).edge(result.edgeId);
  }

  // Returns the point on given result edge that is closest to "point".
  S2Point project(in S2Point point, in Result result) const {
    if (result.edgeId < 0) return point;
    auto edge = getEdge(result);
    return .project(point, edge.v0, edge.v1);
  }

 private:
  Options _options;
  Base _base;
}
