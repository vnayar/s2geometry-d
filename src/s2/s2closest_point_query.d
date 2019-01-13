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
//
// See S2ClosestPointQuery (defined below) for an overview.

module s2.s2closest_point_query;

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cell;
import s2.s2closest_point_query_base;
import s2.s2min_distance_targets;
import s2.s2point;
import s2.s2point_index;
import s2.s2shape_index;

/**
 * Options that control the set of points returned.  Note that by default
 * *all* points are returned, so you will always want to set either the
 * max_points() option or the max_distance() option (or both).
 *
 * This class is also available as S2ClosestPointQuery<Data>::Options.
 * (It is defined here to avoid depending on the "Data" template argument.)
 */
class S2ClosestPointQueryOptions : S2ClosestPointQueryBaseOptions!S2MinDistance {
public:
  alias Distance = S2MinDistance;
  alias Base = S2ClosestPointQueryBaseOptions!Distance;

  // See S2ClosestPointQueryBaseOptions for the full set of options.

  /**
   * Specifies that only points whose distance to the target is less than
   * "max_distance" should be returned.
   *
   * Note that points whose distance is exactly equal to "max_distance" are
   * not returned.  Normally this doesn't matter, because distances are not
   * computed exactly in the first place, but if such points are needed then
   * see set_inclusive_max_distance() below.
   *
   * DEFAULT: Distance::Infinity()
   */
  void setMaxDistance(S1ChordAngle max_distance) {
    super.setMaxDistance(Distance(max_distance));
  }

  // Avoid hiding the base class method setMaxDistance(S1ChordAngle).
  alias setMaxDistance = Base.setMaxDistance;

  /**
   * Like set_max_distance(), except that points whose distance is exactly
   * equal to "max_distance" are also returned.  Equivalent to calling
   * set_max_distance(max_distance.Successor()).
   */
  void setInclusiveMaxDistance(S1ChordAngle max_distance) {
    setMaxDistance(max_distance.successor());
  }

  /**
   * Like set_inclusive_max_distance(), except that "max_distance" is also
   * increased by the maximum error in the distance calculation.  This ensures
   * that all points whose true distance is less than or equal to
   * "max_distance" will be returned (along with some points whose true
   * distance is slightly greater).
   *
   * Algorithms that need to do exact distance comparisons can use this
   * option to find a set of candidate points that can then be filtered
   * further (e.g., using s2pred::CompareDistance).
   */
  void setConservativeMaxDistance(S1ChordAngle max_distance) {
    import s2.s2edge_distances : getUpdateMinDistanceMaxError;
    setMaxDistance(
        Distance(max_distance.plusError(getUpdateMinDistanceMaxError(max_distance)).successor()));
  }

  /**
   * Versions of set_max_distance that take an S1Angle argument.  (Note that
   * these functions require a conversion, and that the S1ChordAngle versions
   * are preferred.)
   */
  void setMaxDistance(S1Angle max_distance) {
    super.setMaxDistance(Distance(S1ChordAngle(max_distance)));
  }

  void setInclusiveMaxDistance(S1Angle max_distance) {
    setInclusiveMaxDistance(S1ChordAngle(max_distance));
  }

  void setConservativeMaxDistance(S1Angle max_distance) {
    setConservativeMaxDistance(S1ChordAngle(max_distance));
  }

  /// See S2ClosestPointQueryBaseOptions for documentation.
  // using Base::set_max_error;              // S1Chordangle version
  void set_max_error(S1Angle max_error) {
    super.setMaxError(S1ChordAngle(max_error));
  }

  /// Inherited options (see s2closest_point_query_base.h for details):
  // using Base::set_max_points;
  // using Base::set_region;
  // using Base::set_use_brute_force;
}

/**
 * S2ClosestPointQueryTarget represents the geometry to which the distance is
 * measured.  There are subtypes for measuring the distance to a point, an
 * edge, an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
 */
alias S2ClosestPointQueryTarget = S2MinDistanceTarget;

/**
 * Target subtype that computes the closest distance to a point.
 *
 * This class is also available as S2ClosestPointQuery<Data>::PointTarget.
 * (It is defined here to avoid depending on the "Data" template argument.)
 */
final class S2ClosestPointQueryPointTarget : S2MinDistancePointTarget {
public:
  this(in S2Point point) {
    super(point);
  }

  override
  int maxBruteForceIndexSize() const {
    // Using BM_FindClosest (which finds the single closest point), the
    // break-even points are approximately X, Y, and Z points for grid,
    // fractal, and regular loop geometry respectively.
    //
    // TODO(ericv): Adjust using benchmarks.
    return 150;
  }
}

/**
 * Target subtype that computes the closest distance to an edge.
 *
 * This class is also available as S2ClosestPointQuery<Data>::EdgeTarget.
 * (It is defined here to avoid depending on the "Data" template argument.)
 */
final class S2ClosestPointQueryEdgeTarget : S2MinDistanceEdgeTarget {
public:
  this(in S2Point a, in S2Point b) {
    super(a, b);
  }

  override
  int maxBruteForceIndexSize() const {
    // Using BM_FindClosestToEdge (which finds the single closest point), the
    // break-even points are approximately X, Y, and Z points for grid,
    // fractal, and regular loop geometry respectively.
    //
    // TODO(ericv): Adjust using benchmarks.
    return 100;
  }
}

/**
 * Target subtype that computes the closest distance to an S2Cell
 * (including the interior of the cell).
 *
 * This class is also available as S2ClosestPointQuery<Data>::CellTarget.
 * (It is defined here to avoid depending on the "Data" template argument.)
 */
final class S2ClosestPointQueryCellTarget : S2MinDistanceCellTarget {
 public:
  this(in S2Cell cell) {
    super(cell);
  }

  override
  int maxBruteForceIndexSize() const {
    // Using BM_FindClosestToCell (which finds the single closest point), the
    // break-even points are approximately X, Y, and Z points for grid,
    // fractal, and regular loop geometry respectively.
    //
    // TODO(ericv): Adjust using benchmarks.
    return 50;
  }
}

/**
 * Target subtype that computes the closest distance to an S2ShapeIndex
 * (an arbitrary collection of points, polylines, and/or polygons).
 *
 * By default, distances are measured to the boundary and interior of
 * polygons in the S2ShapeIndex rather than to polygon boundaries only.
 * If you wish to change this behavior, you may call
 *
 *   target.set_include_interiors(false);
 *
 * (see S2MinDistanceShapeIndexTarget for details).
 *
 * This class is also available as S2ClosestPointQuery<Data>::ShapeIndexTarget.
 * (It is defined here to avoid depending on the "Data" template argument.)
 */
final class S2ClosestPointQueryShapeIndexTarget : S2MinDistanceShapeIndexTarget {
public:
  this(S2ShapeIndex index) {
    super(index);
  }

  override
  int maxBruteForceIndexSize() const {
    // For BM_FindClosestToSameSizeAbuttingIndex (which uses a nearby
    // S2ShapeIndex target of similar complexity), the break-even points are
    // approximately X, Y, and Z points for grid, fractal, and regular loop
    // geometry respectively.
    //
    // TODO(ericv): Adjust using benchmarks.
    return 30;
  }
}

/**
 * Given a set of points stored in an S2PointIndex, S2ClosestPointQuery
 * provides methods that find the closest point(s) to a given query point
 * or query edge.  Example usage:
 *
 * void Test(const vector<S2Point>& index_points,
 *           const vector<S2Point>& target_points) {
 *   // The template argument allows auxiliary data to be attached to each
 *   // point (in this case, the array index).
 *   S2PointIndex<int> index;
 *   for (const S2Point& point : index_points) {
 *     index.Add(point, i);
 *   }
 *   S2ClosestPointQuery<int> query(&index);
 *   query.mutable_options()->set_max_points(5);
 *   for (const S2Point& target_point : target_points) {
 *     S2ClosestPointQueryPointTarget target(target_point);
 *     for (const auto& result : query.FindClosestPoints(&target)) {
 *       // The Result class contains the following methods:
 *       //   distance() is the distance to the target.
 *       //   point() is the indexed point.
 *       //   data() is the auxiliary data.
 *       DoSomething(target_point, result);
 *     }
 *   }
 * }
 *
 * You can find either the k closest points, or all points within a given
 * radius, or both (i.e., the k closest points up to a given maximum radius).
 * E.g. to find all the points within 5 kilometers, call
 *
 *   query.mutable_options()->set_max_distance(
 *       S2Earth::ToAngle(util::units::Kilometers(5)));
 *
 * By default *all* points are returned, so you should always specify either
 * max_points() or max_distance() or both.  There is also a FindClosestPoint()
 * convenience method that returns only the closest point.
 *
 * You can restrict the results to an arbitrary S2Region, for example:
 *
 *   S2LatLngRect rect(...);
 *   query.set_region(&rect);  // Does *not* take ownership.
 *
 * To find the closest points to a query edge rather than a point, use:
 *
 *   S2ClosestPointQueryEdgeTarget target(v0, v1);
 *   query.FindClosestPoints(&target);
 *
 * The implementation is designed to be fast for both small and large
 * point sets.
 */
class S2ClosestPointQuery(Data) {
public:
  // See S2ClosestPointQueryBase for full documentation.

  alias Index = S2PointIndex!Data;
  alias PointData = Index.PointData;

  // S2MinDistance is a thin wrapper around S1ChordAngle that implements the
  // Distance concept required by S2ClosestPointQueryBase.
  alias Distance = S2MinDistance;
  alias Base = S2ClosestPointQueryBase!(Distance, Data);

  /**
   * Each "Result" object represents a closest point.  Here are its main
   * methods (see S2ClosestPointQueryBase::Result for details):
   *
   *   // The distance from the target to this point.
   *   S1ChordAngle distance() const;
   *
   *   // The point itself.
   *   const S2Point& point() const;
   *
   *   // The client-specified data associated with this point.
   *   const Data& data() const;
   */
  alias Result = Base.Result;

  alias Options = S2ClosestPointQueryOptions;

  // The available target types (see definitions above).
  alias Target = S2ClosestPointQueryTarget;
  alias PointTarget = S2ClosestPointQueryPointTarget;
  alias EdgeTarget = S2ClosestPointQueryEdgeTarget;
  alias CellTarget = S2ClosestPointQueryCellTarget;
  alias ShapeIndexTarget = S2ClosestPointQueryShapeIndexTarget;

  /**
   * Convenience constructor that calls Init().  Options may be specified here
   * or changed at any time using the mutable_options() accessor method.
   */
  this(Index index, in Options options) {
    if (options is null)
      options = new Options();
    initialize(index, options);
  }

  /// Default constructor; requires initialize() to be called.
  this() {}

  /**
   * Initializes the query.  Options may be specified here or changed at any
   * time using the mutable_options() accessor method.
   *
   * REQUIRES: "index" must persist for the lifetime of this object.
   * REQUIRES: ReInit() must be called if "index" is modified.
   */
  void initialize(Index index, Options options) {
    if (options is null)
      options = new Options();

    _options = options;
    _base.initialize(index);
  }

  /**
   * Reinitializes the query.  This method must be called whenever the
   * underlying index is modified.
   */
  void reInitialize() {
    _base.reInitialize();
  }

  /// Returns a reference to the underlying S2PointIndex.
  const(Index) index() const {
    return _base.index();
  }

  /// Returns the query options.  Options can be modifed between queries.
  const(Options) options() const {
    return _options;
  }

  Options mutableOptions() {
    return _options;
  }

  /**
   * Returns the closest points to the given target that satisfy the given
   * options.  This method may be called multiple times.
   */
  Result[] findClosestPoints(Target target) {
    return _base.findClosestPoints(target, _options);
  }

  /**
   * This version can be more efficient when this method is called many times,
   * since it does not require allocating a new vector on each call.
   */
  void findClosestPoints(Target target, ref Result[] results) {
    _base.findClosestPoints(target, _options, results);
  }

  //////////////////////// Convenience Methods ////////////////////////

  /**
   * Returns the closest point to the target.  If no point satisfies the search
   * criteria, then a Result object with distance() == Infinity() and
   * is_empty() == true is returned.
   */
  Result findClosestPoint(Target target) {
    static assert(__traits(classInstanceSize, Options) <= 32, "Consider not copying Options here");
    Options tmp_options = _options.dup;
    tmp_options.setMmaxPoints(1);
    return _base.findClosestPoint(target, tmp_options);
  }

  /**
   * Returns the minimum distance to the target.  If the index or target is
   * empty, returns S1ChordAngle::Infinity().
   *
   * Use IsDistanceLess() if you only want to compare the distance against a
   * threshold value, since it is often much faster.
   */
  S1ChordAngle getDistance(Target target) {
    return findClosestPoint(target).distance.s1ChordAngle;
  }

  /**
   * Returns true if the distance to "target" is less than "limit".
   *
   * This method is usually much faster than GetDistance(), since it is much
   * less work to determine whether the minimum distance is above or below a
   * threshold than it is to calculate the actual minimum distance.
   */
  bool isDistanceLess(Target target, S1ChordAngle limit) {
    static assert(__traits(classInstanceSize, Options) <= 32, "Consider not copying Options here");
    Options tmp_options = _options.dup;
    tmp_options.setMaxPoints(1);
    tmp_options.setMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight());
    return !_base.findClosestPoint(target, tmp_options).isEmpty();
  }

  // TODO: Resume here.

  /**
   * Like IsDistanceLess(), but also returns true if the distance to "target"
   * is exactly equal to "limit".
   */
  bool isDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    static assert(__traits(classInstanceSize, Options) <= 32, "Consider not copying Options here");
    Options tmp_options = _options.dup;
    tmp_options.setMaxPoints(1);
    tmp_options.setInclusiveMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight());
    return !_base.findClosestPoint(target, tmp_options).isEmpty();
  }

  /**
   * Like IsDistanceLessOrEqual(), except that "limit" is increased by the
   * maximum error in the distance calculation.  This ensures that this
   * function returns true whenever the true, exact distance is less than
   * or equal to "limit".
   *
   * For example, suppose that we want to test whether two geometries might
   * intersect each other after they are snapped together using S2Builder
   * (using the IdentitySnapFunction with a given "snap_radius").  Since
   * S2Builder uses exact distance predicates (s2predicates.h), we need to
   * measure the distance between the two geometries conservatively.  If the
   * distance is definitely greater than "snap_radius", then the geometries
   * are guaranteed to not intersect after snapping.
   */
  bool isConservativeDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    static_assert(sizeof(Options) <= 32, "Consider not copying Options here");
    Options tmp_options = _options;
    tmp_options.setMaxPoints(1);
    tmp_options.setConservativeMaxDistance(limit);
    tmp_options.setMaxError(S1ChordAngle.straight());
    return !_base.findClosestPoint(target, tmp_options).isEmpty();
  }

  ///////////////////////// Deprecated Methods //////////////////////////

private:
  Options _options;
  Base _base;

  // Deprecated methods that return results using the result interface require
  // keeping a copy of the result vector.
  Result[] _results;
}
