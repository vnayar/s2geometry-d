// Copyright 2017 Google Inc. All Rights Reserved.
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
// See S2ClosestPointQueryBase (defined below) for an overview.

module s2.s2closest_point_query_base;

import s2.logger;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2distance_target;
import s2.s2edge_distances;
import s2.s2point;
import s2.s2point_index;
import s2.s2region;
import s2.s2region_coverer;

import std.algorithm : isSorted, reverse, sort, uniq;
import std.container.binaryheap;
import std.exception;
import std.range;
import std.typecons : Rebindable;


/**
 * Options that control the set of points returned.  Note that by default
 * *all* points are returned, so you will always want to set either the
 * max_points() option or the max_distance() option (or both).
 *
 * This class is also available as S2ClosestPointQueryBase<Data>::Options.
 * (It is defined here to avoid depending on the "Data" template argument.)
 *
 * The Distance template argument is described below.
 */
class S2ClosestPointQueryBaseOptions(Distance) {
public:
  alias Delta = Distance.Delta;

  this() {}

  /**
   * Specifies that at most "max_points" points should be returned.
   *
   * REQUIRES: max_points >= 1
   * DEFAULT: numeric_limits<int>::max()
   */
  int maxPoints() const {
    return _maxPoints;
  }

  void setMaxPoints(int max_points)
  in {
    assert(max_points >= 1);
  } body {
    _maxPoints = max_points;
  }

  enum int MAX_MAX_POINTS = int.max;

  /**
   * Specifies that only points whose distance to the target is less than
   * "max_distance" should be returned.
   *
   * Note that points whose distance is exactly equal to "max_distance" are
   * not returned.  In most cases this doesn't matter (since distances are
   * not computed exactly in the first place), but if such points are needed
   * then you can retrieve them by specifying "max_distance" as the next
   * largest representable Distance.  For example, if Distance is an
   * S1ChordAngle then you can specify max_distance.Successor().
   *
   * DEFAULT: Distance::Infinity()
   */
  Distance maxDistance() const {
    return _maxDistance;
  }

  void setMaxDistance(Distance max_distance) {
    _maxDistance = max_distance;
  }

  /**
   * Specifies that points up to max_error() further away than the true
   * closest points may be substituted in the result set, as long as such
   * points satisfy all the remaining search criteria (such as max_distance).
   * This option only has an effect if max_points() is also specified;
   * otherwise all points closer than max_distance() will always be returned.
   *
   * Note that this does not affect how the distance between points is
   * computed; it simply gives the algorithm permission to stop the search
   * early as soon as the best possible improvement drops below max_error().
   *
   * This can be used to implement distance predicates efficiently.  For
   * example, to determine whether the minimum distance is less than D, the
   * IsDistanceLess() method sets max_points() == 1 and max_distance() ==
   * max_error() == D.  This causes the algorithm to terminate as soon as it
   * finds any point whose distance is less than D, rather than continuing to
   * search for a point that is even closer.
   *
   * DEFAULT: Distance::Delta::Zero()
   */
  Delta maxError() const {
    return _maxError;
  }

  void setMaxError(Delta max_error) {
    _maxError = max_error;
  }

  /**
   * Specifies that points must be contained by the given S2Region.  "region"
   * is owned by the caller and must persist during the lifetime of this
   * object.  The value may be changed between calls to FindClosestPoints(),
   * or reset by calling set_region(nullptr).
   *
   * Note that if you want to set the region to a disc around the target
   * point, it is faster to use set_max_distance() instead.  You can also call
   * both methods, e.g. to set a maximum distance and also require that points
   * lie within a given rectangle.
   */
  inout(S2Region) region() inout {
    return _region;
  }

  void setRegion(S2Region region) {
    _region = region;
  }

  /**
   * Specifies that distances should be computed by examining every point
   * rather than using the S2ShapeIndex.  This is useful for testing,
   * benchmarking, and debugging.
   *
   * DEFAULT: false
   */
  bool useBruteForce() const {
    return _useBruteForce;
  }

  void setUseBruteForce(bool use_brute_force) {
    _useBruteForce = use_brute_force;
  }

  ThisT dup(this ThisT)() {
    ThisT d = new ThisT();
    d._maxDistance = _maxDistance;
    d._maxError = _maxError;
    d._region = _region;
    d._maxPoints = _maxPoints;
    d._useBruteForce = _useBruteForce;
    return d;
  }

  override
  string toString() const {
    import std.conv;
    return "S2ClosestPointQueryBaseOptions"
        ~ "[ maxDistance=" ~ _maxDistance.to!string
        ~ ", maxError=" ~ _maxError.to!string
        ~ ", region=" ~ _region.to!string
        ~ ", maxPoints=" ~ _maxPoints.to!string
        ~ ", useBruteForce=" ~ _useBruteForce.to!string ~ " ]";
  }

protected:
  Distance _maxDistance = Distance.infinity();
  Delta _maxError = Delta.zero();
  S2Region _region = null;
  int _maxPoints = MAX_MAX_POINTS;
  bool _useBruteForce = false;
}

/**
 * S2ClosestPointQueryBase is a templatized class for finding the closest
 * point(s) to a given target.  It is not intended to be used directly, but
 * rather to serve as the implementation of various specialized classes with
 * more convenient APIs (such as S2ClosestPointQuery).  It is flexible enough
 * so that it can be adapted to compute maximum distances and even potentially
 * Hausdorff distances.
 *
 * By using the appropriate options, this class can answer questions such as:
 *
 *  - Find the minimum distance between a point collection A and a target B.
 *  - Find all points in collection A that are within a distance D of target B.
 *  - Find the k points of collection A that are closest to a given point P.
 *
 * The target is any class that implements the S2DistanceTarget interface.
 * There are predefined targets for points, edges, S2Cells, and S2ShapeIndexes
 * (arbitrary collctions of points, polylines, and polygons).
 *
 * The Distance template argument is used to represent distances.  Usually
 * this type is a thin wrapper around S1ChordAngle, but another distance type
 * may be substituted as long as it implements the API below.  This can be
 * used to change the comparison function (e.g., to find the furthest edges
 * from the target), or to get more accuracy if desired.
 *
 * The Distance concept is as follows:
 *
 * class Distance {
 *  public:
 *   // Default and copy constructors, assignment operator:
 *   Distance();
 *   Distance(const Distance&);
 *   Distance& operator=(const Distance&);
 *
 *   // Factory methods:
 *   static Distance Zero();      // Returns a zero distance.
 *   static Distance Infinity();  // Larger than any valid distance.
 *   static Distance Negative();  // Smaller than any valid distance.
 *
 *   // Comparison operators:
 *   friend bool operator==(Distance x, Distance y);
 *   friend bool operator<(Distance x, Distance y);
 *
 *   // Delta represents the positive difference between two distances.
 *   // It is used together with operator-() to implement Options::max_error().
 *   // Typically Distance::Delta is simply S1ChordAngle.
 *   class Delta {
 *    public:
 *     Delta();
 *     Delta(const Delta&);
 *     Delta& operator=(const Delta&);
 *     friend bool operator==(Delta x, Delta y);
 *     static Delta Zero();
 *   };
 *
 *   // Subtraction operator.  Note that the second argument represents a
 *   // delta between two distances.  This distinction is important for
 *   // classes that compute maximum distances (e.g., S2FurthestEdgeQuery).
 *   friend Distance operator-(Distance x, Delta delta);
 *
 *   // Method that returns an upper bound on the S1ChordAngle corresponding
 *   // to this Distance (needed to implement Options::max_distance
 *   // efficiently).  For example, if Distance measures WGS84 ellipsoid
 *   // distance then the corresponding angle needs to be 0.56% larger.
 *   S1ChordAngle GetChordAngleBound() const;
 * };
 */
class S2ClosestPointQueryBase(Distance, Data) {
public:
  alias Delta = Distance.Delta;
  alias Index = S2PointIndex!Data;
  alias PointData = Index.PointData;
  alias Options = S2ClosestPointQueryBaseOptions!Distance;

  /**
   * The Target class represents the geometry to which the distance is
   * measured.  For example, there can be subtypes for measuring the distance
   * to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
   * geometry).
   *
   * Implementations do *not* need to be thread-safe.  They may cache data or
   * allocate temporary data structures in order to improve performance.
   */
  alias Target = S2DistanceTarget!Distance;

  /// Each "Result" object represents a closest point.
  struct Result {
  public:
    /// Constructs a Result object for the given point.
    this(Distance distance, PointData point_data) {
      _distance = distance;
      _pointData = point_data;
    }

    /**
     * Returns true if this Result object does not refer to any data point.
     * (The only case where an empty Result is returned is when the
     * FindClosestPoint() method does not find any points that meet the
     * specified criteria.)
     */
    bool isEmpty() const {
      return _pointData is null;
    }

    /// The distance from the target to this point.
    Distance distance() const {
      return _distance;
    }

    /// The point itself.
    const(S2Point) point() const {
      return _pointData.point();
    }

    /// The client-specified data associated with this point.
    const(Data) data() const {
      return _pointData.data();
    }

    /// Returns true if two Result objects are identical.
    bool opEquals(in Result y) const {
      return _distance == y._distance && _pointData == y._pointData;
    }

    /// Compares two Result objects first by distance, then by point_data().
    int opCmp(in Result y) const {
      if (_distance < y._distance) return -1;
      if (_distance > y._distance) return 1;
      if (_pointData < y._pointData) return -1;
      if (_pointData > y._pointData) return 1;
      return 0;
    }

  private:
    Distance _distance = Distance.infinity();
    PointData _pointData = null;
  }

  /**
   * The minimum number of points that a cell must contain to enqueue it
   * rather than processing its contents immediately.
   */
  enum int MIN_POINTS_TO_ENQUEUE = 13;

  /// Default constructor; requires initialize() to be called.
  this() {
    _queue = CellQueue(new QueueEntry[0]);
    _resultSet = BinaryHeap!(Result[])(new Result[0]);
  }

  /// Convenience constructor that calls Init().
  this(Index index) {
    this();
    initialize(index);
  }

  /**
   * Initializes the query.
   * REQUIRES: ReInit() must be called if "index" is modified.
   */
  void initialize(S2PointIndex!Data index) {
    _index = index;
    reInitialize();
  }

  /**
   * Reinitializes the query.  This method must be called whenever the
   * underlying index is modified.
   */
  void reInitialize() {
    _iter.initialize(_index);
    _indexCovering.length = 0;
  }

  /// Return a reference to the underlying S2PointIndex.
  const(Index) index() const {
    return _index;
  }

  /**
   * Returns the closest points to the given target that satisfy the given
   * options.  This method may be called multiple times.
   */
  Result[] findClosestPoints(Target target, Options options) {
    Result[] results;
    findClosestPoints(target, options, results);
    return results;
  }

  /**
   * This version can be more efficient when this method is called many times,
   * since it does not require allocating a new vector on each call.
   */
  void findClosestPoints(Target target, Options options, ref Result[] results) {
    findClosestPointsInternal(target, options);
    results.length = 0;
    if (options.maxPoints() == 1) {
      if (!_resultSingleton.isEmpty()) {
        results ~= _resultSingleton;
      }
    } else if (options.maxPoints() == Options.MAX_MAX_POINTS) {
      sort(_resultVector);
      results = _resultVector.uniq.array;
      _resultVector.length = 0;
    } else {
      results.reserve(_resultSet.length());
      for (; !_resultSet.empty(); _resultSet.popFront()) {
        results ~= _resultSet.front();  // The highest-priority result.
      }
      // The priority queue returns the largest elements first.
      reverse(results);
      enforce(isSorted(results));
    }
  }

  /**
   * Convenience method that returns exactly one point.  If no points satisfy
   * the given search criteria, then a Result with distance() == Infinity()
   * and is_empty() == true is returned.
   *
   * REQUIRES: options.max_points() == 1
   */
  Result findClosestPoint(Target target, Options options)
  in {
    assert(options.maxPoints() == 1);
  } body {
    findClosestPointsInternal(target, options);
    return _resultSingleton;
  }

private:
  alias Iterator = Index.Iterator;

  inout(Options) options() inout {
    return _options;
  }

  void findClosestPointsInternal(Target target, Options options) {
    _target = target;
    _options = options;

    _distanceLimit = options.maxDistance();
    _resultSingleton = Result();
    enforce(_resultVector.empty());
    enforce(_resultSet.empty());
    enforce(target.maxBruteForceIndexSize() >= 0);
    if (_distanceLimit == Distance.zero()) return;

    if (options.maxPoints() == Options.MAX_MAX_POINTS
        && options.maxDistance() == Distance.infinity()
        && options.region() is null) {
      logger.logWarn("Returning all points (max_points/max_distance/region not set)");
    }

    // Note that given point is processed only once (unlike S2ClosestEdgeQuery),
    // and therefore we don't need to worry about the possibility of having
    // duplicate points in the results.
    if (options.maxError() != Delta.zero()) {
      _target.setMaxError(options.maxError());
    }
    if (options.useBruteForce() || _index.numPoints() <= _target.maxBruteForceIndexSize()) {
      findClosestPointsBruteForce();
    } else {
      findClosestPointsOptimized();
    }
  }

  void findClosestPointsBruteForce() {
    for (_iter.begin(); !_iter.done(); _iter.next()) {
      maybeAddResult(_iter.pointData());
    }
  }

  void findClosestPointsOptimized() {
    initQueue();
    while (!_queue.empty()) {
      // We need to copy the top entry before removing it, and we need to remove
      // it before adding any new entries to the queue.
      QueueEntry entry = _queue.front();
      _queue.popFront();
      // Work around weird parse error in gcc 4.9 by using a local variable for
      // entry.distance.
      Distance distance = entry.distance;
      if (!(distance < _distanceLimit)) {
        _queue.clear();  // Clear any remaining entries.
        break;
      }
      S2CellId child = entry.id.childBegin();
      // We already know that it has too many points, so process its children.
      // Each child may either be processed directly or enqueued again.  The
      // loop is optimized so that we don't seek unnecessarily.
      bool seek = true;
      for (int i = 0; i < 4; ++i, child = child.next()) {
        seek = enqueueCell(child, _iter, seek);
      }
    }
  }

  void initQueue()
  in {
    assert(_queue.empty());
  } body {
    // Optimization: rather than starting with the entire index, see if we can
    // limit the search region to a small disc.  Then we can find a covering for
    // that disc and intersect it with the covering for the index.  This can
    // save a lot of work when the search region is small.
    S2Cap cap = _target.getCapBound();
    if (options().maxPoints() == 1) {
      // If the user is searching for just the closest point, we can compute an
      // upper bound on search radius by seeking to the target point in the
      // index and looking at the adjacent index points (in S2CellId order).
      // The minimum distance to either of these points is an upper bound on the
      // search radius.
      //
      // TODO(ericv): The same strategy would also work for small values of
      // max_points() > 1, e.g. max_points() == 20, except that we would need to
      // examine more neighbors (at least 20, and preferably 20 in each
      // direction).  It's not clear whether this is a common case, though, and
      // also this would require extending MaybeAddResult() so that it can
      // remove duplicate entries.  (The points added here may be re-added by
      // EnqueueCell(), but this is okay when max_points() == 1.)
      _iter.seek(S2CellId(cap.center()));
      if (!_iter.done()) {
        maybeAddResult(_iter.pointData());
      }
      if (_iter.prev()) {
        maybeAddResult(_iter.pointData());
      }
      // Skip the rest of the algorithm if we found a matching point.
      if (_distanceLimit == Distance.zero()) return;
    }
    // We start with a covering of the set of indexed points, then intersect it
    // with the given region (if any) and maximum search radius disc (if any).
    if (_indexCovering.empty()) initCovering();
    const(S2CellId)[] initial_cells = _indexCovering;
    if (options().region()) {
      auto coverer = new S2RegionCoverer();
      coverer.mutableOptions().setMaxCells(4);
      coverer.getCovering(options().region(), _regionCovering);
      S2CellUnion.getIntersection(_indexCovering, _regionCovering, _intersectionWithRegion);
      initial_cells = _intersectionWithRegion;
    }
    if (_distanceLimit < Distance.infinity()) {
      auto coverer = new S2RegionCoverer();
      coverer.mutableOptions().setMaxCells(4);
      S1ChordAngle radius = cap.radius() + _distanceLimit.getChordAngleBound();
      auto search_cap = new S2Cap(cap.center(), radius);
      coverer.getFastCovering(search_cap, _maxDistanceCovering);
      S2CellUnion.getIntersection(
          initial_cells, _maxDistanceCovering, _intersectionWithMaxDistance);
      initial_cells = _intersectionWithMaxDistance;
    }
    _iter.begin();
    for (int i = 0; i < initial_cells.length && !_iter.done(); ++i) {
      S2CellId id = initial_cells[i];
      enqueueCell(id, _iter, id.rangeMin() > _iter.id() /*seek*/);
    }
  }

  void initCovering() {
    // Compute the "index covering", which is a small number of S2CellIds that
    // cover the indexed points.  There are two cases:
    //
    //  - If the index spans more than one face, then there is one covering cell
    // per spanned face, just big enough to cover the index cells on that face.
    //
    //  - If the index spans only one face, then we find the smallest cell "C"
    // that covers the index cells on that face (just like the case above).
    // Then for each of the 4 children of "C", if the child contains any index
    // cells then we create a covering cell that is big enough to just fit
    // those index cells (i.e., shrinking the child as much as possible to fit
    // its contents).  This essentially replicates what would happen if we
    // started with "C" as the covering cell, since "C" would immediately be
    // split, except that we take the time to prune the children further since
    // this will save work on every subsequent query.
    _indexCovering.reserve(6);
    _iter.finish();
    if (!_iter.prev()) return;  // Empty index.
    S2CellId index_last_id = _iter.id();
    _iter.begin();
    if (_iter.id() != index_last_id) {
      // The index has at least two cells.  Choose a level such that the entire
      // index can be spanned with at most 6 cells (if the index spans multiple
      // faces) or 4 cells (it the index spans a single face).
      int level = _iter.id().getCommonAncestorLevel(index_last_id) + 1;

      // Visit each potential covering cell except the last (handled below).
      S2CellId last_id = index_last_id.parent(level);
      for (S2CellId id = _iter.id().parent(level); id != last_id; id = id.next()) {
        // Skip any covering cells that don't contain any index cells.
        if (id.rangeMax() < _iter.id()) continue;

        // Find the range of index cells contained by this covering cell and
        // then shrink the cell if necessary so that it just covers them.
        S2CellId cell_first_id = _iter.id();
        _iter.seek(id.rangeMax().next());
        _iter.prev();
        S2CellId cell_last_id = _iter.id();
        _iter.next();
        addInitialRange(cell_first_id, cell_last_id);
      }
    }
    addInitialRange(_iter.id(), index_last_id);
  }

  /**
   * Adds a cell to index_covering_ that covers the given inclusive range.
   *
   * REQUIRES: "first" and "last" have a common ancestor.
   */
  void addInitialRange(S2CellId first_id, S2CellId last_id) {
    // Add the lowest common ancestor of the given range.
    int level = first_id.getCommonAncestorLevel(last_id);
    debug enforce(level >= 0);
    _indexCovering ~= first_id.parent(level);
  }

  void maybeAddResult(PointData point_data) {
    Distance distance = _distanceLimit;
    if (!_target.updateMinDistance(point_data.point(), distance)) return;
    S2Region region = options().region();
    if (region && !region.contains(point_data.point())) return;

    auto result = Result(distance, point_data);
    if (options().maxPoints() == 1) {
      // Optimization for the common case where only the closest point is wanted.
      _resultSingleton = result;
      _distanceLimit = result.distance() - options().maxError();
    } else if (options().maxPoints() == Options.MAX_MAX_POINTS) {
      _resultVector ~= result;  // Sort/unique at end.
    } else {
      // Add this point to result_set_.  Note that with the current algorithm
      // each candidate point is considered at most once (except for one special
      // case where max_points() == 1, see InitQueue for details), so we don't
      // need to worry about possibly adding a duplicate entry here.
      if (_resultSet.length() >= options().maxPoints()) {
        _resultSet.popFront();  // Replace the furthest result point.
      }
      _resultSet.insert(result);
      if (_resultSet.length() >= options().maxPoints()) {
        _distanceLimit = _resultSet.front().distance() - options().maxError();
      }
    }
  }

  /**
   * Either process the contents of the given cell immediately, or add it to the
   * queue to be subdivided.  If "seek" is false, then "iter" must already be
   * positioned at the first indexed point within this cell.
   *
   * Returns "true" if the cell was added to the queue, and "false" if it was
   * processed immediately, in which case "iter" is left positioned at the next
   * cell in S2CellId order.
   */
  bool enqueueCell(S2CellId id, ref Iterator iter, bool seek) {
    if (seek) iter.seek(id.rangeMin());
    if (id.isLeaf()) {
      // Leaf cells can't be subdivided.
      for (; !iter.done() && iter.id() == id; iter.next()) {
        maybeAddResult(iter.pointData());
      }
      return false;  // No need to seek to next child.
    }
    S2CellId last = id.rangeMax();
    int num_points = 0;
    for (; !iter.done() && iter.id() <= last; iter.next()) {
      if (num_points == MIN_POINTS_TO_ENQUEUE - 1) {
        // This cell has too many points (including this one), so enqueue it.
        auto cell = new S2Cell(id);
        Distance distance = _distanceLimit;
        // We check "region_" second because it may be relatively expensive.
        if (_target.updateMinDistance(cell, distance)
            && (!options().region() || options().region().mayIntersect(cell))) {
          _queue.insert(QueueEntry(distance, id));
        }
        return true;  // Seek to next child.
      }
      _tmpPointData[num_points++] = iter.pointData();
    }
    // There were few enough points that we might as well process them now.
    for (int i = 0; i < num_points; ++i) {
      maybeAddResult(_tmpPointData[i]);
    }
    return false;  // No need to seek to next child.
  }

  /// The maximum number of points to process without subdividing further.
  enum int MAX_LEAF_POINTS = 12;

  Index _index;
  Options _options;
  Target _target;

  /**
   * For the optimized algorihm we precompute the top-level S2CellIds that
   * will be added to the priority queue.  There can be at most 6 of these
   * cells.  Essentially this is just a covering of the indexed points.
   */
  S2CellId[] _indexCovering;

  /**
   * The distance beyond which we can safely ignore further candidate points.
   * (Candidates that are exactly at the limit are ignored; this is more
   * efficient for UpdateMinDistance() and should not affect clients since
   * distance measurements have a small amount of error anyway.)
   *
   * Initially this is the same as the maximum distance specified by the user,
   * but it can also be updated by the algorithm (see MaybeAddResult).
   */
  Distance _distanceLimit;

  /**
   * The current result set is stored in one of three ways:
   *
   *  - If max_points() == 1, the best result is kept in result_singleton_.
   *
   *  - If max_points() == "infinity", results are appended to result_vector_
   *    and sorted/uniqued at the end.
   *
   *  - Otherwise results are kept in a priority queue so that we can
   *    progressively reduce the distance limit once max_points() results have
   *    been found.
   */
  Result _resultSingleton;
  Result[] _resultVector;

  /// Used as a priority queue for the results.
  BinaryHeap!(Result[]) _resultSet;

  /**
   * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
   * in increasing order of distance from the target point.
   */
  struct QueueEntry {
    /// A lower bound on the distance from the target point to any point
    /// within "id".  This is the key of the priority queue.
    Distance distance;

    /// The cell being queued.
    S2CellId id;

    /// The priority queue returns the largest elements first, so we want the
    /// "largest" entry to have the smallest distance.
    int opCmp(in QueueEntry other) const {
      if (distance > other.distance)
        return -1;
      else if (distance < other.distance)
        return 1;
      return 0;
    }
  }

  alias CellQueue = BinaryHeap!(QueueEntry[]);
  CellQueue _queue;

  // Temporaries, defined here to avoid multiple allocations / initializations.
  Iterator _iter;
  S2CellId[] _regionCovering;
  S2CellId[] _maxDistanceCovering;
  S2CellId[] _intersectionWithRegion;
  S2CellId[] _intersectionWithMaxDistance;
  PointData[MIN_POINTS_TO_ENQUEUE - 1] _tmpPointData;
}
