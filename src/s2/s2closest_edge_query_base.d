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

module s2.s2closest_edge_query_base;

import s2.logger;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2distance_target;
import s2.s2point;
import s2.s2region_coverer;
import s2.s2shape;
import s2.s2shape_index;
import s2.shapeutil.count_edges : countEdgesUpTo;
import s2.shapeutil.shape_edge_id;
import s2.util.container.btree;

import std.algorithm;
import std.array;
import std.container.binaryheap;
import std.exception : enforce;
import std.range : retro;


// S2ClosestEdgeQueryBase is a templatized class for finding the closest
// edge(s) between two geometries.  It is not intended to be used directly,
// but rather to serve as the implementation of various specialized classes
// with more convenient APIs (such as S2ClosestEdgeQuery).  It is flexible
// enough so that it can be adapted to compute maximum distances and even
// potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between two geometries A and B.
//  - Find all edges of geometry A that are within a distance D of geometry B.
//  - Find the k edges of geometry A that are closest to a given point P.
//
// You can also specify whether polygons should include their interiors (i.e.,
// if a point is contained by a polygon, should the distance be zero or should
// it be measured to the polygon boundary?)
//
// The input geometries may consist of any number of points, polylines, and
// polygons (collectively referred to as "shapes").  Shapes do not need to be
// disjoint; they may overlap or intersect arbitrarily.  The implementation is
// designed to be fast for both simple and complex geometries.
//
// The Distance template argument is used to represent distances.  Usually
// this type is a thin wrapper around S1ChordAngle, but another distance type
// may be substituted as long as it implements the API below.  This can be
// used to change the comparison function (e.g., to find the furthest edges
// from the target), or to get more accuracy if desired.
//
// The Distance concept is as follows:
//
// class Distance {
//  public:
//   // Default and copy constructors, assignment operator:
//   Distance();
//   Distance(const Distance&);
//   Distance& operator=(const Distance&);
//
//   // Factory methods:
//   static Distance Zero();      // Returns a zero distance.
//   static Distance Infinity();  // Larger than any valid distance.
//   static Distance Negative();  // Smaller than any valid distance.
//
//   // Comparison operators:
//   friend bool operator==(Distance x, Distance y);
//   friend bool operator<(Distance x, Distance y);
//
//   // Delta represents the positive difference between two distances.
//   // It is used together with operator-() to implement Options::max_error().
//   // Typically Distance::Delta is simply S1ChordAngle.
//   class Delta {
//    public:
//     Delta();
//     Delta(const Delta&);
//     Delta& operator=(const Delta&);
//     friend bool operator==(Delta x, Delta y);
//     static Delta Zero();
//   };
//
//   // Subtraction operator.  Note that the second argument represents a
//   // delta between two distances.  This distinction is important for
//   // classes that compute maximum distances (e.g., S2FurthestEdgeQuery).
//   friend Distance operator-(Distance x, Delta delta);
//
//   // Method that returns an upper bound on the S1ChordAngle corresponding
//   // to this Distance (needed to implement Options::max_distance
//   // efficiently).  For example, if Distance measures WGS84 ellipsoid
//   // distance then the corresponding angle needs to be 0.56% larger.
//   S1ChordAngle GetChordAngleBound() const;
// };
class S2ClosestEdgeQueryBase(DistanceT) {
public:
  alias Delta = DistanceT.Delta;

  // Options that control the set of edges returned.  Note that by default
  // *all* edges are returned, so you will always want to set either the
  // max_edges() option or the max_distance() option (or both).
  static class Options {
  public:
    this() { }

    this(Options options) {
      _maxDistance = options._maxDistance;
      _maxError = options._maxError;
      _maxEdges = options._maxEdges;
      _includeInteriors = options._includeInteriors;
      _useBruteForce = options._useBruteForce;
    }

    // Specifies that at most "max_edges" edges should be returned.
    //
    // REQUIRES: max_edges >= 1
    // DEFAULT: numeric_limits<int>::max()
    int maxEdges() const {
      return _maxEdges;
    }

    void setMaxEdges(int max_edges)
    in {
      assert(max_edges >= 1);
    } body {
      _maxEdges = max_edges;
    }

    enum int MAX_MAX_EDGES = int.max;

    // Specifies that only edges whose distance to the target is less than
    // "max_distance" should be returned.
    //
    // Note that edges whose distance is exactly equal to "max_distance" are
    // not returned.  In most cases this doesn't matter (since distances are
    // not computed exactly in the first place), but if such edges are needed
    // then you can retrieve them by specifying "max_distance" as the next
    // largest representable DistanceT.  For example, if DistanceT is an
    // S1ChordAngle then you can specify max_distance.Successor().
    //
    // DEFAULT: DistanceT::Infinity()
    DistanceT maxDistance() const {
      return _maxDistance;
    }

    void setMaxDistance(DistanceT max_distance) {
      _maxDistance = max_distance;
    }

    // Specifies that edges up to max_error() further away than the true
    // closest edges may be substituted in the result set, as long as such
    // edges satisfy all the remaining search criteria (such as max_distance).
    // This option only has an effect if max_edges() is also specified;
    // otherwise all edges closer than max_distance() will always be returned.
    //
    // Note that this does not affect how the distance between edges is
    // computed; it simply gives the algorithm permission to stop the search
    // early as soon as the best possible improvement drops below max_error().
    //
    // This can be used to implement distance predicates efficiently.  For
    // example, to determine whether the minimum distance is less than D, set
    // max_edges() == 1 and max_distance() == max_error() == D.  This causes
    // the algorithm to terminate as soon as it finds any edge whose distance
    // is less than D, rather than continuing to search for an edge that is
    // even closer.
    //
    // DEFAULT: DistanceT::Delta::Zero()
    Delta maxError() const {
      return _maxError;
    }

    void setMaxError(Delta max_error) {
      _maxError = max_error;
    }

    // Specifies that polygon interiors should be included when measuring
    // distances.  In other words, polygons that contain the target should
    // have a distance of zero.  (For targets consisting of multiple connected
    // components, the distance is zero if any component is contained.)  This
    // is indicated in the results by returning a (shape_id, edge_id) pair
    // with edge_id == -1, i.e. this value denotes the polygons's interior.
    //
    // Note that for efficiency, any polygon that intersects the target may or
    // may not have an (edge_id == -1) result.  Such results are optional
    // because in that case the distance to the polygon is already zero.
    //
    // DEFAULT: true
    bool includeInteriors() const {
      return _includeInteriors;
    }

    void setIncludeInteriors(bool include_interiors) {
      _includeInteriors = include_interiors;
    }

    // Specifies that distances should be computed by examining every edge
    // rather than using the S2ShapeIndex.  This is useful for testing,
    // benchmarking, and debugging.
    //
    // DEFAULT: false
    bool useBruteForce() const {
      return _useBruteForce;
    }

    void setUseBruteForce(bool use_brute_force) {
      _useBruteForce = use_brute_force;
    }

  private:
    DistanceT _maxDistance = DistanceT.infinity;
    Delta _maxError = Delta.zero;
    int _maxEdges = MAX_MAX_EDGES;
    bool _includeInteriors = true;
    bool _useBruteForce = false;
  }

  // The Target class represents the geometry to which the distance is
  // measured.  For example, there can be subtypes for measuring the distance
  // to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
  // geometry).
  //
  // Implementations do *not* need to be thread-safe.  They may cache data or
  // allocate temporary data structures in order to improve performance.
  alias Target = S2DistanceTarget!DistanceT;

  // Each "Result" object represents a closest edge.  Note the following
  // special cases:
  //
  //  - (shape_id >= 0) && (edge_id < 0) represents the interior of a shape.
  //    Such results may be returned when options.include_interiors() is true.
  //
  //  - (shape_id < 0) && (edge_id < 0) is returned by `FindClosestEdge` to
  //    indicate that no edge satisfies the requested query options.
  //
  // TODO(ericv): Convert to a class with accessor methods.
  struct Result {
    DistanceT distance = DistanceT.infinity;  // The distance from the target to this edge.
    int shapeId = -1;     // Identifies an indexed shape.
    int edgeId = -1;      // Identifies an edge within the shape.

    // Compares edges first by distance, then by (shape_id, edge_id).
    int opCmp(in Result o) const {
      if (distance != o.distance) return distance.opCmp(o.distance);
      if (shapeId != o.shapeId) return shapeId - o.shapeId;
      if (edgeId != o.edgeId) return edgeId - o.edgeId;
      return 0;
    }
  }

  // Default constructor; requires Init() to be called.
  this() {
    _resultSet = new BTree!Result();
    _queue = BinaryHeap!(QueueEntry[])([]);
  }

  // Convenience constructor that calls Init().
  this(S2ShapeIndex index) {
    initialize(index);
  }

  // Initializes the query.
  // REQUIRES: ReInit() must be called if "index" is modified.
  void initialize(S2ShapeIndex index) {
    _index = index;
    reInitialize();
  }

  // Reinitializes the query.  This method must be called whenever the
  // underlying index is modified.
  void reInitialize() {
    _indexNumEdges = 0;
    _indexNumEdgesLimit = 0;
    _indexCovering.length = 0;
    _indexCells.length = 0;
    // We don't initialize iter_ here to make queries on small indexes a bit
    // faster (i.e., where brute force is used).
  }

  // Returns a reference to the underlying S2ShapeIndex.
  const(S2ShapeIndex) index() const {
    return _index;
  }

  // Returns the closest edges to the given target that satisfy the given
  // options.  This method may be called multiple times.
  //
  // Note that if options().include_interiors() is true, the result vector may
  // include some entries with edge_id == -1.  This indicates that the target
  // intersects the indexed polygon with the given shape_id.
  Result[] findClosestEdges(Target target, Options options) {
    import std.stdio;
    Result[] results;
    findClosestEdges(target, options, results);
    return results;
  }

  // This version can be more efficient when this method is called many times,
  // since it does not require allocating a new vector on each call.
  void findClosestEdges(Target target, Options options, out Result[] results) {
    findClosestEdgesInternal(target, options);
    if (options.maxEdges() == 1) {
      if (_resultSingleton.shapeId >= 0) {
        results ~= _resultSingleton;
      }
    } else if (options.maxEdges() == Options.MAX_MAX_EDGES) {
      results ~= _resultVector;
      _resultVector.length = 0;
    } else {
      results = _resultSet[].array;
      _resultSet.clear();
    }
  }

  // Convenience method that returns exactly one edge.  If no edges satisfy
  // the given search criteria, then a Result with distance == Infinity() and
  // shape_id == edge_id == -1 is returned.
  //
  // Note that if options.include_interiors() is true, edge_id == -1 is also
  // used to indicate that the target intersects an indexed polygon (but in
  // that case distance == Zero() and shape_id >= 0).
  //
  // REQUIRES: options.max_edges() == 1
  Result findClosestEdge(Target target, Options options)
  in {
    assert(options.maxEdges() == 1);
  } body {
    findClosestEdgesInternal(target, options);
    return _resultSingleton;
  }

 private:

  const(Options) options() const {
    return _options;
  }

  void findClosestEdgesInternal(Target target, Options options)
  in {
    assert(_resultVector.empty());
    assert(_resultSet.empty());
    assert(target.maxBruteForceIndexSize() >= 0);
  } body {
    _target = target;
    _options = options;

    _testedEdges.clear();
    _distanceLimit = options.maxDistance();
    _resultSingleton = Result();

    if (_distanceLimit == DistanceT.zero()) return;

    if (options.maxEdges() == Options.MAX_MAX_EDGES
        && options.maxDistance() == DistanceT.infinity()) {
      logger.logWarn("Returning all edges (max_edges/max_distance not set)");
    }

    if (options.includeInteriors()) {
      auto shape_ids = new BTree!int();
      target.visitContainingShapes(
          _index,
          (in S2Shape containing_shape, in S2Point target_point) {
            shape_ids.insert(containing_shape.id());
            return shape_ids.length < options.maxEdges();
          });
      foreach (int shape_id; shape_ids) {
        addResult(Result(DistanceT.zero(), shape_id, -1));
      }
      if (_distanceLimit == DistanceT.zero()) return;
    }

    // If max_error() > 0 and the target takes advantage of this, then we may
    // need to adjust the distance estimates to S2ShapeIndex cells to ensure
    // that they are always a lower bound on the true distance.  For example,
    // suppose max_distance == 100, max_error == 30, and we compute the distance
    // to the target from some cell C0 as d(C0) == 80.  Then because the target
    // takes advantage of max_error(), the true distance could be as low as 50.
    // In order not to miss edges contained by such cells, we need to subtract
    // max_error() from the distance estimates.  This behavior is controlled by
    // the use_conservative_cell_distance_ flag.
    //
    // However there is one important case where this adjustment is not
    // necessary, namely when max_distance() < max_error().  This is because
    // max_error() only affects the algorithm once at least max_edges() edges
    // have been found that satisfy the given distance limit.  At that point,
    // max_error() is subtracted from distance_limit_ in order to ensure that
    // any further matches are closer by at least that amount.  But when
    // max_distance() < max_error(), this reduces the distance limit to 0,
    // i.e. all remaining candidate cells and edges can safely be discarded.
    // (Note that this is how IsDistanceLess() and friends are implemented.)
    //
    // Note that DistanceT::Delta only supports operator==.
    bool target_uses_max_error = (!(options.maxError() == Delta.zero())
        && _target.setMaxError(options.maxError()));

    // Note that we can't compare max_error() and distance_limit_ directly
    // because one is a Delta and one is a DistanceT.  Instead we subtract them.
    _useConservativeCellDistance = target_uses_max_error
        && (_distanceLimit == DistanceT.infinity()
            || DistanceT.zero() < _distanceLimit - options.maxError());

    // Use the brute force algorithm if the index is small enough.  To avoid
    // spending too much time counting edges when there are many shapes, we stop
    // counting once there are too many edges.  We may need to recount the edges
    // if we later see a target with a larger brute force edge threshold.
    int min_optimized_edges = _target.maxBruteForceIndexSize() + 1;
    if (min_optimized_edges > _indexNumEdgesLimit && _indexNumEdges >= _indexNumEdgesLimit) {
      _indexNumEdges = countEdgesUpTo(_index, min_optimized_edges);
      _indexNumEdgesLimit = min_optimized_edges;
    }

    if (options.useBruteForce() || _indexNumEdges < min_optimized_edges) {
      // The brute force algorithm considers each edge exactly once.
      _avoidDuplicates = false;
      findClosestEdgesBruteForce();
    } else {
      // If the target takes advantage of max_error() then we need to avoid
      // duplicate edges explicitly.  (Otherwise it happens automatically.)
      _avoidDuplicates = (target_uses_max_error && options.maxEdges() > 1);
      findClosestEdgesOptimized();
    }
  }

  void findClosestEdgesBruteForce() {
    int num_shape_ids = _index.numShapeIds();
    for (int id = 0; id < num_shape_ids; ++id) {
      const(S2Shape) shape = _index.shape(id);
      if (shape is null) continue;
      int num_edges = shape.numEdges();
      for (int e = 0; e < num_edges; ++e) {
        maybeAddResult(shape, e);
      }
    }
  }

  void findClosestEdgesOptimized() {
    initQueue();
    // Repeatedly find the closest S2Cell to "target" and either split it into
    // its four children or process all of its edges.
    while (!_queue.empty()) {
      // We need to copy the top entry before removing it, and we need to
      // remove it before adding any new entries to the queue.
      QueueEntry entry = _queue.front();
      _queue.popFront();
      // Work around weird parse error in gcc 4.9 by using a local variable for
      // entry.distance.
      DistanceT distance = entry.distance;
      if (!(distance < _distanceLimit)) {
        _queue.clear();  // Clear any remaining entries.
        break;
      }
      // If this is already known to be an index cell, just process it.
      if (entry.indexCell !is null) {
        processEdges(entry);
        continue;
      }
      // Otherwise split the cell into its four children.  Before adding a
      // child back to the queue, we first check whether it is empty.  We do
      // this in two seek operations rather than four by seeking to the key
      // between children 0 and 1 and to the key between children 2 and 3.
      S2CellId id = entry.id;
      _iter.seek(id.child(1).rangeMin());
      if (!_iter.done() && _iter.id() <= id.child(1).rangeMax()) {
        enqueueCurrentCell(id.child(1));
      }
      if (_iter.prev() && _iter.id() >= id.rangeMin()) {
        enqueueCurrentCell(id.child(0));
      }
      _iter.seek(id.child(3).rangeMin());
      if (!_iter.done() && _iter.id() <= id.rangeMax()) {
        enqueueCurrentCell(id.child(3));
      }
      if (_iter.prev() && _iter.id() >= id.child(2).rangeMin()) {
        enqueueCurrentCell(id.child(2));
      }
    }
  }

  void initQueue()
  in {
    assert(_queue.empty());
  } body {
    if (_indexCovering.empty()) {
      // We delay iterator initialization until now to make queries on very
      // small indexes a bit faster (i.e., where brute force is used).
      _iter.initialize(_index, S2ShapeIndex.InitialPosition.UNPOSITIONED);
    }

    // Optimization: if the user is searching for just the closest edge, and the
    // target happens to intersect an index cell, then we try to limit the search
    // region to a small disc by first processing the edges in that cell.  This
    // sets distance_limit_ based on the closest edge in that cell, which we can
    // then use to limit the search area.  This means that the cell containing
    // "target" will be processed twice, but in general this is still faster.
    S2Cap cap = _target.getCapBound();
    if (options().maxEdges() == 1 && _iter.locate(cap.center())) {
      processEdges(new QueueEntry(DistanceT.zero(), _iter.id(), _iter.cell()));
      // Skip the rest of the algorithm if we found an intersecting edge.
      if (_distanceLimit == DistanceT.zero()) return;
    }
    if (_indexCovering.empty()) initCovering();
    if (_distanceLimit == DistanceT.infinity()) {
      // Start with the precomputed index covering.
      for (int i = 0; i < _indexCovering.length; ++i) {
        enqueueCell(_indexCovering[i], _indexCells[i]);
      }
    } else {
      // Compute a covering of the search disc and intersect it with the
      // precomputed index covering.
      auto coverer = new S2RegionCoverer();
      coverer.mutableOptions().setMaxCells(4);
      S1ChordAngle radius = cap.radius() + _distanceLimit.getChordAngleBound();
      auto search_cap = new S2Cap(cap.center(), radius);
      coverer.getFastCovering(search_cap, _maxDistanceCovering);
      S2CellUnion.getIntersection(_indexCovering, _maxDistanceCovering, _initialCells);

      // Now we need to clean up the initial cells to ensure that they all
      // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
      // the index at all, while other may be descendants of an index cell.)
      for (int i = 0, j = 0; i < _initialCells.length; ) {
        S2CellId id_i = _initialCells[i];
        // Find the top-level cell that contains this initial cell.
        while (_indexCovering[j].rangeMax() < id_i) ++j;
        S2CellId id_j = _indexCovering[j];
        if (id_i == id_j) {
          // This initial cell is one of the top-level cells.  Use the
          // precomputed S2ShapeIndexCell pointer to avoid an index seek.
          enqueueCell(id_j, _indexCells[j]);
          ++i, ++j;
        } else {
          // This initial cell is a proper descendant of a top-level cell.
          // Check how it is related to the cells of the S2ShapeIndex.
          S2ShapeIndex.CellRelation r = _iter.locate(id_i);
          if (r == S2ShapeIndex.CellRelation.INDEXED) {
            // This cell is a descendant of an index cell.  Enqueue it and skip
            // any other initial cells that are also descendants of this cell.
            enqueueCell(_iter.id(), _iter.cell());
            const S2CellId last_id = _iter.id().rangeMax();
            while (++i < _initialCells.length && _initialCells[i] <= last_id)
              continue;
          } else {
            // Enqueue the cell only if it contains at least one index cell.
            if (r == S2ShapeIndex.CellRelation.SUBDIVIDED) enqueueCell(id_i, null);
            ++i;
          }
        }
      }
    }
  }

  void initCovering() {
    // Find the range of S2Cells spanned by the index and choose a level such
    // that the entire index can be covered with just a few cells.  These are
    // the "top-level" cells.  There are two cases:
    //
    //  - If the index spans more than one face, then there is one top-level cell
    // per spanned face, just big enough to cover the index cells on that face.
    //
    //  - If the index spans only one face, then we find the smallest cell "C"
    // that covers the index cells on that face (just like the case above).
    // Then for each of the 4 children of "C", if the child contains any index
    // cells then we create a top-level cell that is big enough to just fit
    // those index cells (i.e., shrinking the child as much as possible to fit
    // its contents).  This essentially replicates what would happen if we
    // started with "C" as the top-level cell, since "C" would immediately be
    // split, except that we take the time to prune the children further since
    // this will save work on every subsequent query.

    // Don't need to reserve index_cells_ since it is an InlinedVector.
    _indexCovering.reserve(6);

    // TODO(ericv): Use a single iterator (iter_) below and save position
    // information using pair<S2CellId, const S2ShapeIndexCell*> type.
    auto next = new S2ShapeIndex.Iterator(_index, S2ShapeIndex.InitialPosition.BEGIN);
    auto last = new S2ShapeIndex.Iterator(_index, S2ShapeIndex.InitialPosition.END);
    last.prev();
    if (next.id() != last.id()) {
      // The index has at least two cells.  Choose a level such that the entire
      // index can be spanned with at most 6 cells (if the index spans multiple
      // faces) or 4 cells (it the index spans a single face).
      int level = next.id().getCommonAncestorLevel(last.id()) + 1;

      // Visit each potential top-level cell except the last (handled below).
      S2CellId last_id = last.id().parent(level);
      for (S2CellId id = next.id().parent(level); id != last_id; id = id.next()) {
        // Skip any top-level cells that don't contain any index cells.
        if (id.rangeMax() < next.id()) continue;

        // Find the range of index cells contained by this top-level cell and
        // then shrink the cell if necessary so that it just covers them.
        S2ShapeIndex.Iterator cell_first = next;
        next.seek(id.rangeMax().next());
        S2ShapeIndex.Iterator cell_last = next;
        cell_last.prev();
        addInitialRange(cell_first, cell_last);
      }
    }
    addInitialRange(next, last);
  }

  /**
   * Add an entry to index_covering_ and index_cells_ that covers the given
   * inclusive range of cells.
   *
   * REQUIRES: "first" and "last" have a common ancestor.
   */
  void addInitialRange(S2ShapeIndex.Iterator first, in S2ShapeIndex.Iterator last) {
    if (first.id() == last.id()) {
      // The range consists of a single index cell.
      _indexCovering ~= first.id();
      _indexCells ~= first.cell();
    } else {
      // Add the lowest common ancestor of the given range.
      int level = first.id().getCommonAncestorLevel(last.id());
      enforce(level >= 0);
      _indexCovering ~= first.id().parent(level);
      _indexCells ~= null;
    }
  }

  void maybeAddResult(in S2Shape shape, int edge_id) {
    auto testEdge = ShapeEdgeId(shape.id(), edge_id);
    if (_avoidDuplicates && testEdge !in _testedEdges) {
      return;
    }
    _testedEdges[testEdge] = true;
    auto edge = shape.edge(edge_id);
    DistanceT distance = _distanceLimit;
    if (_target.updateMinDistance(edge.v0, edge.v1, distance)) {
      addResult(Result(distance, shape.id(), edge_id));
    }
  }

  void addResult(in Result result) {
    if (options().maxEdges() == 1) {
      // Optimization for the common case where only the closest edge is wanted.
      _resultSingleton = result;
      _distanceLimit = result.distance - options().maxError();
    } else if (options().maxEdges() == Options.MAX_MAX_EDGES) {
      _resultVector ~= result;  // Sort/unique at end.
    } else {
      // Add this edge to result_set_.  Note that even if we already have enough
      // edges, we can't erase an element before insertion because the "new"
      // edge might in fact be a duplicate.
      _resultSet.insert(result);
      int size = cast(int) _resultSet.length;
      if (size >= options().maxEdges()) {
        if (size > options().maxEdges()) {
          _resultSet.remove(_resultSet.end().getValue());
        }
        _distanceLimit = _resultSet.end().getValue().distance - options().maxError();
      }
    }
  }

  // Process all the edges of the given index cell.
  void processEdges(in QueueEntry entry) {
    const(S2ShapeIndexCell) index_cell = entry.indexCell;
    for (int s = 0; s < index_cell.numClipped(); ++s) {
      const(S2ClippedShape) clipped = index_cell.clipped(s);
      S2Shape shape = _index.shape(clipped.shapeId());
      for (int j = 0; j < clipped.numEdges(); ++j) {
        maybeAddResult(shape, clipped.edge(j));
      }
    }
  }

  // Add the given cell id to the queue.  "index_cell" is the corresponding
  // S2ShapeIndexCell, or nullptr if "id" is not an index cell.
  void enqueueCell(in S2CellId id, in S2ShapeIndexCell index_cell) {
    if (index_cell) {
      // If this index cell has only a few edges, then it is faster to check
      // them directly rather than computing the minimum distance to the S2Cell
      // and inserting it into the queue.
      enum int kMinEdgesToEnqueue = 10;
      int num_edges = countEdges(index_cell);
      if (num_edges == 0) return;
      if (num_edges < kMinEdgesToEnqueue) {
        // Set "distance" to zero to avoid the expense of computing it.
        processEdges(new QueueEntry(DistanceT.zero(), id, index_cell));
        return;
      }
    }
    // Otherwise compute the minimum distance to any point in the cell and add
    // it to the priority queue.
    auto cell = new S2Cell(id);
    DistanceT distance = _distanceLimit;
    if (!_target.updateMinDistance(cell, distance)) return;
    if (_useConservativeCellDistance) {
      // Ensure that "distance" is a lower bound on the true distance to the cell.
      distance = distance - options().maxError();  // operator-=() not defined.
    }
    _queue.insert(new QueueEntry(distance, id, index_cell));
  }

  // Enqueue the given cell id.
  // REQUIRES: iter_ is positioned at a cell contained by "id".
  void enqueueCurrentCell(S2CellId id)
  in {
    assert(id.contains(_iter.id()));
  } body {
    if (_iter.id() == id) {
      enqueueCell(id, _iter.cell());
    } else {
      enqueueCell(id, null);
    }
  }

  // Return the number of edges in the given index cell.
  static int countEdges(in S2ShapeIndexCell cell) {
    int count = 0;
    for (int s = 0; s < cell.numClipped(); ++s) {
      count += cell.clipped(s).numEdges();
    }
    return count;
  }

  S2ShapeIndex _index;
  Options _options;
  Target _target;

  // True if max_error() must be subtracted from S2ShapeIndex cell distances
  // in order to ensure that such distances are measured conservatively.  This
  // is true only if the target takes advantage of max_error() in order to
  // return faster results, and 0 < max_error() < distance_limit_.
  bool _useConservativeCellDistance;

  // For the optimized algorihm we precompute the top-level S2CellIds that
  // will be added to the priority queue.  There can be at most 6 of these
  // cells.  Essentially this is just a covering of the indexed edges, except
  // that we also store pointers to the corresponding S2ShapeIndexCells to
  // reduce the number of index seeks required.
  //
  // The covering needs to be stored in a std::vector so that we can use
  // S2CellUnion::GetIntersection().
  S2CellId[] _indexCovering;
  S2ShapeIndexCell[] _indexCells;

  // The decision about whether to use the brute force algorithm is based on
  // counting the total number of edges in the index.  However if the index
  // contains a large number of shapes, this in itself might take too long.
  // So instead we only count edges up to (max_brute_force_index_size() + 1)
  // for the current target type (stored as index_num_edges_limit_).
  int _indexNumEdges;
  int _indexNumEdgesLimit;

  // The distance beyond which we can safely ignore further candidate edges.
  // (Candidates that are exactly at the limit are ignored; this is more
  // efficient for UpdateMinDistance() and should not affect clients since
  // distance measurements have a small amount of error anyway.)
  //
  // Initially this is the same as the maximum distance specified by the user,
  // but it can also be updated by the algorithm (see MaybeAddResult).
  DistanceT _distanceLimit;

  // The current result set is stored in one of three ways:
  //
  //  - If max_edges() == 1, the best result is kept in result_singleton_.
  //
  //  - If max_edges() == "infinity", results are appended to result_vector_
  //    and sorted/uniqued at the end.
  //
  //  - Otherwise results are kept in a btree_set so that we can progressively
  //    reduce the distance limit once max_edges() results have been found.
  //    (A priority queue is not sufficient because we need to be able to
  //    check whether a candidate edge is already in the result set.)
  Result _resultSingleton;
  Result[] _resultVector;
  BTree!Result _resultSet;

  // When the result edges are stored in a btree_set (see above), usually
  // duplicates can be removed simply by inserting candidate edges in the
  // current set.  However this is not true if Options::max_error() > 0 and
  // the Target subtype takes advantage of this by returning suboptimal
  // distances.  This is because when UpdateMinDistance() is called with
  // different "min_dist" parameters (i.e., the distance to beat), the
  // implementation may return a different distance for the same edge.  Since
  // the btree_set is keyed by (distance, shape_id, edge_id) this can create
  // duplicate edges in the results.
  //
  // The flag below is true when duplicates must be avoided explicitly.  This
  // is achieved by maintaining a separate set keyed by (shape_id, edge_id)
  // only, and checking whether each edge is in that set before computing the
  // distance to it.
  bool _avoidDuplicates;
  bool[ShapeEdgeId] _testedEdges;

  // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
  // in increasing order of distance from the target point.
  class QueueEntry {
    this(in DistanceT distance, in S2CellId id, in S2ShapeIndexCell indexCell) {
      this.distance = distance;
      this.id = id;
      this.indexCell = indexCell;
    }

    // A lower bound on the distance from the target point to any edge point
    // within "id".  This is the key of the priority queue.
    const(DistanceT) distance;

    // The cell being queued.
    const(S2CellId) id;

    // If "id" belongs to the index, this field stores the corresponding
    // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
    // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
    // field is to avoid an extra Seek() when the queue entry is processed.
    const(S2ShapeIndexCell) indexCell;

    int opCmp(in QueueEntry other) const {
      // The priority queue returns the largest elements first, so we want the
      // "largest" entry to have the smallest distance.
      return distance.opCmp(other.distance);
    }
  }

  alias CellQueue = BinaryHeap!(QueueEntry[]);
  CellQueue _queue;

  // Temporaries, defined here to avoid multiple allocations / initializations.

  S2ShapeIndex.Iterator _iter;
  S2CellId[] _maxDistanceCovering;
  S2CellId[] _initialCells;
}
