// Copyright 2012 Google Inc. All Rights Reserved.
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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.mutable_s2shape_index;

import s2.base.spinlock : SpinLock;
import s2.logger;
import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s2cell_id;
import s2.s2edge_clipping : getClippedEdgeBound, FACE_CLIP_ERROR_UV_COORD, EDGE_CLIP_ERROR_UV_COORD;
import s2.s2edge_crosser;
import s2.s2padded_cell;
import s2.s2point;
import s2.s2pointutil;
import s2.s2shape;
import s2.s2shape_index : S2ShapeIndex, S2ShapeIndexCell, S2ClippedShape;
import s2.s2shapeutil_contains_brute_force : containsBruteForce;
import s2.util.container.btree_map;
import s2coords = s2.s2coords;

import core.atomic : atomicLoad, atomicStore, MemoryOrder;
import core.sync.mutex : Mutex;
import std.algorithm : swap;
import std.array;
import std.exception : enforce;

// The default maximum number of edges per cell (not counting "long" edges).
// If a cell has more than this many edges, and it is not a leaf cell, then it
// is subdivided.  This flag can be overridden via MutableS2ShapeIndex::Options.
// Reasonable values range from 10 to about 50 or so.
enum int S2SHAPE_INDEX_DEFAULT_MAX_EDGES_PER_CELL = 10;

// FLAGS_s2shape_index_tmp_memory_budget_mb
//
// Attempt to limit the amount of temporary memory allocated while building or
// updating a MutableS2ShapeIndex to at most this value.  This is achieved by
// splitting the updates into multiple batches when necessary.  (The memory
// required is proportional to the number of edges being updated at once.)
//
// Note that this limit is not a hard guarantee, for several reasons:
//  (1) the memory estimates are only approximations;
//  (2) all edges in a given shape are added or removed at once, so shapes
//      with huge numbers of edges may exceed the budget;
//  (3) shapes being removed are always processed in a single batch.  (This
//      could be fixed, but it seems better to keep the code simpler for now.)
enum int S2SHAPE_INDEX_TMP_MEMORY_BUDGET_MB = 100;

// FLAGS_s2shape_index_cell_size_to_long_edge_ratio
//
// The cell size relative to the length of an edge at which it is first
// considered to be "long".  Long edges do not contribute toward the decision
// to subdivide a cell further.  For example, a value of 2.0 means that the
// cell must be at least twice the size of the edge in order for that edge to
// be counted.  There are two reasons for not counting long edges: (1) such
// edges typically need to be propagated to several children, which increases
// time and memory costs without much benefit, and (2) in pathological cases,
// many long edges close together could force subdivision to continue all the
// way to the leaf cell level.
enum double S2SHAPE_INDEX_CELL_SIZE_TO_LONG_EDGE_RATIO = 1.0;


// MutableS2ShapeIndex is a class for in-memory indexing of polygonal geometry.
// The objects in the index are known as "shapes", and may consist of points,
// polylines, and/or polygons, possibly overlapping.  The index makes it very
// fast to answer queries such as finding nearby shapes, measuring distances,
// testing for intersection and containment, etc.
//
// MutableS2ShapeIndex allows not only building an index, but also updating it
// incrementally by adding or removing shapes (hence its name).  It is one of
// several implementations of the S2ShapeIndex interface.  MutableS2ShapeIndex
// is designed to be compact; usually it is smaller than the underlying
// geometry being indexed.  It is capable of indexing up to hundreds of
// millions of edges.  The index is also fast to construct.
//
// There are a number of built-in classes that work with S2ShapeIndex objects.
// Generally these classes accept any collection of geometry that can be
// represented by an S2ShapeIndex, i.e. any combination of points, polylines,
// and polygons.  Such classes include:
//
// - S2ContainsPointQuery: returns the shape(s) that contain a given point.
//
// - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
//                       S2CellId, or S2ShapeIndex.
//
// - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
//
// - S2BooleanOperation: computes boolean operations such as union,
//                       and boolean predicates such as containment.
//
// - S2ShapeIndexRegion: computes approximations for a collection of geometry.
//
// - S2ShapeIndexBufferedRegion: computes approximations that have been
//                               expanded by a given radius.
//
// Here is an example showing how to build an index for a set of polygons, and
// then then determine which polygon(s) contain each of a set of query points:
//
//   void TestContainment(const vector<S2Point>& points,
//                        const vector<S2Polygon*>& polygons) {
//     MutableS2ShapeIndex index;
//     for (auto polygon : polygons) {
//       index.Add(absl::make_unique<S2Polygon::Shape>(polygon));
//     }
//     auto query = MakeS2ContainsPointQuery(&index);
//     for (const auto& point : points) {
//       for (S2Shape* shape : query.GetContainingShapes(point)) {
//         S2Polygon* polygon = polygons[shape->id()];
//         ... do something with (point, polygon) ...
//       }
//     }
//   }
//
// This example uses S2Polygon::Shape, which is one example of an S2Shape
// object.  S2Polyline and S2Loop also have nested Shape classes, and there are
// additional S2Shape types defined in *_shape.h.
//
// Internally, MutableS2ShapeIndex is essentially a map from S2CellIds to the
// set of shapes that intersect each S2CellId.  It is adaptively refined to
// ensure that no cell contains more than a small number of edges.
//
// For efficiency, updates are batched together and applied lazily on the
// first subsequent query.  Locking is used to ensure that MutableS2ShapeIndex
// has the same thread-safety properties as "vector": const methods are
// thread-safe, while non-const methods are not thread-safe.  This means that
// if one thread updates the index, you must ensure that no other thread is
// reading or updating the index at the same time.
//
// TODO(ericv): MutableS2ShapeIndex has an Encode() method that allows the
// index to be serialized.  An encoded S2ShapeIndex can be decoded either into
// its original form (MutableS2ShapeIndex) or into an EncodedS2ShapeIndex.
// The key property of EncodedS2ShapeIndex is that it can be constructed
// instantaneously, since the index is kept in its original encoded form.
// Data is decoded only when an operation needs it.  For example, to determine
// which shapes(s) contain a given query point only requires decoding the data
// in the S2ShapeIndexCell that contains that point.
class MutableS2ShapeIndex : S2ShapeIndex {
private:

  alias CellMap = BTreeMap!(S2CellId, S2ShapeIndexCell);

public:
  // Options that affect construction of the MutableS2ShapeIndex.
  static struct Options {
  public:
    // The maximum number of edges per cell.  If a cell has more than this
    // many edges that are not considered "long" relative to the cell size,
    // then it is subdivided.  (Whether an edge is considered "long" is
    // controlled by --s2shape_index_cell_size_to_long_edge_ratio flag.)
    //
    // Values between 10 and 50 represent a reasonable balance between memory
    // usage, construction time, and query time.  Small values make queries
    // faster, while large values make construction faster and use less memory.
    // Values higher than 50 do not save significant additional memory, and
    // query times can increase substantially, especially for algorithms that
    // visit all pairs of potentially intersecting edges (such as polygon
    // validation), since this is quadratic in the number of edges per cell.
    //
    // Note that the *average* number of edges per cell is generally slightly
    // less than half of the maximum value defined here.
    //
    // Defaults to value given by --s2shape_index_default_max_edges_per_cell.
    @property
    int maxEdgesPerCell() const {
      return _maxEdgesPerCell;
    }

    @property
    void maxEdgesPerCell(int maxEdgesPerCell) {
      _maxEdgesPerCell = maxEdgesPerCell;
    }

  private:
    int _maxEdgesPerCell = S2SHAPE_INDEX_DEFAULT_MAX_EDGES_PER_CELL;
  }

  // Creates a MutableS2ShapeIndex that uses the default option settings.
  // Option values may be changed by calling Init().
  this() {
    _cellMap = new CellMap();
    _lock = new SpinLock();
    _indexStatus = IndexStatus.FRESH;
  }

  // Create a MutableS2ShapeIndex with the given options.
  this(Options options) {
    this();
    _options = options;
  }

  ~this() {
    //clear();
  }

  // Initialize a MutableS2ShapeIndex with the given options.  This method may
  // only be called when the index is empty (i.e. newly created or Reset() has
  // just been called).
  void init(Options options)
  in {
    assert(_shapes.length == 0);
  } body {
    _options = options;
  }

  const(Options) options() const {
    return _options;
  }

  // The number of distinct shape ids that have been assigned.  This equals
  // the number of shapes in the index provided that no shapes have ever been
  // removed.  (Shape ids are not reused.)
  override
  int numShapeIds() const {
    return cast(int) _shapes.length;
  }

  // Returns a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  override
  inout(S2Shape) shape(int id) inout {
    return _shapes[id];
  }

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  override
  void minimize() {
    // TODO(ericv): Implement.  In theory we should be able to discard the
    // entire index and rebuild it the next time it is needed.
  }

  final static class Iterator : IteratorBase {
  public:
    // Default constructor; must be followed by a call to Init().
    this() {
      _index = null;
    }

    // Constructs an iterator positioned as specified.  By default iterators
    // are unpositioned, since this avoids an extra seek in this situation
    // where one of the seek methods (such as Locate) is immediately called.
    //
    // If you want to position the iterator at the beginning, e.g. in order to
    // loop through the entire index, do this instead:
    //
    //   for (MutableS2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN);
    //        !it.done(); it.Next()) { ... }
    this(MutableS2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) {
      init(index, pos);
    }


    // Initializes an iterator for the given MutableS2ShapeIndex.  This method
    // may also be called in order to restore an iterator to a valid state
    // after the underlying index has been updated (although it is usually
    // easier just to declare a new iterator whenever required, since iterator
    // construction is cheap).
    void init(MutableS2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) {
      index.maybeApplyUpdates();
      initStale(index, pos);
    }

    // Initialize an iterator for the given MutableS2ShapeIndex without
    // applying any pending updates.  This can be used to observe the actual
    // current state of the index without modifying it in any way.
    void initStale(MutableS2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) {
      _index = index;
      _end = _index._cellMap.end();
      if (pos == InitialPosition.BEGIN) {
        _iter = _index._cellMap.begin();
      } else {
        _iter = _end;
      }
      refresh();
    }

    // Inherited non-virtual methods:
    //   S2CellId id() const;
    //   bool done() const;
    //   S2Point center() const;
    override
    inout(S2ShapeIndexCell) cell() inout {
      // Since MutableS2ShapeIndex always sets the "cell_" field, we can skip the
      // logic in the base class that conditionally calls GetCell().
      return super.cell();
    }

    // IteratorBase API:
    override
    void begin()
    in {
      // Make sure that the index has not been modified since Init() was called.
      assert(_index.isFresh());
    } body {
      _iter = _index._cellMap.begin();
      _end = _index._cellMap.end();
      refresh();
    }

    override
    void finish() {
      _iter = _end;
      refresh();
    }

    override
    void next()
    in {
      assert(!done());
    } body {
      ++_iter;
      refresh();
    }

    override
    bool prev() {
      if (_iter == _index._cellMap.begin()) {
        return false;
      }
      --_iter;
      refresh();
      return true;
    }

    override
    void seek(S2CellId target) {
      if (!_index._cellMap.equalRange(target).empty()) {
        _iter = _index._cellMap.equalRange(target).toIterator();
      } else {
        _iter = _index._cellMap.upperRange(target).toIterator();
      }
      refresh();
    }

    override
    bool locate(in S2Point target) {
      return IteratorBase.locateImpl(target, this);
    }

    override
    CellRelation locate(in S2CellId target) {
      return IteratorBase.locateImpl(target, this);
    }

    override
    void copy(IteratorBase other) {
      Iterator iter = cast(Iterator) other;
      enforce(iter !is null, "Only the same concrete Iterator type may be copied.");
      _index = iter._index;
      _iter = iter._iter;
      _end = iter._end;
    }

  protected:
    override
    const(S2ShapeIndexCell) getCell() const {
      enforce(false, "Should never be called");
      return null;
    }

    override
    IteratorBase clone() {
      auto iterator = new Iterator();
      iterator._index = _index;
      iterator._iter = _iter;
      iterator._end = _end;
      return iterator;
    }

  private:
    // Updates the IteratorBase fields.
    void refresh() {
      if (_iter == _end) {
        setFinished();
      } else {
        setState(_iter.getValue().key, _iter.getValue().value);
      }
    }

    MutableS2ShapeIndex _index;
    CellMap.Iterator _iter;
    CellMap.Iterator _end;
  }

  // Takes ownership of the given shape and adds it to the index.  Also
  // assigns a unique id to the shape (shape->id()) and returns that id.
  // Shape ids are assigned sequentially starting from 0 in the order shapes
  // are added.  Invalidates all iterators and their associated data.
  int add(S2Shape shape) {
    // Additions are processed lazily by ApplyUpdates().
    int id = cast(int) _shapes.length;
    shape._id = id;
    _shapes ~= shape;
    atomicStore!(MemoryOrder.raw)(_indexStatus, IndexStatus.STALE);
    return id;
  }

  // Removes the given shape from the index and return ownership to the caller.
  // Invalidates all iterators and their associated data.
  S2Shape release(int shapeId)
  in {
    assert(_shapes[shapeId] !is null);
  } body {
    // This class updates itself lazily, because it is much more efficient to
    // process additions and removals in batches.  However this means that when
    // a shape is removed, we need to make a copy of all its edges, since the
    // client is free to delete "shape" once this call is finished.

    S2Shape shape = _shapes[shapeId];
    if (shapeId >= _pendingAdditionsBegin) {
      // We are removing a shape that has not yet been added to the index,
      // so there is nothing else to do.
    } else {
      // We build the new RemovedShape in place, since it includes a potentially
      // large vector of edges that might be expensive to copy.
      _pendingRemovals ~= RemovedShape();
      RemovedShape* removed = &_pendingRemovals[$ - 1];
      removed.shapeId = shape.id();
      removed.hasInterior = shape.hasInterior();
      removed.containsTrackerOrigin = containsBruteForce(shape, interiorTrackerOrigin());
      int num_edges = shape.numEdges();
      removed.edges.length = 0;
      for (int e = 0; e < num_edges; ++e) {
        removed.edges ~= shape.edge(e);
      }
    }
    atomicStore!(MemoryOrder.raw)(_indexStatus, IndexStatus.STALE);
    return shape;
  }

  private static S2Point interiorTrackerOrigin() {
    return s2coords.FaceUVtoXYZ(0, -1, -1).normalize();
  }

  // Resets the index to its original state and returns ownership of all
  // shapes to the caller.  This method is much more efficient than removing
  // all shapes one at a time.
  S2Shape[] releaseAll()
  in {
    assert(_updateState is null);
  } body {
    _cellMap.clear();
    _pendingAdditionsBegin = 0;
    _pendingRemovals.length = 0;
    atomicStore!(MemoryOrder.raw)(_indexStatus, IndexStatus.FRESH);
    S2Shape[] result = _shapes[];
    swap(result, _shapes);
    return result;
  }

  // Resets the index to its original state and deletes all shapes.  Any
  // options specified via Init() are preserved.
  void clear() {
    releaseAll();
  }

  // Returns the number of bytes currently occupied by the index (including any
  // unused space at the end of vectors, etc). It has the same thread safety
  // as the other "const" methods (see introduction).
  override
  size_t spaceUsed() const {
    // TODO: Implement correctly when needed.
    return this.sizeof;
    /+
    size_t size = sizeof(*this);
    size += shapes_.capacity() * sizeof(std::unique_ptr<S2Shape>);
    // cell_map_ itself is already included in sizeof(*this).
    size += cell_map_.bytes_used() - sizeof(cell_map_);
    size += cell_map_.size() * sizeof(S2ShapeIndexCell);
    Iterator it;
    for (it.InitStale(this, S2ShapeIndex::BEGIN); !it.done(); it.Next()) {
      const S2ShapeIndexCell& cell = it.cell();
      size += cell.shapes_.capacity() * sizeof(S2ClippedShape);
      for (int s = 0; s < cell.num_clipped(); ++s) {
        const S2ClippedShape& clipped = cell.clipped(s);
        if (!clipped.is_inline()) {
          size += clipped.num_edges() * sizeof(int32);
        }
      }
    }
    if (pending_removals_ != nullptr) {
      size += pending_removals_->capacity() * sizeof(RemovedShape);
    }

    return size;
    +/
  }

  // Calls to Add() and Release() are normally queued and processed on the
  // first subsequent query (in a thread-safe way).  This has many advantages,
  // the most important of which is that sometimes there *is* no subsequent
  // query, which lets us avoid building the index completely.
  //
  // This method forces any pending updates to be applied immediately.
  // Calling this method is rarely a good idea.  (One valid reason is to
  // exclude the cost of building the index from benchmark results.)
  void forceBuild() {
    // No locks required because this is not a const method.  It is the client's
    // responsibility to ensure correct thread synchronization.
    if (atomicLoad!(MemoryOrder.raw)(_indexStatus) != IndexStatus.FRESH) {
      applyUpdatesInternal();
      atomicStore!(MemoryOrder.raw)(_indexStatus, IndexStatus.FRESH);
    }
  }

  // Returns true if there are no pending updates that need to be applied.
  // This can be useful to avoid building the index unnecessarily, or for
  // choosing between two different algorithms depending on whether the index
  // is available.
  //
  // The returned index status may be slightly out of date if the index was
  // built in a different thread.  This is fine for the intended use (as an
  // efficiency hint), but it should not be used by internal methods  (see
  // MaybeApplyUpdates).
  bool isFresh() const {
    return atomicLoad!(MemoryOrder.raw)(_indexStatus) == IndexStatus.FRESH;
  }

 protected:
  override
  IteratorBase newIterator(InitialPosition pos) {
    return new Iterator(this, pos);
  }

 private:
  /**
   * A BatchDescriptor represents a set of pending updates that will be applied
   * at the same time.  The batch consists of all updates with shape ids between
   * the current value of "ShapeIndex._pendingAdditionsBegin" (inclusive) and
   * "additionsEnd" (exclusive).  The first batch to be processed also
   * implicitly includes all shapes being removed.  "numEdges" is the total
   * number of edges that will be added or removed in this batch.
   */
  struct BatchDescriptor {
    int additionsEnd;
    int numEdges;
  }

  /**
   * FaceEdge and ClippedEdge store temporary edge data while the index is being
   * updated.  FaceEdge represents an edge that has been projected onto a given
   * face, while ClippedEdge represents the portion of that edge that has been
   * clipped to a given S2Cell.
   *
   * While it would be possible to combine all the edge information into one
   * structure, there are two good reasons for separating it:
   *
   *  - Memory usage.  Separating the two classes means that we only need to
   *    store one copy of the per-face data no matter how many times an edge is
   *    subdivided, and it also lets us delay computing bounding boxes until
   *    they are needed for processing each face (when the dataset spans
   *    multiple faces).
   *
   *  - Performance.  UpdateEdges is significantly faster on large polygons when
   *    the data is separated, because it often only needs to access the data in
   *    ClippedEdge and this data is cached more successfully.
   */
  struct FaceEdge {
    int shapeId;        // The shape that this edge belongs to
    int edgeId;         // Edge id within that shape
    int maxLevel;       // Not desirable to subdivide this edge beyond this level
    bool hasInterior;   // Belongs to a shape that has an interior
    R2Point a, b;       // The edge endpoints, clipped to a given face
    S2Shape.Edge edge;  // The edge endpoints
  }

  struct ClippedEdge {
    FaceEdge faceEdge;  // The original unclipped edge
    R2Rect bound;       // Bounding box for the clipped portion
  }

  /**
   * Given a set of shapes, InteriorTracker keeps track of which shapes contain
   * a particular point (the "focus").  It provides an efficient way to move the
   * focus from one point to another and incrementally update the set of shapes
   * which contain it.  We use this to compute which shapes contain the center
   * of every S2CellId in the index, by advancing the focus from one cell center
   * to the next.
   *
   * Initially the focus is at the start of the S2CellId space-filling curve.
   * We then visit all the cells that are being added to the MutableS2ShapeIndex
   * in increasing order of S2CellId.  For each cell, we draw two edges: one
   * from the entry vertex to the center, and another from the center to the
   * exit vertex (where "entry" and "exit" refer to the points where the
   * space-filling curve enters and exits the cell).  By counting edge crossings
   * we can incrementally compute which shapes contain the cell center.  Note
   * that the same set of shapes will always contain the exit point of one cell
   * and the entry point of the next cell in the index, because either (a) these
   * two points are actually the same, or (b) the intervening cells in S2CellId
   * order are all empty, and therefore there are no edge crossings if we follow
   * this path from one cell to the other.
   */
  class InteriorTracker {
  public:
    /**
     * Constructs the InteriorTracker.  You must call AddShape() for each shape
     * that will be tracked before calling MoveTo() or DrawTo().
     */
    this() {
      // As shapes are added, we compute which ones contain the start of the
      // S2CellId space-filling curve by drawing an edge from S2::Origin() to this
      // point and counting how many shape edges cross this edge.
      _isActive = false;
      _b = origin();
      _nextS2CellId = S2CellId.begin(S2CellId.MAX_LEVEL);
    }

    /**
     * Returns the initial focus point when the InteriorTracker is constructed
     * (corresponding to the start of the S2CellId space-filling curve).
     */
    static S2Point origin() {
      // The start of the S2CellId space-filling curve.
      return s2coords.FaceUVtoXYZ(0, -1, -1).normalize();
    }

    /// Returns the current focus point (see above).
    const(S2Point) focus() {
      return _b;
    }

    /// Returns true if any shapes are being tracked.
    bool isActive() const {
      return _isActive;
    }

    /**
     * Adds a shape whose interior should be tracked.  "is_inside" indicates
     * whether the current focus point is inside the shape.  Alternatively, if
     * the focus point is in the process of being moved (via MoveTo/DrawTo), you
     * can also specify "is_inside" at the old focus point and call TestEdge()
     * for every edge of the shape that might cross the current DrawTo() line.
     * This updates the state to correspond to the new focus point.
     *
     * REQUIRES: shape->has_interior()
     */
    void addShape(int shape_id, bool contains_focus) {
      _isActive = true;
      if (contains_focus) {
        toggleShape(shape_id);
      }
    }

    /**
     * Moves the focus to the given point.  This method should only be used when
     * it is known that there are no edge crossings between the old and new
     * focus locations; otherwise use DrawTo().
     */
    void moveTo(in S2Point b) {
      _b = b;
    }

    /**
     * Moves the focus to the given point.  After this method is called,
     * TestEdge() should be called with all edges that may cross the line
     * segment between the old and new focus locations.
     */
    void drawTo(in S2Point b) {
      _a = _b;
      _b = b;
      _crosser.init(_a, _b);
    }

    // Indicates that the given edge of the given shape may cross the line
    // segment between the old and new focus locations (see DrawTo).
    // REQUIRES: shape->has_interior()
    void testEdge(int shape_id, in S2Shape.Edge edge) {
      if (_crosser.edgeOrVertexCrossing(edge.v0, edge.v1)) {
        toggleShape(shape_id);
      }
    }


    // The set of shape ids that contain the current focus.
    const(ShapeIdSet) shapeIds() const {
      return _shapeIds;
    }

    // Indicates that the last argument to MoveTo() or DrawTo() was the entry
    // vertex of the given S2CellId, i.e. the tracker is positioned at the start
    // of this cell.  By using this method together with at_cellid(), the caller
    // can avoid calling MoveTo() in cases where the exit vertex of the previous
    // cell is the same as the entry vertex of the current cell.
    void setNextS2CellId(S2CellId next_cellid) {
      _nextS2CellId = next_cellid.rangeMin();
    }

    // Returns true if the focus is already at the entry vertex of the given
    // S2CellId (provided that the caller calls set_next_cellid() as each cell
    // is processed).
    bool atCellId(S2CellId cellid) const {
      return cellid.rangeMin() == _nextS2CellId;
    }

    // Makes an internal copy of the state for shape ids below the given limit,
    // and then clear the state for those shapes.  This is used during
    // incremental updates to track the state of added and removed shapes
    // separately.
    void saveAndClearStateBefore(int limit_shape_id)
    in {
      assert(_savedIds.length == 0);
    } body {
      size_t limit = lowerBound(limit_shape_id);
      _savedIds = _shapeIds[0 .. limit];
      _shapeIds = _shapeIds[limit .. $];
    }

    // Restores the state previously saved by SaveAndClearStateBefore().  This
    // only affects the state for shape_ids below "limit_shape_id".
    void restoreStateBefore(int limit_shape_id) {
      _shapeIds = _shapeIds[lowerBound(limit_shape_id) .. $];
      _shapeIds = _savedIds[] ~ _shapeIds[];
      _savedIds.length = 0;
    }

  private:
    // Removes "shape_id" from _shapeIds if it exists, otherwise insert it.
    void toggleShape(int shape_id) {
      import std.algorithm : remove;
      // Since shape_ids_.size() is typically *very* small (0, 1, or 2), it turns
      // out to be significantly faster to maintain a sorted array rather than
      // using an STL set or btree_set.
      if (_shapeIds.length == 0) {
        _shapeIds ~= shape_id;
      } else if (_shapeIds[0] == shape_id) {
        _shapeIds = _shapeIds.remove(0);
      } else {
        size_t pos = _shapeIds[0];
        while (_shapeIds[pos] < shape_id) {
          if (++pos == _shapeIds.length) {
            _shapeIds ~= shape_id;
            return;
          }
        }
        if (_shapeIds[pos] == shape_id) {
          _shapeIds.remove(pos);
        } else {
          _shapeIds.length++;
          _shapeIds[pos+1 .. $] = _shapeIds[pos .. $-1];
          _shapeIds[pos] = shape_id;
        }
      }
    }

    // Returns a pointer to the first entry "x" where x >= shape_id.
    size_t lowerBound(int shape_id) {
      size_t pos = 0;
      while (pos != _shapeIds.length && _shapeIds[pos] < shape_id) {
        ++pos;
      }
      return pos;
    }

    bool _isActive;
    S2Point _a, _b;
    S2CellId _nextS2CellId;
    S2EdgeCrosser _crosser;
    ShapeIdSet _shapeIds;

    // Shape ids saved by SaveAndClearStateBefore().  The state is never saved
    // recursively so we don't need to worry about maintaining a stack.
    ShapeIdSet _savedIds;
  }


  /**
   * EdgeAllocator provides temporary storage for new ClippedEdges that are
   * created during indexing.  It is essentially a stack model, where edges are
   * allocated as the recursion does down and freed as it comes back up.
   *
   * It also provides a mutable vector of FaceEdges that is used when
   * incrementally updating the index (see AbsorbIndexCell).
   */
  class EdgeAllocator {
  public:

    // Return a pointer to a newly allocated edge.  The EdgeAllocator
    // retains ownership.
    ref ClippedEdge newClippedEdge() {
      if (_size == _clippedEdges.length) {
        _clippedEdges ~= ClippedEdge();
      }
      return _clippedEdges[_size++];
    }

    // Return the number of allocated edges.
    size_t size() const {
      return _size;
    }

    // Reset the allocator to only contain the first "size" allocated edges.
    void reset(size_t size) {
      _size = size;
    }

    FaceEdge[]* mutableFaceEdges() {
      return &_faceEdges;
    }

  private:
    size_t _size = 0;
    ClippedEdge[] _clippedEdges;

    FaceEdge[] _faceEdges;
  }

  alias ShapeIdSet = int[];

  // Return true if this is the first update to the index.
  bool isFirstUpdate() const {
    // Note that it is not sufficient to check whether cell_map_ is empty, since
    // entries are added during the update process.
    return _pendingAdditionsBegin == 0;
  }

  // Given that the given shape is being updated, return true if it is being
  // removed (as opposed to being added).
  bool isShapeBeingRemoved(int shape_id) const {
    // All shape ids being removed are less than all shape ids being added.
    return shape_id < _pendingAdditionsBegin;
  }

  // Ensure that any pending updates have been applied.  This method must be
  // called before accessing the cell_map_ field, even if the index_status_
  // appears to be FRESH, because a memory barrier is required in order to
  // ensure that all the index updates are visible if the updates were done in
  // another thread.
  void maybeApplyUpdates() {
    // To avoid acquiring and releasing the spinlock on every query, we use
    // atomic operations when testing whether the status is FRESH and when
    // updating the status to be FRESH.  This guarantees that any thread that
    // sees a status of FRESH will also see the corresponding index updates.
    if (atomicLoad!(MemoryOrder.acq)(_indexStatus) != IndexStatus.FRESH) {
      this.applyUpdatesThreadSafe();
    }
  }

  // Apply any pending updates in a thread-safe way.
  void applyUpdatesThreadSafe() {
    _lock.lock();
    if (atomicLoad!(MemoryOrder.raw)(_indexStatus) == IndexStatus.FRESH) {
      _lock.unlock();
    } else if (atomicLoad!(MemoryOrder.raw)(_indexStatus) == IndexStatus.UPDATING) {
      // Wait until the updating thread is finished.  We do this by attempting
      // to lock a mutex that is held by the updating thread.  When this mutex
      // is unlocked the index_status_ is guaranteed to be FRESH.
      _updateState.numWaiting++;
      _lock.unlock();
      _updateState.waitMutex.lock();
      _lock.lock();
      _updateState.numWaiting--;
      unlockAndSignal();  // Notify other waiting threads.
    } else {
      enforce(_indexStatus == IndexStatus.STALE);
      atomicStore!(MemoryOrder.raw)(_indexStatus, IndexStatus.UPDATING);
      // Allocate the extra state needed for thread synchronization.  We keep
      // the spinlock held while doing this, because (1) memory allocation is
      // fast, so the chance of a context switch while holding the lock is low;
      // (2) by far the most common situation is that there is no contention,
      // and this saves an extra lock and unlock step; (3) even in the rare case
      // where there is contention, the main side effect is that some other
      // thread will burn a few CPU cycles rather than sleeping.
      _updateState = new UpdateState();
      // lock_.Lock wait_mutex *before* calling Unlock() to ensure that all other
      // threads will block on it.
      _updateState.waitMutex.lock();
      // Release the spinlock before doing any real work.
      _lock.unlock();
      applyUpdatesInternal();
      _lock.lock();
      // index_status_ can be updated to FRESH only while locked *and* using
      // an atomic store operation, so that MaybeApplyUpdates() can check
      // whether the index is FRESH without acquiring the spinlock.
      atomicStore!(MemoryOrder.rel)(_indexStatus, IndexStatus.FRESH);
      unlockAndSignal();  // Notify any waiting threads.
    }
  }

  // This method updates the index by applying all pending additions and
  // removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
  void applyUpdatesInternal() {
    // Check whether we have so many edges to process that we should process
    // them in multiple batches to save memory.  Building the index can use up
    // to 20x as much memory (per edge) as the final index size.
    BatchDescriptor[] batches;
    getUpdateBatches(batches);
    int i = 0;
    foreach (ref BatchDescriptor batch; batches) {
      FaceEdge[][6] all_edges;
      logger.logfDebug("Batch %d: shape_limit=%d, edges=%d\n",
          i++, batch.additionsEnd, batch.numEdges);

      reserveSpace(batch, all_edges);
      InteriorTracker tracker = new InteriorTracker();
      if (_pendingRemovals) {
        // The first batch implicitly includes all shapes being removed.
        foreach (ref pending_removal; _pendingRemovals) {
          removeShape(pending_removal, all_edges, tracker);
        }
        _pendingRemovals.length = 0;
      }
      for (int id = _pendingAdditionsBegin; id < batch.additionsEnd; ++id) {
        addShape(id, all_edges, tracker);
      }
      for (int face = 0; face < 6; ++face) {
        updateFaceEdges(face, all_edges[face], tracker);
        // Save memory by clearing vectors after we are done with them.
        foreach(edges; all_edges) {
          edges.length = 0;
        }
      }
      _pendingAdditionsBegin = batch.additionsEnd;
    }
    // It is the caller's responsibility to update index_status_.
  }

  // Count the number of edges being updated, and break them into several
  // batches if necessary to reduce the amount of memory needed.  (See the
  // documentation for FLAGS_s2shape_index_tmp_memory_budget_mb.)
  void getUpdateBatches(ref BatchDescriptor[] batches) const {
    // Count the edges being removed and added.
    int num_edges_removed = 0;
    if (_pendingRemovals) {
      foreach (pending_removal; _pendingRemovals) {
        num_edges_removed += pending_removal.edges.length;
      }
    }
    int num_edges_added = 0;
    for (int id = _pendingAdditionsBegin; id < _shapes.length; ++id) {
      const(S2Shape) shape = this.shape(id);
      if (shape is null) continue;
      num_edges_added += shape.numEdges();
    }
    int num_edges = num_edges_removed + num_edges_added;

    // The following memory estimates are based on heap profiling.
    //
    // The final size of a MutableS2ShapeIndex depends mainly on how finely the
    // index is subdivided, as controlled by Options::max_edges_per_cell() and
    // --s2shape_index_default_max_edges_per_cell. For realistic values of
    // max_edges_per_cell() and shapes with moderate numbers of edges, it is
    // difficult to get much below 8 bytes per edge.  [The minimum possible size
    // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
    // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
    // shapes_ vector.]
    //
    // The temporary memory consists mainly of the FaceEdge and ClippedEdge
    // structures plus a ClippedEdge pointer for every level of recursive
    // subdivision.  For very large indexes this can be 200 bytes per edge.
    const size_t kFinalBytesPerEdge = 8;
    const size_t kTmpBytesPerEdge = 200;
    const size_t kTmpMemoryBudgetBytes = cast(size_t)(S2SHAPE_INDEX_TMP_MEMORY_BUDGET_MB) << 20;

    // We arbitrarily limit the number of batches just as a safety measure.
    // With the current default memory budget of 100 MB, this limit is not
    // reached even when building an index of 350 million edges.
    const int kMaxUpdateBatches = 100;

    if (num_edges * kTmpBytesPerEdge <= kTmpMemoryBudgetBytes) {
      // We can update all edges at once without exceeding kTmpMemoryBudgetBytes.
      batches ~= BatchDescriptor(cast(int) _shapes.length, num_edges);
      return;
    }
    // Otherwise, break the updates into up to several batches, where the size
    // of each batch is chosen so that all batches use approximately the same
    // high-water memory.  GetBatchSizes() returns the recommended number of
    // edges in each batch.
    int[] batch_sizes;
    getBatchSizes(num_edges, kMaxUpdateBatches, kFinalBytesPerEdge,
        kTmpBytesPerEdge, kTmpMemoryBudgetBytes, batch_sizes);

    // We always process removed edges in a single batch, since (1) they already
    // take up a lot of memory because we have copied all their edges, and (2)
    // AbsorbIndexCell() uses (shapes_[id] == nullptr) to detect when a shape is
    // being removed, so in order to split the removals into batches we would
    // need a different approach (e.g., temporarily add fake entries to shapes_
    // and restore them back to nullptr as shapes are actually removed).
    num_edges = 0;
    if (_pendingRemovals) {
      num_edges += num_edges_removed;
      if (num_edges >= batch_sizes[0]) {
        batches ~= BatchDescriptor(_pendingAdditionsBegin, num_edges);
        num_edges = 0;
      }
    }
    // Keep adding shapes to each batch until the recommended number of edges
    // for that batch is reached, then move on to the next batch.
    for (int id = _pendingAdditionsBegin; id < _shapes.length; ++id) {
      const(S2Shape) shape = this.shape(id);
      if (shape is null) continue;
      num_edges += shape.numEdges();
      if (num_edges >= batch_sizes[batches.length]) {
        batches ~= BatchDescriptor(id + 1, num_edges);
        num_edges = 0;
      }
    }
    // Some shapes have no edges.  If a shape with no edges is the last shape to
    // be added or removed, then the final batch may not include it, so we fix
    // that problem here.
    batches[$-1].additionsEnd = cast(int) _shapes.length;
    enforce(batches.length <= kMaxUpdateBatches);
  }

  // Given "num_items" items, each of which uses "tmp_bytes_per_item" while it
  // is being updated but only "final_bytes_per_item" in the end, divide the
  // items into batches that have approximately the same *total* memory usage
  // consisting of the temporary memory needed for the items in the current
  // batch plus the final size of all the items that have already been
  // processed.  Use the fewest number of batches (but never more than
  // "max_batches") such that the total memory usage does not exceed the
  // combined final size of all the items plus "tmp_memory_budget_bytes".
  static void getBatchSizes(int num_items, int max_batches,
      double final_bytes_per_item,
      double tmp_bytes_per_item,
      double tmp_memory_budget_bytes,
      ref int[] batch_sizes) {
    import std.algorithm : min, max;
    import std.math : pow;
    // This code tries to fit all the data into the same memory space
    // ("total_budget_bytes") at every iteration.  The data consists of some
    // number of processed items (at "final_bytes_per_item" each), plus some
    // number being updated (at "tmp_bytes_per_item" each).  The space occupied
    // by the items being updated is the "free space".  At each iteration, the
    // free space is multiplied by (1 - final_bytes_per_item/tmp_bytes_per_item)
    // as the items are converted into their final form.
    double final_bytes = num_items * final_bytes_per_item;
    double final_bytes_ratio = final_bytes_per_item / tmp_bytes_per_item;
    double free_space_multiplier = 1 - final_bytes_ratio;

    // The total memory budget is the greater of the final size plus the allowed
    // temporary memory, or the minimum amount of memory required to limit the
    // number of batches to "max_batches".
    double total_budget_bytes = max(
        final_bytes + tmp_memory_budget_bytes,
        final_bytes / (1 - pow(free_space_multiplier, max_batches)));

    // "max_batch_items" is the number of items in the current batch.
    double max_batch_items = total_budget_bytes / tmp_bytes_per_item;
    batch_sizes.length = 0;
    for (int i = 0; i + 1 < max_batches && num_items > 0; ++i) {
      int batch_items = min(num_items, cast(int)(max_batch_items + 1));
      batch_sizes ~= batch_items;
      num_items -= batch_items;
      max_batch_items *= free_space_multiplier;
    }
    enforce(batch_sizes.length <= max_batches);
  }

  // Reserve an appropriate amount of space for the top-level face edges in the
  // current batch.  This data structure uses about half of the temporary memory
  // needed during index construction.  Furthermore, if the arrays are grown via
  // push_back() then up to 10% of the total run time consists of copying data
  // as these arrays grow, so it is worthwhile to preallocate space for them.
  void reserveSpace(BatchDescriptor batch, ref FaceEdge[][6] all_edges) const {
    import std.algorithm : max;
    // If the number of edges is relatively small, then the fastest approach is
    // to simply reserve space on every face for the maximum possible number of
    // edges.  We use a different threshold for this calculation than for
    // deciding when to break updates into batches, because the cost/benefit
    // ratio is different.  (Here the only extra expense is that we need to
    // sample the edges to estimate how many edges per face there are.)
    const size_t kMaxCheapBytes = 30 << 20;  // 30 MB
    const int kMaxCheapEdges = kMaxCheapBytes / (6 * FaceEdge.sizeof);
    if (batch.numEdges <= kMaxCheapEdges) {
      for (int face = 0; face < 6; ++face) {
        all_edges[face].reserve(batch.numEdges);
      }
      return;
    }
    // Otherwise we estimate the number of edges on each face by taking a random
    // sample.  The goal is to come up with an estimate that is fast and
    // accurate for non-pathological geometry.  If our estimates happen to be
    // wrong, the vector will still grow automatically - the main side effects
    // are that memory usage will be larger (by up to a factor of 3), and
    // constructing the index will be about 10% slower.
    //
    // Given a desired sample size, we choose equally spaced edges from
    // throughout the entire data set.  We use a Bresenham-type algorithm to
    // choose the samples.
    const int kDesiredSampleSize = 10000;
    const int sample_interval = max(1, batch.numEdges / kDesiredSampleSize);

    // Initialize "edge_id" to be midway through the first sample interval.
    // Because samples are equally spaced the actual sample size may differ
    // slightly from the desired sample size.
    int edge_id = sample_interval / 2;
    const int actual_sample_size = (batch.numEdges + edge_id) / sample_interval;
    int[6] face_count = [ 0, 0, 0, 0, 0, 0 ];
    if (_pendingRemovals) {
      foreach (removed; _pendingRemovals) {
        edge_id += removed.edges.length;
        while (edge_id >= sample_interval) {
          edge_id -= sample_interval;
          face_count[s2coords.GetFace(removed.edges[edge_id].v0)] += 1;
        }
      }
    }
    for (int id = _pendingAdditionsBegin; id < batch.additionsEnd; ++id) {
      const(S2Shape) shape = this.shape(id);
      if (shape is null) continue;
      edge_id += shape.numEdges();
      while (edge_id >= sample_interval) {
        edge_id -= sample_interval;
        // For speed, we only count the face containing one endpoint of the
        // edge.  In general the edge could span all 6 faces (with padding), but
        // it's not worth the expense to compute this more accurately.
        face_count[s2coords.GetFace(shape.edge(edge_id).v0)] += 1;
      }
    }
    // Now given the raw face counts, compute a confidence interval such that we
    // will be unlikely to allocate too little space.  Computing accurate
    // binomial confidence intervals is expensive and not really necessary.
    // Instead we use a simple approximation:
    //  - For any face with at least 1 sample, we use at least a 4-sigma
    //    confidence interval.  (The chosen width is adequate for the worst case
    //    accuracy, which occurs when the face contains approximately 50% of the
    //    edges.)  Assuming that our sample is representative, the probability of
    //    reserving too little space is approximately 1 in 30,000.
    //  - For faces with no samples at all, we don't bother reserving space.
    //    It is quite likely that such faces are truly empty, so we save time
    //    and memory this way.  If the face does contain some edges, there will
    //    only be a few so it is fine to let the vector grow automatically.
    // On average, we reserve 2% extra space for each face that has geometry.

    // kMaxSemiWidth is the maximum semi-width over all probabilities p of a
    // 4-sigma binomial confidence interval with a sample size of 10,000.
    const double kMaxSemiWidth = 0.02;
    const double sample_ratio = 1.0 / actual_sample_size;
    for (int face = 0; face < 6; ++face) {
      if (face_count[face] == 0) continue;
      double fraction = sample_ratio * face_count[face] + kMaxSemiWidth;
      all_edges[face].reserve(cast(int)(fraction * batch.numEdges));
    }
  }

  // Clip all edges of the given shape to the six cube faces, add the clipped
  // edges to "all_edges", and start tracking its interior if necessary.
  void addShape(int id, ref FaceEdge[][6] all_edges, InteriorTracker tracker) const {
    const S2Shape shape = this.shape(id);
    if (shape is null) {
      return;  // This shape has already been removed.
    }
    // Construct a template for the edges to be added.
    FaceEdge edge;
    edge.shapeId = id;
    edge.hasInterior = shape.hasInterior();
    if (edge.hasInterior) {
      tracker.addShape(id, containsBruteForce(shape, tracker.focus()));
    }
    int num_edges = shape.numEdges();
    for (int e = 0; e < num_edges; ++e) {
      edge.edgeId = e;
      edge.edge = shape.edge(e);
      edge.maxLevel = getEdgeMaxLevel(edge.edge);
      addFaceEdge(edge, all_edges);
    }
  }

  void removeShape(
      in RemovedShape removed, ref FaceEdge[][6] all_edges, InteriorTracker tracker) const {
    FaceEdge edge;
    edge.edgeId = -1;  // Not used or needed for removed edges.
    edge.shapeId = removed.shapeId;
    edge.hasInterior = removed.hasInterior;
    if (edge.hasInterior) {
      tracker.addShape(edge.shapeId, removed.containsTrackerOrigin);
    }
    foreach (removed_edge; removed.edges) {
      edge.edge = removed_edge;
      edge.maxLevel = getEdgeMaxLevel(edge.edge);
      addFaceEdge(edge, all_edges);
    }
  }

  void addFaceEdge(ref FaceEdge edge, ref FaceEdge[][6] all_edges) const {
    import s2.s2edge_clipping : clipToPaddedFace;
    import std.math : fabs;
    // Fast path: both endpoints are on the same face, and are far enough from
    // the edge of the face that don't intersect any (padded) adjacent face.
    int a_face = s2coords.GetFace(edge.edge.v0);
    if (a_face == s2coords.GetFace(edge.edge.v1)) {
      s2coords.ValidFaceXYZtoUV(a_face, edge.edge.v0, edge.a);
      s2coords.ValidFaceXYZtoUV(a_face, edge.edge.v1, edge.b);
      const double kMaxUV = 1 - CELL_PADDING;
      if (fabs(edge.a[0]) <= kMaxUV && fabs(edge.a[1]) <= kMaxUV
          && fabs(edge.b[0]) <= kMaxUV && fabs(edge.b[1]) <= kMaxUV) {
        all_edges[a_face] ~= edge;
        return;
      }
    }
    // Otherwise we simply clip the edge to all six faces.
    for (int face = 0; face < 6; ++face) {
      if (clipToPaddedFace(edge.edge.v0, edge.edge.v1, face,
              CELL_PADDING, edge.a, edge.b)) {
        all_edges[face] ~= edge;
      }
    }
  }

  // Given a face and a vector of edges that intersect that face, add or remove
  // all the edges from the index.  (An edge is added if shapes_[id] is not
  // nullptr, and removed otherwise.)
  void updateFaceEdges(int face, in FaceEdge[] face_edges, InteriorTracker tracker) {
    int num_edges = cast(int) face_edges.length;
    if (num_edges == 0 && tracker.shapeIds.length == 0) return;

    // Create the initial ClippedEdge for each FaceEdge.  Additional clipped
    // edges are created when edges are split between child cells.  We create
    // two arrays, one containing the edge data and another containing pointers
    // to those edges, so that during the recursion we only need to copy
    // pointers in order to propagate an edge to the correct child.
    ClippedEdge[] clipped_edge_storage;
    ClippedEdge[] clipped_edges;
    clipped_edge_storage.reserve(num_edges);
    clipped_edges.reserve(num_edges);
    R2Rect bound = R2Rect.empty();
    for (int e = 0; e < num_edges; ++e) {
      ClippedEdge clipped;
      clipped.faceEdge = face_edges[e];
      clipped.bound = R2Rect.fromPointPair(face_edges[e].a, face_edges[e].b);
      clipped_edge_storage ~= clipped;
      clipped_edges ~= clipped_edge_storage[$-1];
      bound.addRect(clipped.bound);
    }
    // Construct the initial face cell containing all the edges, and then update
    // all the edges in the index recursively.
    EdgeAllocator alloc = new EdgeAllocator();
    S2CellId face_id = S2CellId.fromFace(face);
    S2PaddedCell pcell = new S2PaddedCell(face_id, CELL_PADDING);

    // "disjoint_from_index" means that the current cell being processed (and
    // all its descendants) are not already present in the index.
    bool disjoint_from_index = isFirstUpdate();
    if (num_edges > 0) {
      S2CellId shrunk_id = shrinkToFit(pcell, bound);
      if (shrunk_id != pcell.id()) {
        // All the edges are contained by some descendant of the face cell.  We
        // can save a lot of work by starting directly with that cell, but if we
        // are in the interior of at least one shape then we need to create
        // index entries for the cells we are skipping over.
        skipCellRange(face_id.rangeMin(), shrunk_id.rangeMin(),
            tracker, alloc, disjoint_from_index);
        pcell = new S2PaddedCell(shrunk_id, CELL_PADDING);
        updateEdges(pcell, clipped_edges, tracker, alloc, disjoint_from_index);
        skipCellRange(shrunk_id.rangeMax().next(), face_id.rangeMax().next(),
            tracker, alloc, disjoint_from_index);
        return;
      }
    }
    // Otherwise (no edges, or no shrinking is possible), subdivide normally.
    updateEdges(pcell, clipped_edges, tracker, alloc, disjoint_from_index);
  }

  S2CellId shrinkToFit(in S2PaddedCell pcell, in R2Rect bound) {
    S2CellId shrunk_id = pcell.shrinkToFit(bound);
    if (!isFirstUpdate() && shrunk_id != pcell.id()) {
      // Don't shrink any smaller than the existing index cells, since we need
      // to combine the new edges with those cells.
      // Use InitStale() to avoid applying updated recursively.
      Iterator iter = new Iterator();
      iter.initStale(this);
      CellRelation r = iter.locate(shrunk_id);
      if (r == CellRelation.INDEXED) { shrunk_id = iter.id(); }
    }
    return shrunk_id;
  }

  // Skip over the cells in the given range, creating index cells if we are
  // currently in the interior of at least one shape.
  void skipCellRange(
      S2CellId begin, S2CellId end, InteriorTracker tracker, EdgeAllocator alloc,
      bool disjoint_from_index) {
    import s2.s2cell_union;
    // If we aren't in the interior of a shape, then skipping over cells is easy.
    if (tracker.shapeIds().empty()) return;

    // Otherwise generate the list of cell ids that we need to visit, and create
    // an index entry for each one.
    foreach (S2CellId skipped_id; S2CellUnion.fromBeginEnd(begin, end).cellIds()) {
      ClippedEdge[] clipped_edges;
      updateEdges(new S2PaddedCell(skipped_id, CELL_PADDING),
          clipped_edges, tracker, alloc, disjoint_from_index);
    }
  }

  // Given a cell and a set of ClippedEdges whose bounding boxes intersect that
  // cell, add or remove all the edges from the index.  Temporary space for
  // edges that need to be subdivided is allocated from the given EdgeAllocator.
  // "disjoint_from_index" is an optimization hint indicating that cell_map_
  // does not contain any entries that overlap the given cell.
  void updateEdges(S2PaddedCell pcell, ref ClippedEdge[] edges,
      InteriorTracker tracker,
      EdgeAllocator alloc,
      bool disjoint_from_index)
  in {
    // Cases where an index cell is not needed should be detected before this.
    assert(edges.length != 0 || tracker.shapeIds().length != 0);
  } body {
    // This function is recursive with a maximum recursion depth of 30
    // (S2CellId::kMaxLevel).  Note that using an explicit stack does not seem
    // to be any faster based on profiling.

    // Incremental updates are handled as follows.  All edges being added or
    // removed are combined together in "edges", and all shapes with interiors
    // are tracked using "tracker".  We subdivide recursively as usual until we
    // encounter an existing index cell.  At this point we "absorb" the index
    // cell as follows:
    //
    //   - Edges and shapes that are being removed are deleted from "edges" and
    //     "tracker".
    //   - All remaining edges and shapes from the index cell are added to
    //     "edges" and "tracker".
    //   - Continue subdividing recursively, creating new index cells as needed.
    //   - When the recursion gets back to the cell that was absorbed, we
    //     restore "edges" and "tracker" to their previous state.
    //
    // Note that the only reason that we include removed shapes in the recursive
    // subdivision process is so that we can find all of the index cells that
    // contain those shapes efficiently, without maintaining an explicit list of
    // index cells for each shape (which would be expensive in terms of memory).
    bool index_cell_absorbed = false;
    if (!disjoint_from_index) {
      // There may be existing index cells contained inside "pcell".  If we
      // encounter such a cell, we need to combine the edges being updated with
      // the existing cell contents by "absorbing" the cell.
      // Use InitStale() to avoid applying updated recursively.
      Iterator iter = new Iterator();
      iter.initStale(this);
      CellRelation r = iter.locate(pcell.id());
      if (r == CellRelation.DISJOINT) {
        disjoint_from_index = true;
      } else if (r == CellRelation.INDEXED) {
        // Absorb the index cell by transferring its contents to "edges" and
        // deleting it.  We also start tracking the interior of any new shapes.
        absorbIndexCell(pcell, iter, edges, tracker, alloc);
        index_cell_absorbed = true;
        disjoint_from_index = true;
      } else {
        enforce(r == CellRelation.SUBDIVIDED);
      }
    }

    // If there are existing index cells below us, then we need to keep
    // subdividing so that we can merge with those cells.  Otherwise,
    // MakeIndexCell checks if the number of edges is small enough, and creates
    // an index cell if possible (returning true when it does so).
    if (!disjoint_from_index || !makeIndexCell(pcell, edges, tracker)) {
      // Reserve space for the edges that will be passed to each child.  This is
      // important since otherwise the running time is dominated by the time
      // required to grow the vectors.  The amount of memory involved is
      // relatively small, so we simply reserve the maximum space for every child.
      ClippedEdge[][2][2] child_edges;  // [i][j]
      int num_edges = cast(int) edges.length;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          child_edges[i][j].reserve(num_edges);
        }
      }

      // Remember the current size of the EdgeAllocator so that we can free any
      // edges that are allocated during edge splitting.
      size_t alloc_size = alloc.size();

      // Compute the middle of the padded cell, defined as the rectangle in
      // (u,v)-space that belongs to all four (padded) children.  By comparing
      // against the four boundaries of "middle" we can determine which children
      // each edge needs to be propagated to.
      const R2Rect middle = pcell.middle();

      // Build up a vector edges to be passed to each child cell.  The (i,j)
      // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
      // Note that the vast majority of edges are propagated to a single child.
      // This case is very fast, consisting of between 2 and 4 floating-point
      // comparisons and copying one pointer.  (ClipVAxis is inline.)
      for (int e = 0; e < num_edges; ++e) {
        const ClippedEdge edge = edges[e];
        if (edge.bound[0].hi() <= middle[0].lo()) {
          // Edge is entirely contained in the two left children.
          clipVAxis(edge, middle[1], child_edges[0], alloc);
        } else if (edge.bound[0].lo() >= middle[0].hi()) {
          // Edge is entirely contained in the two right children.
          clipVAxis(edge, middle[1], child_edges[1], alloc);
        } else if (edge.bound[1].hi() <= middle[1].lo()) {
          // Edge is entirely contained in the two lower children.
          child_edges[0][0] ~= clipUBound(edge, 1, middle[0].hi(), alloc);
          child_edges[1][0] ~= clipUBound(edge, 0, middle[0].lo(), alloc);
        } else if (edge.bound[1].lo() >= middle[1].hi()) {
          // Edge is entirely contained in the two upper children.
          child_edges[0][1] ~= clipUBound(edge, 1, middle[0].hi(), alloc);
          child_edges[1][1] ~= clipUBound(edge, 0, middle[0].lo(), alloc);
        } else {
          // The edge bound spans all four children.  The edge itself intersects
          // either three or four (padded) children.
          ClippedEdge left = clipUBound(edge, 1, middle[0].hi(), alloc);
          clipVAxis(left, middle[1], child_edges[0], alloc);
          ClippedEdge right = clipUBound(edge, 0, middle[0].lo(), alloc);
          clipVAxis(right, middle[1], child_edges[1], alloc);
        }
      }

      // Now recursively update the edges in each child.  We call the children in
      // increasing order of S2CellId so that when the index is first constructed,
      // all insertions into cell_map_ are at the end (which is much faster).
      for (int pos = 0; pos < 4; ++pos) {
        int i, j;
        pcell.getChildIJ(pos, i, j);
        if (child_edges[i][j].length != 0 || tracker.shapeIds().length != 0) {
          updateEdges(new S2PaddedCell(pcell, i, j), child_edges[i][j],
              tracker, alloc, disjoint_from_index);
        }
      }
      // Free any temporary edges that were allocated during clipping.
      alloc.reset(alloc_size);
    }
    if (index_cell_absorbed) {
      // Restore the state for any edges being removed that we are tracking.
      tracker.restoreStateBefore(_pendingAdditionsBegin);
    }
  }

  // Absorb an index cell by transferring its contents to "edges" and/or
  // "tracker", and then delete this cell from the index.  If "edges" includes
  // any edges that are being removed, this method also updates their
  // InteriorTracker state to correspond to the exit vertex of this cell, and
  // saves the InteriorTracker state by calling SaveAndClearStateBefore().  It
  // is the caller's responsibility to restore this state by calling
  // RestoreStateBefore() when processing of this cell is finished.
  void absorbIndexCell(in S2PaddedCell pcell, in Iterator iter, ref ClippedEdge[] edges,
      InteriorTracker tracker, EdgeAllocator alloc)
  in {
    assert(pcell.id() == iter.id());
  } body {
    import s2.s2edge_clipping : clipToPaddedFace;
    // When we absorb a cell, we erase all the edges that are being removed.
    // However when we are finished with this cell, we want to restore the state
    // of those edges (since that is how we find all the index cells that need
    // to be updated).  The edges themselves are restored automatically when
    // UpdateEdges returns from its recursive call, but the InteriorTracker
    // state needs to be restored explicitly.
    //
    // Here we first update the InteriorTracker state for removed edges to
    // correspond to the exit vertex of this cell, and then save the
    // InteriorTracker state.  This state will be restored by UpdateEdges when
    // it is finished processing the contents of this cell.
    if (tracker.isActive() && edges.length != 0 &&
        isShapeBeingRemoved(edges[0].faceEdge.shapeId)) {
      // We probably need to update the InteriorTracker.  ("Probably" because
      // it's possible that all shapes being removed do not have interiors.)
      if (!tracker.atCellId(pcell.id())) {
        tracker.moveTo(pcell.getEntryVertex());
      }
      tracker.drawTo(pcell.getExitVertex());
      tracker.setNextS2CellId(pcell.id().next());
      foreach (ref const ClippedEdge edge; edges) {
        const FaceEdge face_edge = edge.faceEdge;
        if (!isShapeBeingRemoved(face_edge.shapeId)) {
          break;  // All shapes being removed come first.
        }
        if (face_edge.hasInterior) {
          tracker.testEdge(face_edge.shapeId, face_edge.edge);
        }
      }
    }
    // Save the state of the edges being removed, so that it can be restored
    // when we are finished processing this cell and its children.  We don't
    // need to save the state of the edges being added because they aren't being
    // removed from "edges" and will therefore be updated normally as we visit
    // this cell and its children.
    tracker.saveAndClearStateBefore(_pendingAdditionsBegin);

    // Create a FaceEdge for each edge in this cell that isn't being removed.
    FaceEdge[]* face_edges = alloc.mutableFaceEdges();
    face_edges.length = 0;
    bool tracker_moved = false;
    const S2ShapeIndexCell cell = iter.cell();
    for (int s = 0; s < cell.numClipped(); ++s) {
      const S2ClippedShape clipped = cell.clipped(s);
      int shape_id = clipped.shapeId();
      const S2Shape shape = this.shape(shape_id);
      if (shape is null) continue;  // This shape is being removed.
      int num_edges = clipped.numEdges();

      // If this shape has an interior, start tracking whether we are inside the
      // shape.  UpdateEdges() wants to know whether the entry vertex of this
      // cell is inside the shape, but we only know whether the center of the
      // cell is inside the shape, so we need to test all the edges against the
      // line segment from the cell center to the entry vertex.
      FaceEdge edge;
      edge.shapeId = shape.id();
      edge.hasInterior = shape.hasInterior();
      if (edge.hasInterior) {
        tracker.addShape(shape_id, clipped.containsCenter());
        // There might not be any edges in this entire cell (i.e., it might be
        // in the interior of all shapes), so we delay updating the tracker
        // until we see the first edge.
        if (!tracker_moved && num_edges > 0) {
          tracker.moveTo(pcell.getCenter());
          tracker.drawTo(pcell.getEntryVertex());
          tracker.setNextS2CellId(pcell.id());
          tracker_moved = true;
        }
      }
      for (int i = 0; i < num_edges; ++i) {
        int e = clipped.edge(i);
        edge.edgeId = e;
        edge.edge = shape.edge(e);
        edge.maxLevel = getEdgeMaxLevel(edge.edge);
        if (edge.hasInterior) tracker.testEdge(shape_id, edge.edge);
        if (!clipToPaddedFace(edge.edge.v0, edge.edge.v1,
                pcell.id().face(), CELL_PADDING,
                edge.a, edge.b)) {
          enforce(false, "Invariant failure in MutableS2ShapeIndex");

        }
        *face_edges ~= edge;
      }
    }
    // Now create a ClippedEdge for each FaceEdge, and put them in "new_edges".
    ClippedEdge[] new_edges;
    foreach (const FaceEdge face_edge; *face_edges) {
      ClippedEdge clipped = alloc.newClippedEdge();
      clipped.faceEdge = face_edge;
      clipped.bound = getClippedEdgeBound(face_edge.a, face_edge.b,
          pcell.bound());
      new_edges ~= clipped;
    }
    // Discard any edges from "edges" that are being removed, and append the
    // remainder to "new_edges".  (This keeps the edges sorted by shape id.)
    for (int i = 0; i < edges.length; ++i) {
      ClippedEdge clipped = edges[i];
      if (!isShapeBeingRemoved(clipped.faceEdge.shapeId)) {
        new_edges ~= edges[i .. $];
        break;
      }
    }
    // Update the edge list and delete this cell from the index.
    edges = new_edges;
    _cellMap.remove(pcell.id());
  }

  // Return the first level at which the edge will *not* contribute towards
  // the decision to subdivide.
  package int getEdgeMaxLevel(in S2Shape.Edge edge) const {
    import s2.s2metrics : AVG_EDGE;
    // Compute the maximum cell size for which this edge is considered "long".
    // The calculation does not need to be perfectly accurate, so we use Norm()
    // rather than Angle() for speed.
    double cell_size = (edge.v0 - edge.v1).norm() * S2SHAPE_INDEX_CELL_SIZE_TO_LONG_EDGE_RATIO;
    // Now return the first level encountered during subdivision where the
    // average cell size is at most "cell_size".
    return AVG_EDGE.getLevelForMaxValue(cell_size);
  }

  // Return the number of distinct shapes that are either associated with the
  // given edges, or that are currently stored in the InteriorTracker.
  static int countShapes(in ClippedEdge[] edges, in ShapeIdSet cshape_ids) {
    int count = 0;
    int last_shape_id = -1;
    size_t cnext_pos = 0;  // Next shape
    foreach (ClippedEdge edge; edges) {
      if (edge.faceEdge.shapeId != last_shape_id) {
        ++count;
        last_shape_id = edge.faceEdge.shapeId;
        // Skip over any containing shapes up to and including this one,
        // updating "count" appropriately.
        for (; cnext_pos != cshape_ids.length; ++cnext_pos) {
          if (cshape_ids[cnext_pos] > last_shape_id) break;
          if (cshape_ids[cnext_pos] < last_shape_id) ++count;
        }
      }
    }
    // Count any remaining containing shapes.
    count += (cshape_ids.length - cnext_pos);
    return count;
  }

  // Attempt to build an index cell containing the given edges, and return true
  // if successful.  (Otherwise the edges should be subdivided further.)
  bool makeIndexCell(in S2PaddedCell pcell, in ClippedEdge[] edges, InteriorTracker tracker) {
    if (edges.length == 0 && tracker.shapeIds().length == 0) {
      // No index cell is needed.  (In most cases this situation is detected
      // before we get to this point, but this can happen when all shapes in a
      // cell are removed.)
      return true;
    }

    // Count the number of edges that have not reached their maximum level yet.
    // Return false if there are too many such edges.
    int count = 0;
    foreach (const ClippedEdge edge; edges) {
      count += (pcell.level() < edge.faceEdge.maxLevel);
      if (count > _options.maxEdgesPerCell()) {
        return false;
      }
    }

    // Possible optimization: Continue subdividing as long as exactly one child
    // of "pcell" intersects the given edges.  This can be done by finding the
    // bounding box of all the edges and calling ShrinkToFit():
    //
    // S2CellId cellid = pcell.ShrinkToFit(GetRectBound(edges));
    //
    // Currently this is not beneficial; it slows down construction by 4-25%
    // (mainly computing the union of the bounding rectangles) and also slows
    // down queries (since more recursive clipping is required to get down to
    // the level of a spatial index cell).  But it may be worth trying again
    // once "contains_center" is computed and all algorithms are modified to
    // take advantage of it.

    // We update the InteriorTracker as follows.  For every S2Cell in the index
    // we construct two edges: one edge from entry vertex of the cell to its
    // center, and one from the cell center to its exit vertex.  Here "entry"
    // and "exit" refer the S2CellId ordering, i.e. the order in which points
    // are encountered along the S2 space-filling curve.  The exit vertex then
    // becomes the entry vertex for the next cell in the index, unless there are
    // one or more empty intervening cells, in which case the InteriorTracker
    // state is unchanged because the intervening cells have no edges.

    // Shift the InteriorTracker focus point to the center of the current cell.
    if (tracker.isActive() && edges.length != 0) {
      if (!tracker.atCellId(pcell.id())) {
        tracker.moveTo(pcell.getEntryVertex());
      }
      tracker.drawTo(pcell.getCenter());
      testAllEdges(edges, tracker);
    }
    // Allocate and fill a new index cell.  To get the total number of shapes we
    // need to merge the shapes associated with the intersecting edges together
    // with the shapes that happen to contain the cell center.
    const ShapeIdSet cshape_ids = tracker.shapeIds();
    int num_shapes = countShapes(edges, cshape_ids);
    S2ShapeIndexCell cell = new S2ShapeIndexCell();
    auto shapesAppender = cell.addShapes(num_shapes);

    // To fill the index cell we merge the two sources of shapes: "edge shapes"
    // (those that have at least one edge that intersects this cell), and
    // "containing shapes" (those that contain the cell center).  We keep track
    // of the index of the next intersecting edge and the next containing shape
    // as we go along.  Both sets of shape ids are already sorted.
    int enext = 0;
    int cnext_pos = 0;
    for (int i = 0; i < num_shapes; ++i) {
      S2ClippedShape clipped = S2ClippedShape();
      int eshape_id = numShapeIds(), cshape_id = eshape_id;  // Sentinels
      if (enext != edges.length) {
        eshape_id = edges[enext].faceEdge.shapeId;
      }
      if (cnext_pos != cshape_ids.length) {
        cshape_id = cnext_pos;
      }
      int ebegin = enext;
      if (cshape_id < eshape_id) {
        // The entire cell is in the shape interior.
        clipped.initialize(cshape_id, 0);
        clipped.setContainsCenter(true);
        ++cnext_pos;
      } else {
        // Count the number of edges for this shape and allocate space for them.
        while (enext < edges.length && edges[enext].faceEdge.shapeId == eshape_id) {
          ++enext;
        }
        clipped.initialize(eshape_id, enext - ebegin);
        for (int e = ebegin; e < enext; ++e) {
          clipped.setEdge(e - ebegin, edges[e].faceEdge.edgeId);
        }
        if (cshape_id == eshape_id) {
          clipped.setContainsCenter(true);
          ++cnext_pos;
        }
        shapesAppender ~= clipped;
      }
    }
    // UpdateEdges() visits cells in increasing order of S2CellId, so during
    // initial construction of the index all insertions happen at the end.  It
    // is much faster to give an insertion hint in this case.  Otherwise the
    // hint doesn't do much harm.  With more effort we could provide a hint even
    // during incremental updates, but this is probably not worth the effort.
    _cellMap[pcell.id()] = cell;

    // Shift the InteriorTracker focus point to the exit vertex of this cell.
    if (tracker.isActive() && edges.length != 0) {
      tracker.drawTo(pcell.getExitVertex());
      testAllEdges(edges, tracker);
      tracker.setNextS2CellId(pcell.id().next());
    }
    return true;
  }

  // Call tracker.testEdge() on all edges from shapes that have interiors.
  static void testAllEdges(in ClippedEdge[] edges, InteriorTracker tracker) {
    foreach (ClippedEdge edge; edges) {
      FaceEdge face_edge = edge.faceEdge;
      if (face_edge.hasInterior) {
        tracker.testEdge(face_edge.shapeId, face_edge.edge);
      }
    }
  }

  // Given an edge and two bound endpoints that need to be updated, allocate and
  // return a new edge with the updated bound.
  static const(ClippedEdge) updateBound(
      in ClippedEdge edge,
      int u_end, double u,
      int v_end, double v,
      EdgeAllocator alloc) {
    ClippedEdge clipped = alloc.newClippedEdge();
    clipped.faceEdge = edge.faceEdge;
    clipped.bound[0][u_end] = u;
    clipped.bound[1][v_end] = v;
    clipped.bound[0][1-u_end] = edge.bound[0][1 - u_end];
    clipped.bound[1][1-v_end] = edge.bound[1][1 - v_end];
    enforce(!clipped.bound.isEmpty());
    enforce(edge.bound.contains(clipped.bound));
    return clipped;
  }

  // Given an edge, clip the given endpoint (lo=0, hi=1) of the u-axis so that
  // it does not extend past the given value.
  static const(ClippedEdge) clipUBound(
      in ClippedEdge edge, int u_end, double u, EdgeAllocator alloc) {
    import s2.s2edge_clipping : interpolateDouble;
    // First check whether the edge actually requires any clipping.  (Sometimes
    // this method is called when clipping is not necessary, e.g. when one edge
    // endpoint is in the overlap area between two padded child cells.)
    if (u_end == 0) {
      if (edge.bound[0].lo() >= u) return edge;
    } else {
      if (edge.bound[0].hi() <= u) return edge;
    }
    // We interpolate the new v-value from the endpoints of the original edge.
    // This has two advantages: (1) we don't need to store the clipped endpoints
    // at all, just their bounding box; and (2) it avoids the accumulation of
    // roundoff errors due to repeated interpolations.  The result needs to be
    // clamped to ensure that it is in the appropriate range.
    FaceEdge e = edge.faceEdge;
    double v = edge.bound[1].project(interpolateDouble(u, e.a[0], e.b[0], e.a[1], e.b[1]));

    // Determine which endpoint of the v-axis bound to update.  If the edge
    // slope is positive we update the same endpoint, otherwise we update the
    // opposite endpoint.
    int v_end = u_end ^ ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1]));
    return updateBound(edge, u_end, u, v_end, v, alloc);
  }

  // Given an edge, clip the given endpoint (lo=0, hi=1) of the v-axis so that
  // it does not extend past the given value.
  static ClippedEdge clipVBound(
      in ClippedEdge edge,
      int v_end, double v,
      EdgeAllocator alloc) {
    import s2.s2edge_clipping : interpolateDouble;
    // See comments in ClipUBound.
    if (v_end == 0) {
      if (edge.bound[1].lo() >= v) return edge;
    } else {
      if (edge.bound[1].hi() <= v) return edge;
    }
    const FaceEdge e = edge.faceEdge;
    double u = edge.bound[0].project(
        interpolateDouble(v, e.a[1], e.b[1], e.a[0], e.b[0]));
    int u_end = v_end ^ ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1]));
    return updateBound(edge, u_end, u, v_end, v, alloc);
  }

  // Given an edge and an interval "middle" along the v-axis, clip the edge
  // against the boundaries of "middle" and add the edge to the corresponding
  // children.
  static void clipVAxis(in ClippedEdge edge, in R1Interval middle,
      ref ClippedEdge[][2] child_edges, EdgeAllocator alloc) {
    if (edge.bound[1].hi() <= middle.lo()) {
      // Edge is entirely contained in the lower child.
      child_edges[0] ~= edge;
    } else if (edge.bound[1].lo() >= middle.hi()) {
      // Edge is entirely contained in the upper child.
      child_edges[1] ~= edge;
    } else {
      // The edge bound spans both children.
      child_edges[0] ~= clipVBound(edge, 1, middle.hi(), alloc);
      child_edges[1] ~= clipVBound(edge, 0, middle.lo(), alloc);
    }
  }


  // The amount by which cells are "padded" to compensate for numerical errors
  // when clipping line segments to cell boundaries.
  // The total error when clipping an edge comes from two sources:
  // (1) Clipping the original spherical edge to a cube face (the "face edge").
  //     The maximum error in this step is S2::kFaceClipErrorUVCoord.
  // (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
  //     The maximum error in this step is S2::kEdgeClipErrorUVCoord.
  // Finally, since we encounter the same errors when clipping query edges, we
  // double the total error so that we only need to pad edges during indexing
  // and not at query time.
  package enum double CELL_PADDING = 2 * (FACE_CLIP_ERROR_UV_COORD + EDGE_CLIP_ERROR_UV_COORD);

  // The shapes in the index, accessed by their shape id.  Removed shapes are
  // replaced by nullptr pointers.
  /**/ package /**/ S2Shape[] _shapes;

  // A map from S2CellId to the set of clipped shapes that intersect that
  // cell.  The cell ids cover a set of non-overlapping regions on the
  // sphere.  Note that this field is updated lazily (see below).  Const
  // methods *must* call MaybeApplyUpdates() before accessing this field.
  // (The easiest way to achieve this is simply to use an Iterator.)
  /**/ package /**/ CellMap _cellMap;

  // The options supplied for this index.
  Options _options;

  // The id of the first shape that has been queued for addition but not
  // processed yet.
  int _pendingAdditionsBegin = 0;

  // The representation of an edge that has been queued for removal.
  struct RemovedShape {
    int shapeId;
    bool hasInterior;
    bool containsTrackerOrigin;
    S2Shape.Edge[] edges;
  }

  // The set of shapes that have been queued for removal but not processed
  // yet.  Note that we need to copy the edge data since the caller is free to
  // destroy the shape once Release() has been called.  This field is present
  // only when there are removed shapes to process (to save memory).
  RemovedShape[] _pendingRemovals;

  // Additions and removals are queued and processed on the first subsequent
  // query.  There are several reasons to do this:
  //
  //  - It is significantly more efficient to process updates in batches.
  //  - Often the index will never be queried, in which case we can save both
  //    the time and memory required to build it.  Examples:
  //     + S2Loops that are created simply to pass to an S2Polygon.  (We don't
  //       need the S2Loop index, because S2Polygon builds its own index.)
  //     + Applications that load a database of geometry and then query only
  //       a small fraction of it.
  //     + Applications that only read and write geometry (Decode/Encode).
  //
  // The main drawback is that we need to go to some extra work to ensure that
  // "const" methods are still thread-safe.  Note that the goal is *not* to
  // make this class thread-safe in general, but simply to hide the fact that
  // we defer some of the indexing work until query time.
  //
  // The textbook approach to this problem would be to use a mutex and a
  // condition variable.  Unfortunately pthread mutexes are huge (40 bytes).
  // Instead we use spinlock (which is only 4 bytes) to guard a few small
  // fields representing the current update status, and only create additional
  // state while the update is actually occurring.
  shared(SpinLock) _lock;

  enum IndexStatus {
    STALE,     // There are pending updates.
    UPDATING,  // Updates are currently being applied.
    FRESH     // There are no pending updates.
  }
  // Reads and writes to this field are guarded by "lock_".
  shared IndexStatus _indexStatus;

  // UpdateState holds temporary data related to thread synchronization.  It
  // is only allocated while updates are being applied.
  class UpdateState {
    // This mutex is used as a condition variable.  It is locked by the
    // updating thread for the entire duration of the update; other threads
    // lock it in order to wait until the update is finished.
    Mutex waitMutex;

    // The number of threads currently waiting on "wait_mutex_".  The
    // UpdateState can only be freed when this number reaches zero.
    //
    // Reads and writes to this field are guarded by "lock_".
    int numWaiting;

    this() {
      waitMutex = new Mutex();
    }
  }
  UpdateState _updateState;

  // Documented in the .cc file.
  void unlockAndSignal() {
    _lock.unlock();
    _updateState.waitMutex.unlock();
  }
}
