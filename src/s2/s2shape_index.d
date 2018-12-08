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

// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The main documentation is with the class definition below.
// (Some helper classes are defined first.)

module s2.s2shape_index;

import std.array : appender, RefAppender;
import std.bitmanip : bitfields;
import core.atomic : MemoryOrder, atomicStore, atomicLoad;
import s2.s2shape;
import s2.r1interval;
import s2.s2cell_id;
import s2.s2point;

import std.conv : to;


/**
 * S2ClippedShape represents the part of a shape that intersects an S2Cell.
 * It consists of the set of edge ids that intersect that cell, and a boolean
 * indicating whether the center of the cell is inside the shape (for shapes
 * that have an interior).
 *
 * Note that the edges themselves are not clipped; we always use the original
 * edges for intersection tests so that the results will be the same as the
 * original shape.
 */
struct S2ClippedShape {
public:

  /// The shape id of the clipped shape.
  int shapeId() const {
    return _shapeId;
  }

  /**
   * Returns true if the center of the S2CellId is inside the shape.  Returns
   * false for shapes that do not have an interior.
   */
  bool containsCenter() const {
    return _containsCenter;
  }

  /// The number of edges that intersect the S2CellId.
  int numEdges() const {
    return _numEdges;
  }

  /**
   * Returns the edge id of the given edge in this clipped shape.  Edges are
   * sorted in increasing order of edge id.
   */
  int edge(int i) const
  in {
    assert(0 <= i && i < numEdges());
  } body {
    return isInline() ? _inlineEdges[i] : _edges[i];
  }

  /// Returns true if the clipped shape contains the given edge id.
  bool containsEdge(int id) const {
    // Linear search is fast because the number of edges per shape is typically
    // very small (less than 10).
    for (int e = 0; e < numEdges(); ++e) {
      if (edge(e) == id) return true;
    }
    return false;
  }

package:
  // Initialize an S2ClippedShape to hold the given number of edges.
  void initialize(int shape_id, int num_edges) {
    _shapeId = shape_id;
    _numEdges = cast(uint) num_edges;
    _containsCenter = false;
    if (!isInline()) {
      _edges = new int[_numEdges];
    }
  }

  // This class may be copied by value, but note that it does *not* own its
  // underlying data.  (It is owned by the containing S2ShapeIndexCell.)

  // Free any memory allocated by this S2ClippedShape.  We don't do this in
  // the destructor because S2ClippedShapes are copied by STL code, and we
  // don't want to repeatedly copy and free the edge data.  Instead the data
  // is owned by the containing S2ShapeIndexCell.
  void destruct() {
    if (!isInline()) _edges.length = 0;
  }

  bool isInline() const {
    return _numEdges <= _inlineEdges.length;
  }

  // Set "contains_center_" to indicate whether this clipped shape contains the
  // center of the cell to which it belongs.
  void setContainsCenter(bool contains_center) {
    _containsCenter = contains_center;
  }

  // Set the i-th edge of this clipped shape to be the given edge of the
  // original shape.
  void setEdge(int i, int edge) {
    if (isInline()) {
      _inlineEdges[i] = edge;
    } else {
      _edges[i] = edge;
    }
  }

private:
  // All fields are packed into 16 bytes (assuming 64-bit pointers).  Up to
  // two edge ids are stored inline; this is an important optimization for
  // clients that use S2Shapes consisting of a single edge.
  int _shapeId;
  mixin(bitfields!(
      bool, "_containsCenter", 1,  // shape contains the cell center
      uint, "_numEdges",  31));

  // If there are more than two edges, this field holds a pointer.
  // Otherwise it holds an array of edge ids.
  union {
    int[] _edges;           // This pointer is owned by the containing Cell.
    int[2] _inlineEdges;
  }
}

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
class S2ShapeIndexCell {
public:

  // Returns the number of clipped shapes in this cell.
  int numClipped() const {
    return cast(int) _shapes.length;
  }

  // Returns the clipped shape at the given index.  Shapes are kept sorted in
  // increasing order of shape id.
  //
  // REQUIRES: 0 <= i < num_clipped()
  ref const(S2ClippedShape) clipped(int i) const
  in {
    assert(0 <= i && i < numClipped(),
        "clipped(" ~ i.to!string ~ ") but numClipped()=" ~ numClipped.to!string);
  } body {
    return _shapes[i];
  }

  // Returns a pointer to the clipped shape corresponding to the given shape,
  // or nullptr if the shape does not intersect this cell.
  const(S2ClippedShape)* findClipped(in S2Shape shape) const {
    return findClipped(shape.id());
  }

  const(S2ClippedShape)* findClipped(int shape_id) const {
    // Linear search is fine because the number of shapes per cell is typically
    // very small (most often 1), and is large only for pathological inputs
    // (e.g. very deeply nested loops).
    foreach (ref const s; _shapes) {
      if (s.shapeId() == shape_id) return &s;
    }
    return null;
  }

  // Convenience method that returns the total number of edges in all clipped
  // shapes.
  int numEdges() const {
    int n = 0;
    foreach (i; 0 .. numClipped()) {
      n += clipped(i).numEdges();
    }
    return n;
  }

  override
  string toString() const {
    import std.format : format;
    return format("S2ShapeIndexCell[shapes.length = %d]", _shapes.length);
  }

package:

  // Allocate room for "n" additional clipped shapes in the cell, and return a
  // pointer to the first new clipped shape.  Expects that all new clipped
  // shapes will have a larger shape id than any current shape, and that shapes
  // will be added in increasing shape id order.
  RefAppender!(S2ClippedShape[]) addShapes(int n) {
    _shapes.reserve(_shapes.length + n);
    return appender(&_shapes);
  }

private:
  S2ClippedShape[] _shapes;
}

/**
 * S2ShapeIndex is an abstract base class for indexing polygonal geometry in
 * memory.  The objects in the index are known as "shapes", and may consist of
 * points, polylines, and/or polygons, possibly overlapping.  The index makes
 * it very fast to answer queries such as finding nearby shapes, measuring
 * distances, testing for intersection and containment, etc.
 *
 * Each object in the index implements the S2Shape interface.  An S2Shape is a
 * collection of edges that optionally defines an interior.  The edges do not
 * need to be connected, so for example an S2Shape can represent a polygon
 * with multiple shells and/or holes, or a set of polylines, or a set of
 * points.  All geometry within a single S2Shape must have the same dimension,
 * so for example if you want to create an S2ShapeIndex containing a polyline
 * and 10 points, then you will need at least two different S2Shape objects.
 *
 * The most important type of S2ShapeIndex is MutableS2ShapeIndex, which
 * allows you to build an index incrementally by adding or removing shapes.
 * Soon there will also be an EncodedS2ShapeIndex type that makes it possible
 * to keep the index data in encoded form.  Code that only needs read-only
 * ("const") access to an index should use the S2ShapeIndex base class as the
 * parameter type, so that it will work with any S2ShapeIndex subtype.  For
 * example:
 *
 *   void DoSomething(const S2ShapeIndex& index) {
 *     ... works with MutableS2ShapeIndex or EncodedS2ShapeIndex ...
 *   }
 *
 * There are a number of built-in classes that work with S2ShapeIndex objects.
 * Generally these classes accept any collection of geometry that can be
 * represented by an S2ShapeIndex, i.e. any combination of points, polylines,
 * and polygons.  Such classes include:
 *
 * - S2ContainsPointQuery: returns the shape(s) that contain a given point.
 *
 * - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
 *                       S2CellId, or S2ShapeIndex.
 *
 * - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
 *
 * - S2BooleanOperation: computes boolean operations such as union,
 *                       and boolean predicates such as containment.
 *
 * - S2ShapeIndexRegion: computes approximations for a collection of geometry.
 *
 * - S2ShapeIndexBufferedRegion: computes approximations that have been
 *                               expanded by a given radius.
 *
 * Here is an example showing how to index a set of polygons and then
 * determine which polygon(s) contain each of a set of query points:
 *
 *   void TestContainment(const vector<S2Point>& points,
 *                        const vector<S2Polygon*>& polygons) {
 *     MutableS2ShapeIndex index;
 *     for (auto polygon : polygons) {
 *       index.Add(absl::make_unique<S2Polygon::Shape>(polygon));
 *     }
 *     auto query = MakeS2ContainsPointQuery(&index);
 *     for (const auto& point : points) {
 *       for (S2Shape* shape : query.GetContainingShapes(point)) {
 *         S2Polygon* polygon = polygons[shape->id()];
 *         ... do something with (point, polygon) ...
 *       }
 *     }
 *   }
 *
 * This example uses S2Polygon::Shape, which is one example of an S2Shape
 * object.  S2Polyline and S2Loop also have nested Shape classes, and there are
 * additional S2Shape types defined in *_shape.h.
 *
 * Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set
 * of shapes that intersect each S2CellId.  It is adaptively refined to ensure
 * that no cell contains more than a small number of edges.
 *
 * In addition to implementing a shared set of virtual methods, all
 * S2ShapeIndex subtypes define an Iterator type with the same API.  This
 * makes it easy to convert code that uses a particular S2ShapeIndex subtype
 * to instead use the abstract base class (or vice versa).  You can also
 * choose to avoid the overhead of virtual method calls by making the
 * S2ShapeIndex type a template argument, like this:
 *
 *   template <class IndexType>
 *   void DoSomething(const IndexType& index) {
 *     for (typename IndexType::Iterator it(&index, S2ShapeIndex::BEGIN);
 *          !it.done(); it.Next()) {
 *       ...
 *     }
 *   }
 *
 * Subtypes provided by the S2 library have the same thread-safety properties
 * as std::vector.  That is, const methods may be called concurrently from
 * multiple threads, and non-const methods require exclusive access to the
 * S2ShapeIndex.
 */
abstract class S2ShapeIndex {
public:
  // Returns the number of distinct shape ids in the index.  This is the same
  // as the number of shapes provided that no shapes have ever been removed.
  // (Shape ids are never reused.)
  abstract int numShapeIds() const;

  // Returns a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  abstract inout(S2Shape) shape(int id) inout;

  // Returns the number of bytes currently occupied by the index (including any
  // unused space at the end of vectors, etc).
  abstract size_t spaceUsed() const;

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  abstract void minimize();

  // The possible relationships between a "target" cell and the cells of the
  // S2ShapeIndex.  If the target is an index cell or is contained by an index
  // cell, it is "INDEXED".  If the target is subdivided into one or more
  // index cells, it is "SUBDIVIDED".  Otherwise it is "DISJOINT".
  enum CellRelation {
    INDEXED,       // Target is contained by an index cell
    SUBDIVIDED,    // Target is subdivided into one or more index cells
    DISJOINT       // Target does not intersect any index cells
  }

  // When passed to an Iterator constructor, specifies whether the iterator
  // should be positioned at the beginning of the index (BEGIN), the end of
  // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
  // unpositioned, since this avoids an extra seek in this situation where one
  // of the seek methods (such as Locate) is immediately called.
  enum InitialPosition {
    BEGIN,
    END,
    UNPOSITIONED
  }

  /**
   * A random access iterator that provides low-level access to the cells of
   * the index.  Cells are sorted in increasing order of S2CellId.
   */
  static class Iterator {
  public:
    /// Default constructor; must be followed by a call to Init().
    this() {
      _iter = null;
    }

    /**
     * Constructs an iterator positioned as specified.  By default iterators
     * are unpositioned, since this avoids an extra seek in this situation
     * where one of the seek methods (such as Locate) is immediately called.
     *
     * If you want to position the iterator at the beginning, e.g. in order to
     * loop through the entire index, do this instead:
     *
     *   for (auto it = S2ShapeIndex.Iterator(&index, S2ShapeIndex.InitialPosition.BEGIN);
     *        !it.done(); it.Next()) { ... }
     */
    this(S2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) {
      _iter = index.newIterator(pos);
    }

    /**
     * Initializes an iterator for the given S2ShapeIndex.  This method may
     * also be called in order to restore an iterator to a valid state after
     * the underlying index has been updated (although it is usually easier
     * just to declare a new iterator whenever required, since iterator
     * construction is cheap).
     */
    void initialize(S2ShapeIndex index, InitialPosition pos = InitialPosition.UNPOSITIONED) {
      _iter = index.newIterator(pos);
    }

    /**
     * Returns the S2CellId of the current index cell.  If done() is true,
     * returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
     */
    S2CellId id() const {
      return _iter.id();
    }

    /// Returns the center point of the cell.
    S2Point center() const
    in {
      assert(!done());
    } body {
      return id().toS2Point();
    }

    /// Returns a reference to the contents of the current index cell.
    inout(S2ShapeIndexCell) cell() inout
    in {
      assert(!done());
    } body {
      return _iter.cell();
    }

    /// Returns true if the iterator is positioned past the last index cell.
    bool done() const {
      return _iter.done();
    }

    /// Positions the iterator at the first index cell (if any).
    void begin() {
      _iter.begin();
    }

    /// Positions the iterator past the last index cell.
    void finish() {
      _iter.finish();
    }

    /// Positions the iterator at the next index cell.
    void next()
    in {
      assert(!done());
    } body {
      _iter.next();
    }

    /**
     * If the iterator is already positioned at the beginning, returns false.
     * Otherwise positions the iterator at the previous entry and returns true.
     */
    bool prev() {
      return _iter.prev();
    }

    /**
     * Positions the iterator at the first cell with id() >= target, or at the
     * end of the index if no such cell exists.
     */
    void seek(S2CellId target) {
      _iter.seek(target);
    }

    /**
     * Positions the iterator at the cell containing "target".  If no such cell
     * exists, returns false and leaves the iterator positioned arbitrarily.
     * The returned index cell is guaranteed to contain all edges that might
     * intersect the line segment between "target" and the cell center.
     */
    bool locate(in S2Point target) {
      return IteratorBase.locateImpl(target, this);
    }

    /**
     * Let T be the target S2CellId.  If T is contained by some index cell I
     * (including equality), this method positions the iterator at I and
     * returns INDEXED.  Otherwise if T contains one or more (smaller) index
     * cells, it positions the iterator at the first such cell I and returns
     * SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
     * positioned arbitrarily.
     */
    CellRelation locate(S2CellId target) {
      return IteratorBase.locateImpl(target, this);
    }

    void copy(Iterator other) {
      _iter = other._iter.clone();
    }

  private:
    // Although S2ShapeIndex::Iterator can be used to iterate over any
    // index subtype, it is more efficient to use the subtype's iterator when
    // the subtype is known at compile time.  For example, MutableS2ShapeIndex
    // should use a MutableS2ShapeIndex::Iterator.
    //
    // The following declarations prevent accidental use of
    // S2ShapeIndex::Iterator when the actual subtype is known.  (If you
    // really want to do this, you can down_cast the index argument to
    // S2ShapeIndex.)
    // template <class T>
    // explicit Iterator(const T* index, InitialPosition pos = UNPOSITIONED) {}

    // template <class T>
    // void Init(const T* index, InitialPosition pos = UNPOSITIONED) {}

    IteratorBase _iter;
  }

protected:
  /**
   * Each subtype of S2ShapeIndex should define an Iterator type derived
   * from the following base class.
   */
  abstract static class IteratorBase {
  public:
    this(IteratorBase other) {
      _id = other._id;
      _cell = other.rawCell();
    }

    /**
     * Returns the S2CellId of the current index cell.  If done() is true,
     * returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
     */
    @property
    S2CellId id() const {
      return _id;
    }

    /// Returns the center point of the cell.
    S2Point center() const
    in {
      assert(!done());
    } body {
      return id().toS2Point();
    }

    /// Returns a reference to the contents of the current index cell.
    inout(S2ShapeIndexCell) cell() inout
    in {
      assert(!done());
    } body {
      auto cell = _cell;
      // if (cell == null) {
      //   cell = getCell();
      //   setCell(cell);
      // }
      return cell;
    }

    /// Returns true if the iterator is positioned past the last index cell.
    bool done() const {
      return _id == S2CellId.sentinel();
    }

    /// Positions the iterator at the first index cell (if any).
    abstract void begin();

    /// Positions the iterator past the last index cell.
    abstract void finish();

    /// Positions the iterator at the next index cell.
    abstract void next()
    in {
      assert(!done());
    }

    /**
     * If the iterator is already positioned at the beginning, returns false.
     * Otherwise positions the iterator at the previous entry and returns true.
     */
    abstract bool prev();

    /**
     * Positions the iterator at the first cell with id() >= target, or at the
     * end of the index if no such cell exists.
     */
    abstract void seek(S2CellId target);

    /**
     * Positions the iterator at the cell containing "target".  If no such cell
     * exists, returns false and leaves the iterator positioned arbitrarily.
     * The returned index cell is guaranteed to contain all edges that might
     * intersect the line segment between "target" and the cell center.
     */
    abstract bool locate(in S2Point target);

    /**
     * Let T be the target S2CellId.  If T is contained by some index cell I
     * (including equality), this method positions the iterator at I and
     * returns INDEXED.  Otherwise if T contains one or more (smaller) index
     * cells, it positions the iterator at the first such cell I and returns
     * SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
     * positioned arbitrarily.
     */
    abstract CellRelation locate(in S2CellId target);

    /**
     * Makes a copy of the given source iterator.
     * REQUIRES: "other" has the same concrete type as "this".
     */
    abstract void copy(IteratorBase other);

  protected:
    this() {
      _id = S2CellId.sentinel();
      _cell = null;
    }

    /**
     * Sets the iterator state.  "cell" typically points to the cell contents,
     * but may also be given as "nullptr" in order to implement decoding on
     * demand.  In that situation, the first that the client attempts to
     * access the cell contents, the GetCell() method is called and "cell_" is
     * updated in a thread-safe way.
     */
    void setState(S2CellId id, S2ShapeIndexCell cell)
    in {
      assert(cell !is null);
    } body {
      _id = id;
      setCell(cell);
    }

    /// Sets the iterator state so that done() is true.
    void setFinished() {
      _id = S2CellId.sentinel();
      setCell(null);
    }

    /**
     * Returns the current contents of the "cell_" field, which may be null
     * if the cell contents have not been decoded yet.
     */
    inout(S2ShapeIndexCell) rawCell() inout {
      return _cell;
    }

    /**
     * This method is called to decode the contents of the current cell, if
     * set_state() was previously called with a nullptr "cell" argument.  This
     * allows decoding on demand for subtypes that keep the cell contents in
     * an encoded state.  It does not need to be implemented at all if
     * set_state() is always called with (cell != nullptr).
     *
     * REQUIRES: This method is thread-safe.
     * REQUIRES: Multiple calls to this method return the same value.
     */
    abstract const(S2ShapeIndexCell) getCell() const;

    /// Returns an exact copy of this iterator.
    abstract IteratorBase clone();

    /**
     * The default implementation of Locate(S2Point).  It is instantiated by
     * each subtype in order to (1) minimize the number of virtual method
     * calls (since subtypes are typically "final") and (2) ensure that the
     * correct versions of non-virtual methods such as cell() are called.
     */
    static bool locateImpl(IterT)(in S2Point target_point, IterT it) {
      // Let I = cell_map_->lower_bound(T), where T is the leaf cell containing
      // "target_point".  Then if T is contained by an index cell, then the
      // containing cell is either I or I'.  We test for containment by comparing
      // the ranges of leaf cells spanned by T, I, and I'.

      auto target = S2CellId(target_point);
      it.seek(target);
      if (!it.done() && it.id().rangeMin() <= target) return true;
      if (it.prev() && it.id().rangeMax() >= target) return true;
      return false;
    }


    // The default implementation of Locate(S2CellId) (see comments above).
    static CellRelation locateImpl(IterT)(in S2CellId target, IterT it) {
      // Let T be the target, let I = cell_map_->lower_bound(T.range_min()), and
      // let I' be the predecessor of I.  If T contains any index cells, then T
      // contains I.  Similarly, if T is contained by an index cell, then the
      // containing cell is either I or I'.  We test for containment by comparing
      // the ranges of leaf cells spanned by T, I, and I'.

      it.seek(target.rangeMin());
      if (!it.done()) {
        if (it.id() >= target && it.id().rangeMin() <= target) return CellRelation.INDEXED;
        if (it.id() <= target.rangeMax()) return CellRelation.SUBDIVIDED;
      }
      if (it.prev() && it.id().rangeMax() >= target) return CellRelation.INDEXED;
      return CellRelation.DISJOINT;
    }

  private:

    void setCell(S2ShapeIndexCell cell) {
      _cell = cell;
    }

    S2CellId _id;

    // Must be accessed atomically using atomicLoad and atomicStore.
    S2ShapeIndexCell _cell;
  }

  // Returns a new iterator positioned as specified.
  abstract IteratorBase newIterator(InitialPosition pos);
}
