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

module s2.s2cell_union;

import s2.logger;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2region;

import std.algorithm : countUntil, min, max, isSorted, sort;
import std.exception : enforce;
import std.range;

// The maximum number of cells allowed by S2CellUnion.Decode
enum int S2CELL_UNION_DECODE_MAX_NUM_CELLS = 1_000_000;

/**
 * An S2CellUnion is a region consisting of cells of various sizes.  Typically
 * a cell union is used to approximate some other shape.  There is a tradeoff
 * between the accuracy of the approximation and how many cells are used.
 * Unlike polygons, cells have a fixed hierarchical structure.  This makes
 * them more suitable for optimizations based on preprocessing.
 *
 * An S2CellUnion is represented as a vector of sorted, non-overlapping
 * S2CellIds.  By default the vector is also "normalized", meaning that groups
 * of 4 child cells have been replaced by their parent cell whenever possible.
 * S2CellUnions are not required to be normalized, but certain operations will
 * return different results if they are not (e.g., Contains(S2CellUnion).)
 *
 * S2CellUnion is movable and copyable.
 */
final class S2CellUnion : S2Region {
public:
  /// Creates an empty cell union.
  this() {}

  /**
   * Constructs a cell union with the given S2CellIds, then calls Normalize()
   * to sort them, remove duplicates, and merge cells when possible.  (See
   * FromNormalized if your vector is already normalized.)
   *
   * The argument is passed by value, so if you are passing a named variable
   * and have no further use for it, consider using std::move().
   *
   * A cell union containing a single S2CellId may be constructed like this:
   *
   *     S2CellUnion example({cell_id});
   */
  this(in S2CellId[] cell_ids) {
    _cellIds = cell_ids.dup;
    normalize();
  }

  // Convenience constructor that accepts a vector of uint64.  Note that
  // unlike the constructor above, this one makes a copy of "cell_ids".
  this(in ulong[] cell_ids) {
    _cellIds = toS2CellIds(cell_ids.dup);
    normalize();
  }

  /**
   * Constructs a cell union from S2CellIds that have already been normalized
   * (typically because they were extracted from another S2CellUnion).
   *
   * The argument is passed by value, so if you are passing a named variable
   * and have no further use for it, consider using std::move().
   *
   * REQUIRES: "cell_ids" satisfies the requirements of IsNormalized().
   */
  static S2CellUnion fromNormalized(S2CellId[] cell_ids) {
    auto result = new S2CellUnion(cell_ids, VerbatimFlag.VERBATIM);
    enforce(result.isNormalized());
    return result;
  }

  /**
   * Constructs a cell union from a vector of sorted, non-overlapping
   * S2CellIds.  Unlike the other constructors, FromVerbatim does not require
   * that groups of 4 child cells have been replaced by their parent cell.  In
   * other words, "cell_ids" must satisfy the requirements of IsValid() but
   * not necessarily IsNormalized().
   *
   * Note that if the cell union is not normalized, certain operations may
   * return different results (e.g., Contains(S2CellUnion)).
   *
   * REQUIRES: "cell_ids" satisfies the requirements of IsValid().
   */
  static S2CellUnion fromVerbatim(S2CellId[] cell_ids) {
    auto result = new S2CellUnion(cell_ids, VerbatimFlag.VERBATIM);
    enforce(result.isValid());
    return result;
  }

  /**
   * Constructs a cell union that corresponds to a continuous range of cell
   * ids.  The output is a normalized collection of cell ids that covers the
   * leaf cells between "min_id" and "max_id" inclusive.
   *
   * REQUIRES: min_id.is_leaf(), max_id.is_leaf(), min_id <= max_id.
   */
  static S2CellUnion fromMinMax(S2CellId min_id, S2CellId max_id) {
    auto result = new S2CellUnion();
    result.initFromMinMax(min_id, max_id);
    return result;
  }

  /**
   * Like FromMinMax() except that the union covers the range of leaf cells
   * from "begin" (inclusive) to "end" (exclusive), as with Python ranges or
   * STL iterator ranges.  If (begin == end) the result is empty.
   *
   * REQUIRES: begin.is_leaf(), end.is_leaf(), begin <= end.
   */
  static S2CellUnion fromBeginEnd(S2CellId begin, S2CellId end) {
    auto result = new S2CellUnion();
    result.initFromBeginEnd(begin, end);
    return result;
  }

  /**
   * initialize() methods corresponding to the constructors/factory methods above.
   * TODO(ericv): Consider deprecating these methods in favor of using the
   * constructors and move assignment operator.
   */
  void initialize(S2CellId[] cell_ids) {
    _cellIds = cell_ids;
    normalize();
  }

  void initialize(ulong[] cell_ids) {
    _cellIds = toS2CellIds(cell_ids);
    normalize();
  }

  void initFromMinMax(S2CellId min_id, S2CellId max_id)
  in {
    assert(max_id.isValid());
  } body {
    initFromBeginEnd(min_id, max_id.next());
  }

  void initFromBeginEnd(S2CellId begin, S2CellId end)
  in {
    assert(begin.isLeaf());
    assert(end.isLeaf());
    assert(begin <= end);
  } out {
    // The output is already normalized.
    assert(isNormalized());
  } body {
    // We repeatedly add the largest cell we can.
    _cellIds.length = 0;
    for (S2CellId id = begin.maximumTile(end);
         id != end; id = id.next().maximumTile(end)) {
      _cellIds ~= id;
    }
  }

  // Returns true if the given four cells have a common parent.
  // REQUIRES: The four cells are distinct.
  static bool areSiblings(S2CellId a, S2CellId b, S2CellId c, S2CellId d) {
    // A necessary (but not sufficient) condition is that the XOR of the
    // four cells must be zero.  This is also very fast to test.
    if ((a.id() ^ b.id() ^ c.id()) != d.id()) return false;

    // Now we do a slightly more expensive but exact test.  First, compute a
    // mask that blocks out the two bits that encode the child position of
    // "id" with respect to its parent, then check that the other three
    // children all agree with "mask".
    ulong mask = d.lsb() << 1;
    mask = ~(mask + (mask << 1));
    ulong id_masked = (d.id() & mask);
    return ((a.id() & mask) == id_masked &&
        (b.id() & mask) == id_masked &&
        (c.id() & mask) == id_masked &&
        !d.isFace());
  }

  /// Clears the contents of the cell union and minimizes memory usage.
  void clear() {
    _cellIds.length = 0;
  }

  /**
   * Gives ownership of the vector data to the client without copying, and
   * clears the content of the cell union.  The original data in cell_ids
   * is lost if there was any.
   */
  S2CellId[] release() {
    // vector's rvalue reference constructor does not necessarily leave
    // moved-from value in empty state, so swap instead.
    S2CellId[] cell_ids = _cellIds;
    _cellIds.length = 0;
    return cell_ids;
  }

  /// Convenience methods for accessing the individual cell ids.
  int numCells() const {
    return cast(int)(_cellIds.length);
  }

  S2CellId cellId(int i) const {
    return _cellIds[i];
  }

  /// Vector-like methods for accessing the individual cell ids.
  size_t size() const {
    return _cellIds.length;
  }

  bool empty() const {
    return _cellIds.empty();
  }

  S2CellId opIndex(size_t i) const {
    return _cellIds[i];
  }

  struct Iterator {
    const(S2CellId)[] cellIds;
    size_t pos;

    inout(S2CellId) getValue() inout {
      return cellIds[pos];
    }

    void opUnary(string op)()
    if (op == "++") {
      if (pos < cellIds.length) {
        pos++;
      }
    }

    void opUnary(string op)()
    if (op == "--") {
      if (pos > 0) {
        pos--;
      }
    }
  }

  /**
   * Standard begin/end methods, to allow range-based for loops:
   *
   *  for (S2CellId id : cell_union) { ... }
   */
  Iterator begin() const {
    return Iterator(_cellIds, 0);
  }

  Iterator end() const {
    return Iterator(_cellIds, _cellIds.length);
  }

  /// Direct access to the underlying vector for STL algorithms.
  inout(S2CellId[]) cellIds() inout {
    return _cellIds;
  }

  /**
   * Returns true if the cell union is valid, meaning that the S2CellIds are
   * valid, non-overlapping, and sorted in increasing order.
   */
  bool isValid() const {
    if (numCells() > 0 && !cellId(0).isValid()) {
      return false;
    }
    for (int i = 1; i < numCells(); ++i) {
      if (!cellId(i).isValid()) {
        return false;
      }
      if (cellId(i - 1).rangeMax() >= cellId(i).rangeMin()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if the cell union is normalized, meaning that it is
   * satisfies IsValid() and that no four cells have a common parent.
   * Certain operations such as Contains(S2CellUnion) will return a different
   * result if the cell union is not normalized.
   */
  bool isNormalized() const {
    if (numCells() > 0 && !cellId(0).isValid()) {
      return false;
    }
    for (int i = 1; i < numCells(); ++i) {
      if (!cellId(i).isValid()) {
        return false;
      }
      if (cellId(i - 1).rangeMax() >= cellId(i).rangeMin()) {
        return false;
      }
      if (i >= 3 && areSiblings(cellId(i - 3), cellId(i - 2), cellId(i - 1), cellId(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Normalizes the cell union by discarding cells that are contained by other
   * cells, replacing groups of 4 child cells by their parent cell whenever
   * possible, and sorting all the cell ids in increasing order.
   *
   * Returns true if the number of cells was reduced.
   * TODO(ericv): Change this method to return void.
   */
  bool normalize() {
    return normalize(_cellIds);
  }

  /**
   * Replaces "output" with an expanded version of the cell union where any
   * cells whose level is less than "min_level" or where (level - min_level)
   * is not a multiple of "level_mod" are replaced by their children, until
   * either both of these conditions are satisfied or the maximum level is
   * reached.
   *
   * This method allows a covering generated by S2RegionCoverer using
   * min_level() or level_mod() constraints to be stored as a normalized cell
   * union (which allows various geometric computations to be done) and then
   * converted back to the original list of cell ids that satisfies the
   * desired constraints.
   */
  void denormalize(int min_level, int level_mod, ref S2CellId[] output) const {
    denormalize(_cellIds, min_level, level_mod, output);
  }

  /**
   * If there are more than "excess" elements of the cell_ids() vector that
   * are allocated but unused, reallocates the array to eliminate the excess
   * space.  This reduces memory usage when many cell unions need to be held
   * in memory at once.
   */
  void pack(int excess = 0) {
    // TODO: For now, let the D runtime manage this.
    return;
  }

  /**
   * Returns true if the cell union contains the given cell id.  Containment
   * is defined with respect to regions, e.g. a cell contains its 4 children.
   * This is a fast operation (logarithmic in the size of the cell union).
   *
   * CAVEAT: If you have constructed a non-normalized S2CellUnion using
   * FromVerbatim, note that groups of 4 child cells are *not* considered to
   * contain their parent cell.  To get this behavior you must use one of the
   * other constructors or call Normalize() explicitly.
   */
  bool contains(S2CellId id) const {
    // This is an exact test.  Each cell occupies a linear span of the S2
    // space-filling curve, and the cell id is simply the position at the center
    // of this span.  The cell union ids are sorted in increasing order along
    // the space-filling curve.  So we simply find the pair of cell ids that
    // surround the given cell id (using binary search).  There is containment
    // if and only if one of these two cell ids contains this cell.

    // Split the cellIds into a ranges [0] = less, [1] = equal, [2] = greater
    auto trisectRanges = assumeSorted(_cellIds[]).trisect(id);
    return !trisectRanges[1].empty
        || !trisectRanges[2].empty && trisectRanges[2].front.rangeMin() <= id
        || !trisectRanges[0].empty && trisectRanges[0].back.rangeMax() >= id;
  }

  /**
   * Returns true if the cell union intersects the given cell id.
   * This is a fast operation (logarithmic in the size of the cell union).
   */
  bool intersects(S2CellId id) const {
    // This is an exact test; see the comments for Contains() above.
    auto trisectRanges = assumeSorted(_cellIds[]).trisect(id);
    return !trisectRanges[1].empty
        || !trisectRanges[2].empty && trisectRanges[2].front.rangeMin() <= id.rangeMax()
        || !trisectRanges[0].empty && trisectRanges[0].back.rangeMax() >= id.rangeMin();
  }

  // Returns true if this cell union contains the given other cell union.
  //
  // CAVEAT: If you have constructed a non-normalized S2CellUnion using
  // FromVerbatim, note that groups of 4 child cells are *not* considered to
  // contain their parent cell.  To get this behavior you must use one of the
  // other constructors or call Normalize() explicitly.
  bool contains(in S2CellUnion y) const {
    // TODO(ericv): A divide-and-conquer or alternating-skip-search
    // approach may be sigificantly faster in both the average and worst case.

    foreach (S2CellId y_id; y._cellIds) {
      if (!contains(y_id)) return false;
    }
    return true;
  }

  // Returns true if this cell union intersects the given other cell union.
  bool intersects(in S2CellUnion y) const {
    // TODO(ericv): A divide-and-conquer or alternating-skip-search
    // approach may be sigificantly faster in both the average and worst case.

    foreach (S2CellId y_id; y._cellIds) {
      if (intersects(y_id)) return true;
    }
    return false;
  }

  // Returns the union of the two given cell unions.
  S2CellUnion unite(in S2CellUnion y) const {
    return new S2CellUnion(array(chain(_cellIds, y._cellIds)));
  }

  // Returns the intersection of the two given cell unions.
  S2CellUnion intersect(in S2CellUnion y) const
  out (result) {
    // The output is normalized as long as at least one input is normalized.
    assert(result.isNormalized() || (!isNormalized() && !y.isNormalized()));
  } body {
    auto result = new S2CellUnion();
    getIntersection(_cellIds, y._cellIds, result._cellIds);
    return result;
  }

  // Specialized version of GetIntersection() that returns the intersection of
  // a cell union with an S2CellId.  This can be useful for splitting a cell
  // union into pieces.
  S2CellUnion intersect(S2CellId id) const
  out (result) {
    assert(result.isNormalized() || !isNormalized());
  } body {
    auto result = new S2CellUnion();
    if (contains(id)) {
      result._cellIds ~= id;
    } else {
      auto ranges = assumeSorted(_cellIds).trisect(id.rangeMin());
      // Get the range of values >= id.rangeMin().
      auto geCells = chain(ranges[1], ranges[2]);
      S2CellId id_max = id.rangeMax();
      while (!geCells.empty() && geCells.front() <= id_max) {
        result._cellIds ~= geCells.front();
        geCells.popFront();
      }
    }
    return result;
  }

  /// Returns the difference of the two given cell unions.
  S2CellUnion difference(in S2CellUnion y) const
  out (result) {
    // The output is normalized as long as the first argument is normalized.
    assert(result.isNormalized() || !isNormalized());
  } body {
    // TODO(ericv): this is approximately O(N*log(N)), but could probably
    // use similar techniques as GetIntersection() to be more efficient.

    auto result = new S2CellUnion();
    for (auto iter = begin(); iter != end(); iter++) {
      getDifferenceInternal(iter.getValue(), y, result._cellIds);
    }
    return result;
  }

  private static void getDifferenceInternal(
      S2CellId cell, in S2CellUnion y, ref S2CellId[] cell_ids) {
    // Add the difference between cell and y to cell_ids.
    // If they intersect but the difference is non-empty, divide and conquer.
    if (!y.intersects(cell)) {
      cell_ids ~= cell;
    } else if (!y.contains(cell)) {
      S2CellId child = cell.childBegin();
      for (int i = 0; ; ++i) {
        getDifferenceInternal(child, y, cell_ids);
        if (i == 3) break;  // Avoid unnecessary next() computation.
        child = child.next();
      }
    }
  }

  /**
   * Expands the cell union by adding a buffer of cells at "expand_level"
   * around the union boundary.
   *
   * For each cell "c" in the union, we add all neighboring cells at level
   * "expand_level" that are adjacent to "c".  Note that there can be many
   * such cells if "c" is large compared to "expand_level".  If "c" is smaller
   * than "expand_level", we first add the parent of "c" at "expand_level" and
   * then add all the neighbors of that cell.
   *
   * Note that the size of the output is exponential in "expand_level".  For
   * example, if expand_level == 20 and the input has a cell at level 10,
   * there will be on the order of 4000 adjacent cells in the output.  For
   * most applications the Expand(min_radius, max_level_diff) method below is
   * easier to use.
   */
  void expand(int expand_level) {
    S2CellId[] output;
    ulong level_lsb = S2CellId.lsbForLevel(expand_level);
    for (int i = numCells(); --i >= 0; ) {
      S2CellId id = cellId(i);
      if (id.lsb() < level_lsb) {
        id = id.parent(expand_level);
        // Optimization: skip over any cells contained by this one.  This is
        // especially important when very small regions are being expanded.
        while (i > 0 && id.contains(cellId(i - 1))) --i;
      }
      output ~= id;
      id.appendAllNeighbors(expand_level, output);
    }
    initialize(output);
  }

  /**
   * Expands the cell union such that it contains all points whose distance to
   * the cell union is at most "min_radius", but do not use cells that are
   * more than "max_level_diff" levels higher than the largest cell in the
   * input.  The second parameter controls the tradeoff between accuracy and
   * output size when a large region is being expanded by a small amount
   * (e.g. expanding Canada by 1km).  For example, if max_level_diff == 4 the
   * region will always be expanded by approximately 1/16 the width of its
   * largest cell.  Note that in the worst case, the number of cells in the
   * output can be up to 4 * (1 + 2 ** max_level_diff) times larger than the
   * number of cells in the input.
   */
  void expand(S1Angle min_radius, int max_level_diff) {
    import s2.s2metrics : MIN_WIDTH;

    int min_level = S2CellId.MAX_LEVEL;
    foreach (S2CellId id; _cellIds) {
      min_level = min(min_level, id.level());
    }
    // Find the maximum level such that all cells are at least "min_radius" wide.
    int radius_level = MIN_WIDTH.getLevelForMinValue(min_radius.radians());
    if (radius_level == 0 && min_radius.radians() > MIN_WIDTH.getValue(0)) {
      // The requested expansion is greater than the width of a face cell.
      // The easiest way to handle this is to expand twice.
      expand(0);
    }
    expand(min(min_level + max_level_diff, radius_level));
  }

  /**
   * The number of leaf cells covered by the union.
   * This will be no more than 6*2^60 for the whole sphere.
   */
  ulong leafCellsCovered() const {
    ulong num_leaves = 0;
    foreach (S2CellId id; _cellIds) {
      int inverted_level = S2CellId.MAX_LEVEL - id.level();
      num_leaves += (1uL << (inverted_level << 1));
    }
    return num_leaves;
  }

  /**
   * Approximates this cell union's area in steradians by summing the average
   * area of each contained cell's average area, using the AverageArea method
   * from the S2Cell class.  This is equivalent to the number of leaves covered,
   * multiplied by the average area of a leaf.  Note that AverageArea does not
   * take into account distortion of cell, and thus may be off by up to a
   * factor of up to 1.7.
   *
   * NOTE: Since this is proportional to LeafCellsCovered(), it is
   * always better to use that function if all you care about is
   * the relative average area between objects.
   */
  double averageBasedArea() const {
    return S2Cell.averageArea(S2CellId.MAX_LEVEL) * leafCellsCovered();
  }

  // Calculates this cell union's area in steradians by summing the approximate
  // area for each contained cell, using the ApproxArea method from the S2Cell
  // class.
  double approxArea() const {
    double area = 0;
    foreach (S2CellId id; _cellIds) {
      area += (new S2Cell(id)).approxArea();
    }
    return area;
  }

  /**
   * Calculates this cell union's area in steradians by summing the exact area
   * for each contained cell, using the Exact method from the S2Cell class.
   */
  double ExactArea() const {
    double area = 0;
    foreach (S2CellId id; _cellIds) {
      area += (new S2Cell(id)).exactArea();
    }
    return area;
  }

  /// Return true if two cell unions are identical.
  override
  bool opEquals(Object o) {
    S2CellUnion y = cast(S2CellUnion) o;
    if (y is null) return false;
    return cellIds() == y.cellIds();
  }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  override
  S2Region clone() {
    return new S2CellUnion(_cellIds, VerbatimFlag.VERBATIM);
  }

  override
  S2Cap getCapBound() {
    // Compute the approximate centroid of the region.  This won't produce the
    // bounding cap of minimal area, but it should be close enough.
    if (_cellIds.length == 0) {
      return S2Cap.empty();
    }

    auto centroid = S2Point(0, 0, 0);
    foreach (S2CellId id; _cellIds) {
      double area = S2Cell.averageArea(id.level());
      centroid += area * id.toS2Point();
    }
    if (centroid == S2Point(0, 0, 0)) {
      centroid = S2Point(1, 0, 0);
    } else {
      centroid = centroid.normalize();
    }

    // Use the centroid as the cap axis, and expand the cap angle so that it
    // contains the bounding caps of all the individual cells.  Note that it is
    // *not* sufficient to just bound all the cell vertices because the bounding
    // cap may be concave (i.e. cover more than one hemisphere).
    S2Cap cap = S2Cap.fromPoint(centroid);
    for (Iterator iter = begin(); iter != end(); iter++) {
      S2CellId id = iter.getValue();
      cap.addCap((new S2Cell(id)).getCapBound());
    }
    return cap;
  }

  override
  S2LatLngRect getRectBound() {
    S2LatLngRect bound = S2LatLngRect.empty();
    for (Iterator iter = begin(); iter != end(); iter++) {
      S2CellId id = iter.getValue();
      bound = bound.unite((new S2Cell(id)).getRectBound());
    }
    return bound;
  }

  override
  void getCellUnionBound(out S2CellId[] cellIds) {
    cellIds = _cellIds.dup;
  }

  /// This is a fast operation (logarithmic in the size of the cell union).
  override
  bool contains(in S2Cell cell) const {
    return contains(cell.id());
  }

  /// This is a fast operation (logarithmic in the size of the cell union).
  override
  bool mayIntersect(in S2Cell cell) const {
    return intersects(cell.id());
  }

  /**
   * The point 'p' does not need to be normalized.
   * This is a fast operation (logarithmic in the size of the cell union).
   */
  override
  bool contains(in S2Point p) const {
    return contains(S2CellId(p));
  }

  // TODO: Uncomment when support for encode/decode is added.
  /**
   * Appends a serialized representation of the S2CellUnion to "encoder".
   *
   * REQUIRES: "encoder" uses the default constructor, so that its buffer
   *           can be enlarged as necessary by calling Ensure(int).
   */
  //void Encode(Encoder* const encoder) const;

  /// Decodes an S2CellUnion encoded with Encode().  Returns true on success.
  //bool Decode(Decoder* const decoder);

  ////////////////////////////////////////////////////////////////////////
  // Static methods intended for high-performance clients that prefer to
  // manage their own storage.

  /**
   * Like Normalize(), but works with a vector of S2CellIds.
   * Equivalent to:
   *   *cell_ids = S2CellUnion(std::move(*cell_ids)).Release();
   */
  static bool normalize(ref S2CellId[] ids) {
    // Optimize the representation by discarding cells contained by other cells,
    // and looking for cases where all subcells of a parent cell are present.
    ids = array(sort(ids));
    int output = 0;
    foreach (id; ids) {
      // Check whether this cell is contained by the previous cell.
      if (output > 0 && ids[output - 1].contains(id)) continue;

      // Discard any previous cells contained by this cell.
      while (output > 0 && id.contains(ids[output - 1])) --output;

      // Check whether the last 3 elements plus "id" can be collapsed into a
      // single parent cell.
      while (output >= 3 && areSiblings(ids[output - 3], ids[output - 2], ids[output - 1], id)) {
        // Replace four children by their parent cell.
        id = id.parent();
        output -= 3;
      }
      ids[output++] = id;
    }
    if (ids.length == output) return false;
    ids.length = output;
    return true;
  }


  // Like Denormalize(), but works with a vector of S2CellIds.
  // REQUIRES: out != &in
  static void denormalize(
      in S2CellId[] input, int min_level, int level_mod, out S2CellId[] output)
  in {
    assert(min_level >= 0);
    assert(min_level <=  S2CellId.MAX_LEVEL);
    assert(level_mod >= 1);
    assert(level_mod <= 3);
    assert(input.length == 0 || output !is input);
  } body {
    output.reserve(input.length);
    foreach (S2CellId id; input) {
      int level = id.level();
      int new_level = max(min_level, level);
      if (level_mod > 1) {
        // Round up so that (new_level - min_level) is a multiple of level_mod.
        // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
        new_level += (S2CellId.MAX_LEVEL - (new_level - min_level)) % level_mod;
        new_level = min(S2CellId.MAX_LEVEL, new_level);
      }
      if (new_level == level) {
        output ~= id;
      } else {
        S2CellId end = id.childEnd(new_level);
        for (id = id.childBegin(new_level); id != end; id = id.next()) {
          output ~= id;
        }
      }
    }
  }

  /**
   * Like GetIntersection(), but works directly with vectors of S2CellIds,
   * Equivalent to:
   *
   *    *out = S2CellUnion(x).Intersection(S2CellUnion(y)).Release()
   *
   * except that this method has slightly more relaxed normalization
   * requirements: the input vectors may contain groups of 4 child cells that
   * all have the same parent.  (In a normalized S2CellUnion, such groups are
   * always replaced by the parent cell.)
   */
  static void getIntersection(in S2CellId[] x, in S2CellId[] y, out S2CellId[] output)
  in {
    import std.conv : to;
    assert(output.length == 0 || output !is x);
    assert(output.length == 0 || output !is y);
    assert(isSorted(x), "x is unsorted: x=" ~ to!string(x));
    assert(isSorted(y), "y is unsorted: y=" ~ to!string(y));
  } out {
    // The output is generated in sorted order.
    assert(isSorted(output));
  } body {
    // This is a fairly efficient calculation that uses binary search to skip
    // over sections of both input vectors.  It takes logarithmic time if all the
    // cells of "x" come before or after all the cells of "y" in S2CellId order.
    size_t xPos = 0;
    size_t yPos = 0;
    while (xPos < x.length && yPos < y.length) {
      S2CellId xPosMin = x[xPos].rangeMin();
      S2CellId yPosMin = y[yPos].rangeMin();
      if (xPosMin > yPosMin) {
        // Either x[xPos].contains(*i) or the two cells are disjoint.
        if (x[xPos] <= y[yPos].rangeMax()) {
          output ~= x[xPos++];
        } else {
          // Advance "j" to the first cell possibly contained by *i.
          yPos++;
          long skipAmount = countUntil!(a => a >= xPosMin)(y[yPos .. $]);
          if (skipAmount == -1) {
            yPos = y.length;
          } else {
            yPos += skipAmount;
          }
          // The previous cell *(j-1) may now contain *i.
          if (x[xPos] <= y[yPos - 1].rangeMax()) {
            yPos--;
          }
        }
      } else if (yPosMin > xPosMin) {
        // Identical to the code above with "i" and "j" reversed.
        if (y[yPos] <= x[xPos].rangeMax()) {
          output ~= y[yPos++];
        } else {
          xPos++;
          long skipAmount = countUntil!(a => a >= yPosMin)(x[xPos .. $]);
          if (skipAmount == -1) {
            xPos = x.length;
          } else {
            xPos += skipAmount;
          }
          if (y[yPos] <= x[xPos - 1].rangeMax()) {
            --xPos;
          }
        }
      } else {
        // "i" and "j" have the same range_min(), so one contains the other.
        if (x[xPos] < y[yPos])
          output ~= x[xPos++];
        else
          output ~= y[yPos++];
      }
    }
  }

package:
  // Internal constructor that does not check "cell_ids" for validity.
  this(in S2CellId[] cell_ids, VerbatimFlag verbatim) {
    _cellIds = cell_ids.dup;
  }

  // Internal constructor that does not check "cell_ids" for validity.
  enum VerbatimFlag {
    VERBATIM
  }

private:
  // Converts a vector of uint64 to a vector of S2CellIds.
  static S2CellId[] toS2CellIds(in ulong[] ids) {
    S2CellId[] cell_ids;
    cell_ids.reserve(ids.length);
    foreach (id; ids) cell_ids ~= S2CellId(id);
    return cell_ids;
  }

  S2CellId[] _cellIds;
}
