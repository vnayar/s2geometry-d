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

module s2.s2region_coverer;

import s2.logger;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2point;
import s2.s2region;

import std.algorithm : isSorted, min, max;
import std.array : array;
import std.container : BinaryHeap, heapify;
import std.range : assumeSorted, chain, SortedRange, takeOne;

/**
 * An S2RegionCoverer is a class that allows arbitrary regions to be
 * approximated as unions of cells (S2CellUnion).  This is useful for
 * implementing various sorts of search and precomputation operations.
 *
 * Typical usage:
 *
 * S2RegionCoverer::Options options;
 * options.set_max_cells(5);
 * S2RegionCoverer coverer(options);
 * S2Cap cap(center, radius);
 * S2CellUnion covering = coverer.GetCovering(cap);
 *
 * This yields a vector of at most 5 cells that is guaranteed to cover the
 * given cap (a disc-shaped region on the sphere).
 *
 * The approximation algorithm is not optimal but does a pretty good job in
 * practice.  The output does not always use the maximum number of cells
 * allowed, both because this would not always yield a better approximation,
 * and because max_cells() is a limit on how much work is done exploring the
 * possible covering as well as a limit on the final output size.
 *
 * Because it is an approximation algorithm, one should not rely on the
 * stability of the output.  In particular, the output of the covering algorithm
 * may change across different versions of the library.
 *
 * One can also generate interior coverings, which are sets of cells which
 * are entirely contained within a region.  Interior coverings can be
 * empty, even for non-empty regions, if there are no cells that satisfy
 * the provided constraints and are contained by the region.  Note that for
 * performance reasons, it is wise to specify a max_level when computing
 * interior coverings - otherwise for regions with small or zero area, the
 * algorithm may spend a lot of time subdividing cells all the way to leaf
 * level to try to find contained cells.
 */
class S2RegionCoverer {
public:

  static class Options {
  public:
    /**
     * Sets the desired maximum number of cells in the approximation.  Note
     * the following:
     *
     *  - For any setting of max_cells(), up to 6 cells may be returned if
     *    that is the minimum number required (e.g. if the region intersects
     *    all six cube faces).  Even for very tiny regions, up to 3 cells may
     *    be returned if they happen to be located at the intersection of
     *    three cube faces.
     *
     *  - min_level() takes priority over max_cells(), i.e. cells below the
     *    given level will never be used even if this causes a large number of
     *    cells to be returned.
     *
     *  - If max_cells() is less than 4, the area of the covering may be
     *    arbitrarily large compared to the area of the original region even
     *    if the region is convex (e.g. an S2Cap or S2LatLngRect).
     *
     * Accuracy is measured by dividing the area of the covering by the area
     * of the original region.  The following table shows the median and worst
     * case values for this area ratio on a test case consisting of 100,000
     * spherical caps of random size (generated using s2region_coverer_test):
     *
     *   max_cells:        3      4     5     6     8    12    20   100   1000
     *   median ratio:  5.33   3.32  2.73  2.34  1.98  1.66  1.42  1.11   1.01
     *   worst case:  215518  14.41  9.72  5.26  3.91  2.75  1.92  1.20   1.02
     *
     * The default value of 8 gives a reasonable tradeoff between the number
     * of cells used and the accuracy of the approximation.
     *
     * DEFAULT: kDefaultMaxCells
     */
    enum int DEFAULT_MAX_CELLS = 8;

    int maxCells() const {
      return _maxCells;
    }

    void setMaxCells(int max_cells) {
        _maxCells = max_cells;
    }

    /**
     * Sets the minimum and maximum cell levels to be used.  The default is to
     * use all cell levels.
     *
     * To find the cell level corresponding to a given physical distance, use
     * the S2Cell metrics defined in s2metrics.h.  For example, to find the
     * cell level that corresponds to an average edge length of 10km, use:
     *
     *   int level =
     *       S2::kAvgEdge.GetClosestLevel(S2Earth::KmToRadians(length_km));
     *
     * Note that min_level() takes priority over max_cells(), i.e. cells below
     * the given level will never be used even if this causes a large number
     * of cells to be returned.  (This doesn't apply to interior coverings,
     * since interior coverings make no completeness guarantees -- the result
     * is simply a set of cells that covers as much of the interior as
     * possible while satisfying the given restrictions.)
     *
     * REQUIRES: max_level() >= min_level()
     * DEFAULT: 0
     */
    int minLevel() const {
      return _minLevel;
    }

    void setMinLevel(int min_level)
    in {
      assert(min_level >=  0);
      assert(min_level <= S2CellId.MAX_LEVEL);
    } do {
      _minLevel = max(0, min(S2CellId.MAX_LEVEL, min_level));
    }

    // DEFAULT: S2CellId::kMaxLevel
    int maxLevel() const {
      return _maxLevel;
    }

    void setMaxLevel(int max_level)
    in {
      assert(max_level >= 0);
      assert(max_level <= S2CellId.MAX_LEVEL);
    } do {
      _maxLevel = max(0, min(S2CellId.MAX_LEVEL, max_level));
    }

    /// Convenience function that sets both the maximum and minimum cell levels.
    void setFixedLevel(int level) {
      setMinLevel(level);
      setMaxLevel(level);
    }

    /**
     * If specified, then only cells where (level - min_level) is a multiple
     * of "level_mod" will be used (default 1).  This effectively allows the
     * branching factor of the S2CellId hierarchy to be increased.  Currently
     * the only parameter values allowed are 1, 2, or 3, corresponding to
     * branching factors of 4, 16, and 64 respectively.
     *
     * DEFAULT: 1
     */
    int levelMod() const {
      return _levelMod;
    }

    void setLevelMod(int level_mod)
    in {
      assert(level_mod >= 1);
      assert(level_mod <= 3);
    } do {
      _levelMod = max(1, min(3, level_mod));
    }

    /**
     * Convenience function that returns the maximum level such that
     *
     *   (level <= max_level()) && (level - min_level()) % level_mod() == 0.
     *
     * This is the maximum level that will actually be used in coverings.
     */
    int trueMaxLevel() const {
      if (_levelMod == 1) return _maxLevel;
      return _maxLevel - (_maxLevel - _minLevel) % _levelMod;
    }

    override
    string toString() const {
      import std.format : format;
      return format(
          "[maxCells=%d, minLevel=%d, maxLevel=%d, levelMod=%d, trueMaxLevel=%d]",
          _maxCells, _minLevel, _maxLevel, _levelMod, trueMaxLevel());
    }

  protected:
    int _maxCells = DEFAULT_MAX_CELLS;
    int _minLevel = 0;
    int _maxLevel = S2CellId.MAX_LEVEL;
    int _levelMod = 1;
  }

  this() {
    QueueEntry[] queueEntries;
    _pq = heapify(queueEntries);
    _options = new Options();
  }

  // Constructs an S2RegionCoverer with the given options.
  this(Options options) {
    this();
    _options = options;
  }

  // Returns the current options.  Options can be modifed between calls.
  const(Options) options() const {
    return _options;
  }

  Options mutableOptions() {
    return _options;
  }

  /**
   * Returns an S2CellUnion that covers (GetCovering) or is contained within
   * (GetInteriorCovering) the given region and satisfies the current options.
   *
   * Note that if options().min_level() > 0 or options().level_mod() > 1, the
   * by definition the S2CellUnion may not be normalized, i.e. there may be
   * groups of four child cells that can be replaced by their parent cell.
   */
  S2CellUnion getCovering(S2Region region) {
    _interiorCovering = false;
    getCoveringInternal(region);
    auto r = S2CellUnion.fromVerbatim(_result);
    _result = null;
    return r;
  }

  S2CellUnion getInteriorCovering(S2Region region) {
    _interiorCovering = true;
    getCoveringInternal(region);
    auto r = S2CellUnion.fromVerbatim(_result);
    _result = null;
    return r;
  }

  /**
   * Like the methods above, but works directly with a vector of S2CellIds.
   * This version can be more efficient when this method is called many times,
   * since it does not require allocating a new vector on each call.
   */
  void getCovering(S2Region region, ref S2CellId[] covering) {
    _interiorCovering = false;
    getCoveringInternal(region);
    covering = _result;
    _result = null;
  }

  void getInteriorCovering(S2Region region, ref S2CellId[] interior) {
    _interiorCovering = true;
    getCoveringInternal(region);
    interior = _result;
    _result = null;
  }

  /**
   * Like GetCovering(), except that this method is much faster and the
   * coverings are not as tight.  All of the usual parameters are respected
   * (max_cells, min_level, max_level, and level_mod), except that the
   * implementation makes no attempt to take advantage of large values of
   * max_cells().  (A small number of cells will always be returned.)
   *
   * This function is useful as a starting point for algorithms that
   * recursively subdivide cells.
   */
  void getFastCovering(S2Region region, ref S2CellId[] covering) {
    region.getCellUnionBound(covering);
    canonicalizeCovering(covering);
  }

  /**
   * Given a connected region and a starting point, return a set of cells at
   * the given level that cover the region.
   *
   * Note that this method is *not* faster than the regular GetCovering()
   * method for most region types, such as S2Cap or S2Polygon, and in fact it
   * can be much slower when the output consists of a large number of cells.
   * Currently it can be faster at generating coverings of long narrow regions
   * such as polylines, but this may change in the future, in which case this
   * method will most likely be removed.
   */
  static void getSimpleCovering(
      S2Region region, in S2Point start, int level, ref S2CellId[] output) {
    return floodFill(region, S2CellId(start).parent(level), output);
  }

  /**
   * Given a region and a starting cell, returns the set of all the
   * edge-connected cells at the same level that intersect "region".
   * The output cells are returned in arbitrary order.
   */
  static void floodFill(S2Region region, S2CellId start, out S2CellId[] output) {
    bool[S2CellId] all;
    S2CellId[] frontier;
    all[start] = true;
    frontier ~= start;
    while (frontier.length != 0) {
      S2CellId id = frontier[$-1];
      frontier.length--;
      if (!region.mayIntersect(new S2Cell(id))) continue;
      output ~= id;

      S2CellId[4] neighbors;
      id.getEdgeNeighbors(neighbors);
      for (int edge = 0; edge < 4; ++edge) {
        S2CellId nbr = neighbors[edge];
        if (nbr !in all) {
          all[nbr] = true;
          frontier ~= nbr;
        }
      }
    }
  }

  /**
   * Returns true if the given S2CellId vector represents a valid covering
   * that conforms to the current covering parameters.  In particular:
   *
   *  - All S2CellIds must be valid.
   *
   *  - S2CellIds must be sorted and non-overlapping.
   *
   *  - S2CellId levels must satisfy min_level(), max_level(), and level_mod().
   *
   *  - If covering.size() > max_cells(), there must be no two cells with
   *    a common ancestor at min_level() or higher.
   *
   *  - There must be no sequence of cells that could be replaced by an
   *    ancestor (i.e. with level_mod() == 1, the 4 child cells of a parent).
   */
  bool isCanonical(in S2CellUnion covering) const {
    return isCanonical(covering.cellIds());
  }

  bool isCanonical(in S2CellId[] covering) const {
    const int min_level = _options.minLevel();
    const int max_level = _options.trueMaxLevel();
    const int level_mod = _options.levelMod();
    const bool too_many_cells = covering.length > _options.maxCells();
    int same_parent_count = 1;
    S2CellId prev_id = S2CellId.none();
    foreach (const S2CellId id; covering) {
      if (!id.isValid()) {
        return false;
      }

      // Check that the S2CellId level is acceptable.
      const int level = id.level();
      if (level < min_level || level > max_level) {
        return false;
      }
      if (level_mod > 1 && (level - min_level) % level_mod != 0) {
        return false;
      }

      if (prev_id != S2CellId.none()) {
        // Check that cells are sorted and non-overlapping.
        if (prev_id.rangeMax() >= id.rangeMin()) {
          return false;
        }

        // If there are too many cells, check that no pair of adjacent cells
        // could be replaced by an ancestor.
        if (too_many_cells && id.getCommonAncestorLevel(prev_id) >= min_level) {
          return false;
        }

        // Check that there are no sequences of (4 ** level_mod) cells that all
        // have the same parent (considering only multiples of "level_mod").
        int plevel = level - level_mod;
        if (plevel < min_level || level != prev_id.level()
            || id.parent(plevel) != prev_id.parent(plevel)) {
          same_parent_count = 1;
        } else if (++same_parent_count == (1 << (2 * level_mod))) {
          return false;
        }
      }
      prev_id = id;
    }
    return true;
  }

  /**
   * Modify "covering" if necessary so that it conforms to the current
   * covering parameters (max_cells, min_level, max_level, and level_mod).
   * There are no restrictions on the input S2CellIds (they may be unsorted,
   * overlapping, etc).
   */
  S2CellUnion canonicalizeCovering(in S2CellUnion covering) {
    auto ids = covering.cellIds().dup;
    canonicalizeCovering(ids);
    return new S2CellUnion(ids);
  }

  void canonicalizeCovering(ref S2CellId[] covering)
  out {
    assert(isCanonical(covering));
  } do {
    // Note that when the covering parameters have their default values, almost
    // all of the code in this function is skipped.

    // If any cells are too small, or don't satisfy level_mod(), then replace
    // them with ancestors.
    if (_options.maxLevel() < S2CellId.MAX_LEVEL || _options.levelMod() > 1) {
      for (int i = 0; i < covering.length; ++i) {
        S2CellId id = covering[i];
        int level = id.level();
        int new_level = adjustLevel(min(level, _options.maxLevel()));
        if (new_level != level) {
          covering[i] = id.parent(new_level);
        }
      }
    }

    // Sort the cells and simplify them.
    S2CellUnion.normalize(covering);

    // Make sure that the covering satisfies min_level() and level_mod(),
    // possibly at the expense of satisfying max_cells().
    if (_options.minLevel() > 0 || _options.levelMod() > 1) {
      S2CellUnion.denormalize(covering, _options.minLevel(), _options.levelMod(), _result);
      covering = _result;
      _result = null;
    }

    // If there are too many cells and the covering is very large, use the
    // S2RegionCoverer to compute a new covering.  (This avoids possible O(n^2)
    // behavior of the simpler algorithm below.)
    long excess = covering.length - _options.maxCells();
    if (excess <= 0 || isCanonical(covering)) {
      return;
    }
    if (excess * covering.length > 10000) {
      getCovering(new S2CellUnion(covering), covering);
    } else {
      // Repeatedly replace two adjacent cells in S2CellId order by their lowest
      // common ancestor until the number of cells is acceptable.
      while (covering.length > _options.maxCells()) {
        int best_index = -1, best_level = -1;
        for (int i = 0; i + 1 < covering.length; ++i) {
          int level = covering[i].getCommonAncestorLevel(covering[i+1]);
          level = adjustLevel(level);
          if (level > best_level) {
            best_level = level;
            best_index = i;
          }
        }
        if (best_level < _options.minLevel()) break;

        // Replace all cells contained by the new ancestor cell.
        S2CellId id = covering[best_index].parent(best_level);
        replaceCellsWithAncestor(covering, id);

        // Now repeatedly check whether all children of the parent cell are
        // present, in which case we can replace those cells with their parent.
        while (best_level > _options.minLevel()) {
          best_level -= _options.levelMod();
          id = id.parent(best_level);
          if (!containsAllChildren(covering, id)) break;
          replaceCellsWithAncestor(covering, id);
        }
      }
    }
  }


 private:
  class Candidate {
    S2Cell cell;
    bool isTerminal;        // Cell should not be expanded further.
    int numChildren;        // Number of children that intersect the region.
    Candidate[] children;  // Actual size may be 0, 4, 16, or 64 elements.

    override
    string toString() const {
      import std.format : format;
      return format("Candidate[cell=%s, isTerminal=%s, numChildren=%d]",
          cell, isTerminal, numChildren);
    }
  }

  /**
   * If the cell intersects the given region, return a new candidate with no
   * children, otherwise return nullptr.  Also marks the candidate as "terminal"
   * if it should not be expanded further.
   */
  Candidate newCandidate(S2Cell cell) {
    if (!_region.mayIntersect(cell)) return null;

    bool is_terminal = false;
    if (cell.level() >= _options.minLevel()) {
      if (_interiorCovering) {
        if (_region.contains(cell)) {
          is_terminal = true;
        } else if (cell.level() + _options.levelMod() > _options.maxLevel()) {
          return null;
        }
      } else {
        if (cell.level() + _options.levelMod() > _options.maxLevel() || _region.contains(cell)) {
          is_terminal = true;
        }
      }
    }
    size_t children_size = 0;
    if (!is_terminal) {
      children_size = Candidate.sizeof << maxChildrenShift();
    }
    Candidate candidate = new Candidate();
    candidate.cell = cell;
    candidate.isTerminal = is_terminal;
    candidate.numChildren = 0;
    ++_candidatesCreatedCounter;
    return candidate;
  }

  /// Return the log base 2 of the maximum number of children of a candidate.
  int maxChildrenShift() const {
    return 2 * options().levelMod();
  }

  // Process a candidate by either adding it to the result_ vector or
  // expanding its children and inserting it into the priority queue.
  // Passing an argument of nullptr does nothing.
  void addCandidate(Candidate candidate) {
    if (candidate is null) return;

    if (candidate.isTerminal) {
      _result ~= candidate.cell.id();
      return;
    }
    assert(candidate.numChildren == 0);

    // Expand one level at a time until we hit min_level() to ensure that we
    // don't skip over it.
    int num_levels = (candidate.cell.level() < _options.minLevel()) ? 1 : _options.levelMod();
    int num_terminals = expandChildren(candidate, candidate.cell, num_levels);

    if (candidate.numChildren != 0 && !_interiorCovering
        && num_terminals == 1 << maxChildrenShift()
        && candidate.cell.level() >= _options.minLevel()) {
      // Optimization: add the parent cell rather than all of its children.
      // We can't do this for interior coverings, since the children just
      // intersect the region, but may not be contained by it - we need to
      // subdivide them further.
      candidate.isTerminal = true;
      addCandidate(candidate);
    } else {
      // We negate the priority so that smaller absolute priorities are returned
      // first.  The heuristic is designed to refine the largest cells first,
      // since those are where we have the largest potential gain.  Among cells
      // of the same size, we prefer the cells with the fewest children.
      // Finally, among cells with equal numbers of children we prefer those
      // with the smallest number of children that cannot be refined further.
      int priority = -((((candidate.cell.level() << maxChildrenShift())
                  + candidate.numChildren) << maxChildrenShift())
          + num_terminals);
      _pq.insert(QueueEntry(priority, candidate));
      logger.logDebug("Push: ", candidate.cell.id(), " (", priority, ") ");
    }
  }

  /**
   * Populate the children of "candidate" by expanding the given number of
   * levels from the given cell.  Returns the number of children that were
   * marked "terminal".
   */
  int expandChildren(Candidate candidate, const S2Cell cell, int num_levels) {
    num_levels--;
    S2Cell[4] child_cells;
    cell.subdivide(child_cells);
    int num_terminals = 0;
    foreach (int i; 0 .. 4) {
      if (num_levels > 0) {
        if (_region.mayIntersect(child_cells[i])) {
          num_terminals += expandChildren(candidate, child_cells[i], num_levels);
        }
        continue;
      }
      Candidate child = newCandidate(child_cells[i]);
      if (child) {
        candidate.children ~= child;
        candidate.numChildren++;
        if (child.isTerminal) ++num_terminals;
      }
    }
    return num_terminals;
  }

  // Computes a set of initial candidates that cover the given region.
  void getInitialCandidates() {
    // Optimization: start with a small (usually 4 cell) covering of the
    // region's bounding cap.
    auto tmp_coverer = new S2RegionCoverer();
    tmp_coverer.mutableOptions().setMaxCells(min(4, _options.maxCells()));
    tmp_coverer.mutableOptions().setMaxLevel(_options.maxLevel());
    S2CellId[] cells;
    tmp_coverer.getFastCovering(_region, cells);
    adjustCellLevels(cells);
    foreach (S2CellId cell_id; cells) {
      addCandidate(newCandidate(new S2Cell(cell_id)));
    }
  }

  /// Generates a covering and stores it in result_.
  void getCoveringInternal(S2Region region)
  in {
    assert(_pq.length == 0);
    assert(_result.length == 0);
  } out {
    assert(isCanonical(_result));
  } do {
    // Strategy: Start with the 6 faces of the cube.  Discard any
    // that do not intersect the shape.  Then repeatedly choose the
    // largest cell that intersects the shape and subdivide it.
    //
    // result_ contains the cells that will be part of the output, while pq_
    // contains cells that we may still subdivide further.  Cells that are
    // entirely contained within the region are immediately added to the output,
    // while cells that do not intersect the region are immediately discarded.
    // Therefore pq_ only contains cells that partially intersect the region.
    // Candidates are prioritized first according to cell size (larger cells
    // first), then by the number of intersecting children they have (fewest
    // children first), and then by the number of fully contained children
    // (fewest children first).

    _region = region;
    _candidatesCreatedCounter = 0;

    getInitialCandidates();
    while (!_pq.empty()
        && (!_interiorCovering || _result.length < _options.maxCells())) {
      Candidate candidate = _pq.front().candidate;
      _pq.popFront();
      logger.logDebug("Pop: ", candidate.cell.id());
      // For interior coverings we keep subdividing no matter how many children
      // the candidate has.  If we reach max_cells() before expanding all
      // children, we will just use some of them.  For exterior coverings we
      // cannot do this, because the result has to cover the whole region, so
      // all children have to be used.  The (candidate->num_children == 1) case
      // takes care of the situation when we already have more than max_cells()
      // in results (min_level is too high).  Subdividing the candidate with one
      // child does no harm in this case.
      if (_interiorCovering
          || candidate.cell.level() < _options.minLevel()
          || candidate.numChildren == 1
          || (_result.length + _pq.length + candidate.numChildren <= _options.maxCells())) {
        // Expand this candidate into its children.
        for (int i = 0; i < candidate.numChildren; ++i) {
          if (!_interiorCovering || _result.length < _options.maxCells()) {
            addCandidate(candidate.children[i]);
          }
        }
      } else {
        candidate.isTerminal = true;
        addCandidate(candidate);
      }
    }
    logger.logDebug("Created ", _result.length, " cells, ",
        _candidatesCreatedCounter, " candidates created, ",
        _pq.length, " left");
    while (!_pq.empty()) {
      _pq.popFront();
    }
    _region = null;

    // Rather than just returning the raw list of cell ids, we construct a cell
    // union and then denormalize it.  This has the effect of replacing four
    // child cells with their parent whenever this does not violate the covering
    // parameters specified (min_level, level_mod, etc).  This significantly
    // reduces the number of cells returned in many cases, and it is cheap
    // compared to computing the covering in the first place.
    S2CellUnion.normalize(_result);
    if (_options.minLevel() > 0 || _options.levelMod() > 1) {
      auto result_copy = _result;
      S2CellUnion.denormalize(result_copy, _options.minLevel(), _options.levelMod(), _result);
    }
  }

  /**
   * If level > min_level(), then reduce "level" if necessary so that it also
   * satisfies level_mod().  Levels smaller than min_level() are not affected
   * (since cells at these levels are eventually expanded).
   */
  int adjustLevel(int level) const {
    if (_options.levelMod() > 1 && level > _options.minLevel()) {
      level -= (level - _options.minLevel()) % _options.levelMod();
    }
    return level;
  }

  // Ensure that all cells with level > min_level() also satisfy level_mod(),
  // by replacing them with an ancestor if necessary.  Cell levels smaller
  // than min_level() are not modified (see AdjustLevel).  The output is
  // then normalized to ensure that no redundant cells are present.
  void adjustCellLevels(ref S2CellId[] cells) const
  in {
    assert(isSorted(cells));
  } do {
    if (_options.levelMod() == 1) return;

    int output = 0;
    foreach (S2CellId id; cells) {
      int level = id.level();
      int new_level = adjustLevel(level);
      if (new_level != level) id = id.parent(new_level);
      if (output > 0 && cells[output-1].contains(id)) continue;
      while (output > 0 && id.contains(cells[output-1])) --output;
      cells[output++] = id;
    }
    cells.length = output;
  }

  // Returns true if "covering" contains all children of "id" at level
  // (id.level() + options_.level_mod()).
  bool containsAllChildren(in S2CellId[] covering, S2CellId id) const {
    auto ranges = assumeSorted(covering).trisect(id.rangeMin());
    auto geRange = chain(ranges[1], ranges[2]);
    int level = id.level() + _options.levelMod();
    for (S2CellId child = id.childBegin(level);
         child != id.childEnd(level);
         geRange.popFront(), child = child.next()) {
      if (geRange.empty() || geRange.front != child) return false;
    }
    return true;
  }

  // Replaces all descendants of "id" in "covering" with "id".
  // REQUIRES: "covering" contains at least one descendant of "id".
  void replaceCellsWithAncestor(ref S2CellId[] covering, S2CellId id) const {
    auto ltRange = assumeSorted(covering).lowerBound(id.rangeMin());
    auto gtRange = assumeSorted(covering).upperBound(id.rangeMax());
    auto cutRange = chain(ltRange, [id], gtRange);
    covering = array(cutRange);
  }

  Options _options;

  // We save a temporary copy of the pointer passed to GetCovering() in order
  // to avoid passing this parameter around internally.  It is only used (and
  // only valid) for the duration of a single GetCovering() call.
  S2Region _region = null;

  // The set of S2CellIds that have been added to the covering so far.
  S2CellId[] _result;

  // We keep the candidates in a priority queue.  We specify a vector to hold
  // the queue entries since for some reason priority_queue<> uses a deque by
  // default.  We define our own own comparison function on QueueEntries in
  // order to make the results deterministic.  (Using the default
  // less<QueueEntry>, entries of equal priority would be sorted according to
  // the memory address of the candidate.)

  struct QueueEntry {
    int priority;
    Candidate candidate;

    int opCmp(ref in QueueEntry b) {
      return priority - b.priority;
    }
  }
  alias CandidateQueue = BinaryHeap!(QueueEntry[]);
  CandidateQueue _pq;

  // True if we're computing an interior covering.
  bool _interiorCovering;

  // Counter of number of candidates created, for performance evaluation.
  int _candidatesCreatedCounter;
}
