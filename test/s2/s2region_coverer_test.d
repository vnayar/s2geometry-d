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

module s2.s2region_coverer_test;

import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2latlng;
import s2.s2region;
import s2.s2region_coverer;
import s2.s2testing;
import s2.util.math.s2const;

import fluent.asserts;

import math = std.math;
import std.algorithm : min, max, sort;
import std.array;
import std.container.binaryheap;
import std.range;
import std.stdio;

// List of values to use for 'max_cells'.
enum int[] max_cells = [4, 8];

// Number of random caps to try for each max_cells value.
enum int iters = 1000;

@("RandomCells") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(1);
  auto coverer = new S2RegionCoverer(options);
  // Test random cell ids at all levels.
  for (int i = 0; i < 10_000; ++i) {
    S2CellId id = S2Testing.getRandomCellId();
    //SCOPED_TRACE(StrCat("Iteration ", i, ", cell ID token ", id.ToToken()));
    S2CellId[] covering = coverer.getCovering(new S2Cell(id)).release();
    Assert.equal(covering.length, 1);
    Assert.equal(covering[0], id);
  }
}

void checkCovering(in S2RegionCoverer.Options options,
    S2Region region, in S2CellId[] covering, bool interior) {
  // Keep track of how many cells have the same options.min_level() ancestor.
  int[S2CellId] min_level_cells;
  foreach (S2CellId cell_id; covering) {
    int level = cell_id.level();
    Assert.notGreaterThan(options.minLevel(), level);
    Assert.notLessThan(options.maxLevel(), level);
    Assert.equal((level - options.minLevel()) % options.levelMod(), 0);
    min_level_cells[cell_id.parent(options.minLevel())] += 1;
  }
  if (covering.length > options.maxCells()) {
    // If the covering has more than the requested number of cells, then check
    // that the cell count cannot be reduced by using the parent of some cell.
    foreach (S2CellId id, int count; min_level_cells) {
      Assert.equal(count, 1);
    }
  }
  if (interior) {
    foreach (S2CellId cell_id; covering) {
      Assert.equal(region.contains(new S2Cell(cell_id)), true);
    }
  } else {
    auto cell_union = new S2CellUnion(covering);
    S2Testing.checkCovering(region, cell_union, true);
  }
}

@("RandomCaps")
unittest {
  static const int kMaxLevel = S2CellId.MAX_LEVEL;
  auto options = new S2RegionCoverer.Options();
  for (int i = 0; i < 1000; ++i) {
    do {
      options.setMinLevel(S2Testing.rnd.uniform(kMaxLevel + 1));
      options.setMaxLevel(S2Testing.rnd.uniform(kMaxLevel + 1));
    } while (options.minLevel() > options.maxLevel());
    options.setMaxCells(S2Testing.rnd.skewed(10));
    options.setLevelMod(1 + S2Testing.rnd.uniform(3));
    double max_area =  min(
        4 * M_PI,
        (3 * options.maxCells() + 1) * S2Cell.averageArea(options.minLevel()));
    S2Cap cap = S2Testing.getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), max_area);
    auto coverer = new S2RegionCoverer(options);
    S2CellId[] covering, interior;
    coverer.getCovering(cap, covering);
    checkCovering(options, cap, covering, false);
    coverer.getInteriorCovering(cap, interior);
    checkCovering(options, cap, interior, true);

    // Check that GetCovering is deterministic.
    S2CellId[] covering2;
    coverer.getCovering(cap, covering2);
    Assert.equal(covering, covering2);

    // Also check S2CellUnion::Denormalize().  The denormalized covering
    // may still be different and smaller than "covering" because
    // S2RegionCoverer does not guarantee that it will not output all four
    // children of the same parent.
    auto cells = new S2CellUnion(covering);
    S2CellId[] denormalized;
    cells.denormalize(options.minLevel(), options.levelMod(), denormalized);
    checkCovering(options, cap, denormalized, false);
  }
}

@("SimpleCoverings") unittest {
  static const int kMaxLevel = S2CellId.MAX_LEVEL;
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(int.max);
  foreach (i; 0 .. 1000) {
    int level = S2Testing.rnd.uniform(kMaxLevel + 1);
    options.setMinLevel(level);
    options.setMaxLevel(level);
    double max_area =  min(4 * M_PI, 1000 * S2Cell.averageArea(level));
    S2Cap cap = S2Testing.getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), max_area);
    S2CellId[] covering;
    S2RegionCoverer.getSimpleCovering(cap, cap.center(), level, covering);
    checkCovering(options, cap, covering, false);
  }
}

// We keep a priority queue of the caps that had the worst approximation
// ratios so that we can print them at the end.
struct WorstCap {
  double ratio;
  S2Cap cap;
  int num_cells;

  double opCmp(in WorstCap o) const {
    return ratio - o.ratio;
  }
}

static void checkAccuracy(int max_cells) {
  // SCOPED_TRACE(StrCat(max_cells, " cells"));

  enum int kNumMethods = 1;
  // This code is designed to evaluate several approximation algorithms and
  // figure out which one works better.  The way to do this is to hack the
  // S2RegionCoverer interface to add a global variable to control which
  // algorithm (or variant of an algorithm) is selected, and then assign to
  // this variable in the "method" loop below.  The code below will then
  // collect statistics on all methods, including how often each one wins in
  // terms of cell count and approximation area.

  auto coverer = new S2RegionCoverer();
  coverer.mutableOptions().setMaxCells(max_cells);

  double[kNumMethods] ratio_total = [0];
  double[kNumMethods] min_ratio;  // initialized in loop below
  double[kNumMethods] max_ratio = [0];
  double[][kNumMethods] ratios;
  int[kNumMethods] cell_total = [0];
  int[kNumMethods] area_winner_tally = [0];
  int[kNumMethods] cell_winner_tally = [0];
  static const int kMaxWorstCaps = 10;
  BinaryHeap!(WorstCap[])[kNumMethods] worst_caps =
      array(generate(() => heapify(new WorstCap[0])).take(kNumMethods));

  for (int method = 0; method < kNumMethods; ++method) {
    min_ratio[method] = 1e20;
  }
  for (int i = 0; i < iters; ++i) {
    // Choose the log of the cap area to be uniformly distributed over
    // the allowable range.  Don't try to approximate regions that are so
    // small they can't use the given maximum number of cells efficiently.
    const double min_cap_area = S2Cell.averageArea(S2CellId.MAX_LEVEL) * max_cells * max_cells;
    // Coverings for huge caps are not interesting, so limit the max area too.
    S2Cap cap = S2Testing.getRandomCap(min_cap_area, 0.1 * M_PI);
    double cap_area = cap.getArea();

    double min_area = 1e30;
    int min_cells = 1 << 30;
    double[kNumMethods] area;
    int[kNumMethods] cells;
    for (int method = 0; method < kNumMethods; ++method) {
      // If you want to play with different methods, do this:
      // S2RegionCoverer::method_number = method;

      S2CellId[] covering;
      coverer.getCovering(cap, covering);

      double union_area = 0;
      foreach (S2CellId cell_id; covering) {
        union_area += (new S2Cell(cell_id)).exactArea();
      }
      cells[method] = cast(int) covering.length;
      min_cells = min(cells[method], min_cells);
      area[method] = union_area;
      min_area = min(area[method], min_area);
      cell_total[method] += cells[method];
      double ratio = area[method] / cap_area;
      ratio_total[method] += ratio;
      min_ratio[method] = min(ratio, min_ratio[method]);
      max_ratio[method] = max(ratio, max_ratio[method]);
      ratios[method] ~= ratio;
      if (worst_caps[method].length < kMaxWorstCaps) {
        worst_caps[method].insert(WorstCap(ratio, cap, cells[method]));
      } else if (ratio > worst_caps[method].front().ratio) {
        worst_caps[method].popFront();
        worst_caps[method].insert(WorstCap(ratio, cap, cells[method]));
      }
    }
    for (int method = 0; method < kNumMethods; ++method) {
      if (area[method] == min_area) ++area_winner_tally[method];
      if (cells[method] == min_cells) ++cell_winner_tally[method];
    }
  }
  for (int method = 0; method < kNumMethods; ++method) {
    writefln("\n  Max cells %d, method %d:", max_cells, method);
    writefln("  Average cells: %.4f", cell_total[method] / cast(double)(iters));
    writefln("  Average area ratio: %.4f", ratio_total[method] / iters);
    double[] mratios = ratios[method];
    mratios = array(sort(mratios));
    writefln("  Median ratio: %.4f", mratios[mratios.length / 2]);
    writefln("  Max ratio: %.4f", max_ratio[method]);
    writefln("  Min ratio: %.4f", min_ratio[method]);
    if (kNumMethods > 1) {
      writefln("  Cell winner probability: %.4f", cell_winner_tally[method] / cast(double)(iters));
      writefln("  Area winner probability: %.4f", area_winner_tally[method] / cast(double)(iters));
    }
    writefln("  Caps with the worst approximation ratios:");
    for (; !worst_caps[method].empty(); worst_caps[method].popFront()) {
      const(WorstCap) w = worst_caps[method].front();
      S2LatLng ll = S2LatLng(w.cap.center());
      writefln("    Ratio %.4f, Cells %d, Center (%.8f, %.8f), Km %.6f",
             w.ratio, w.num_cells,
             ll.lat().degrees(), ll.lng().degrees(),
             w.cap.getRadius().radians() * 6367.0);
    }
  }
}

@("Accuracy") unittest {
  for (int i = 0; i < max_cells.length; ++i) {
    checkAccuracy(max_cells[i]);
  }
}

@("InteriorCovering") unittest {
  // We construct the region the following way. Start with S2 cell of level l.
  // Remove from it one of its grandchildren (level l+2). If we then set
  //   min_level < l + 1
  //   max_level > l + 2
  //   max_cells = 3
  // the best interior covering should contain 3 children of the initial cell,
  // that were not effected by removal of a grandchild.
  enum int level = 12;
  S2CellId small_cell = S2CellId(S2Testing.randomPoint()).parent(level + 2);
  S2CellId large_cell = small_cell.parent(level);
  S2CellUnion diff = (new S2CellUnion([large_cell])).difference(new S2CellUnion([small_cell]));
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(3);
  options.setMaxLevel(level + 3);
  options.setMinLevel(level);
  auto coverer = new S2RegionCoverer(options);
  S2CellId[] interior;
  coverer.getInteriorCovering(diff, interior);
  Assert.equal(interior.length, 3);
  for (int i = 0; i < 3; ++i) {
    Assert.equal(interior[i].level(), level + 1);
  }
}

@("HugeFixedLevelCovering") unittest {
  // Test a "fast covering" with a huge number of cells due to min_level().
  auto options = new S2RegionCoverer.Options();
  options.setMinLevel(10);
  auto coverer = new S2RegionCoverer(options);
  S2CellId[] covering;
  auto region = new S2Cell(S2CellId.fromDebugString("1/23"));
  coverer.getFastCovering(region, covering);
  Assert.notLessThan(covering.length, 1 << 16);
}

bool isCanonical(in string[] input_str, S2RegionCoverer.Options options) {
  S2CellId[] input;
  foreach (str; input_str) {
    input ~= S2CellId.fromDebugString(str);
  }
  auto coverer = new S2RegionCoverer(options);
  return coverer.isCanonical(input);
}

@("IsCanonical.InvalidS2CellId") unittest {
  Assert.equal(isCanonical(["1/"], new S2RegionCoverer.Options()), true);
  Assert.equal(isCanonical(["invalid"], new S2RegionCoverer.Options()), false);
}

@("IsCanonical.Unsorted") unittest {
  Assert.equal(isCanonical(["1/1", "1/3"], new S2RegionCoverer.Options()), true);
  Assert.equal(isCanonical(["1/3", "1/1"], new S2RegionCoverer.Options()), false);
}

@("IsCanonical.Overlapping") unittest {
  Assert.equal(isCanonical(["1/2", "1/33"], new S2RegionCoverer.Options()), true);
  Assert.equal(isCanonical(["1/3", "1/33"], new S2RegionCoverer.Options()), false);
}

@("IsCanonical.MinLevel") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMinLevel(2);
  Assert.equal(isCanonical(["1/31"], options), true);
  Assert.equal(isCanonical(["1/3"], options), false);
}

@("IsCanonical.MaxLevel") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMaxLevel(2);
  Assert.equal(isCanonical(["1/31"], options), true);
  Assert.equal(isCanonical(["1/312"], options), false);
}

@("IsCanonical.LevelMod") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setLevelMod(2);
  Assert.equal(isCanonical(["1/31"], options), true);
  Assert.equal(isCanonical(["1/312"], options), false);
}

@("IsCanonical.MaxCells") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(2);
  Assert.equal(isCanonical(["1/1", "1/3"], options), true);
  Assert.equal(isCanonical(["1/1", "1/3", "2/"], options), false);
  Assert.equal(isCanonical(["1/123", "2/1", "3/0122"], options), true);
}

@("IsCanonical.Normalized") unittest {
  // Test that no sequence of cells could be replaced by an ancestor.
  auto options = new S2RegionCoverer.Options();
  Assert.equal(isCanonical(["1/01", "1/02", "1/03", "1/10", "1/11"], options), true);
  Assert.equal(isCanonical(["1/00", "1/01", "1/02", "1/03", "1/10"], options), false);

  Assert.equal(isCanonical(["0/22", "1/01", "1/02", "1/03", "1/10"], options), true);
  Assert.equal(isCanonical(["0/22", "1/00", "1/01", "1/02", "1/03"], options), false);

  options.setMaxCells(20);
  options.setLevelMod(2);
  Assert.equal(isCanonical([
      "1/1101", "1/1102", "1/1103", "1/1110",
      "1/1111", "1/1112", "1/1113", "1/1120",
      "1/1121", "1/1122", "1/1123", "1/1130",
      "1/1131", "1/1132", "1/1133", "1/1200"], options), true);
  Assert.equal(isCanonical([
      "1/1100", "1/1101", "1/1102", "1/1103",
      "1/1110", "1/1111", "1/1112", "1/1113",
      "1/1120", "1/1121", "1/1122", "1/1123",
      "1/1130", "1/1131", "1/1132", "1/1133"], options), false);
}

void checkCanonicalizeCovering(
    string[] input_str, string[] expected_str, S2RegionCoverer.Options options) {
  S2CellId[] actual, expected;
  foreach (str; input_str) {
    actual ~= S2CellId.fromDebugString(str);
  }
  foreach (str; expected_str) {
    expected ~= S2CellId.fromDebugString(str);
  }
  auto coverer = new S2RegionCoverer(options);
  Assert.equal(coverer.isCanonical(actual), false);
  coverer.canonicalizeCovering(actual);
  Assert.equal(coverer.isCanonical(actual), true);
  string[] actual_str;
  Assert.equal(expected, actual);
}

@("CanonicalizeCovering.UnsortedDuplicateCells") unittest {
  auto options = new S2RegionCoverer.Options();
  checkCanonicalizeCovering(
      ["1/200", "1/13122", "1/20", "1/131", "1/13100"],
      ["1/131", "1/20"], options);
}

@("CanonicalizeCovering.MaxLevelExceeded") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMaxLevel(2);
  checkCanonicalizeCovering(
      ["0/3001", "0/3002", "4/012301230123"],
      ["0/30", "4/01"], options);
}

@("CanonicalizeCovering.WrongLevelMod") unittest {
  auto options = new S2RegionCoverer.Options();
  options.setMinLevel(1);
  options.setLevelMod(3);
  checkCanonicalizeCovering(
      ["0/0", "1/11", "2/222", "3/3333"],
      ["0/0", "1/1", "2/2", "3/3333"], options);
}

@("CanonicalizeCovering.ReplacedByParent") unittest {
  // Test that 16 children are replaced by their parent when level_mod == 2.
  auto options = new S2RegionCoverer.Options();
  options.setLevelMod(2);
  checkCanonicalizeCovering(
      ["0/00", "0/01", "0/02", "0/03", "0/10", "0/11", "0/12", "0/13",
       "0/20", "0/21", "0/22", "0/23", "0/30", "0/31", "0/32", "0/33"],
      ["0/"], options);
}

@("CanonicalizeCovering.DenormalizedCellUnion") unittest {
  // Test that all 4 children of a cell may be used when this is necessary to
  // satisfy min_level() or level_mod();
  auto options = new S2RegionCoverer.Options();
  options.setMinLevel(1);
  options.setLevelMod(2);
  checkCanonicalizeCovering(
      ["0/", "1/130", "1/131", "1/132", "1/133"],
      ["0/0", "0/1", "0/2", "0/3", "1/130", "1/131", "1/132", "1/133"],
      options);
}

@("CanonicalizeCovering.MaxCellsMergesSmallest") unittest {
  // When there are too many cells, the smallest cells should be merged first.
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(3);
  checkCanonicalizeCovering(
      ["0/", "1/0", "1/1", "2/01300", "2/0131313"],
      ["0/", "1/", "2/013"], options);
}

@("CanonicalizeCovering.MaxCellsMergesRepeatedly") unittest {
  // Check that when merging creates a cell when all 4 children are present,
  // those cells are merged into their parent (repeatedly if necessary).
  auto options = new S2RegionCoverer.Options();
  options.setMaxCells(8);
  checkCanonicalizeCovering(
      ["0/0121", "0/0123", "1/0", "1/1", "1/2", "1/30", "1/32", "1/33",
       "1/311", "1/312", "1/313", "1/3100", "1/3101", "1/3103",
       "1/31021", "1/31023"],
      ["0/0121", "0/0123", "1/"], options);
}
