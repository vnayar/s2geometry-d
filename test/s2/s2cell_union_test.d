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

module s2.s2cell_union_test;

import s2.s2cell_union;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2edge_distances;
import s2.s2metrics;
import s2.s2point;
import s2.s2testing;
import s2.s2region_coverer;

import fluent.asserts;
import std.stdio;
import std.range;
import std.exception : enforce;
import std.algorithm;
import math = std.math;

Random rnd;

static this() {
  rnd = S2Testing.rnd;
}

class S2CellUnionTestPeer {
public:
  // Creates a possibly invalid S2CellUnion without any checks.
  static S2CellUnion fromVerbatimNoChecks(S2CellId[] cell_ids) {
    return new S2CellUnion(cell_ids, S2CellUnion.VerbatimFlag.VERBATIM);
  }
}

@("DefaultConstructor")
unittest {
  S2CellId[] ids;
  auto empty = new S2CellUnion(ids);
  Assert.equal(empty.empty(), true);
}

@("S2CellIdConstructor")
unittest {
  auto face1_id = S2CellId.fromFace(1);
  auto face1_union = new S2CellUnion([face1_id]);
  Assert.equal(face1_union.numCells(), 1);
  Assert.equal(face1_union.cellId(0), face1_id);
}

@("DuplicateCellsNotValid")
unittest {
  S2CellId id = S2CellId(S2Point(1, 0, 0));
  auto cell_union = S2CellUnionTestPeer.fromVerbatimNoChecks([id, id]);
  Assert.equal(cell_union.isValid(), false);
}

@("UnsortedCellsNotValid")
unittest {
  S2CellId id = S2CellId(S2Point(1, 0, 0)).parent(10);
  auto cell_union = S2CellUnionTestPeer.fromVerbatimNoChecks([id, id.prev()]);
  Assert.equal(cell_union.isValid(), false);
}

@("InvalidCellIdNotValid")
unittest {
  Assert.equal(S2CellId.none().isValid(), false);
  auto cell_union = S2CellUnionTestPeer.fromVerbatimNoChecks([S2CellId.none()]);
  Assert.equal(cell_union.isValid(), false);
}

@("IsNormalized")
unittest {
  auto id = S2CellId(S2Point(1, 0, 0)).parent(10);
  auto cell_union = S2CellUnion.fromVerbatim([id.child(0), id.child(1), id.child(2), id.child(3)]);
  Assert.equal(cell_union.isValid(), true);
  Assert.equal(cell_union.isNormalized(), false);
}

static void addCells(
    S2CellId id, bool selected, ref S2CellId[] input, ref S2CellId[] expected) {
  // Decides whether to add "id" and/or some of its descendants to the
  // test case.  If "selected" is true, then the region covered by "id"
  // *must* be added to the test case (either by adding "id" itself, or
  // some combination of its descendants, or both).  If cell ids are to
  // the test case "input", then the corresponding expected result after
  // simplification is added to "expected".

  if (id == S2CellId.none()) {
    // Initial call: decide whether to add cell(s) from each face.
    for (int face = 0; face < 6; ++face) {
      addCells(S2CellId.fromFace(face), false, input, expected);
    }
    return;
  }
  if (id.isLeaf()) {
    // The rnd.OneIn() call below ensures that the parent of a leaf cell
    // will always be selected (if we make it that far down the hierarchy).
    enforce(selected);
    input ~= id;
    return;
  }
  // The following code ensures that the probability of selecting a cell
  // at each level is approximately the same, i.e. we test normalization
  // of cells at all levels.
  if (!selected && rnd.oneIn(S2CellId.MAX_LEVEL - id.level())) {
    // Once a cell has been selected, the expected output is predetermined.
    // We then make sure that cells are selected that will normalize to
    // the desired output.
    expected ~= id;
    selected = true;
  }

  // With the rnd.OneIn() constants below, this function adds an average
  // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
  // level at which the cell was first selected (level 15 on average).
  // Therefore the average number of input cells in a test case is about
  // (5/6 * 15 * 6) = 75.  The average number of output cells is about 6.

  // If a cell is selected, we add it to "input" with probability 5/6.
  bool added = false;
  if (selected && !rnd.oneIn(6)) {
    input ~= id;
    added = true;
  }
  int num_children = 0;
  S2CellId child = id.childBegin();
  for (int pos = 0; pos < 4; ++pos, child = child.next()) {
    // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
    // This intentionally may result in a cell and some of its children
    // being included in the test case.
    //
    // If the cell is not selected, on average we recurse on one child.
    // We also make sure that we do not recurse on all 4 children, since
    // then we might include all 4 children in the input case by accident
    // (in which case the expected output would not be correct).
    if (rnd.oneIn(selected ? 12 : 4) && num_children < 3) {
      addCells(child, selected, input, expected);
      ++num_children;
    }
    // If this cell was selected but the cell itself was not added, we
    // must ensure that all 4 children (or some combination of their
    // descendants) are added.
    if (selected && !added) addCells(child, selected, input, expected);
  }
}

@("Normalize")
unittest {
  // Try a bunch of random test cases, and keep track of average
  // statistics for normalization (to see if they agree with the
  // analysis above).
  double in_sum = 0, out_sum = 0;
  enum int kIters = 2000;
  for (int i = 0; i < kIters; ++i) {
    S2CellId[] input, expected;
    addCells(S2CellId.none(), false, input, expected);
    in_sum += input.length;
    out_sum += expected.length;
    auto cellunion = new S2CellUnion(input);
    Assert.equal(cellunion.size(), expected.length);
    for (int j = 0; j < expected.length; ++j) {
      Assert.equal(cellunion[j], expected[j]);
    }

    // Test GetCapBound().
    S2Cap cap = cellunion.getCapBound();
    foreach (S2CellId id; cellunion.cellIds()) {
      Assert.equal(cap.contains(new S2Cell(id)), true);
    }

    // Test Contains(S2CellId) and Intersects(S2CellId).
    foreach (S2CellId input_id; input) {
      Assert.equal(cellunion.contains(input_id), true);
      Assert.equal(cellunion.contains(input_id.toS2Point()), true);
      Assert.equal(cellunion.intersects(input_id), true);
      if (!input_id.isFace()) {
        Assert.equal(cellunion.intersects(input_id.parent()), true);
        if (input_id.level() > 1) {
          Assert.equal(cellunion.intersects(input_id.parent().parent()), true);
          Assert.equal(cellunion.intersects(input_id.parent(0)), true);
        }
      }
      if (!input_id.isLeaf()) {
        Assert.equal(cellunion.contains(input_id.childBegin()), true);
        Assert.equal(cellunion.intersects(input_id.childBegin()), true);
        Assert.equal(cellunion.contains(input_id.childEnd().prev()), true);
        Assert.equal(cellunion.intersects(input_id.childEnd().prev()), true);
        Assert.equal(cellunion.contains(input_id.childBegin(S2CellId.MAX_LEVEL)), true);
        Assert.equal(cellunion.intersects(input_id.childBegin(S2CellId.MAX_LEVEL)), true);
      }
    }
    foreach (S2CellId expected_id; expected) {
      if (!expected_id.isFace()) {
        Assert.equal(!cellunion.contains(expected_id.parent()), true);
        Assert.equal(!cellunion.contains(expected_id.parent(0)), true);
      }
    }

    // Test Contains(S2CellUnion), Intersects(S2CellUnion), Union(),
    // Intersection(), and Difference().
    S2CellId[] x, y, x_or_y, x_and_y;
    foreach (S2CellId input_id; input) {
      bool in_x = rnd.oneIn(2);
      bool in_y = rnd.oneIn(2);
      if (in_x) x ~= input_id;
      if (in_y) y ~= input_id;
      if (in_x || in_y) x_or_y ~= input_id;
    }
    auto xcells = new S2CellUnion(x);
    auto ycells = new S2CellUnion(y);
    auto x_or_y_expected = new S2CellUnion(x_or_y);
    auto x_or_y_cells = xcells.unite(ycells);
    Assert.equal(x_or_y_cells, x_or_y_expected);

    // Compute the intersection of "x" with each cell of "y",
    // check that this intersection is correct, and append the
    // results to x_and_y_expected.
    foreach (S2CellId yid; ycells.cellIds()) {
      S2CellUnion ucells = xcells.intersect(yid);
      foreach (S2CellId xid; xcells.cellIds()) {
        if (xid.contains(yid)) {
          Assert.equal(ucells.size(), 1);
          Assert.equal(ucells[0], yid);
        } else if (yid.contains(xid)) {
          Assert.equal(ucells.contains(xid), true);
        }
      }
      foreach (S2CellId uid; ucells.cellIds()) {
        Assert.equal(xcells.contains(uid), true);
        Assert.equal(yid.contains(uid), true);
      }
      x_and_y ~= ucells.cellIds();
    }
    auto x_and_y_expected = new S2CellUnion(x_and_y);
    auto x_and_y_cells = xcells.intersect(ycells);
    Assert.equal(x_and_y_cells, x_and_y_expected);

    S2CellUnion x_minus_y_cells = xcells.difference(ycells);
    S2CellUnion y_minus_x_cells = ycells.difference(xcells);
    Assert.equal(xcells.contains(x_minus_y_cells), true);
    Assert.equal(!x_minus_y_cells.intersects(ycells), true);
    Assert.equal(ycells.contains(y_minus_x_cells), true);
    Assert.equal(!y_minus_x_cells.intersects(xcells), true);
    Assert.equal(!x_minus_y_cells.intersects(y_minus_x_cells), true);

    S2CellUnion diff_intersection_union =
        x_minus_y_cells.unite(y_minus_x_cells).unite(x_and_y_cells);
    Assert.equal(diff_intersection_union, x_or_y_cells);

    S2CellId[] test, dummy;
    addCells(S2CellId.none(), false, test, dummy);
    foreach (S2CellId test_id; test) {
      bool contains = false, intersects = false;
      foreach (S2CellId expected_id; expected) {
        if (expected_id.contains(test_id)) contains = true;
        if (expected_id.intersects(test_id)) intersects = true;
      }
      Assert.equal(cellunion.contains(test_id), contains);
      Assert.equal(cellunion.intersects(test_id), intersects);
    }
  }
  writefln("avg in %.2f, avg out %.2f\n", in_sum / kIters, out_sum / kIters);
}

// Return the maximum geodesic distance from "axis" to any point of
// "covering".
static double getRadius(in S2CellUnion covering, in S2Point axis) {
  double max_dist = 0;
  foreach (S2CellId id; covering.cellIds()) {
    auto cell = new S2Cell(id);
    for (int j = 0; j < 4; ++j) {
      S2Point a = cell.getVertex(j);
      S2Point b = cell.getVertex(j + 1);
      double dist;
      // The maximum distance is not always attained at a cell vertex: if at
      // least one vertex is in the opposite hemisphere from "axis" then the
      // maximum may be attained along an edge.  We solve this by computing
      // the minimum distance from the edge to (-axis) instead.  We can't
      // simply do this all the time because S2::GetDistance() has
      // poor accuracy when the result is close to Pi.
      //
      // TODO(ericv): Improve S2::GetDistance() accuracy near Pi.
      if (a.angle(axis) > math.PI_2 || b.angle(axis) > math.PI_2) {
        dist = math.PI - getDistance(-axis, a, b).radians();
      } else {
        dist = a.angle(axis);
      }
      max_dist = max(max_dist, dist);
    }
  }
  return max_dist;
}

@("S2CellUnion.Expand") unittest {
  // This test generates coverings for caps of random sizes, expands
  // the coverings by a random radius, and then make sure that the new
  // covering covers the expanded cap.  It also makes sure that the
  // new covering is not too much larger than expected.

  auto coverer = new S2RegionCoverer();
  for (int i = 0; i < 1000; ++i) {
    //SCOPED_TRACE(StrCat("Iteration ", i));
    S2Cap cap = S2Testing.getRandomCap(
        S2Cell.averageArea(S2CellId.MAX_LEVEL), 4 * math.PI);

    // Expand the cap area by a random factor whose log is uniformly
    // distributed between 0 and log(1e2).
    S2Cap expanded_cap = S2Cap.fromCenterHeight(
        cap.center(), min(2.0, math.pow(1e2, rnd.randDouble()) * cap.height()));

    double radius = (expanded_cap.getRadius() - cap.getRadius()).radians();
    int max_level_diff = rnd.uniform(8);

    // Generate a covering for the original cap, and measure the maximum
    // distance from the cap center to any point in the covering.
    coverer.mutableOptions().setMaxCells(1 + rnd.skewed(10));
    S2CellUnion covering = coverer.getCovering(cap);
    S2Testing.checkCovering(cap, covering, true);
    double covering_radius = getRadius(covering, cap.center());

    // This code duplicates the logic in Expand(min_radius, max_level_diff)
    // that figures out an appropriate cell level to use for the expansion.
    int min_level = S2CellId.MAX_LEVEL;
    foreach (S2CellId id; covering.cellIds()) {
      min_level = min(min_level, id.level());
    }
    int expand_level = min(min_level + max_level_diff,
        MIN_WIDTH.getLevelForMinValue(radius));

    // Generate a covering for the expanded cap, and measure the new maximum
    // distance from the cap center to any point in the covering.
    covering.expand(S1Angle.fromRadians(radius), max_level_diff);
    S2Testing.checkCovering(expanded_cap, covering, false);
    double expanded_covering_radius = getRadius(covering, cap.center());

    // If the covering includes a tiny cell along the boundary, in theory the
    // maximum angle of the covering from the cap center can increase by up to
    // twice the maximum length of a cell diagonal.
    Assert.notGreaterThan(expanded_covering_radius - covering_radius,
        2 * MAX_DIAG.getValue(expand_level));
  }
}

/+ TODO: Add when encode/decode is implemented.
TEST(S2CellUnion, EncodeDecode) {
  vector<S2CellId> cell_ids = {S2CellId(0x33),
                               S2CellId(0x8e3748fab),
                               S2CellId(0x91230abcdef83427)};
  auto cell_union = S2CellUnion::FromVerbatim(std::move(cell_ids));

  Encoder encoder;
  cell_union.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2CellUnion decoded_cell_union;
  EXPECT_TRUE(decoded_cell_union.Decode(&decoder));
  EXPECT_EQ(cell_union, decoded_cell_union);
}

TEST(S2CellUnion, EncodeDecodeEmpty) {
  S2CellUnion empty_cell_union;

  Encoder encoder;
  empty_cell_union.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2CellUnion decoded_cell_union;
  EXPECT_TRUE(decoded_cell_union.Decode(&decoder));
  EXPECT_EQ(empty_cell_union, decoded_cell_union);
}
+/


static void testFromMinMax(S2CellId min_id, S2CellId max_id) {
  auto cell_union = S2CellUnion.fromMinMax(min_id, max_id);
  S2CellId[] cell_ids = cell_union.cellIds();

  Assert.greaterThan(cell_ids.length, 0);
  Assert.equal(cell_ids.front().rangeMin(), min_id);
  Assert.equal(cell_ids.back().rangeMax(), max_id);
  for (int i = 1; i < cell_ids.length; ++i) {
    Assert.equal(cell_ids[i].rangeMin(), cell_ids[i-1].rangeMax().next());
  }
  Assert.equal(cell_union.isNormalized(), true);
}

@("FromMinMax")
unittest {
  // Check the very first leaf cell and face cell.
  S2CellId face1_id = S2CellId.fromFace(0);
  testFromMinMax(face1_id.rangeMin(), face1_id.rangeMin());
  testFromMinMax(face1_id.rangeMin(), face1_id.rangeMax());

  // Check the very last leaf cell and face cell.
  S2CellId face5_id = S2CellId.fromFace(5);
  testFromMinMax(face5_id.rangeMin(), face5_id.rangeMax());
  testFromMinMax(face5_id.rangeMax(), face5_id.rangeMax());

  // Check random ranges of leaf cells.
  for (int iter = 0; iter < 100; ++iter) {
    S2CellId x = S2Testing.getRandomCellId(S2CellId.MAX_LEVEL);
    S2CellId y = S2Testing.getRandomCellId(S2CellId.MAX_LEVEL);
    if (x > y) swap(x, y);
    testFromMinMax(x, y);
  }
}

@("FromBeginEnd")
unittest {
  // Since FromMinMax() is implemented in terms of FromBeginEnd(), we
  // focus on test cases that generate an empty range.
  S2CellId initial_id = S2CellId.fromFace(3);

  // Test an empty range before the minimum S2CellId.
  auto cell_union = new S2CellUnion([initial_id]);
  S2CellId id_begin = S2CellId.begin(S2CellId.MAX_LEVEL);
  cell_union.initFromBeginEnd(id_begin, id_begin);
  Assert.equal(cell_union.empty(), true);

  // Test an empty range after the maximum S2CellId.
  cell_union.init([initial_id]);
  S2CellId id_end = S2CellId.end(S2CellId.MAX_LEVEL);
  cell_union.initFromBeginEnd(id_end, id_end);
  Assert.equal(cell_union.empty(), true);

  // Test the full sphere.
  cell_union = S2CellUnion.fromBeginEnd(id_begin, id_end);
  Assert.equal(cell_union.numCells(), 6);
  foreach (S2CellId id; cell_union.cellIds()) {
    Assert.equal(id.isFace(), true);
  }
}

@("Empty")
unittest {
  auto empty_cell_union = new S2CellUnion();
  auto face1_id = S2CellId.fromFace(1);

  // Normalize()
  empty_cell_union.normalize();
  Assert.equal(empty_cell_union.empty(), true);

  // Denormalize(...)
  S2CellId[] output;
  empty_cell_union.denormalize(0, 2, output);
  Assert.equal(empty_cell_union.empty(), true);

  // Pack(...)
  empty_cell_union.pack();

  // Contains(...)
  Assert.equal(empty_cell_union.contains(face1_id), false);
  Assert.equal(empty_cell_union.contains(empty_cell_union), true);

  // Intersects(...)
  Assert.equal(empty_cell_union.intersects(face1_id), false);
  Assert.equal(empty_cell_union.intersects(empty_cell_union), false);

  // Union(...)
  S2CellUnion cell_union = empty_cell_union.unite(empty_cell_union);
  Assert.equal(cell_union.empty(), true);

  // Intersection(...)
  S2CellUnion intersection = empty_cell_union.intersect(face1_id);
  Assert.equal(intersection.empty(), true);
  intersection = empty_cell_union.intersect(empty_cell_union);
  Assert.equal(intersection.empty(), true);

  // Difference(...)
  S2CellUnion difference = empty_cell_union.difference(empty_cell_union);
  Assert.equal(difference.numCells(), 0);

  // Expand(...)
  empty_cell_union.expand(S1Angle.fromRadians(1), 20);
  Assert.equal(empty_cell_union.empty(), true);
  empty_cell_union.expand(10);
  Assert.equal(empty_cell_union.empty(), true);
}

@("Clear")
unittest {
  auto face1_id = S2CellId.fromFace(1);
  auto face1_union = new S2CellUnion([face1_id]);

  Assert.equal(face1_union.numCells(), 1);
  Assert.equal(face1_union.cellIds().length, 1);

  face1_union.clear();
  Assert.equal(face1_union.numCells(), 0);
  Assert.equal(face1_union.cellIds().length, 0);
  Assert.equal(face1_union.cellIds().capacity(), 0);
}

/+ TODO: Add when decode is added.
TEST(S2CellUnion, RefuseToDecode) {
  std::vector<S2CellId> cellids;
  S2CellId id = S2CellId::Begin(S2CellId::kMaxLevel);
  for (int i = 0; i <= FLAGS_s2cell_union_decode_max_num_cells; ++i) {
    cellids.push_back(id);
    id = id.next();
  }
  S2CellUnion cell_union = S2CellUnion::FromVerbatim(cellids);
  Encoder encoder;
  cell_union.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2CellUnion decoded_cell_union;
  EXPECT_FALSE(decoded_cell_union.Decode(&decoder));
}
+/

@("Release")
unittest {
  auto face1_id = S2CellId.fromFace(1);
  auto face1_union = new S2CellUnion([face1_id]);
  Assert.equal(face1_union.numCells(), 1);
  Assert.equal(face1_union.cellId(0), face1_id);

  S2CellId[] released = face1_union.release();
  Assert.equal(released.length, 1);
  Assert.equal(released[0], face1_id);
  Assert.equal(face1_union.numCells(), 0);
}

@("LeafCellsCovered")
unittest {
  auto cell_union = new S2CellUnion();
  Assert.equal(cell_union.leafCellsCovered(), 0);

  S2CellId[] ids;
  // One leaf cell on face 0.
  ids ~= S2CellId.fromFace(0).childBegin(S2CellId.MAX_LEVEL);
  cell_union.init(ids);
  Assert.equal(cell_union.leafCellsCovered(), 1uL);

  // Face 0 itself (which includes the previous leaf cell).
  ids ~= S2CellId.fromFace(0);
  cell_union.init(ids);
  Assert.equal(cell_union.leafCellsCovered(), 1uL << 60);
  // Five faces.
  cell_union.expand(0);
  Assert.equal(cell_union.leafCellsCovered(), 5uL << 60);
  // Whole world.
  cell_union.expand(0);
  Assert.equal(cell_union.leafCellsCovered(), 6uL << 60);

  // Add some disjoint cells.
  ids ~= S2CellId.fromFace(1).childBegin(1);
  ids ~= S2CellId.fromFace(2).childBegin(2);
  ids ~= S2CellId.fromFace(2).childEnd(2).prev();
  ids ~= S2CellId.fromFace(3).childBegin(14);
  ids ~= S2CellId.fromFace(4).childBegin(27);
  ids ~= S2CellId.fromFace(4).childEnd(15).prev();
  ids ~= S2CellId.fromFace(5).childBegin(30);
  cell_union.init(ids);
  ulong expected = 1uL + (1uL << 6) + (1uL << 30) + (1uL << 32)
      + (2uL << 56) + (1uL << 58) + (1uL << 60);
  Assert.equal(cell_union.leafCellsCovered(), expected);
}

@("WorksInContainers")
unittest {
  S2CellId[] ids = [S2CellId.fromFace(1)];
  auto cell_union0 = new S2CellUnion(ids);

  // This gives a compilation error if the S2CellUnion is neither movable nor
  // copyable.
  S2CellUnion[] union_vector;
  union_vector ~= cell_union0;

  Assert.equal(union_vector.back().cellIds(), ids);
}
