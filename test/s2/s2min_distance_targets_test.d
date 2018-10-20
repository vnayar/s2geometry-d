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

module s2.s2min_distance_targets_test;

import s2.s2min_distance_targets;

import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2edge_distances;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2text_format : parsePointsOrDie, makePointOrDie;
import s2.util.container.btree;

import fluent.asserts;
import std.array;

@("PointTarget.UpdateMinDistanceToEdgeWhenEqual") unittest {
  // Verifies that UpdateMinDistance only returns true when the new distance
  // is less than the old distance (not less than or equal to).
  auto target = new S2MinDistancePointTarget(makePointOrDie("1:0"));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto edge = parsePointsOrDie("0:-1, 0:1");
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), true);
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), false);
}

@("PointTarget.UpdateMinDistanceToCellWhenEqual") unittest {
  // Verifies that UpdateMinDistance only returns true when the new distance
  // is less than the old distance (not less than or equal to).
  auto target = new S2MinDistancePointTarget(makePointOrDie("1:0"));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto cell = new S2Cell(S2CellId(makePointOrDie("0:0")));
  Assert.equal(target.updateMinDistance(cell, dist), true);
  Assert.equal(target.updateMinDistance(cell, dist), false);
}

@("EdgeTarget.UpdateMinDistanceToEdgeWhenEqual") unittest {
  auto target = new S2MinDistanceEdgeTarget(makePointOrDie("1:0"), makePointOrDie("1:1"));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto edge = parsePointsOrDie("0:-1, 0:1");
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), true);
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), false);
}

@("EdgeTarget.UpdateMinDistanceToCellWhenEqual") unittest {
  auto target = new S2MinDistanceEdgeTarget(makePointOrDie("1:0"), makePointOrDie("1:1"));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto cell = new S2Cell(S2CellId(makePointOrDie("0:0")));
  Assert.equal(target.updateMinDistance(cell, dist), true);
  Assert.equal(target.updateMinDistance(cell, dist), false);
}

@("CellTarget.UpdateMinDistanceToEdgeWhenEqual") unittest {
  auto target = new S2MinDistanceCellTarget(new S2Cell(S2CellId(makePointOrDie("0:1"))));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto edge = parsePointsOrDie("0:-1, 0:1");
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), true);
  Assert.equal(target.updateMinDistance(edge[0], edge[1], dist), false);
}

@("CellTarget.UpdateMinDistanceToCellWhenEqual") unittest {
  auto target = new S2MinDistanceCellTarget(new S2Cell(S2CellId(makePointOrDie("0:1"))));
  auto dist = S2MinDistance(S1ChordAngle.infinity());
  auto cell = new S2Cell(S2CellId(makePointOrDie("0:0")));
  Assert.equal(target.updateMinDistance(cell, dist), true);
  Assert.equal(target.updateMinDistance(cell, dist), false);
}

/+ TODO: Add when s2min_distance_targets.S2MinDistanceShapeIndexTarget is added.
@("ShapeIndexTarget.UpdateMinDistanceToEdgeWhenEqual") unittest {
  auto target_index = makeIndexOrDie("1:0 # #");
  auto target = new S2MinDistanceShapeIndexTarget(target_index.get());
  S2MinDistance dist(S1ChordAngle::Infinity());
  auto edge = ParsePointsOrDie("0:-1, 0:1");
  EXPECT_TRUE(target.UpdateMinDistance(edge[0], edge[1], &dist));
  EXPECT_FALSE(target.UpdateMinDistance(edge[0], edge[1], &dist));
}

TEST(ShapeIndexTarget, UpdateMinDistanceToCellWhenEqual) {
  auto target_index = MakeIndexOrDie("1:0 # #");
  S2MinDistanceShapeIndexTarget target(target_index.get());
  S2MinDistance dist(S1ChordAngle::Infinity());
  S2Cell cell{S2CellId(MakePointOrDie("0:0"))};
  EXPECT_TRUE(target.UpdateMinDistance(cell, &dist));
  EXPECT_FALSE(target.UpdateMinDistance(cell, &dist));
}
+/

int[] getContainingShapes(S2MinDistanceTarget target, S2ShapeIndex index, int max_shapes) {
  auto shape_ids = new BTree!int();
  target.visitContainingShapes(
      index, (in S2Shape containing_shape, in S2Point target_point) {
        shape_ids.insert(containing_shape.id());
        return shape_ids.length < max_shapes;
      });
  return shape_ids[].array;
}

/+ TODO: Resume when makeIndexOrDie is added.
@("PointTarget.VisitContainingShapes") unittest {
  // Only shapes 2 and 4 should contain the target point.
  auto index = makeIndexOrDie(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
  auto target = new S2MinDistancePointTarget(makePointOrDie("1:1"));
  Assert.equal(getContainingShapes(target, index, 1), [2]);
  Assert.equal(getContainingShapes(target, index, 5), [2, 4]);
}

TEST(EdgeTarget, VisitContainingShapes) {
  // Only shapes 2 and 4 should contain the target point.
  auto index = MakeIndexOrDie(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
  S2MinDistanceEdgeTarget target(MakePointOrDie("1:2"), MakePointOrDie("2:1"));
  EXPECT_EQ((vector<int>{2}), GetContainingShapes(&target, *index, 1));
  EXPECT_EQ((vector<int>{2, 4}), GetContainingShapes(&target, *index, 5));
}

TEST(CellTarget, VisitContainingShapes) {
  auto index = MakeIndexOrDie(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
  // Only shapes 2 and 4 should contain a very small cell near 1:1.
  S2CellId cellid1(MakePointOrDie("1:1"));
  S2MinDistanceCellTarget target1{S2Cell(cellid1)};
  EXPECT_EQ((vector<int>{2}), GetContainingShapes(&target1, *index, 1));
  EXPECT_EQ((vector<int>{2, 4}), GetContainingShapes(&target1, *index, 5));

  // For a larger cell that properly contains one or more index cells, all
  // shapes that intersect the first such cell in S2CellId order are returned.
  // In the test below, this happens to again be the 1st and 3rd polygons
  // (whose shape_ids are 2 and 4).
  S2CellId cellid2 = cellid1.parent(5);
  S2MinDistanceCellTarget target2{S2Cell(cellid2)};
  EXPECT_EQ((vector<int>{2, 4}), GetContainingShapes(&target2, *index, 5));
}

TEST(ShapeIndexTarget, VisitContainingShapes) {
  // Create an index containing a repeated grouping of one point, one
  // polyline, and one polygon.
  auto index = MakeIndexOrDie(
      "1:1 | 4:4 | 7:7 | 10:10 # "
      "1:1, 1:2 | 4:4, 4:5 | 7:7, 7:8 | 10:10, 10:11 # "
      "0:0, 0:3, 3:0 | 3:3, 3:6, 6:3 | 6:6, 6:9, 9:6 | 9:9, 9:12, 12:9");

  // Construct a target consisting of one point, one polyline, and one polygon
  // with two loops where only the second loop is contained by a polygon in
  // the index above.
  auto target_index = MakeIndexOrDie(
      "1:1 # 4:5, 5:4 # 20:20, 20:21, 21:20; 10:10, 10:11, 11:10");

  S2MinDistanceShapeIndexTarget target(target_index.get());
  // These are the shape_ids of the 1st, 2nd, and 4th polygons of "index"
  // (noting that the 4 points are represented by one S2PointVectorShape).
  EXPECT_EQ((vector<int>{5, 6, 8}), GetContainingShapes(&target, *index, 5));
}

TEST(ShapeIndexTarget, VisitContainingShapesEmptyAndFull) {
  // Verify that VisitContainingShapes never returns empty polygons and always
  // returns full polygons (i.e., those containing the entire sphere).

  // Creating an index containing one empty and one full polygon.
  auto index = MakeIndexOrDie("# # empty | full");

  // Check only the full polygon is returned for a point target.
  auto point_index = MakeIndexOrDie("1:1 # #");
  S2MinDistanceShapeIndexTarget point_target(point_index.get());
  EXPECT_EQ((vector<int>{1}), GetContainingShapes(&point_target, *index, 5));

  // Check only the full polygon is returned for a full polygon target.
  auto full_polygon_index = MakeIndexOrDie("# # full");
  S2MinDistanceShapeIndexTarget full_target(full_polygon_index.get());
  EXPECT_EQ((vector<int>{1}), GetContainingShapes(&full_target, *index, 5));

  // Check that nothing is returned for an empty polygon target.  (An empty
  // polygon has no connected components and does not intersect anything, so
  // according to the API of GetContainingShapes nothing should be returned.)
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  S2MinDistanceShapeIndexTarget empty_target(empty_polygon_index.get());
  EXPECT_EQ((vector<int>{}), GetContainingShapes(&empty_target, *index, 5));
}
+/
