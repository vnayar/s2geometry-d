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

module s2.shapeutil.get_reference_point_test;

import s2.shapeutil.get_reference_point;

import s2.s2lax_polygon_shape;
import s2.s2polygon;
import s2.shapeutil.contains_brute_force;
import s2.s2testing;
import s2.s2text_format;

import fluent.asserts;

@("GetReferencePoint.EmptyPolygon") unittest {
  auto shape = new S2LaxPolygonShape(new S2Polygon());
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("GetReferencePoint.FullPolygon") unittest {
  auto shape = new S2LaxPolygonShape(new S2Polygon(makeLoopOrDie("full")));
  Assert.equal(shape.getReferencePoint().contained, true);
}

/+

TEST(GetReferencePoint, DegenerateLoops) {
  vector<S2LaxPolygonShape::Loop> loops = {
    s2textformat::ParsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
    s2textformat::ParsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
    s2textformat::ParsePoints("5:5, 6:6")
  };
  S2LaxPolygonShape shape(loops);
  EXPECT_FALSE(shape.GetReferencePoint().contained);
}

TEST(GetReferencePoint, InvertedLoops) {
  vector<S2LaxPolygonShape::Loop> loops = {
    s2textformat::ParsePoints("1:2, 1:1, 2:2"),
    s2textformat::ParsePoints("3:4, 3:3, 4:4")
  };
  S2LaxPolygonShape shape(loops);
  EXPECT_TRUE(s2shapeutil::ContainsBruteForce(shape, S2::Origin()));
}

TEST(GetReferencePoint, PartiallyDegenerateLoops) {
  for (int iter = 0; iter < 100; ++iter) {
    // First we construct a long convoluted edge chain that follows the
    // S2CellId Hilbert curve.  At some random point along the curve, we
    // insert a small triangular loop.
    vector<S2LaxPolygonShape::Loop> loops(1);
    S2LaxPolygonShape::Loop* loop = &loops[0];
    const int num_vertices = 100;
    S2CellId start = S2Testing::GetRandomCellId(S2CellId::kMaxLevel - 1);
    S2CellId end = start.advance_wrap(num_vertices);
    S2CellId loop_cellid = start.advance_wrap(
        S2Testing::rnd.Uniform(num_vertices - 2) + 1);
    vector<S2Point> triangle;
    for (S2CellId cellid = start; cellid != end; cellid = cellid.next_wrap()) {
      if (cellid == loop_cellid) {
        // Insert a small triangular loop.  We save the loop so that we can
        // test whether it contains the origin later.
        triangle.push_back(cellid.child(0).ToPoint());
        triangle.push_back(cellid.child(1).ToPoint());
        triangle.push_back(cellid.child(2).ToPoint());
        loop->insert(loop->end(), triangle.begin(), triangle.end());
        loop->push_back(cellid.child(0).ToPoint());
      } else {
        loop->push_back(cellid.ToPoint());
      }
    }
    // Now we retrace our steps, except that we skip the three edges that form
    // the triangular loop above.
    for (S2CellId cellid = end; cellid != start; cellid = cellid.prev_wrap()) {
      if (cellid == loop_cellid) {
        loop->push_back(cellid.child(0).ToPoint());
      } else {
        loop->push_back(cellid.ToPoint());
      }
    }
    S2LaxPolygonShape shape(loops);
    S2Loop triangle_loop(triangle);
    auto ref = shape.GetReferencePoint();
    EXPECT_EQ(triangle_loop.Contains(ref.point), ref.contained);
  }
}
+/
