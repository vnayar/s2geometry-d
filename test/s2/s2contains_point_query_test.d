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

module s2.s2contains_point_query_test;

import s2.s2contains_point_query;

import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s2cap;
import s2.s2loop;
import s2.s2point;
import s2.s2shape;
import s2.s2testing;
import s2.s2text_format : makePointOrDie, makeIndexOrDie;
import s2.shapeutil.shape_edge;
import s2.shapeutil.shape_edge_id;

import fluent.asserts;

@("S2ContainsPointQuery.VertexModelOpen") unittest {
  auto index = makeIndexOrDie("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  auto options = new S2ContainsPointQueryOptions(S2VertexModel.OPEN);
  auto q = makeS2ContainsPointQuery(index, options);
  Assert.equal(q.contains(makePointOrDie("0:0")), false);
  Assert.equal(q.contains(makePointOrDie("0:1")), false);
  Assert.equal(q.contains(makePointOrDie("0:2")), false);
  Assert.equal(q.contains(makePointOrDie("0:5")), false);
  Assert.equal(q.contains(makePointOrDie("0:7")), false);
  Assert.equal(q.contains(makePointOrDie("2:6")), false);
  Assert.equal(q.contains(makePointOrDie("1:6")), true);
  Assert.equal(q.contains(makePointOrDie("10:10")), false);

  // Test the last few cases using the Init() method instead.
  auto q2 = new S2ContainsPointQuery!MutableS2ShapeIndex();
  q2.initialize(index, options);
  Assert.equal(q2.shapeContains(index.shape(1), makePointOrDie("1:6")), false);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("1:6")), true);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:5")), false);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:7")), false);
}

@("S2ContainsPointQuery.VertexModelSemiOpen") unittest {
  auto index = makeIndexOrDie("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  auto options = new S2ContainsPointQueryOptions(S2VertexModel.SEMI_OPEN);
  auto q = makeS2ContainsPointQuery(index, options);
  Assert.equal(q.contains(makePointOrDie("0:0")), false);
  Assert.equal(q.contains(makePointOrDie("0:1")), false);
  Assert.equal(q.contains(makePointOrDie("0:2")), false);
  Assert.equal(q.contains(makePointOrDie("0:5")), false);
  Assert.equal(q.contains(makePointOrDie("0:7")), true);  // Contained vertex.
  Assert.equal(q.contains(makePointOrDie("2:6")), false);
  Assert.equal(q.contains(makePointOrDie("1:6")), true);
  Assert.equal(q.contains(makePointOrDie("10:10")), false);

  // Test the last few cases using the Init() method instead.
  auto q2 = new S2ContainsPointQuery!MutableS2ShapeIndex();
  q2.initialize(index, options);
  Assert.equal(q2.shapeContains(index.shape(1), makePointOrDie("1:6")), false);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("1:6")), true);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:5")), false);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:7")), true);
}

@("S2ContainsPointQuery.VertexModelClosed") unittest {
  auto index = makeIndexOrDie("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  auto options = new S2ContainsPointQueryOptions(S2VertexModel.CLOSED);
  auto q = makeS2ContainsPointQuery(index, options);
  Assert.equal(q.contains(makePointOrDie("0:0")), true);
  Assert.equal(q.contains(makePointOrDie("0:1")), true);
  Assert.equal(q.contains(makePointOrDie("0:2")), true);
  Assert.equal(q.contains(makePointOrDie("0:5")), true);
  Assert.equal(q.contains(makePointOrDie("0:7")), true);
  Assert.equal(q.contains(makePointOrDie("2:6")), true);
  Assert.equal(q.contains(makePointOrDie("1:6")), true);
  Assert.equal(q.contains(makePointOrDie("10:10")), false);

  // Test the last few cases using the Init() method instead.
  auto q2 = new S2ContainsPointQuery!MutableS2ShapeIndex();
  q2.initialize(index, options);
  Assert.equal(q2.shapeContains(index.shape(1), makePointOrDie("1:6")), false);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("1:6")), true);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:5")), true);
  Assert.equal(q2.shapeContains(index.shape(2), makePointOrDie("0:7")), true);
}

@("S2ContainsPointQuery.GetContainingShapes") unittest {
  // Also tests ShapeContains().
  const int kNumVerticesPerLoop = 10;
  const S1Angle kMaxLoopRadius = S2Testing.kmToAngle(10);
  const S2Cap center_cap = new S2Cap(S2Testing.randomPoint(), kMaxLoopRadius);
  auto index = new MutableS2ShapeIndex();
  for (int i = 0; i < 100; ++i) {
    S2Loop loop = S2Loop.makeRegularLoop(
        S2Testing.samplePoint(center_cap),
        S2Testing.rnd.randDouble() * kMaxLoopRadius, kNumVerticesPerLoop);
    index.add(new S2Loop.Shape(loop));
  }
  auto query = makeS2ContainsPointQuery(index);
  for (int i = 0; i < 100; ++i) {
    S2Point p = S2Testing.samplePoint(center_cap);
    S2Shape[] expected;
    for (int j = 0; j < index.numShapeIds(); ++j) {
      S2Shape shape = index.shape(j);
      S2Loop loop = (cast(S2Loop.Shape) shape).loop();
      if (loop.contains(p)) {
        Assert.equal(query.shapeContains(shape, p), true);
        expected ~= shape;
      } else {
        Assert.equal(query.shapeContains(shape, p), false);
      }
    }
    S2Shape[] actual = query.getContainingShapes(p);
    Assert.equal(actual, expected);
  }
}

alias EdgeIdVector = ShapeEdgeId[];

void expectIncidentEdgeIds(
    in EdgeIdVector expected, MutableS2ShapeIndex index, in S2Point p) {
  EdgeIdVector actual;
  auto q = makeS2ContainsPointQuery(index);
  Assert.equal(
      q.visitIncidentEdges(
          p,
          (in ShapeEdge e) {
            actual ~= e.id();
            return true;
          }),
      true);
  Assert.equal(expected, actual);
}

@("S2ContainsPointQuery.VisitIncidentEdges") unittest {
  auto index = makeIndexOrDie("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2");
  expectIncidentEdgeIds([ShapeEdgeId(0, 0)], index, makePointOrDie("0:0"));
  expectIncidentEdgeIds([ShapeEdgeId(0, 1), ShapeEdgeId(1, 0)], index, makePointOrDie("1:1"));
  expectIncidentEdgeIds(
      [ShapeEdgeId(1, 0), ShapeEdgeId(2, 0), ShapeEdgeId(2, 2)],
      index, makePointOrDie("1:2"));
  expectIncidentEdgeIds([ShapeEdgeId(2, 0), ShapeEdgeId(2, 1)], index, makePointOrDie("1:3"));
  expectIncidentEdgeIds([ShapeEdgeId(2, 1), ShapeEdgeId(2, 2)], index, makePointOrDie("2:2"));
}
