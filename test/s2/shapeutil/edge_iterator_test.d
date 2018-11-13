// Copyright 2018 Google Inc. All Rights Reserved.
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

// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.shapeutil.edge_iterator_test;

import s2.shapeutil.edge_iterator;

import s2.s2shape;
import s2.s2shape_index;
import s2.mutable_s2shape_index;
import s2.s2text_format;

import fluent.asserts;

// Returns the full list of edges in g.
// The edges are collected from points, lines, and polygons in that order.
S2Shape.Edge[] getEdges(S2ShapeIndex index) {
  S2Shape.Edge[] result;
  for (int i = 0; i < index.numShapeIds(); i++) {
    if (index.shape(i) is null) {
      continue;
    }
    for (int j = 0; j < index.shape(i).numEdges(); j++) {
      result ~= index.shape(i).edge(j);
    }
  }
  return result;
}

// Verifies that the edges produced by an EdgeIterator matches GetEdges.
void verify(S2ShapeIndex index) {
  S2Shape.Edge[] expected = getEdges(index);

  int i = 0;
  for (auto it = EdgeIterator(index); !it.done(); it.next(), ++i) {
    Assert.equal(i < expected.length, true);
    Assert.equal(expected[i], it.edge());
  }
}

@("S2ShapeutilEdgeIteratorTest.Empty") unittest {
  auto index = makeIndexOrDie("##");
  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.Points") unittest {
  auto index = makeIndexOrDie("0:0|1:1##");
  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.Lines") unittest {
  auto index = makeIndexOrDie("#0:0,10:10|5:5,5:10|1:2,2:1#");
  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.Polygons") unittest {
  auto index = makeIndexOrDie("##10:10,10:0,0:0|-10:-10,-10:0,0:0,0:-10");
  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.Collection") unittest {
  auto index = makeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      ~ "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.Remove") unittest {
  auto index = makeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      ~ "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
  index.release(0);

  verify(index);
}

@("S2ShapeutilEdgeIteratorTest.AssignmentAndEquality") unittest {
  auto index1 = makeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      ~ "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

  auto index2 = makeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      ~ "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

  auto it1 = EdgeIterator(index1);
  auto it2 = EdgeIterator(index2);

  // Different indices.
  Assert.notEqual(it1, it2);

  it1 = it2;
  Assert.equal(it1, it2);

  it1.next();
  Assert.notEqual(it1, it2);

  it2.next();
  Assert.equal(it1, it2);
}
