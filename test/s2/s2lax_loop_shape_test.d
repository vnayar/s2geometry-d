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

module s2.s2lax_loop_shape_test;

import s2.s2lax_loop_shape;

import s2.s2loop;
import s2.s2point;
import s2.s2text_format;
import s2.shapeutil.contains_brute_force;

import fluent.asserts;

@("S2LaxLoopShape.EmptyLoop") unittest {
  // Test S2Loop constructor.
  auto shape = new S2LaxLoopShape();
  shape.initialize(new S2Loop(S2Loop.empty()));
  Assert.equal(shape.numVertices(), 0);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LaxLoopShape.NonEmptyLoop") unittest {
  // Test vector<S2Point> constructor.
  S2Point[] vertices = parsePointsOrDie("0:0, 0:1, 1:1, 1:0");
  auto shape = new S2LaxLoopShape(vertices);
  Assert.equal(shape.numVertices(), cast(int) vertices.length);
  Assert.equal(shape.numEdges(), cast(int) vertices.length);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, cast(int) vertices.length);
  for (int i = 0; i < vertices.length; ++i) {
    Assert.equal(shape.vertex(i), vertices[i]);
    auto edge = shape.edge(i);
    Assert.equal(edge.v0, vertices[i]);
    Assert.equal(edge.v1, vertices[(i + 1) % vertices.length]);
  }
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LaxClosedPolylineShape.NoInterior") unittest {
  S2Point[] vertices = parsePointsOrDie("0:0, 0:1, 1:1, 1:0");
  auto shape = new S2LaxClosedPolylineShape(vertices);
  Assert.equal(shape.dimension(), 1);
  Assert.equal(shape.hasInterior(), false);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2VertexIdLaxLoopShape.EmptyLoop") unittest {
  auto shape = new S2VertexIdLaxLoopShape(new int[0], null);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numVertices(), 0);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2VertexIdLaxLoopShape.InvertedLoop") unittest {
  import s2.s2pointutil : origin;

  S2Point[] vertex_array = parsePointsOrDie("0:0, 0:1, 1:1, 1:0");
  int[] vertex_ids = [ 0, 3, 2, 1 ];  // Inverted.
  auto shape = new S2VertexIdLaxLoopShape(vertex_ids, vertex_array);
  Assert.equal(shape.numEdges(), 4);
  Assert.equal(shape.numVertices(), 4);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, 4);
  Assert.equal(shape.vertex(0), vertex_array[0]);
  Assert.equal(shape.vertex(1), vertex_array[3]);
  Assert.equal(shape.vertex(2), vertex_array[2]);
  Assert.equal(shape.vertex(3), vertex_array[1]);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(containsBruteForce(shape, origin()), true);
}
