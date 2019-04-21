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

module s2.s2lax_polygon_shape_test;

import s2.s2lax_polygon_shape;

//import s2.s2lax_loop_shape;
import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s2cap;
import s2.s2contains_point_query;
import s2.s2latlng;
import s2.s2loop;
import s2.s2point;
import s2.s2polygon;
import s2.s2shape;
import s2.s2testing;
import s2.s2text_format;
import s2.shapeutil.contains_brute_force;

import fluent.asserts;

import std.stdio;

@("S2LaxPolygonShape.EmptyPolygon") unittest {
  auto shape = new S2LaxPolygonShape(new S2Polygon());
  Assert.equal(shape.numLoops(), 0);
  Assert.equal(shape.numVertices(), 0);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LaxPolygonShape.FullPolygon") unittest {
  auto shape = new S2LaxPolygonShape(new S2Polygon(makeLoopOrDie("full")));
  Assert.equal(shape.numLoops(), 1);
  Assert.equal(shape.numVertices(), 0);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, true);
}

@("S2LaxPolygonShape.SingleVertexPolygon") unittest {
  // S2Polygon doesn't support single-vertex loops, so we need to construct
  // the S2LaxPolygonShape directly.
  S2Point[][] loops;
  loops ~= parsePoints("0:0");
  auto shape = new S2LaxPolygonShape(loops);
  Assert.equal(shape.numLoops(), 1);
  Assert.equal(shape.numVertices(), 1);
  Assert.equal(shape.numEdges(), 1);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, 1);
  auto edge = shape.edge(0);
  Assert.equal(loops[0][0], edge.v0);
  Assert.equal(loops[0][0], edge.v1);
  Assert.equal(edge == shape.chainEdge(0, 0), true);
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LaxPolygonShape.SingleLoopPolygon") unittest {
  // Test S2Polygon constructor.
  S2Point[] vertices = parsePoints("0:0, 0:1, 1:1, 1:0");
  auto shape = new S2LaxPolygonShape(new S2Polygon(new S2Loop(vertices)));
  Assert.equal(shape.numLoops(), 1);
  Assert.equal(shape.numVertices(), cast(int) vertices.length);
  Assert.equal(shape.numLoopVertices(0), cast(int) vertices.length);
  Assert.equal(shape.numEdges(), cast(int) vertices.length);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, cast(int) vertices.length);
  for (int i = 0; i < vertices.length; ++i) {
    Assert.equal(shape.loopVertex(0, i), vertices[i]);
    auto edge = shape.edge(i);
    Assert.equal(edge.v0, vertices[i]);
    Assert.equal(edge.v1, vertices[(i + 1) % vertices.length]);
    Assert.equal(edge.v0, shape.chainEdge(0, i).v0);
    Assert.equal(edge.v1, shape.chainEdge(0, i).v1);
  }
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(containsBruteForce(shape, origin()), false);
}

@("S2LaxPolygonShape.MultiLoopPolygon") unittest {
  // Test vector<vector<S2Point>> constructor.  Make sure that the loops are
  // oriented so that the interior of the polygon is always on the left.
  S2LaxPolygonShape.Loop[] loops = [
    parsePoints("0:0, 0:3, 3:3"),  // CCW
    parsePoints("1:1, 2:2, 1:2")   // CW
  ];
  auto shape = new S2LaxPolygonShape(loops);

  Assert.equal(loops.length, shape.numLoops());
  int num_vertices = 0;
  Assert.equal(loops.length, shape.numChains());
  for (int i = 0; i < loops.length; ++i) {
    Assert.equal(loops[i].length, shape.numLoopVertices(i));
    Assert.equal(num_vertices, shape.chain(i).start);
    Assert.equal(loops[i].length, shape.chain(i).length);
    for (int j = 0; j < loops[i].length; ++j) {
      Assert.equal(loops[i][j], shape.loopVertex(i, j));
      auto edge = shape.edge(num_vertices + j);
      Assert.equal(loops[i][j], edge.v0);
      Assert.equal(loops[i][(j + 1) % loops[i].length], edge.v1);
    }
    num_vertices += loops[i].length;
  }
  Assert.equal(num_vertices, shape.numVertices());
  Assert.equal(num_vertices, shape.numEdges());
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(containsBruteForce(shape, origin()), false);
}

@("S2LaxPolygonShape.ManyLoopPolygon") unittest {
  // Test a polygon enough loops so that cumulative_vertices_ is used.
  S2Point[][] loops;
  for (int i = 0; i < 100; ++i) {
    auto center = S2LatLng.fromDegrees(0, i).toS2Point();
    loops ~= S2Testing.makeRegularPoints(
        center, S1Angle.fromDegrees(0.1), i % 3 /*S2Testing.rnd.uniform(3)*/);
  }
  auto shape = new S2LaxPolygonShape(loops);

  Assert.equal(loops.length, shape.numLoops());
  int num_vertices = 0;
  Assert.equal(loops.length, shape.numChains());
  for (int i = 0; i < loops.length; ++i) {
    writeln("Test 0: loops[i].length=", loops[i].length);
    Assert.equal(loops[i].length, shape.numLoopVertices(i));
    Assert.equal(num_vertices, shape.chain(i).start);
    Assert.equal(loops[i].length, shape.chain(i).length);
    for (int j = 0; j < loops[i].length; ++j) {
      Assert.equal(loops[i][j], shape.loopVertex(i, j));
      auto edge = shape.edge(num_vertices + j);
      writeln("Test 1: edge=", edge);
      Assert.equal(loops[i][j], edge.v0);
      Assert.equal(loops[i][(j + 1) % loops[i].length], edge.v1);
    }
    num_vertices += loops[i].length;
  }
  Assert.equal(num_vertices, shape.numVertices());
  Assert.equal(num_vertices, shape.numEdges());
}

@("S2LaxPolygonShape.DegenerateLoops") unittest {
  S2LaxPolygonShape.Loop[] loops = [
    parsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
    parsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
    parsePoints("5:5, 6:6")
  ];
  auto shape = new S2LaxPolygonShape(loops);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LaxPolygonShape.InvertedLoops") unittest {
  S2LaxPolygonShape.Loop[] loops = [
    parsePoints("1:2, 1:1, 2:2"),
    parsePoints("3:4, 3:3, 4:4")
  ];
  auto shape = new S2LaxPolygonShape(loops);
  Assert.equal(containsBruteForce(shape, origin()), true);
}

void compareS2LoopToShape(S2Loop loop, S2Shape shape) {
  auto index = new MutableS2ShapeIndex();
  index.add(shape);
  S2Cap cap = loop.getCapBound();
  auto query = makeS2ContainsPointQuery(index);
  for (int iter = 0; iter < 100; ++iter) {
    S2Point point = S2Testing.samplePoint(cap);
    Assert.equal(loop.contains(point), query.shapeContains(index.shape(0), point));
  }
}

/+ TODO: Implement when S2LaxLoopShape is complete.
@("S2LaxPolygonShape.CompareToS2Loop") unittest {
  for (int iter = 0; iter < 100; ++iter) {
    auto fractal = new S2Testing.Fractal();
    fractal.setMaxLevel(S2Testing.rnd.uniform(5));
    fractal.setFractalDimension(1 + S2Testing.rnd.randDouble());
    S2Point center = S2Testing.randomPoint();
    auto loop = new S2Loop(fractal.makeLoop(
            S2Testing.getRandomFrameAt(center), S1Angle.fromDegrees(5)));

    // Compare S2Loop to S2LaxLoopShape.
    compareS2LoopToShape(loop, new S2LaxLoopShape(loop));

    // Compare S2Loop to S2LaxPolygonShape.
    S2LaxPolygonShape.Loop[] loops = [loop];
    compareS2LoopToShape(loop, new S2LaxPolygonShape(loops));
  }
}
+/
