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

module s2.s2lax_polyline_shape_test;

import s2.s2lax_polyline_shape;
import s2.s2text_format : parsePointsOrDie;
import s2.s2point;

import fluent.asserts;

@("S2LaxPolylineShape.NoVertices") unittest {
  S2Point[] vertices;
  auto shape = new S2LaxPolylineShape(vertices);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
  Assert.equal(shape.dimension(), 1);
}

@("S2LaxPolylineShape.OneVertex") unittest {
  S2Point[] vertices = [ S2Point(1, 0, 0) ];
  auto shape = new S2LaxPolylineShape(vertices);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
  Assert.equal(shape.dimension(), 1);
}

@("S2LaxPolylineShape.EdgeAccess") unittest {
  S2Point[] vertices = parsePointsOrDie("0:0, 0:1, 1:1");
  auto shape = new S2LaxPolylineShape(vertices);
  Assert.equal(shape.numEdges(), 2);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, 2);
  Assert.equal(shape.dimension(), 1);
  auto edge0 = shape.edge(0);
  Assert.equal(edge0.v0, vertices[0]);
  Assert.equal(edge0.v1, vertices[1]);
  auto edge1 = shape.edge(1);
  Assert.equal(edge1.v0, vertices[1]);
  Assert.equal(edge1.v1, vertices[2]);
}
