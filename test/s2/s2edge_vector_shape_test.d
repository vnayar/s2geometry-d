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

module s2.s2edge_vector_shape_test;

import s2.s2point;
import s2.s2edge_vector_shape;
import s2.s2testing;

import fluent.asserts;

@("S2EdgeVectorShape.EdgeAccess")
unittest {
  auto shape = new S2EdgeVectorShape();
  S2Testing.rnd.reset(s2RandomSeed);
  enum int kNumEdges = 100;
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point a = S2Testing.randomPoint();  // Control the evaluation order
    shape.add(a, S2Testing.randomPoint());
  }
  Assert.equal(kNumEdges, shape.numEdges());
  Assert.equal(kNumEdges, shape.numChains());
  Assert.equal(1, shape.dimension());
  S2Testing.rnd.reset(s2RandomSeed);
  for (int i = 0; i < kNumEdges; ++i) {
    Assert.equal(i, shape.chain(i).start);
    Assert.equal(1, shape.chain(i).length);
    auto edge = shape.edge(i);
    Assert.equal(S2Testing.randomPoint(), edge.v0);
    Assert.equal(S2Testing.randomPoint(), edge.v1);
  }
}

@("S2EdgeVectorShape.SingletonConstructor")
unittest {
  auto a = S2Point(1, 0, 0);
  auto b = S2Point(0, 1, 0);
  auto shape = new S2EdgeVectorShape(a, b);
  Assert.equal(1, shape.numEdges());
  Assert.equal(1, shape.numChains());
  auto edge = shape.edge(0);
  Assert.equal(a, edge.v0);
  Assert.equal(b, edge.v1);
}
