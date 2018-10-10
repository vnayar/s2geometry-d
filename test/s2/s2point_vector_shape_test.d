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

module s2.s2point_vector_shape_test;

import s2.s2point_vector_shape;
import s2.s2testing;
import s2.s2point;

import fluent.asserts;

@("S2PointVectorShape.ConstructionAndAccess") unittest {
  S2Point[] points;
  S2Testing.rnd.reset(s2RandomSeed);
  const int kNumPoints = 100;
  for (int i = 0; i < kNumPoints; ++i) {
    points ~= S2Testing.randomPoint();
  }
  auto shape = new S2PointVectorShape(points);

  Assert.equal(shape.numEdges(), kNumPoints);
  Assert.equal(shape.numChains(), kNumPoints);
  Assert.equal(shape.dimension(), 0);
  S2Testing.rnd.reset(s2RandomSeed);
  for (int i = 0; i < kNumPoints; ++i) {
    Assert.equal(shape.chain(i).start, i);
    Assert.equal(shape.chain(i).length, 1);
    auto edge = shape.edge(i);
    S2Point pt = S2Testing.randomPoint();
    Assert.equal(edge.v0, pt);
    Assert.equal(edge.v1, pt);
  }
}
