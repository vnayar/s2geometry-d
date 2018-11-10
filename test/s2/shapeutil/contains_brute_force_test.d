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

module s2.shapeutil.contains_brute_force_test;

import s2.s2lax_polygon_shape;
import s2.s2lax_polyline_shape;
import s2.s2text_format;
import s2.shapeutil.contains_brute_force : containsBruteForce;

import fluent.asserts;

@("ContainsBruteForce.NoInterior") unittest {
  // Defines a polyline that almost entirely encloses the point 0:0.
  auto polyline = makeLaxPolylineOrDie("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
  Assert.equal(containsBruteForce(polyline, makePointOrDie("0:0")), false);
}

@("ContainsBruteForce.ContainsReferencePoint") unittest {
  // Checks that ContainsBruteForce agrees with GetReferencePoint.
  auto polygon = makeLaxPolygonOrDie("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
  auto refPoint = polygon.getReferencePoint();
  Assert.equal(refPoint.contained, containsBruteForce(polygon, refPoint.point));
}

/+ TODO: Resume when S2Loop is implemented.
@("ContainsBruteForce.ConsistentWithS2Loop") unittest {
  // Checks that ContainsBruteForce agrees with S2Loop::Contains().
  auto loop = S2Loop::MakeRegularLoop(MakePoint("89:-179"),
                                      S1Angle::Degrees(10), 100);
  S2Loop::Shape shape(loop.get());
  for (int i = 0; i < loop->num_vertices(); ++i) {
    EXPECT_EQ(loop->Contains(loop->vertex(i)),
              ContainsBruteForce(shape, loop->vertex(i)));
  }
}
+/
