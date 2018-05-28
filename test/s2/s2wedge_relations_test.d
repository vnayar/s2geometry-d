module s2.s2wedge_relations_test;

// Copyright 2005 Google Inc. All Rights Reserved.
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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

import fluent.asserts;
import s2.s2wedge_relations;
import s2.s2point;

void testWedge(
    S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2,
    bool contains, bool intersects, WedgeRelation wedge_relation) {
  a0 = a0.normalize();
  ab1 = ab1.normalize();
  a2 = a2.normalize();
  b0 = b0.normalize();
  b2 = b2.normalize();
  Assert.equal(wedgeContains(a0, ab1, a2, b0, b2), contains);
  Assert.equal(wedgeIntersects(a0, ab1, a2, b0, b2), intersects);
  Assert.equal(getWedgeRelation(a0, ab1, a2, b0, b2), wedge_relation);
}

@("Wedges")
unittest {
  // For simplicity, all of these tests use an origin of (0, 0, 1).
  // This shouldn't matter as long as the lower-level primitives are
  // implemented correctly.

  // Intersection in one wedge.
  testWedge(S2Point(-1, 0, 10), S2Point(0, 0, 1), S2Point(1, 2, 10),
            S2Point(0, 1, 10), S2Point(1, -2, 10),
            false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);
  // Intersection in two wedges.
  testWedge(S2Point(-1, -1, 10), S2Point(0, 0, 1), S2Point(1, -1, 10),
            S2Point(1, 0, 10), S2Point(-1, 1, 10),
            false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);

  // Normal containment.
  testWedge(S2Point(-1, -1, 10), S2Point(0, 0, 1), S2Point(1, -1, 10),
            S2Point(-1, 0, 10), S2Point(1, 0, 10),
            true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
  // Containment with equality on one side.
  testWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(-1, -1, 10),
            S2Point(2, 1, 10), S2Point(1, -5, 10),
            true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
  // Containment with equality on the other side.
  testWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(-1, -1, 10),
            S2Point(1, -2, 10), S2Point(-1, -1, 10),
            true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);

  // Containment with equality on both sides.
  testWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(-2, 3, 10), S2Point(4, -5, 10),
            true, true, WedgeRelation.WEDGE_EQUALS);

  // Disjoint with equality on one side.
  testWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(4, -5, 10), S2Point(-2, -3, 10),
            false, false, WedgeRelation.WEDGE_IS_DISJOINT);
  // Disjoint with equality on the other side.
  testWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(0, 5, 10),
            S2Point(4, -5, 10), S2Point(-2, 3, 10),
            false, false, WedgeRelation.WEDGE_IS_DISJOINT);
  // Disjoint with equality on both sides.
  testWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(4, -5, 10), S2Point(-2, 3, 10),
            false, false, WedgeRelation.WEDGE_IS_DISJOINT);

  // B contains A with equality on one side.
  testWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(1, -5, 10),
            S2Point(2, 1, 10), S2Point(-1, -1, 10),
            false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
  // B contains A with equality on the other side.
  testWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(1, -5, 10),
            S2Point(-2, 1, 10), S2Point(1, -5, 10),
            false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
}
