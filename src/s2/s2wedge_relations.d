module s2.s2wedge_relations;

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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

// Defines functions for determining the relationship between two angles
// ("wedges") that share a common vertex.

import s2.s2point;
import s2.s2predicates : orderedCCW;

// Given an edge chain (x0, x1, x2), the wedge at x1 is the region to the
// left of the edges.  More precisely, it is the set of all rays from x1x0
// (inclusive) to x1x2 (exclusive) in the *clockwise* direction.
//
// The following functions compare two *non-empty* wedges that share the
// same middle vertex: A=(a0, ab1, a2) and B=(b0, ab1, b2).

// Detailed relation from one wedge A to another wedge B.
enum WedgeRelation {
  WEDGE_EQUALS,                 // A and B are equal.
  WEDGE_PROPERLY_CONTAINS,      // A is a strict superset of B.
  WEDGE_IS_PROPERLY_CONTAINED,  // A is a strict subset of B.
  WEDGE_PROPERLY_OVERLAPS,      // A-B, B-A, and A intersect B are non-empty.
  WEDGE_IS_DISJOINT,            // A and B are disjoint.
}

// Returns the relation from wedge A to B.
// REQUIRES: A and B are non-empty.
WedgeRelation getWedgeRelation(
    in S2Point a0, in S2Point ab1, in S2Point a2, in S2Point b0, in S2Point b2) {
  // There are 6 possible edge orderings at a shared vertex (all
  // of these orderings are circular, i.e. abcd == bcda):
  //
  //  (1) a2 b2 b0 a0: A contains B
  //  (2) a2 a0 b0 b2: B contains A
  //  (3) a2 a0 b2 b0: A and B are disjoint
  //  (4) a2 b0 a0 b2: A and B intersect in one wedge
  //  (5) a2 b2 a0 b0: A and B intersect in one wedge
  //  (6) a2 b0 b2 a0: A and B intersect in two wedges
  //
  // We do not distinguish between 4, 5, and 6.
  // We pay extra attention when some of the edges overlap.  When edges
  // overlap, several of these orderings can be satisfied, and we take
  // the most specific.
  if (a0 == b0 && a2 == b2) return WedgeRelation.WEDGE_EQUALS;

  if (orderedCCW(a0, a2, b2, ab1)) {
    // The cases with this vertex ordering are 1, 5, and 6,
    // although case 2 is also possible if a2 == b2.
    if (orderedCCW(b2, b0, a0, ab1)) return WedgeRelation.WEDGE_PROPERLY_CONTAINS;

    // We are in case 5 or 6, or case 2 if a2 == b2.
    return (a2 == b2)
        ? WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED : WedgeRelation.WEDGE_PROPERLY_OVERLAPS;
  }

  // We are in case 2, 3, or 4.
  if (orderedCCW(a0, b0, b2, ab1)) return WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED;
  return orderedCCW(a0, b0, a2, ab1)
      ? WedgeRelation.WEDGE_IS_DISJOINT : WedgeRelation.WEDGE_PROPERLY_OVERLAPS;
}

// Returns true if wedge A contains wedge B.  Equivalent to but faster than
// GetWedgeRelation() == WEDGE_PROPERLY_CONTAINS || WEDGE_EQUALS.
// REQUIRES: A and B are non-empty.
bool wedgeContains(
    in S2Point a0, in S2Point ab1, in S2Point a2, in S2Point b0, in S2Point b2) {
  // For A to contain B (where each loop interior is defined to be its left
  // side), the CCW edge order around ab1 must be a2 b2 b0 a0.  We split
  // this test into two parts that test three vertices each.
  return (orderedCCW(a2, b2, b0, ab1) && orderedCCW(b0, a0, a2, ab1));
}

// Returns true if wedge A intersects wedge B.  Equivalent to but faster
// than GetWedgeRelation() != WEDGE_IS_DISJOINT.
// REQUIRES: A and B are non-empty.
bool wedgeIntersects(
    in S2Point a0, in S2Point ab1, in S2Point a2, in S2Point b0, in S2Point b2) {
  // For A not to intersect B (where each loop interior is defined to be
  // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
  // Note that it's important to write these conditions as negatives
  // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
  // results when two vertices are the same.
  return !(orderedCCW(a0, b2, b0, ab1) && orderedCCW(b0, a2, a0, ab1));
}
