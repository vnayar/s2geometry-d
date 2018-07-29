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

module s2.s2shapeutil_contains_brute_force;

import s2.s2shape;
import s2.s2point;
import s2.s2shape_index;
import s2.s2edge_crosser : S2CopyingEdgeCrosser;

// Returns true if the given shape contains the given point.  Most clients
// should not use this method, since its running time is linear in the number
// of shape edges.  Instead clients should create an S2ShapeIndex and use
// S2ContainsPointQuery, since this strategy is much more efficient when many
// points need to be tested.
//
// Polygon boundaries are treated as being semi-open (see S2ContainsPointQuery
// and S2VertexModel for other options).
//
// CAVEAT: Typically this method is only used internally.  Its running time is
//         linear in the number of shape edges.
bool containsBruteForce(in S2Shape shape, in S2Point point) {
  if (!shape.hasInterior()) return false;

  S2Shape.ReferencePoint ref_point = shape.getReferencePoint();
  if (ref_point.point == point) {
    return ref_point.contained;
  }

  auto crosser = new S2CopyingEdgeCrosser(ref_point.point, point);
  bool inside = ref_point.contained;
  for (int e = 0; e < shape.numEdges(); ++e) {
    auto edge = shape.edge(e);
    inside ^= crosser.edgeOrVertexCrossing(edge.v0, edge.v1);
  }
  return inside;
}
