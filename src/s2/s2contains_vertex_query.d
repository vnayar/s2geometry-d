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

module s2.s2contains_vertex_query;

import s2.s2point;
import s2.s2pointutil : ortho;
import s2.s2predicates : orderedCCW;

import std.exception : enforce;
import std.math : abs;
import std.typecons : tuple;

// This class determines whether a polygon contains one of its vertices given
// the edges incident to that vertex.  The result is +1 if the vertex is
// contained, -1 if it is not contained, and 0 if the incident edges consist
// of matched sibling pairs (in which case the result cannot be determined
// locally).
//
// Point containment is defined according to the "semi-open" boundary model
// (see S2VertexModel), which means that if several polygons tile the region
// around a vertex, then exactly one of those polygons contains that vertex.
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
class S2ContainsVertexQuery {
public:
  // "target" is the vertex whose containment will be determined.
  this(in S2Point target) {
    _target = target;
  }

  // Indicates that the polygon has an edge between "target" and "v" in the
  // given direction (+1 = outgoing, -1 = incoming, 0 = degenerate).
  void addEdge(in S2Point v, int direction) {
    _edgeMap[v] += direction;
  }

  // Returns +1 if the vertex is contained, -1 if it is not contained, and 0
  // if the incident edges consisted of matched sibling pairs.
  int containsSign() {
    // Find the unmatched edge that is immediately clockwise from S2::Ortho(P).
    S2Point reference_dir = ortho(_target);
    auto best = tuple(reference_dir, 0);
    foreach (S2Point point, int dir; _edgeMap) {
      enforce(abs(dir) <= 1);
      if (dir == 0) continue;  // This is a "matched" edge.
      if (orderedCCW(reference_dir, best[0], point, _target)) {
        best = tuple(point, dir);
      }
    }
    return best[1];
  }

private:
  S2Point _target;
  int[S2Point] _edgeMap;
}
