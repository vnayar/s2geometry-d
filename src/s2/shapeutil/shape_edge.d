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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.shapeutil.shape_edge;

import s2.s2point;
import s2.s2shape;
import s2.shapeutil.shape_edge_id;


// A class representing a ShapeEdgeId together with the two endpoints of that
// edge.  It should be passed by reference.
class ShapeEdge {
public:
  this() {}
  this(in S2Shape shape, int edge_id) {
    this(shape.id(), edge_id, shape.edge(edge_id));
  }

  this(int shape_id, int edge_id, in S2Shape.Edge edge) {
    _id = ShapeEdgeId(shape_id, edge_id);
    _edge = edge;
  }

  ShapeEdgeId id() const {
    return _id;
  }

  ref const(S2Point) v0() const {
    return _edge.v0;
  }

  ref const(S2Point) v1() const {
    return _edge.v1;
  }

private:
  ShapeEdgeId _id;
  S2Shape.Edge _edge;
}
