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

module s2.s2edge_vector_shape;

import s2.s2shape;
import s2.s2point;

/**
 * S2EdgeVectorShape is an S2Shape representing an arbitrary set of edges.  It
 * is mainly used for testing, but it can also be useful if you have, say, a
 * collection of polylines and don't care about memory efficiency (since this
 * class would store most of the vertices twice).
 *
 * Note that if you already have data stored in an S2Loop, S2Polyline, or
 * S2Polygon, then you would be better off using the "Shape" class defined
 * within those classes (e.g., S2Loop::Shape).  Similarly, if the vertex data
 * is stored in your own data structures, you can easily write your own
 * subclass of S2Shape that points to the existing vertex data rather than
 * copying it.
 */
class S2EdgeVectorShape : S2Shape {
public:
  // Constructs an empty edge vector.
  this() {}

  // Constructs an S2EdgeVectorShape from a vector of edges.
  this(S2Point[2][] edges) {
    _edges = edges;
  }

  // Creates an S2EdgeVectorShape containing a single edge.
  this(in S2Point a, in S2Point b) {
    _edges ~= [a, b];
  }

  // Adds an edge to the vector.
  //
  // IMPORTANT: This method should only be called *before* adding the
  // S2EdgeVectorShape to an S2ShapeIndex.  S2Shapes can only be modified by
  // removing them from the index, making changes, and adding them back again.
  void add(in S2Point a, in S2Point b) {
    _edges ~= [a, b];
  }

  // S2Shape interface:
  override
  int numEdges() const {
    return cast(int) _edges.length;
  }

  override
  S2Shape.Edge edge(int e) const {
    return S2Shape.Edge(_edges[e][0], _edges[e][1]);
  }

  override
  int dimension() const {
    return 1;
  }

  override
  S2Shape.ReferencePoint getReferencePoint() const {
    return ReferencePoint(false);
  }

  override
  int numChains() const {
    return cast(int) _edges.length;
  }

  override
  S2Shape.Chain chain(int i) const {
    return Chain(i, 1);
  }

  override
  S2Shape.Edge chainEdge(int i, int j) const
  in {
    assert(j == 0);
  } do {
    return S2Shape.Edge(_edges[i][0], _edges[i][1]);
  }

  override
  S2Shape.ChainPosition chainPosition(int e) const {
    return S2Shape.ChainPosition(e, 0);
  }

private:
  S2Point[2][] _edges;
}
