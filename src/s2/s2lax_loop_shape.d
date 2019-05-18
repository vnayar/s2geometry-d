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
//
// This file defines various S2Shape types representing loops:
//
// S2LaxLoopShape
//   - like S2Loop::Shape but allows duplicate vertices & edges, more compact
//     representation, and faster to initialize.
//
// S2LaxClosedPolylineShape
//   - like S2LaxLoopShape, but defines a loop that does not have an interior
//     (a closed polyline).
//
// S2VertexIdLaxLoopShape
//   - like S2LaxLoopShape, but vertices are specified as indices into an
//     existing vertex array.

module s2.s2lax_loop_shape;

import s2.s2loop;
import s2.s2shape;
import s2.shapeutil.get_reference_point;
import s2.s2point;

import std.algorithm;
import std.typecons;

/**
 * S2LaxLoopShape represents a closed loop of edges surrounding an interior
 * region.  It is similar to S2Loop::Shape except that this class allows
 * duplicate vertices and edges.  Loops may have any number of vertices,
 * including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
 * consisting of a single point.)
 *
 * Note that S2LaxLoopShape is faster to initialize and more compact than
 * S2Loop::Shape, but does not support the same operations as S2Loop.
 */
class S2LaxLoopShape : S2Shape {
public:
  /// Constructs an empty loop.
  this() { }

  /// Constructs an S2LaxLoopShape with the given vertices.
  this(in S2Point[] vertices) {
    initialize(vertices);
  }

  /// Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
  this(in S2Loop loop) {
    initialize(loop);
  }

  /// Initializes an S2LaxLoopShape with the given vertices.
  void initialize(in S2Point[] vertices) {
    _vertices = new S2Point[vertices.length];
    copy(vertices, _vertices);
  }

  /**
   * Initializes an S2LaxLoopShape from the given S2Loop, by copying its data.
   *
   * REQUIRES: !loop->is_full()
   *           [Use S2LaxPolygonShape if you need to represent a full loop.]
   */
  void initialize(in S2Loop loop)
  in {
    assert(!loop.isFull(), "Full loops not supported; use S2LaxPolygonShape");
  } do {
    if (loop.isEmpty()) {
      _vertices.length = 0;
    } else {
      _vertices = new S2Point[loop.numVertices()];
      copy(loop.vertices(), _vertices);
    }
  }

  int numVertices() const {
    return cast(int) _vertices.length;
  }

  const(S2Point) vertex(int i) const {
    return _vertices[i];
  }

  /// S2Shape interface:
  final override
  int numEdges() const {
    return numVertices();
  }

  final override
  Edge edge(int e0) const
  in {
    assert(e0 < numEdges());
  } do {
    int e1 = e0 + 1;
    if (e1 == numVertices()) e1 = 0;
    return Edge(_vertices[e0], _vertices[e1]);
  }

  /// Not final; overridden by S2LaxClosedPolylineShape.
  override
  int dimension() const {
    return 2;
  }

  /// Not final; overridden by S2LaxClosedPolylineShape.
  override
  ReferencePoint getReferencePoint() const {
    // GetReferencePoint interprets a loop with no vertices as "full".
    if (numVertices() == 0) return ReferencePoint(false);
    return .getReferencePoint(this);
  }

  final override
  int numChains() const {
    return min(1, cast(int) _vertices.length);
  }

  final override
  Chain chain(int i) const {
    return Chain(0, cast(int) _vertices.length);
  }

  final override
  Edge chainEdge(int i, int j) const
  in {
    assert(i == 0);
    assert(j < numEdges());
  } do {
    int k = (j + 1 == numVertices()) ? 0 : j + 1;
    return Edge(_vertices[j], _vertices[k]);
  }

  final override
  ChainPosition chainPosition(int e) const {
    return ChainPosition(0, e);
  }

private:
  // For clients that have many small loops, we save some memory by
  // representing the vertices as an array rather than using std::vector.
  int _numVertices;
  S2Point[] _vertices;
}

/**
 * S2LaxClosedPolylineShape is like S2LaxPolylineShape except that the last
 * vertex is implicitly joined to the first.  It is also like S2LaxLoopShape
 * except that it does not have an interior (which makes it more efficient to
 * index).
 */
class S2LaxClosedPolylineShape : S2LaxLoopShape {
public:

  this() {
    super();
  }

  this(in S2Point[] vertices) {
    super(vertices);
  }

  this(in S2Loop loop) {
    super(loop);
  }

  final override
  int dimension() const {
    return 1;
  }

  final override
  ReferencePoint getReferencePoint() const {
    return ReferencePoint(false);
  }
}

/**
 * S2VertexIdLaxLoopShape is just like S2LaxLoopShape, except that vertices are
 * specified as indices into a vertex array.  This representation can be more
 * compact when many loops are arranged in a mesh structure.
 */
class S2VertexIdLaxLoopShape : S2Shape {
public:
  // Constructs an empty loop.
  this() { }

  // Constructs the shape from the given vertex array and indices.
  // "vertex_ids" is a vector of indices into "vertex_array".
  //
  // ENSURES:  loop->vertex(i) == (*vertex_array)[vertex_ids[i]]
  // REQUIRES: "vertex_array" persists for the lifetime of this object.
  this(in int[] vertex_ids, in S2Point[] vertex_array) {
      initialize(vertex_ids, vertex_array);
  }

  // Initializes the shape from the given vertex array and indices.
  // "vertex_ids" is a vector of indices into "vertex_array".
  void initialize(in int[] vertex_ids, in S2Point[] vertex_array) {
    _vertexIds = new int[vertex_ids.length];
    copy(vertex_ids, _vertexIds);
    _vertexArray = vertex_array;
  }

  /// Returns the number of vertices in the loop.
  int numVertices() const {
    return cast(int) _vertexIds.length;
  }

  int vertexId(int i) const {
    return _vertexIds[i];
  }

  const(S2Point) vertex(int i) const {
    return _vertexArray[vertexId(i)];
  }

  /// S2Shape interface:
  final override
  int numEdges() const {
    return numVertices();
  }

  final override
  Edge edge(int e0) const
  in {
    assert(e0 < numEdges());
  } do {
    int e1 = e0 + 1;
    if (e1 == numVertices()) e1 = 0;
    return Edge(vertex(e0), vertex(e1));
  }

  final override
  int dimension() const {
    return 2;
  }

  final override
  ReferencePoint getReferencePoint() const {
    // GetReferencePoint interprets a loop with no vertices as "full".
    if (numVertices() == 0) return ReferencePoint(false);
    return .getReferencePoint(this);
  }

  final override
  int numChains() const {
    return 1;
  }

  final override
  Chain chain(int i) const {
    return Chain(0, cast(int) _vertexIds.length);
  }

  final override
  Edge chainEdge(int i, int j) const
  in {
    assert(i == 0);
    assert(j < numEdges());
  } do {
    int k = (j + 1 == numVertices()) ? 0 : j + 1;
    return Edge(vertex(j), vertex(k));
  }

  final override
  ChainPosition chainPosition(int e) const {
    return ChainPosition(0, e);
  }

 private:
  int[] _vertexIds;
  Rebindable!(const(S2Point[])) _vertexArray;
}
