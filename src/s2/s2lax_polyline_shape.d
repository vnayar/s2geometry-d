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

module s2.s2lax_polyline_shape;

import s2.s2polyline;
import s2.s2shape;
import s2.logger;
import s2.s2point;

import std.algorithm;

/**
 * S2LaxPolylineShape represents a polyline.  It is similar to
 * S2Polyline::Shape except that duplicate vertices are allowed, and the
 * representation is slightly more compact.
 *
 * Polylines may have any number of vertices, but note that polylines with
 * fewer than 2 vertices do not define any edges.  (To create a polyline
 * consisting of a single degenerate edge, either repeat the same vertex twice
 * or use S2LaxClosedPolylineShape defined in s2_lax_loop_shape.h.)
 */
class S2LaxPolylineShape : S2Shape {
public:
  // Constructs an empty polyline.
  this() { }

  // Constructs an S2LaxPolylineShape with the given vertices.
  this(S2Point[] vertices) {
    init(vertices);
  }

  // Constructs an S2LaxPolylineShape from the given S2Polyline, by copying
  // its data.
  this(in S2Polyline polyline) {
    init(polyline);
  }

  // Initializes an S2LaxPolylineShape with the given vertices.
  void init(S2Point[] vertices) {
    if (vertices.length == 1)
      logger.logWarn("s2shapeutil::S2LaxPolylineShape with one vertex has no edges");
    _vertices = vertices;
  }

  // Initializes an S2LaxPolylineShape from the given S2Polyline, by copying
  // its data.
  void init(in S2Polyline polyline) {
    if (polyline.vertices().length == 1)
        logger.logWarn("s2shapeutil::S2LaxPolylineShape with one vertex has no edges");
    _vertices = polyline.vertices().dup;
  }

  int numVertices() const {
    return cast(int) _vertices.length;
  }

  const(S2Point) vertex(int i) const {
    return _vertices[i];
  }

  const(S2Point[]) vertices() const {
    return _vertices;
  }

  // S2Shape interface:
  final override
  int numEdges() const {
    return max(0, numVertices() - 1);
  }

  final override
  Edge edge(int e) const
  in {
    assert(e < numEdges());
  } do {
    return Edge(_vertices[e], _vertices[e + 1]);
  }

  final override
  int dimension() const {
    return 1;
  }

  final override
  ReferencePoint getReferencePoint() const {
    return ReferencePoint(false);
  }

  final override
  int numChains() const {
    return min(1, numEdges());
  }

  final override
  Chain chain(int i) const {
    return Chain(0, numEdges());
  }

  final override
  Edge chainEdge(int i, int j) const
  in {
    assert(i == 0);
    assert(j < numEdges());
  } do {
    return Edge(_vertices[j], _vertices[j + 1]);
  }

  final override
  ChainPosition chainPosition(int e) const {
    return ChainPosition(0, e);
  }

 private:
  // For clients that have many small polylines, we save some memory by
  // representing the vertices as an array rather than using std::vector.
  S2Point[] _vertices;
}
