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

module s2.s2point_vector_shape;

import s2.s2shape;
import s2.s2point;

// S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is reprsented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.
class S2PointVectorShape : S2Shape {
public:
  // Constructs an empty point vector.
  this() {}

  // Constructs an S2PointVectorShape from a vector of points.
  this(S2Point[] points) {
    _points = points;
  }

  int numPoints() const {
    return cast(int)(_points.length);
  }

  S2Point point(int i) const {
    return _points[i];
  }

  // S2Shape interface:
  final override
  int numEdges() const {
    return cast(int)(_points.length);
  }

  final override
  Edge edge(int e) const {
    return Edge(_points[e], _points[e]);
  }

  final override
  int dimension() const {
    return 0;
  }

  final override
  ReferencePoint getReferencePoint() const {
    return ReferencePoint(false);
  }

  final override
  int numChains() const {
    return cast(int)(_points.length);
  }

  final override
  Chain chain(int i) const {
    return Chain(i, 1);
  }

  final override
  Edge chainEdge(int i, int j) const
  in {
    assert(j == 0);
  } do {
    return Edge(_points[i], _points[i]);
  }

  final override
  ChainPosition chainPosition(int e) const {
    return ChainPosition(e, 0);
  }

private:
  S2Point[] _points;
}
