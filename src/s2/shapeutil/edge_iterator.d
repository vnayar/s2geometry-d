// Copyright 2018 Google Inc. All Rights Reserved.
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

// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.shapeutil.edge_iterator;

import s2.s2shape;
import s2.s2shape_index;
import s2.shapeutil.shape_edge_id;

import std.format;

/**
 * An iterator that advances through all edges in an S2ShapeIndex.
 *
 * Example usage:
 *
 * for (EdgeIterator it(index); !it.Done(); it.Next()) {
 *   auto edge = it.edge();
 *   //...
 * }
 */
struct EdgeIterator {
public:
  this(S2ShapeIndex index) {
    _index = index;
    _shapeId = -1;
    _numEdges = 0;
    _edgeId = -1;
    next();
  }

  // Returns the current shape id.
  int shapeId() const {
    return _shapeId;
  }

  // Returns the current edge id.
  int edgeId() const {
    return _edgeId;
  }

  // Returns the current (shape_id, edge_id).
  ShapeEdgeId shapeEdgeId() const {
    return ShapeEdgeId(_shapeId, _edgeId);
  }

  // Returns the current edge.
  S2Shape.Edge edge() const
  in {
    assert(!done());
  } do {
    return _index.shape(_shapeId).edge(_edgeId);
  }

  // Returns true if there are no more edges in the index.
  bool done() const {
    return shapeId() >= _index.numShapeIds();
  }

  // Advances to the next edge.
  void next() {
    while (++_edgeId >= _numEdges) {
      if (++_shapeId >= _index.numShapeIds()) break;
      S2Shape shape = _index.shape(_shapeId);
      _numEdges = (shape is null) ? 0 : shape.numEdges();
      _edgeId = -1;
    }
  }

  bool opEquals(in EdgeIterator v) const {
    return _index is v._index && _shapeId == v._shapeId && _edgeId == v._edgeId;
  }

  string debugString() const {
    return format("(shape=%d, edge=%d)", _shapeId, _edgeId);
  }

 private:
  S2ShapeIndex _index;
  int _shapeId;
  int _numEdges;
  int _edgeId;
}
