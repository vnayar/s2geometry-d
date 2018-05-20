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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar

module s2.shapeutil.shape_edge_id;

import std.format : format;

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32, int32> except that it has named fields.
// It should be passed and returned by value.
struct ShapeEdgeId {
public:
  int shapeId = -1;
  int edgeId = -1;

  int opCmp(ref const ShapeEdgeId other) const {
    if (shapeId > other.shapeId) return 1;
    if (shapeId < other.shapeId) return -1;
    if (edgeId > other.edgeId) return 1;
    if (edgeId < other.edgeId) return -1;
    return 0;
  }

  string toString() const {
    return format("%d:%d", shapeId, edgeId);
  }

  size_t toHash() const nothrow @trusted {
    // The following preserves all bits even when edgeId < 0.
    return (cast(size_t) shapeId << 32) | edgeId;
  }
}
