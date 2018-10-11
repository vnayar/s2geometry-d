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

module s2.shapeutil.count_edges;

import s2.s2shape_index;
import s2.s2shape;

// Returns the total number of edges in all indexed shapes.  This method takes
// time linear in the number of shapes.
int countEdges(S2ShapeIndexT)(in S2ShapeIndexT index) {
  return countEdgesUpTo(index, int.max);
}

// Like CountEdges(), but stops once "max_edges" edges have been found (in
// which case the current running total is returned).
int countEdgesUpTo(S2ShapeIndexT)(in S2ShapeIndexT index, int max_edges) {
  const int num_shape_ids = index.numShapeIds();
  int num_edges = 0;
  for (int s = 0; s < num_shape_ids; ++s) {
    S2Shape shape = index.shape(s);
    if (shape is null) continue;
    num_edges += shape.numEdges();
    if (num_edges >= max_edges) break;
  }
  return num_edges;
}
