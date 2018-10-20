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

module s2.shapeutil.count_edges_test;

import s2.shapeutil.count_edges : countEdges, countEdgesUpTo;
import s2.mutable_s2shape_index;
//import s2.s2text_format : makeIndex;

import fluent.asserts : Assert;

// @("CountEdgesUpTo.StopsEarly") unittest {
//   auto index = makeIndex("0:0 | 0:1 | 0:2 | 0:3 | 0:4 # 1:0, 1:1 | 1:2, 1:3 | 1:4, 1:5, 1:6 #");
//   // Verify the test parameters.
//   Assert.equal(index.numShapeIds(), 4);
//   Assert.equal(index.shape(0).numEdges(), 5);
//   Assert.equal(index.shape(1).numEdges(), 1);
//   Assert.equal(index.shape(2).numEdges(), 1);
//   Assert.equal(index.shape(3).numEdges(), 2);

//   Assert.equal(countEdges(index), 9);
//   Assert.equal(countEdgesUpTo(index, 1), 5);
//   Assert.equal(countEdgesUpTo(index, 5), 5);
//   Assert.equal(countEdgesUpTo(index, 6), 6);
//   Assert.equal(countEdgesUpTo(index, 8), 9);
// }
