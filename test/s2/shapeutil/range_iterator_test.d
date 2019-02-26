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

module s2.shapeutil.range_iterator_test;

import s2.shapeutil.range_iterator;
import s2.mutable_s2shape_index;
import s2.s2text_format;
import s2.s2cell_id;

import fluent.asserts;

@("RangeIterator.Next") unittest {
  // Create an index with one point each on S2CellId faces 0, 1, and 2.
  auto index = makeIndexOrDie("0:0 | 0:90 | 90:0 # #");
  auto it = new RangeIterator(index);
  Assert.equal(it.id().face(), 0);
  it.next();
  Assert.equal(it.id().face(), 1);
  it.next();
  Assert.equal(it.id().face(), 2);
  it.next();
  Assert.equal(it.id(), S2CellId.sentinel());
  Assert.equal(it.done(), true);
}

@("RangeIterator.EmptyIndex") unittest {
  auto empty = makeIndexOrDie("# #");
  auto non_empty = makeIndexOrDie("0:0 # #");
  auto empty_it = new RangeIterator(empty);
  auto non_empty_it = new RangeIterator(non_empty);
  Assert.equal(non_empty_it.done(), false);
  Assert.equal(empty_it.done(), true);

  empty_it.seekTo(non_empty_it);
  Assert.equal(empty_it.done(), true);

  empty_it.seekBeyond(non_empty_it);
  Assert.equal(empty_it.done(), true);

  empty_it.seekTo(empty_it);
  Assert.equal(empty_it.done(), true);

  empty_it.seekBeyond(empty_it);
  Assert.equal(empty_it.done(), true);
}
