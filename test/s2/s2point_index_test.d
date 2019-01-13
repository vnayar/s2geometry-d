// Copyright 2015 Google Inc. All Rights Reserved.
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

module s2.s2point_index_test;

import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2point;
import s2.s2point_index;
import s2.s2testing;
import s2.util.container.btree;

import fluent.asserts;

import std.stdio;


class S2PointIndexTest {
private:
  alias Index = S2PointIndex!int;
  alias PointData = Index.PointData;
  alias Contents = BTree!PointData;
  Index _index;
  Contents _contents;

public:
  this() {
    _index = new S2PointIndex!int();
    _contents = new Contents();
  }

  void add(in S2Point point, int data) {
    _index.add(point, data);
    _contents.insert(new PointData(point, data));
  }

  void remove(in S2Point point, int data) {
    _index.remove(point, data);
    // If there are multiple copies, remove only one.
    _contents.remove(new PointData(point, data));
  }

  void verify() {
    Contents remaining = _contents;
    for (auto it = Index.Iterator(_index); !it.done(); it.next()) {
      Contents.Iterator element = remaining.equalRange(it.pointData()).toIterator();
      Assert.notEqual(element, remaining.end());
      remaining.remove(element.getValue());
    }
    Assert.equal(remaining.empty(), true);
  }

  void checkIteratorMethods() {
    auto it = Index.Iterator(_index);
    Assert.equal(it.prev(), false);
    it.finish();
    Assert.equal(it.done(), true);

    // Iterate through all the cells in the index.
    S2CellId prev_cellid = S2CellId.none();
    S2CellId min_cellid = S2CellId.begin(S2CellId.MAX_LEVEL);
    for (it.begin(); !it.done(); it.next()) {
      S2CellId cellid = it.id();
      Assert.equal(cellid, S2CellId(it.point()));

      auto it2 = Index.Iterator(_index);
      if (cellid == prev_cellid) {
        it2.seek(cellid);
      }

      // Generate a cellunion that covers the range of empty leaf cells between
      // the last cell and this one.  Then make sure that seeking to any of
      // those cells takes us to the immediately following cell.
      auto skipped = S2CellUnion.fromBeginEnd(min_cellid, cellid.rangeMin());
      foreach (S2CellId skipped_id; skipped.cellIds()) {
        it2.seek(skipped_id);
        Assert.equal(cellid, it2.id());
      }
      // Test Prev(), Next(), and Seek().
      if (prev_cellid.isValid()) {
        it2 = it;
        Assert.equal(it2.prev(), true);
        Assert.equal(prev_cellid, it2.id());
        it2.next();
        Assert.equal(cellid, it2.id());
        it2.seek(prev_cellid);
        Assert.equal(prev_cellid, it2.id());
      }
      prev_cellid = cellid;
      min_cellid = cellid.rangeMax().next();
    }
  }
}

@("S2PointIndexTest.NoPoints") unittest {
  auto t = new S2PointIndexTest();
  t.checkIteratorMethods();
}

@("S2PointIndexTest.RandomPoints") unittest {
  auto t = new S2PointIndexTest();
  foreach (int i; 0..1000) {
    t.add(S2Testing.randomPoint(), S2Testing.rnd.uniform(100));
  }
  t.verify();
  t.checkIteratorMethods();
}
