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

module s2.shapeutil.range_iterator;

import s2.s2cell_id;
import s2.s2shape_index;

// RangeIterator is a wrapper over S2ShapeIndex::Iterator with extra methods
// that are useful for merging the contents of two or more S2ShapeIndexes.
class RangeIterator {
private:
  // Updates internal state after the iterator has been repositioned.
  void refresh() {
    _rangeMin = id().rangeMin();
    _rangeMax = id().rangeMax();
  }

  S2ShapeIndex.Iterator _it;
  S2CellId _rangeMin, _rangeMax;

public:
  // Construct a new RangeIterator positioned at the first cell of the index.
  this(in S2ShapeIndex index) {
    refresh();
  }

  // The current S2CellId and cell contents.
  S2CellId id() const {
    return _it.id();
  }

  const(S2ShapeIndexCell) cell() const {
    return _it.cell();
  }

  // The min and max leaf cell ids covered by the current cell.  If done() is
  // true, these methods return a value larger than any valid cell id.
  S2CellId rangeMin() const {
    return _rangeMin;
  }
  S2CellId rangeMax() const {
    return _rangeMax;
  }

  void next() {
    _it.next();
    refresh();
  }

  bool done() {
    return _it.done();
  }

  // Position the iterator at the first cell that overlaps or follows
  // "target", i.e. such that range_max() >= target.range_min().
  void seekTo(in RangeIterator target) {
    _it.seek(target.rangeMin());
    // If the current cell does not overlap "target", it is possible that the
    // previous cell is the one we are looking for.  This can only happen when
    // the previous cell contains "target" but has a smaller S2CellId.
    if (_it.done() || _it.id().rangeMin() > target.rangeMax()) {
      if (_it.prev() && _it.id().rangeMax() < target.id()) _it.next();
    }
    refresh();
  }

  // Position the iterator at the first cell that follows "target", i.e. the
  // first cell such that range_min() > target.range_max().
  void seekBeyond(in RangeIterator target) {
    _it.seek(target.rangeMax().next());
    if (!_it.done() && _it.id().rangeMin() <= target.rangeMax()) {
      _it.next();
    }
    refresh();
  }

}
