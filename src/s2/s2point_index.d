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

module s2.s2point_index;

import s2.s2cell_id;
import s2.util.container.btree_map;


/**
 * S2PointIndex maintains an index of points sorted by leaf S2CellId.  Each
 * point has some associated client-supplied data, such as an integer or
 * pointer.  This can be used to map results back to client data structures.
 *
 * The class supports adding or removing points dynamically, and provides a
 * seekable iterator interface for navigating the index.
 *
 * You can use this class in conjuction with S2ClosestPointQuery to find the
 * closest index points to a given query point.  For example:
 *
 * void Test(const vector<S2Point>& points, const S2Point& target) {
 *   // The template argument allows auxiliary data to be attached to each
 *   // point (in this case, the array index).
 *   S2PointIndex<int> index;
 *   for (int i = 0; i < points.size(); ++i) {
 *     index.Add(points[i], i);
 *   }
 *   S2ClosestPointQuery<int> query(&index);
 *   query.FindClosestPoint(target);
 *   if (query.num_points() > 0) {
 *     // query.point(0) is the closest point (result 0).
 *     // query.distance(0) is the distance to the target.
 *     // query.data(0) is the auxiliary data (the array index set above).
 *     DoSomething(query.point(0), query.data(0), query.distance(0));
 *   }
 * }
 *
 * Alternatively, you can access the index directly using the iterator
 * interface.  For example, here is how to iterate through all the points in a
 * given S2CellId "target_id":
 *
 *   S2PointIndex<int>::Iterator it(&index);
 *   it.Seek(target_id.range_min());
 *   for (; !it.done() && it.id() <= target_id.range_max(); it.Next()) {
 *     DoSomething(it.id(), it.point(), it.data());
 *   }
 *
 * Points can be added or removed from the index at any time by calling Add()
 * or Remove().  However when the index is modified, you must call Init() on
 * each iterator before using it again (or simply create a new iterator).
 *
 *   index.Add(new_point, 123456);
 *   it.Init(&index);
 *   it.Seek(target.range_min());
 *
 * TODO(ericv): Make this a subtype of S2Region, so that it can also be used
 * to efficiently compute coverings of a collection of S2Points.
 *
 * REQUIRES: "Data" has default and copy constructors.
 * REQUIRES: "Data" has operator== and operator<.
 */
class S2PointIndex(DataT) {
public:
  // PointData is essentially std::pair with named fields.  It stores an
  // S2Point and its associated client data.
  static struct PointData {
  public:

    this() {}  // Needed by STL

    this(in S2Point point, in DataT data) {
      _point = point;
      _data = data;
    }

    S2Point point() const {
      return _point;
    }

    const(DataT) data() const {
      return _data;
    }

    bool opEquals(PointData other) const {
      if (other is null) return false;
      return _point == other._point && _data == other._data;
    }

    // Not required by S2PointIndex but useful for tests.
    int opCmp(PointData other) const {
      if (_point < other._point) return -1;
      if (_point != other._point) return 1;
      if (_data < other._data) return -1;
      if (_data != other._data) return 1;
      return 0;
    }

  private:
    S2Point _point;
    DataT _data;
  }

  // Default constructor.
  this() { }

  // Returns the number of points in the index.
  int numPoints() const {
    return _map.length;
  }

  // Adds the given point to the index.  Invalidates all iterators.
  void add(in S2Point point, in DataT data) {
    add(PointData(point, data));
  }

  void add(in PointData point_data) {
    auto id = S2CellId(point_data.point());
    _map.insert(id, point_data);
  }

  // Removes the given point from the index.  Both the "point" and "data"
  // fields must match the point to be removed.  Returns false if the given
  // point was not present.  Invalidates all iterators.
  bool remove(in S2Point point, in DataT data) {
    return remove(PointData(point, data));
  }

  bool remove(in PointData point_data) {
    auto id = S2CellId(point_data.point());
    foreach (Map.Pair entry; _map.equalRange(id)) {
      if (entry.value == point_data) {
        _map.remove(point_data);
        return true;
      }
    }
    return false;
  }

  // Resets the index to its original empty state.  Invalidates all iterators.
  void clear() {
    _map.clear();
  }


private:
  // Defined here because the Iterator class below uses it.
  alias Map = BTreeMap!(S2CellId, PointData);

public:
  static class Iterator {
  public:
    // Default constructor; must be followed by a call to Init().
    this() {
      _map = null;
    }

    // Convenience constructor that calls Init().
    this(in S2PointIndex index) {
      init(index);
    }

    // Initializes an iterator for the given S2PointIndex.  If the index is
    // non-empty, the iterator is positioned at the first cell.
    //
    // This method may be called multiple times, e.g. to make an iterator
    // valid again after the index is modified.
    void initialize(in S2PointIndex index) {
      _map = index._map;
      _iter = _map.begin();
      _end = _map.end();
    }

    // The S2CellId for the current index entry.
    // REQUIRES: !done()
    S2CellId id() const
    in {
      assert(!done());
    } body {
      return _iter.getKey();
    }

    // The point associated with the current index entry.
    // REQUIRES: !done()
    const(S2Point) point() const
    in {
      assert(!done());
    } body {
      return _iter.getValue().point();
    }

    // The client-supplied data associated with the current index entry.
    // REQUIRES: !done()
    const(DataT) data() const
    in {
      assert(!done());
    } body {
      return _iter.getValue().data();
    }


    // The (S2Point, data) pair associated with the current index entry.
    const(PointData) pointData() const
    in {
      assert(!done());
    } body {
      return _iter.getValue();
    }

    // Returns true if the iterator is positioned past the last index entry.
    bool done() const {
      return _iter == _end;
    }

    // Positions the iterator at the first index entry (if any).
    void begin() {
      _iter = _map.begin();
    }

    // Positions the iterator so that done() is true.
    void finish() {
      _iter = _end;
    }

    // Advances the iterator to the next index entry.
    // REQUIRES: !done()
    void next()
    in {
      assert(!done());
    } body {
      ++_iter;
    }

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    bool prev() {
      if (_iter == _map.begin()) return false;
      --_iter;
      return true;
    }

    // Positions the iterator at the first entry with id() >= target, or at the
    // end of the index if no such entry exists.
    void seek(S2CellId target) {
      _iter = _map.equalRange().toIterator();
    }

  private:
    const(Map) _map;
    Map.Iterator _iter, _end;
  }

 private:
  Map _map;
}
