module s2.util.container.btree;

import std.algorithm : max;
import std.functional : binaryFun;
import std.traits : ReturnType;
import std.format : format;


/**
 * A B-Tree implementation based upon "Introduction to Algorithms" by Cormen, Leiserson, Rivest,
 * and Stein.
 *
 * B-Trees are both smaller and faster than most implementations of set/map. The red-black tree
 * implementation of a set/map typically has an overhead of 3 pointers (left, right, and parent)
 * plus the node color information for each stored value.  This B-Tree implementation stores
 * multiple values on fixed size nodes (usually 256 bytes) and does not store child pointers for
 * leaf nodes.
 *
 * The packing of multiple values into each node of a B-Tree also improves cache locality which
 * translates into faster operations.
 *
 * Because nodes contain both node pointers and values, a BTree will not have good performance
 * when used with large types that are passed by value, such as large structs.
 *
 * As always, using the BTree with class objects, which are passed by reference, also avoids
 * storage in the BTree itself.
 *
 * Params:
 *   T = The element type to be organized by the BTree.
 *   NodeSizeV = The size in bytes of a node in the BinaryTree. Values are chosen to optimize
 *     performance depending on usage. For example, a value of 4096 bytes may be useful when
 *     retrieving items from disk, or the value 256 bytes to assure good use of CPU caches.
 *     The default value is 256 bytes.
 *   ValueLessF = A less-than comparison function for values, either a function or a string
 *     representing a function as described by
 *     $(LINK2 https://dlang.org/phobos/std_functional.html#binaryFun, std.functional.binaryFun).
 */
final class BTree(ValueT, size_t NodeSizeV = 256, alias ValueLessF = "a < b")
if (is(typeof(binaryFun!ValueLessF(ValueT.init, ValueT.init)) : bool)) {

private:
  Node _root;
  size_t _length;

  // TODO: Determine how to enforce const functions are passed in, which is
  // currently preventing the use of const with reference types.
  alias _isValueLess = binaryFun!ValueLessF;

  static bool _isValueEqual(in ValueT v1, in ValueT v2) {
    return !_isValueLess(v1, v2) && !_isValueLess(v2, v1);
  }

  // The B-Tree is defined by a minimum degree "t".
  // Each node other than the root must have at least t-1 values.
  // Each node may have up to 2t - 1 values, thus an internal node may have up to 2*t children.
  //
  // Thus, the maximum size of a non-leaf node is:
  //   [values]    (2t-1) * (size of a value)
  //   [numValues] + 1 * (size of size_t)
  //   [isLeaf]    + 1 * (size of a bool)
  //   [children]  + 2*t * (size of pointer to a node)
  //   = NodeSize
  static size_t getMinDegree() @nogc @safe pure nothrow {
    size_t nodeSizeBase = bool.sizeof + Node.sizeof + size_t.sizeof + size_t.sizeof;
    size_t nodeExtraSizePerChild = ValueT.sizeof + Node.sizeof;

    if (NodeSizeV < nodeSizeBase) {
      return 1;
    }
    size_t maxChildren = (NodeSizeV - nodeSizeBase) / nodeExtraSizePerChild;
    if (maxChildren < 2) {
      return 1;
    } else {
      return maxChildren / 2;
    }
  }

public:

  enum MIN_DEGREE = getMinDegree();
  enum MAX_DEGREE = MIN_DEGREE * 2;

  this() {
    // TODO: ALLOCATE-NODE() which allocates a disk page in O(1).
    _root = new Node();
    _length = 0;
    // TODO: DISK-WRITE(_root)
  }

  /**
   * The root Node of the B-Tree.
   */
  @property
  inout(Node) root() inout {
    return _root;
  }

  @property
  size_t length() const {
    return _length;
  }

  ////
  // Iterator Style Methods
  ////

  Iterator begin() {
    return Iterator(_root.leftmost(), 0);
  }

  Iterator end() {
    Node right = _root.rightmost();
    return Iterator(right, right.numValues());
  }

  ////
  // Range Style Methods
  ////

  /**
   * A full slice for expressions of the form: myBTree[]
   */
  Range opIndex() {
    return Range(begin(), end());
  }

  /**
   * Defines the structure of a slice for slice expression of the form: myBTree[a..b]
   */
  ValueT[2] opSlice(size_t dim)(ValueT a, ValueT b) {
    return [a, b];
  }

  Range opIndex(ValueT[2] slice) {
    return Range(_root.findFirstGE(slice[0]), _root.findFirstGE(slice[1]));
  }

  /**
   * Find a range of elements whose value is equal to 'v'.
   *
   * If no elements are found, the result will be an empty range. If there is more than one
   * element with the same value, the range will have more than 1 element.
   */
  Range equalRange(ValueT v) {
    return Range(_root.findFirstGE(v), _root.findFirstGT(v));
  }

  /**
   * Returns a range of all elements strictly less than 'v'.
   */
  Range lowerRange(ValueT v) {
    return Range(begin(), _root.findFirstGE(v));
  }

  /**
   * Returns a range of all elements strictly greater than 'v'.
   */
  Range upperRange(ValueT v) {
    return Range(_root.findFirstGT(v), end());
  }

  ////
  // Core BTree algorithms
  ////

  void clear() {
    _root = new Node();
    _length = 0;
  }

  /**
   * Inserts a new value into the BTree.
   */
  void insert(ValueT v) {
    _length++;
    Node curRoot = _root;
    if (curRoot.isFull()) {
      Node newRoot = new Node();
      newRoot._isLeaf = false;
      newRoot._numValues = 0;
      newRoot._parent = null;
      newRoot.setChild(0, _root);
      _root = newRoot;
      _root.splitChild(0);
    }
    _root.insertNonFull(v);
  }

  /**
   * Membership
   */
  bool opBinaryRight(string op)(ValueT value)
  if (op == "in") {
    return !equalRange(value).empty();
  }

  /**
   * Removes a value and value from the BTree if it exists.
   */
  void remove(ValueT v) {
    bool isFound = _root.remove(v);
    if (isFound) {
      _length--;
    }
    if (!_root._isLeaf && _root._numValues == 0) {
      _root = _root._children[0];
      _root._parent = null;
    }
  }

  /**
   * Lower level iterators on the tree.
   *
   * While not in the D-style, these methods may be used for easier compatibility when working
   * with code bases not in the D Language which are based on iterators.
   */
  static struct Iterator {
  private:
    Node _node;
    size_t _position;
  public:

    bool isFound() {
      return _node !is null;
    }

    inout(ValueT) getValue() inout {
      return _node.getValue(_position);
    }

    void increment() {
      if (_node.isLeaf() && ++_position < _node.numValues()) {
        return;
      }
      // We've gone past the last position.
      if (_node.isLeaf()) {
        Iterator save = this;
        while (_position == _node.numValues() && !_node.isRoot()) {
          _position = _node._position;
          _node = _node._parent;
        }
        if (_position == _node.numValues()) {
          this = save;
        }
      } else {
        _node = _node._children[_position + 1];
        while (!_node.isLeaf()) {
          _node = _node._children[0];
        }
        _position = 0;
      }
    }

    void opUnary(string op)()
    if (op == "++") {
      increment();
    }

    void decrement() {
      if (_node.isLeaf() && _position > 0) {
        _position--;
        return;
      }
      if (_node.isLeaf()) {
        Iterator save = this;
        while (_position == 0 && !_node.isRoot()) {
          _position = _node._position;
          _node = _node._parent;
        }
        _position--;
        if (_position < 0) {
          this = save;
        }
      } else {
        _node = _node._children[_position];
        while (!_node.isLeaf()) {
          _node = _node._children[_node._numValues];
        }
        _position = _node._numValues - 1;
      }
    }

    void opUnary(string op)()
    if (op == "--") {
      decrement();
    }
  }

  /**
   * D-style Ranges
   */
  static struct Range {
  private:
    Iterator _begin;
    Iterator _end;

  public:
    bool empty() {
      return _begin == _end;
    }

    inout(ValueT) front() inout {
      return _begin.getValue();
    }

    void popFront()
    in {
      assert(!empty());
    } body {
      _begin.increment();
    }

    Iterator toIterator() {
      return _begin;
    }
  }

  /**
   * A node in the BTree.
   *
   * A notable feature is that the values that are inserted are stored inside the BTree itself
   * in the case of value data types, such as integers, floats, static arrays, or structs. In
   * other cases, such as dynamic arrays and classes, on the reference is stored.
   */
  static class Node {
  package:
    /// Indicates if the node has reached the size limit specified in $(D_INLINECODE NodeSizeV).
    bool isFull() {
      return _numValues >= MAX_DEGREE - 1;
    }

    void setChild(size_t i, Node child)
    in {
      assert(!_isLeaf);
      assert(i >= 0);
    } body {
      _children[i] = child;
      _children[i]._position = i;
      _children[i]._parent = this;
    }

    Node leftmost() {
      Node seek = this;
      while (!seek._isLeaf) {
        seek = seek._children[0];
      }
      return seek;
    }

    Node rightmost() {
      Node seek = this;
      while (!seek._isLeaf) {
        seek = seek._children[seek._numValues];
      }
      return seek;
    }

  package:
    // A flag indicating whether the node is a leaf node or not.
    bool _isLeaf = true;

    // A reference to the node's parent.
    Node _parent = null;

    // The position of the node in the node's parent.
    size_t _position = 0;

    // A count of the number of values in the node.
    size_t _numValues = 0;
    ValueT[MAX_DEGREE - 1] _values;

    // Only non-leaf (internal) nodes have children.
    Node[MAX_DEGREE] _children;

    invariant {
      // Assure that the number of values stays within the bounds of the MAX_DEGREE.
      assert(_numValues <= MAX_DEGREE - 1);

      // Assure that values remain in a consistent order.
      foreach (i; 1 .. _numValues) {
        assert(!_isValueLess(_values[i], _values[i - 1]));
      }

      // Assure that no child that should be present is null.
      if (!_isLeaf) {
        foreach (i; 0 .. _numValues + 1) {
          assert(_children[i] !is null);
          assert(_children[i]._parent is this);
          assert(
              _children[i]._position == i,
              format("Found position %d, expected %d in child with values %s. _numValues=%d",
                  _children[i]._position, i, _children[i]._values, _numValues));
        }
      }
    }

  public:
    /**
     * Given that this node is non-full, but a child node that is, split the child node into two
     * separate nodes that are half full, and insert a value into this node between them.
     */
    void splitChild(size_t i)
    in {
      assert(!isFull(), "this = " ~ toString());
      assert(_children[i].isFull(), format("_children[%d] = %s", i, _children[i].toString()));
    } out {
      assert(_children[i]._numValues == MIN_DEGREE - 1);
      assert(_children[i+1]._numValues == MIN_DEGREE - 1);
    } body {
      Node toSplitNode = _children[i];

      // Prepare a new node containing the right half of the node at _children[i].
      // TODO: ALLOCATE-NODE(newNode)
      Node newNode = new Node();
      newNode._parent = this;
      newNode._isLeaf = toSplitNode._isLeaf;
      newNode._numValues = MIN_DEGREE - 1;
      foreach (j; 0 .. MIN_DEGREE - 1) {
        newNode._values[j] = toSplitNode._values[j + MIN_DEGREE];
      }
      if (!toSplitNode._isLeaf) {
        foreach (j; 0 .. MIN_DEGREE) {
          newNode.setChild(j, toSplitNode._children[j + MIN_DEGREE]);
        }
      }

      // Make way for the new value added at position i.
      for (auto j = _numValues; j > i; j--) {
        _values[j] = _values[j - 1];
      }
      _values[i] = toSplitNode._values[MIN_DEGREE - 1];
      _numValues++;

      // Make way for a newNode to be added at position i + 1.
      // There can be up to _numValues + 1 children.
      for (auto j = _numValues; j > i + 1; j--) {
        setChild(j, _children[j - 1]);
      }
      setChild(i + 1, newNode);

      // Reduce the size of values in the node being split, cutting it in half.
      toSplitNode._numValues = MIN_DEGREE - 1;

      // TODO: DISK-WRITE(toSplitNode)
      // TODO: DISK-WRITE(newNode)
      // TODO: DISK-WRITE(this)
    }

    /**
     * Inserts a new value into the BTree, provided that this node is not already full.
     */
    void insertNonFull(ValueT v)
    in {
      assert(!isFull());
    } body {
      int i = cast(int) _numValues - 1;
      if (_isLeaf) {
        // Shift over the values to make room.
        for (; i >= 0 && _isValueLess(v, _values[i]); i--) {
          _values[i + 1] = _values[i];
        }
        _values[i + 1] = v;
        _numValues++;
        // TODO: DISK-WRITE(this)
      } else {
        // Find the index of the first value less than the inserted value.
        for (; i >= 0 && _isValueLess(v, _values[i]); i--) {}
        // The child to insert into is one past the value's index, which may be -1.
        //     v[0]    v[1]    v[2]
        // c[0]    c[1]    c[2]    c[3]
        i++;
        // TODO: DISK-READ(_children[i])
        if (_children[i].isFull()) {
          splitChild(i);
          // After splitting, the median value of the child is not in this node at index i.
          if (_isValueLess(_values[i], v)) {
            i++;
          }
        }
        _children[i].insertNonFull(v);
      }
    }

    static if (MIN_DEGREE < 16) {
      // If degree is small, use a simple linear search.
      size_t findFirstGEIndex(ValueT v) const {
        size_t i = 0;
        while (i < _numValues && _isValueLess(_values[i], v)) {
          i++;
        }
        return i;
      }

      size_t findFirstGTIndex(ValueT v) const {
        size_t i = 0;
        while (i < _numValues && !_isValueLess(v, _values[i])) {
          i++;
        }
        return i;
      }
    } else {
      // If degree is higher, use a binary search.
      size_t findFirstGEIndex(ValueT v) const {
        size_t i = 0;
        size_t j = _numValues;
        while (i < j) {
          size_t mid = (i + j) / 2;
          const(ValueT) midValue = _values[mid];
          if ((mid == 0 || _isValueLess(_values[mid - 1], v)) && !_isValueLess(midValue, v)) {
            return mid;
          } else if (!_isValueLess(midValue, v)) {
            j = mid;
          } else {
            i = mid + 1;
          }
        }
        return _numValues;
      }

      size_t findFirstGTIndex(ValueT v) const {
        size_t i = 0;
        size_t j = _numValues;
        while (i < j) {
          size_t mid = (i + j) / 2;
          const(ValueT) midValue = _values[mid];
          if ((mid == 0 || !_isValueLess(v, _values[mid - 1])) && _isValueLess(v, midValue)) {
            return mid;
          } else if (_isValueLess(v, midValue)) {
            j = mid;
          } else {
            i = mid + 1;
          }
        }
        return _numValues;
      }
    }

    bool isRoot() {
      return _parent is null;
    }

    bool isLeaf() {
      return _isLeaf;
    }

    Node getChild(size_t i)
    in {
      assert(!_isLeaf);
      assert(i >= 0);
      assert(i <= _numValues);
    } body {
      return _children[i];
    }

    size_t numChildren() const {
      return _isLeaf ? 0 : _numValues + 1;
    }

    ValueT[] getValues() {
      return _values[0 .. _numValues];
    }

    bool remove(ValueT v) {
      size_t i = findFirstGEIndex(v);
      // 1. If the value v is in this node and this is a leaf, delete the value from this.
      if (_isLeaf) {
        if (i != _numValues && _isValueEqual(v, _values[i])) {
          foreach (j; i .. _numValues - 1) {
            _values[j] = _values[j+1];
          }
          _numValues--;
          return true;
        }
        // Otherwise the value was not in the BTree.
      }
      // 2. If the value v is in this node and this is an internal node:
      else if (i != _numValues && _isValueEqual(v, _values[i])) {
        // 2a. If child y that precedes v in this node has at least t values, then find the
        // predecessor v' of v in the subtree rooted at y. Recursively delete v' and replace
        // v by v' in this.
        if (_children[i]._numValues >= MIN_DEGREE) {
          Node predecessorNode = _children[i].rightmost();
          ValueT predecessor = predecessorNode._values[predecessorNode._numValues - 1];
          _values[i] = predecessor;
          _children[i].remove(predecessor);
        }
        // 2b. If child z that follows v in this node has at least t values, then find the
        // successor v' of v in the subtree rooted at z.  Recursively delete v', and replace
        // v by v' in this.
        else if (_children[i+1]._numValues >= MIN_DEGREE) {
          ValueT successor = _children[i+1].leftmost()._values[0];
          _values[i] = successor;
          _children[i+1].remove(successor);
        }
        // 2c. Otherwise, if both y and z have only t-1 values, merge v and all of z into y,
        // so that x loses both v and the pointer to z, and y now contains 2t-1 values.
        // Then, free z and recursively delete v from y.
        else {
          Node y = _children[i];
          Node z = _children[i + 1];
          // Add v and z to y.
          y._values[y._numValues++] = _values[i];
          foreach (j; 0 .. z._numValues) {
            y._values[y._numValues + j] = z._values[j];
          }
          // Note: All leaves have the same depth in a BTree.
          assert(y._isLeaf == z._isLeaf);
          if (!y._isLeaf) {
            foreach (j; 0 .. z._numValues + 1) {
              y.setChild(y._numValues + j, z._children[j]);
            }
          }
          y._numValues += z._numValues;

          // Remove v and z from this.
          foreach (j; i .. _numValues - 1) {
            _values[j] = _values[j + 1];
          }
          foreach (j; i + 1 .. _numValues) {
            setChild(j, _children[j + 1]);
          }
          _numValues--;

          // Recursively remove v from y.
          y.remove(v);
        }
        return true;
      }
      // 3. The internal node does not have the value, but maybe a child does.
      else if (!_isLeaf) {
        // Assure that the child to descend to has at least MIN_DEGREE values.
        // If not, it needs to be adjusted.
        if (_children[i]._numValues == MIN_DEGREE - 1) {
          // 3a1. If _children[i] has a sibling with at least MIN_DEGREE values, move one value
          // into this, and one one of this's values into _children.
          // Handle if the right sibling has an extra value.
          if (i < _numValues && _children[i + 1]._numValues >= MIN_DEGREE) {
            Node y = _children[i];
            Node z = _children[i + 1];
            y._values[y._numValues] = _values[i];
            y._numValues++;

            _values[i] = z._values[0];
            foreach (j; 0 .. z._numValues - 1) {
              z._values[j] = z._values[j + 1];
            }

            if (!y._isLeaf) {
              y.setChild(y._numValues, z._children[0]);
              foreach (j; 0 .. z._numValues) {
                z.setChild(j, z._children[j + 1]);
              }
            }
            z._numValues--;
          }
          // 3a2. Handle if the left sibling has an extra value.
          else if (i > 0 && _children[i - 1]._numValues >= MIN_DEGREE) {
            Node x = _children[i - 1];
            Node y = _children[i];

            // Make room for 1 more value at the start of y.
            for (auto j = y._numValues; j > 0; j--) {
              y._values[j] = y._values[j - 1];
            }
            y._values[0] = _values[i - 1];
            _values[i - 1] = x._values[x._numValues - 1];

            if (!x._isLeaf) {
              for (auto j = y._numValues + 1; j > 0; j--) {
                y.setChild(j, y._children[j - 1]);
              }
              y.setChild(0, x._children[x._numValues]);
            }

            y._numValues++;
            x._numValues--;
          }
          // 3b. If both siblings of the child that may have the value are of size
          // MIN_DEGREE - 1, then merge the child with one of those siblings, merging
          // a value from this which becomes the median node.
          else if ((i == _numValues || _children[i + 1]._numValues == MIN_DEGREE - 1)
              && (_children[i]._numValues == MIN_DEGREE - 1)) {
            // 3b1. First try to merge with the right node.
            if (i < _numValues) {
              Node y = _children[i];
              Node z = _children[i + 1];
              // Add the value from this into y.
              y._values[y._numValues] = _values[i];
              y._numValues++;
              // Merge the values and values from z into y also.
              foreach (j; 0 .. z._numValues) {
                y._values[y._numValues + j] = z._values[j];
              }
              if (!y._isLeaf) {
                foreach (j; 0 .. z._numValues + 1) {
                  y.setChild(y._numValues + j, z._children[j]);
                }
              }
              y._numValues += z._numValues;

              // Now remove _values[i] and _children[i+1] from this.
              foreach (j; i .. _numValues - 1) {
                _values[j] = _values[j + 1];
              }
              foreach (j; i + 1 .. _numValues) {
                setChild(j, _children[j + 1]);
              }
              _numValues--;
            }
            // 3b2. Otherwise merge with the left node.
            else {
              Node x = _children[i - 1];
              Node y = _children[i];

              // Move the value to the left of y into x.
              x._values[x._numValues] = _values[i - 1];
              x._numValues++;
              // Now merge y into x.
              foreach (j; 0 .. y._numValues) {
                x._values[x._numValues + j] = y._values[j];
              }
              if (!x._isLeaf) {
                foreach (j; 0 .. y._numValues + 1) {
                  x.setChild(x._numValues + j, y._children[j]);
                }
              }
              x._numValues += y._numValues;

              // Erase the value from this that was merged into x.
              foreach (j; i - 1 .. _numValues - 1) {
                _values[j] = _values[j + 1];
              }
              foreach (j; i .. _numValues) {
                setChild(j, _children[j + 1]);
              }
              _numValues--;
              i--;
            }
          }
        }
        return _children[i].remove(v);
      }
      return false;
    }

    override
    string toString() {
      return format("[isLeaf=%d, numValues=%d, position=%d]", _isLeaf, _numValues, _position);
    }

    /// Indicates how many values are in this node.
    size_t numValues() const {
      return _numValues;
    }

    /// Retrieves a values stored in this node.
    inout(ValueT) getValue(size_t i) inout
    in {
      assert(i >= 0);
      assert(i < _numValues, format("Value at %d is beyond numValues %d.", i, _numValues));
    } body {
      return _values[i];
    }

    /**
     * Recursively search for a given value and produce a result indicating whether a match
     * was found and what it's value is.
     */
    Iterator findFirstGE(ValueT v) {
      size_t i = findFirstGEIndex(v);
      if (i < _numValues && _isValueEqual(_values[i], v)) {
        return Iterator(this, i);
      }
      if (_isLeaf) {
        // The index 'i' is one past the end.
        return Iterator(this, i);
      } else {
        return _children[i].findFirstGE(v);
      }
    }

    /**
     * Recursively search for a given value and produce a result indicating whether a match
     * was found and what it's value is.
     */
    Iterator findFirstGT(ValueT v) {
      size_t i = findFirstGTIndex(v);

      // As a leaf, no further processing is possible.
      if (_isLeaf) {
        return Iterator(this, i);
      }

      // We have an index with the first greater-than value, but was there a lower value in
      // one of the node's children?
      Iterator childIterator = _children[i].findFirstGT(v);
      if (i < _numValues &&
          childIterator._position == childIterator._node.numValues()) {
        // The value found in this node was less than what could be found in a child node.
        return Iterator(this, i);
      } else {
        return childIterator;
      }
    }
  }
}

/// Simple use case with primitive types.
unittest {
  auto btree = new BTree!int();
  btree.insert(2);
  btree.insert(3);
  btree.insert(4);
  btree.insert(1);
  btree.insert(4);

  // equalRange() provides a range used to get the tree values equal to the given search value.
  auto r = btree.equalRange(2);
  assert(!r.empty());
  assert(r.front() == 2);

  // There can be more than one match.
  r = btree.equalRange(4);
  assert(!r.empty());
  assert(r.front() == 4);
  r.popFront();
  assert(!r.empty());
  assert(r.front() == 4);

  // If there is not match, an empty range is returned.
  r = btree.equalRange(-3);
  assert(r.empty());
}

/// Use case using comparison by a specific value.
unittest {
  struct Structo {
    int _id;
    string _data;
  }

  // Organize the BTree using the _id field as the value used for comparison.
  auto btree = new BTree!(Structo, 1024, "a._id < b._id");
  btree.insert(Structo(1, "Good Day"));
  btree.insert(Structo(2, "Guten Tag"));
  btree.insert(Structo(3, "G'Day Mate"));
  assert(!btree.equalRange(Structo(2)).empty());
  assert(btree.equalRange(Structo(4)).empty());
  assert(btree.equalRange(Structo(3)).front()._data == "G'Day Mate");
}

/// Use case using comparison by a specific value.
unittest {
  struct Structo {
    int _id;
    string _data;
  }

  // This time use the string _data field, but only the first two characters.
  auto btree2 = new BTree!(Structo, 1024, "a._data[0..2] > b._data[0..2]");

  // Lambdas may also be used.
  auto btree3 = new BTree!(
      Structo,  // The type of thing being stored in the BTree.
      1024,  // The size of a node in bytes.
      (a, b) => a._data[0..2] > b._data[0..2]);  // Determine if one value is less than another.

  btree3.insert(Structo(1, "RW-Fish"));
  btree3.insert(Structo(2, "LG-Sheep"));
  btree3.insert(Structo(3, "BK-Bunny"));
  assert(!btree3.equalRange(Structo(0, "RW")).empty());
  assert(btree3.equalRange(Structo(0, "ZM")).empty());
  assert(btree3.equalRange(Structo(0, "BK")).front()._data == "BK-Bunny");
}

/// Use case compatible with class comparing operator overrides:
unittest {
  class Thingy {
    int _a, _b;

    this(int a, int b) {
      _a = a;
      _b = b;
    }

    int opCmp(in Thingy o) const {
      return _a < o._a ? -1
          : _a > o._a ? 1
          : _b < o._b ? -1
          : _b > o._b ? 1
          : 0;
    }

    override
    bool opEquals(in Object o) const {
      Thingy other = cast(Thingy) o;
      if (other is null) return false;
      return _a == other._a && _b == other._b;
    }
  }

  auto btree = new BTree!Thingy();
  btree.insert(new Thingy(1, 2));
  btree.insert(new Thingy(1, 3));
  btree.insert(new Thingy(2, 1));
  btree.insert(new Thingy(1, 2));

  assert(!btree.equalRange(new Thingy(2, 1)).empty());
  assert(btree.equalRange(new Thingy(2, 2)).empty());
}
