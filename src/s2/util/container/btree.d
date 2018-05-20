module s2.util.container.btree;

import std.algorithm : max;
import std.functional : unaryFun, binaryFun;
import std.traits : ReturnType;

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
 * This implementation chooses to keep the data inside the B-Tree itself rather than only storing
 * references to the data. However, if a solution that only stores a struct id is desired,
 * simply use the BTree with the object id (and an appropriate KeyF if the key is not the object
 * id).
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
 *   KeyT = The key used to compare elements. Defaults to T.
 *   KeyF = An element type to key conversion function, either a function or a string
 *     representing a function as described by
 *     $(LINK2 https://dlang.org/phobos/std_functional.html#unaryFun, std.functional.unaryFun).
 *     When working with class types, this method must be $(D_INLINECODE const) or
 *     $(D_INLINECODE inout).
 *   KeyLessF = A less-than comparison function for keys, either a function or a string
 *     representing a function as described by
 *     $(LINK2 https://dlang.org/phobos/std_functional.html#binaryFun, std.functional.binaryFun).
 *   KeyEqualF = An equality function.
 */
final class BTree(
    ValueT, size_t NodeSizeV = 256, KeyT = ValueT,
    alias KeyF = "a", alias KeyLessF = "a < b", alias KeyEqualF = "a == b")
if (is(typeof(binaryFun!KeyLessF(KeyT.init, KeyT.init)) : bool)
    && is(typeof(binaryFun!KeyEqualF(KeyT.init, KeyT.init)) : bool)
    && is(typeof(unaryFun!KeyF(ValueT.init)) : KeyT)) {

private:
  Node _root;

  // TODO: Determine how to enforce const functions are passed in, which is
  // currently preventing the use of const with reference types.
  alias _getKey = unaryFun!KeyF;
  alias _keyLess = binaryFun!KeyLessF;
  alias _keyEqual = binaryFun!KeyEqualF;

  // The B-Tree is defined by a minimum degree "t".
  // Each node other than the root must have at least t-1 keys.
  // Each node may have up to 2t - 1 keys, and thus an internal node may have up to 2*t children.
  //
  // Thus, the maximum size of a non-leaf node is:
  //   [values]    (2t-1) * (size of a value)
  //   [numValues] + 1 * (size of size_t)
  //   [isLeaf]    + 1 * (size of a bool)
  //   [children]  + 2*t * (size of pointer to a node)
  //   = NodeSize
  static size_t getMinDegree() @nogc @safe pure nothrow {
    return (NodeSizeV - bool.sizeof - size_t.sizeof + ValueT.sizeof)
        / (2 * (ValueT.sizeof + (Node*).sizeof));
  }

public:

  enum MIN_DEGREE = getMinDegree();
  enum MAX_DEGREE = MIN_DEGREE * 2;

  struct Result {
  private:
    Node* _node;
    size_t _position;
  public:

    bool isFound() {
      return _node !is null;
    }

    inout(KeyT) getKey() inout {
      return _node.getKey(_position);
    }

    inout(ValueT) getValue() inout {
      return _node.getValue(_position);
    }
  }

  struct Node {
  private:
    // Dynamic arrays have a length field of type size_t.
    // The key is extracted from the value using KeyF.
    ValueT[MAX_DEGREE - 1] _values;
    size_t _numValues;
    bool _isLeaf;
    // Only non-leaf (internal) nodes have children.
    Node*[MAX_DEGREE] _children;

    void splitChild(size_t i)
    in {
      assert(!isFull());
      assert(_children[i].isFull());
    } body {
      Node* toSplitNode = _children[i];

      // Prepare a new node containig the right half of the node at _children[i].
      // TODO: ALLOCATE-NODE(newNode)
      Node newNode = Node();
      newNode._isLeaf = toSplitNode._isLeaf;
      newNode._numValues = MIN_DEGREE - 1;
      foreach (j; 0 .. MIN_DEGREE - 1) {
        newNode._values[j] = toSplitNode._values[j + MIN_DEGREE];
      }
      if (!toSplitNode._isLeaf) {
        foreach (j; 1 .. MIN_DEGREE) {
          newNode._children[j] = toSplitNode._children[j + MIN_DEGREE];
        }
      }

      // Make way for a key value to be added from _children[i].
      for (auto j = _numValues; j >= i; j--) {
        _children[j + 1] = _children[j];
      }
      _children[i] = &newNode;

      for (auto j = _numValues - 1; j >= i; j--) {
        _values[j+1] = _values[j];
      }
      _values[i] = toSplitNode._values[MIN_DEGREE];
      _numValues++;

      // Finally reduce the size of key/values in the node being split, cutting it in half.
      toSplitNode._numValues = MIN_DEGREE - 1;

      // TODO: DISK-WRITE(toSplitNode)
      // TODO: DISK-WRITE(newNode)
      // TODO: DISK-WRITE(this)
    }

    void insertNonFull(ValueT v) {
      int i = cast(int) _numValues - 1;
      KeyT k = _getKey(v);
      if (_isLeaf) {
        for (; i >= 0 && _keyLess(k, getKey(i)); i--) {
          _values[i+1] = _values[i];
        }
        _values[i+1] = v;
        _numValues++;
        // TODO: DISK-WRITE(this)
      } else {
        // Find the index of the first key less than the inserted value.
        for (; i >= 0 && _keyLess(k, getKey(i)); i--) {}
        // The child to insert into is one past the key's index, which may be -1.
        //     k[0]    k[1]    k[2]
        // c[0]    c[1]    c[2]    c[3]
        i++;
        // TODO: DISK-READ(_children[i])
        if (_children[i].isFull()) {
          splitChild(i);
          // After splitting, the median value of the child is not in this node at index i.
          if (k > getKey(i)) {
            i++;
          }
        }
        _children[i].insertNonFull(v);
      }
    }

  public:
    inout(KeyT) getKey(size_t i) inout
    in {
      assert(i >= 0);
      assert(i < MAX_DEGREE - 1);
    } body {
      return _getKey(getValue(i));
    }

    size_t numKeys() const {
      return _numValues;
    }

    bool isFull() {
      return _numValues >= MAX_DEGREE - 1;
    }

    inout(ValueT) getValue(size_t i) inout
    in {
      assert(i >= 0);
      assert(i < MAX_DEGREE);
    } body {
      return _values[i];
    }

    size_t numValues() const {
      return _numValues;
    }

    inout(Result) search(KeyT k) inout {
      size_t i = 0;
      while (i < numKeys() && _keyLess(getKey(i), k)) {
        i++;
      }
      if (i < numKeys() && _keyEqual(getKey(i), k)) {
        return inout(Result)(&this, i);
      }
      if (_isLeaf) {
        return inout(Result)(null, -1);
      } else {
        return _children[i].search(k);
      }
    }
  }

  this() {
    // TODO: ALLOCATE-NODE() which allocates a disk page in O(1).
    _root = Node();
    _root._isLeaf = true;
    _root._numValues = 0;
    // TODO: DISK-WRITE(_root)
  }

  @property
  inout(Node) root() inout {
    return _root;
  }

  inout(Result) search(KeyT k) inout {
    return _root.search(k);
  }

  void insert(ValueT v) {
    Node* r = &_root;
    if (r.isFull()) {
      Node newRoot = Node();
      _root = newRoot;
      newRoot._isLeaf = false;
      newRoot._numValues = 0;
      newRoot._children[0] = r;
      newRoot.splitChild(0);
      newRoot.insertNonFull(v);
    } else {
      r.insertNonFull(v);
    }
  }
}

unittest {
  assert(BTree!(int, 256).MIN_DEGREE == 10);
  assert(BTree!(int, 256).MAX_DEGREE == 20);

  assert(BTree!(int, 4096).MIN_DEGREE == 170);
  assert(BTree!(int, 4096).MAX_DEGREE == 340);

  // Structs are passed by value, so each value needs 6*4 = 24 bytes.
  struct S {
    int a, b, c, d, e, f;
  }
  assert(BTree!(S, 1024, int, "a.a").MIN_DEGREE == 16);
  assert(BTree!(S, 1024, int, "a.a").MAX_DEGREE == 32);

  // Classes are passed by reference, and thus the value needs only 8 bytes.
  class C {
    int a, b, c, d, e, f;
    int getVal() const { return d; }
  }
  assert(BTree!(C, 1024, int, "a.getVal()").MIN_DEGREE == 31);
  assert(BTree!(C, 1024, int, "a.getVal()").MAX_DEGREE == 62);
}

/// Simple use case with primitive types.
unittest {
  auto btree = new BTree!int();
  btree.insert(2);
  btree.insert(3);
  btree.insert(4);
  btree.insert(1);

  auto r = btree.search(2);
  assert(r.isFound());
  assert(r.getValue() == 2);

  r = btree.search(-3);
  assert(!r.isFound());
}

/// Use case using comparison by a specific value.
unittest {
  struct Structo {
    int _id;
    string _data;
  }

  // Organize the BTree using the _id field as the key used for comparison.
  auto btree = new BTree!(Structo, 1024, int, "a._id");

  // This time use the string _data field, but only the first two characters.
  auto btree2 = new BTree!(
      Structo, 1024, string, "a._data", "a[0..2] > b[0..2]", "a[0..2] == b[0..2]");
}

/// Use case compatible with class comparing operator overrides:
unittest {
  class Thingy {
    int _a, _b;

    this(int a, int b) {
      _a = a;
      _b = b;
    }

    int opCmp(Thingy o) const {
      return _a < o._a ? -1
          : _a > o._a ? 1
          : _b < o._b ? -1
          : _b > o._b ? 1
          : 0;
    }

    override
    bool opEquals(Object o) const {
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

  assert(btree.search(new Thingy(2, 1)).isFound());
  assert(!btree.search(new Thingy(2, 2)).isFound());
}

