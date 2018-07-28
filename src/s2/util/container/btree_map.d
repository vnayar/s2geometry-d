module s2.util.container.btree_map;

import s2.util.container.btree;

import std.algorithm : map;
import std.functional : binaryFun;
import std.meta : allSatisfy;
import std.range : ElementType, isInputRange;
import std.traits : isDynamicArray, isImplicitlyConvertible;

/**
 * An associative-array or map implementation based upon a B-Tree.
 */
final class BTreeMap(KeyT, ValueT, size_t NodeSizeV = 256, alias KeyLessF = "a < b") {
public:
  static struct Pair {
    KeyT key;
    ValueT value;
  }

  alias _isKeyLess = binaryFun!KeyLessF;

  alias BTreeT = BTree!(Pair, NodeSizeV, (pair1, pair2) => _isKeyLess(pair1.key, pair2.key));

  BTreeT bTree;

  // Forward compatible methods like: empty(), length(), opSlice(), etc.
  alias bTree this;

  this() {
    bTree = new BTreeT();
  }

  /// Insertion
  void insert(K, V)(K key, V value)
  if (isImplicitlyConvertible!(K, KeyT) && isImplicitlyConvertible!(V, ValueT)) {
    return bTree.insert(Pair(key, value));
  }

  ValueT opIndexAssign(ValueT value, KeyT key) {
    bTree.insert(Pair(key, value));
    return value;
  }

  /// Membership
  bool opBinaryRight(string op)(KeyT key) const
  if (op == "in") {
    return !bTree.equalRange(Pair(key)).empty();
  }

  /// Removal
  void remove(KeyT key) {
    bTree.remove(Pair(key));
  }

  /// Ranges
  BTreeT.Range upperRange(KeyT key) {
    return bTree.upperRange(Pair(key));
  }

  BTreeT.Range lowerRange(KeyT key) {
    return bTree.lowerRange(Pair(key));
  }

  BTreeT.Range equalRange(KeyT key) {
    return bTree.equalRange(Pair(key));
  }
}
