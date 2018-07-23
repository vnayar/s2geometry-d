module s2.util.container.btree_map;

import s2.util.container.btree;

import std.functional : binaryFun;

/**
 * An associative-array or map implementation based upon a B-Tree.
 */
final class BTreeMap(KeyT, ValueT, alias KeyLessF = "a < b") {
public:
  static struct Pair {
    KeyT key;
    ValueT value;
  }

  alias _isKeyLess = binaryFun!KeyLessF;

  alias BTreeT = BTree!(Pair, (pair1, pair2) => _isKeyLess(pair1.key, pair2.key));

  BTreeT bTree;

  // Forward compatible methods like: empty(), length(), opSlice(), etc.
  alias bTree this;

  this() {
    bTree = new BTree();
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
  size_t remove(K...)(K keys)
  if (allSatisfy!(isImplicitlyConvertibleToKey, K)) {
    KeyT[K.length] toRemove = [keys];
    return remove(toRemove[]);
  }

  //Helper for removeKey.
  private template isImplicitlyConvertibleToKey(K)
  {
    enum isImplicitlyConvertibleToKey = isImplicitlyConvertible!(K, KeyT);
  }

  size_t remove(K)(K[] keys)
  if (isImplicitlyConvertible!(K, KeyT)) {
    auto keyPairs = keys.map!(key => Pair(key));
    return bTree.remove(keyPairs);
  }

  size_t remove(KeyRange)(KeyRange keyRange)
  if (isInputRange!KeyRange
      && isImplicitlyConvertible!(ElementType!KeyRange, KeyT)
      && !isDynamicArray!KeyRange) {
    auto keyPairs = keys.map(key => Pair(key));
    return bTree.remove(keyPairs);
  }

  /// Ranges
  BTreeT.Range upperRange(KeyT key) {
    return bTree.upperBound(Pair(key));
  }

  BTreeT.Range lowerRange(KeyT key) {
    return bTree.lowerRange(Pair(key));
  }

  BTreeT.Range equalRange(KeyT key) {
    return bTree.equalRange(Pair(key));
  }
}
