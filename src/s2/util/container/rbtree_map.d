module s2.util.container.rbtree_map;

import std.container.rbtree;
import std.algorithm : map;
import std.functional : binaryFun;
import std.traits : isImplicitlyConvertible;

/**
 * A dictionary or associative array backed by a Red-Black tree.
 */
final class RBTreeMap(KeyT, ValueT, alias KeyLessF = "a < b", bool allowDuplicates = false) {
private:
  static struct Pair {
    KeyT key;
    ValueT value;
  }
  alias keyLess = binaryFun!KeyLessF;

  alias RedBlackTreeT =
      RedBlackTree!(Pair, (pair1, pair2) => keyLess(pair1.key, pair2.key), allowDuplicates);

public:
  RedBlackTreeT rbtree;

  // Forward compatible methods like: empty(), length(), opSlice(), etc.
  alias rbtree this;

  this() {
    rbtree = new RedBlackTreeT();
  }

  /// Insertion
  size_t stableInsert(K, V)(K key, V value)
  if (isImplicitlyConvertible!(K, KeyT) && isImplicitlyConvertible!(V, ValueT)) {
    return rbtree.stableInsert(Pair(key, value));
  }
  alias insert = stableInsert;

  ValueT opIndexAssign(ValueT value, KeyT key) {
    rbtree.stableInsert(Pair(key, value));
    return value;
  }

  /// Membership
  bool opBinaryRight(string op)(KeyT key) const
  if (op == "in") {
    return Pair(key) in rbtree;
  }

  /// Removal
  size_t removeKey(K...)(K keys)
  if (allSatisfy!(isImplicitlyConvertibleToKey, K)) {
    KeyT[U.length] toRemove = [keys];
    return removeKey(toRemove[]);
  }

  //Helper for removeKey.
  private template isImplicitlyConvertibleToKey(K)
  {
    enum isImplicitlyConvertibleToKey = isImplicitlyConvertible!(K, KeyT);
  }

  size_t removeKey(K)(K[] keys)
  if (isImplicitlyConvertible!(K, KeyT)) {
    auto keyPairs = keys.map(key => Pair(key));
    return rbtree.removeKey(keyPairs);
  }

  size_t removeKey(KeyRange)(KeyRange keyRange)
  if (isInputRange!KeyRange
      && isImplicitlyConvertible!(ElementType!KeyRange, KeyT)
      && !isDynamicArray!KeyRange) {
    auto keyPairs = keys.map(key => Pair(key));
    return rbtree.removeKey(keyPairs);
  }

  /// Ranges
  auto /*RedBlackTreeT.Range*/ upperBound(KeyT key) {
    return rbtree.upperBound(Pair(key));
  }

  auto /*RedBlackTree.ConstRange*/ upperBound(KeyT key) const {
    return rbtree.upperBound(Pair(key));
  }

  auto /*RedBlackTree.ImmutableRange*/ upperBound(KeyT key) immutable {
    return rbtree.upperBound(Pair(key));
  }

  auto /*RedBlackTreeT.Range*/ lowerBound(KeyT key) {
    return rbtree.lowerBound(Pair(key));
  }

  auto /*RedBlackTreeT.ConstRange*/ lowerBound(KeyT key) const {
    return rbtree.lowerBound(Pair(key));
  }

  auto /*RedBlackTreeT.ImmutableRange*/ lowerBound(KeyT key) immutable {
    return rbtree.lowerBound(Pair(key));
  }

  auto equalRange(this This)(Key key) {
    return rbtree.equalRange(key);
  }

}
