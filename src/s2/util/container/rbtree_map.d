module s2.util.container.rbtree_map;

import std.container.rbtree;
import std.algorithm : map;
import std.functional : binaryFun;
import std.meta : allSatisfy;
import std.range : ElementType, isInputRange;
import std.traits : isDynamicArray, isImplicitlyConvertible;

/**
 * A dictionary or associative array backed by a Red-Black tree.
 */
final class RBTreeMap(KeyT, ValueT, alias KeyLessF = "a < b", bool allowDuplicates = false) {
public:
  static struct Pair {
    KeyT key;
    ValueT value;
  }

  alias keyLess = binaryFun!KeyLessF;

  alias RedBlackTreeT =
      RedBlackTree!(Pair, (pair1, pair2) => keyLess(pair1.key, pair2.key), allowDuplicates);

  RedBlackTreeT rbTree;

  // Forward compatible methods like: empty(), length(), opSlice(), etc.
  alias rbTree this;

  this() {
    rbTree = new RedBlackTreeT();
  }

  this(Pair[] elems...) {
    rbTree = new RedBlackTreeT(elems);
  }

  this(PairRange)(PairRange pairRange)
  if (isInputRange!PairRange && isImplicitlyConvertible!(ElementType!PairRange, Pair)) {
    rbTree = new RedBlackTreeT(pairRange);
  }

  override
  bool opEquals(Object rhs) {
    RBTreeMap that = cast(RBTreeMap) rhs;
    if (that is null) return false;

    return rbTree == that.rbTree;
  }

  /// Insertion
  size_t stableInsert(K, V)(K key, V value)
  if (isImplicitlyConvertible!(K, KeyT) && isImplicitlyConvertible!(V, ValueT)) {
    return rbTree.stableInsert(Pair(key, value));
  }
  alias insert = stableInsert;

  ValueT opIndexAssign(ValueT value, KeyT key) {
    rbTree.stableInsert(Pair(key, value));
    return value;
  }

  /// Membership
  bool opBinaryRight(string op)(KeyT key) const
  if (op == "in") {
    return Pair(key) in rbTree;
  }

  /// Removal
  size_t removeKey(K...)(K keys)
  if (allSatisfy!(isImplicitlyConvertibleToKey, K)) {
    KeyT[K.length] toRemove = [keys];
    return removeKey(toRemove[]);
  }

  //Helper for removeKey.
  private template isImplicitlyConvertibleToKey(K)
  {
    enum isImplicitlyConvertibleToKey = isImplicitlyConvertible!(K, KeyT);
  }

  size_t removeKey(K)(K[] keys)
  if (isImplicitlyConvertible!(K, KeyT)) {
    auto keyPairs = keys.map!(key => Pair(key));
    return rbTree.removeKey(keyPairs);
  }

  size_t removeKey(KeyRange)(KeyRange keyRange)
  if (isInputRange!KeyRange
      && isImplicitlyConvertible!(ElementType!KeyRange, KeyT)
      && !isDynamicArray!KeyRange) {
    auto keyPairs = keys.map(key => Pair(key));
    return rbTree.removeKey(keyPairs);
  }

  /// Ranges
  RedBlackTreeT.Range upperBound(KeyT key) {
    return rbTree.upperBound(Pair(key));
  }

  RedBlackTreeT.ConstRange upperBound(KeyT key) const {
    return rbTree.upperBound(Pair(key));
  }

  RedBlackTreeT.ImmutableRange upperBound(KeyT key) immutable {
    return rbTree.upperBound(Pair(key));
  }

  RedBlackTreeT.Range lowerBound(KeyT key) {
    return rbTree.lowerBound(Pair(key));
  }

  RedBlackTreeT.ConstRange lowerBound(KeyT key) const {
    return rbTree.lowerBound(Pair(key));
  }

  RedBlackTreeT.ImmutableRange lowerBound(KeyT key) immutable {
    return rbTree.lowerBound(Pair(key));
  }

  auto equalRange(KeyT key) {
    return rbTree.equalRange(Pair(key));
  }

}
