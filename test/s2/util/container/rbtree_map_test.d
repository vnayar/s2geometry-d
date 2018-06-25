module s2.util.container.rbtree_map_test;

import s2.util.container.rbtree_map;
import fluent.asserts;

@("insert")
unittest {
  auto rbTreeMap = new RBTreeMap!(string, int)();
  rbTreeMap.insert("a", 4);
  rbTreeMap.insert("b", 2);
  rbTreeMap.insert("c", 3);
  rbTreeMap.insert("d", 1);

  Assert.equal(rbTreeMap.length(), 4);

  auto rbTreeMap2 = new RBTreeMap!(string, int)();
  rbTreeMap2["a"] = 4;
  rbTreeMap2["b"] = 2;
  rbTreeMap2["c"] = 3;
  rbTreeMap2["d"] = 1;

  Assert.equal(rbTreeMap2.length(), 4);

  Assert.equal(rbTreeMap, rbTreeMap2);
}

@("membership")
unittest {
  alias RBTreeMapT = RBTreeMap!(string, int);
  auto rbTreeMap = new RBTreeMapT([
    RBTreeMapT.Pair("a", 4),
    RBTreeMapT.Pair("b", 2),
    RBTreeMapT.Pair("c", 3),
    RBTreeMapT.Pair("d", 1)
  ]);

  assert("b" in rbTreeMap);
  assert("d" in rbTreeMap);
  assert("e" !in rbTreeMap);
}

@("removal")
unittest {
  auto rbTreeMap = new RBTreeMap!(string, int)();
  rbTreeMap["a"] = 4;
  rbTreeMap["b"] = 2;
  rbTreeMap["c"] = 3;
  rbTreeMap["d"] = 1;

  rbTreeMap.removeKey("b");
  rbTreeMap.removeKey("c");

  Assert.equal(rbTreeMap.length(), 2);
}

@("lowerBound")
unittest {
  auto rbTreeMap = new RBTreeMap!(string, int)();
  rbTreeMap["a"] = 4;
  rbTreeMap["b"] = 2;
  rbTreeMap["c"] = 3;
  rbTreeMap["d"] = 1;

  auto pairRange = rbTreeMap.lowerBound("c");
  Assert.equal(pairRange.front(), rbTreeMap.Pair("a", 4));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), false);
  Assert.equal(pairRange.front(), rbTreeMap.Pair("b", 2));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), true);
}

@("equalRange")
unittest {
  auto rbTreeMap = new RBTreeMap!(string, int, "a < b", true)();
  rbTreeMap["a"] = 4;
  rbTreeMap["b"] = 2;
  rbTreeMap["c"] = 3;
  rbTreeMap["d"] = 1;
  rbTreeMap["b"] = 5;

  auto pairRange = rbTreeMap.equalRange("b");
  Assert.equal(pairRange.front(), rbTreeMap.Pair("b", 2));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), false);
  Assert.equal(pairRange.front(), rbTreeMap.Pair("b", 5));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), true);
}

@("upperBound")
unittest {
  auto rbTreeMap = new RBTreeMap!(string, int)();
  rbTreeMap["a"] = 4;
  rbTreeMap["b"] = 2;
  rbTreeMap["c"] = 3;
  rbTreeMap["d"] = 1;

  auto pairRange = rbTreeMap.upperBound("b");
  Assert.equal(pairRange.front(), rbTreeMap.Pair("c", 3));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), false);
  Assert.equal(pairRange.front(), rbTreeMap.Pair("d", 1));
  pairRange.popFront();
  Assert.equal(pairRange.empty(), true);
}
