module s2.util.container.btree_map_test;

import fluent.asserts;
import s2.util.container.btree_map;

@("BTreeMap.insert")
unittest {
  auto btreeMap = new BTreeMap!(string, int)();
  btreeMap.insert("bob", 42);
  btreeMap.insert("sam", 56);
  Assert.equal(btreeMap.length, 2);
}

@("BTreeMap.opIndexAssign")
unittest {
  auto btreeMap = new BTreeMap!(string, int)();
  btreeMap["bob"] = 42;
  btreeMap["sam"] = 56;
  Assert.equal(btreeMap.length, 2);
}

@("BTreeMap.opIndexAssign")
unittest {
  auto btreeMap = new BTreeMap!(string, int)();
  btreeMap["bob"] = 42;
  btreeMap["sam"] = 56;
  Assert.equal(btreeMap.length, 2);
}

auto createTestBTreeMap() {
  auto btreeMap = new BTreeMap!(string, int)();
  btreeMap["Alpha"] = 14;
  btreeMap["Bravo"] = 67;
  btreeMap["Charlie"] = 34;
  btreeMap["Delta"] = 89;
  btreeMap["Echo"] = 23;
  btreeMap["Foxtrot"] = 76;
  btreeMap["Golf"] = 31;
  btreeMap["Hotel"] = 76;
  Assert.equal(btreeMap.length, 8);
  return btreeMap;
}

@("BTreeMap.remove")
unittest {
  auto btreeMap = createTestBTreeMap();
  btreeMap.remove("Echo");
  Assert.equal(btreeMap.length, 7);
  btreeMap.remove("Alpha");
  Assert.equal(btreeMap.length, 6);
}

@("BTreeMap.in")
unittest {
  auto btreeMap = createTestBTreeMap();
  Assert.equal("Echo" in btreeMap, true);
  Assert.equal("Elephant" !in btreeMap, true);
}

@("BTreeMap.upperRange")
unittest {
  import std.array;
  auto btreeMap = createTestBTreeMap();
  Assert.equal(array(btreeMap.upperRange("Echo")).length, 3);
}

@("BTreeMap.equalRange")
unittest {
  import std.array;
  auto btreeMap = createTestBTreeMap();
  Assert.equal(array(btreeMap.equalRange("Echo")).length, 1);
}

@("BTreeMap.lowerRange")
unittest {
  import std.array;
  auto btreeMap = createTestBTreeMap();
  Assert.equal(array(btreeMap.lowerRange("Echo")).length, 4);
}
