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
}
