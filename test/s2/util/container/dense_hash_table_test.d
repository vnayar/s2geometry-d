module s2.util.container.dense_hash_table_test;

import s2.util.container.dense_hash_table;
import fluent.asserts;

@("insert-simple") unittest {
  auto ht = new DenseHashTable!(int, int, (int key) => cast(size_t) key, (int key) => key)();
}
