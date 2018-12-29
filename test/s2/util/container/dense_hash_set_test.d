module s2.util.container.dense_hash_set_test;

import s2.util.container.dense_hash_set;
import fluent.asserts;

@("int-set-insert") unittest {
  auto intSet = denseHashSet!int();

  intSet.setEmptyKey(-1);
  intSet.insert(4);
  intSet.insert(7);
  intSet.insert(5);
  Assert.equal(intSet.size(), 3);
}

@("char-set-insert") unittest {
  auto charSet = denseHashSet!char();

  charSet.setEmptyKey(cast(char) 255);
  charSet.insert('c');
  charSet.insert('o');
  charSet.insert('w');
  Assert.equal(charSet.size(), 3);
}

@("int-set-find") unittest {
  auto hs = denseHashSet!int();
  hs.setEmptyKey(-1);
  hs.insert(3);
  hs.insert(5);
  hs.insert(8);
  auto it = hs.find(5);  // Iterator is returned to find an element.
  Assert.equal(*it, 5);  // Iterator's value matches the search target.
  it = hs.find(10);      // Iterators are re-assignable with no error.
  Assert.equal(it, hs.end());  // Find for non-existent elements returns end().
}

@("int-set-erase") unittest {
  auto hs = denseHashSet!int();
  hs.setEmptyKey(-1);
  hs.setDeletedKey(-2);
  hs.insert(3);
  hs.insert(5);
  hs.insert(8);

  Assert.equal(hs.erase(5), 1);  // Erase returns the number of items removed.
  Assert.equal(hs.erase(5), 0);  // Test duplicate erase.
  Assert.equal(hs.erase(7), 0);  // Test erase of non-existent element.
  Assert.equal(hs.size(), 2);    // Verify the new size of the DenseHashTable.
}

