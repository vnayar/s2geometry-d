module s2.util.container.dense_hash_table_test;

import s2.util.container.dense_hash_table;
import fluent.asserts;

auto makeBasicDenseHashTable(T)() {
  return denseHashTable!(T, T)(
      0,
      (in T key) => cast(size_t) key,
      (in T key) => key,
      (ref T val, T key) => val = key,
      (in T k1, in T k2) => k1 == k2);
}

@("basic-insert") unittest {
  auto ht = makeBasicDenseHashTable!int();
  ht.setEmptyKey(-1);
  ht.insert(3);
  ht.insert(5);
  ht.insert(5); // Duplicate
  ht.insert(8);
  Assert.equal(ht.size(), 3);
}

@("basic-find") unittest {
  auto ht = makeBasicDenseHashTable!int();
  ht.setEmptyKey(-1);
  ht.insert(3);
  ht.insert(5);
  ht.insert(8);
  auto it = ht.find(5);  // Iterator is returned to find an element.
  Assert.equal(*it, 5);  // Iterator's value matches the search target.
  it = ht.find(10);      // Iterators are re-assignable with no error.
  Assert.equal(it, ht.end());  // Find for non-existent elements returns end().
}

@("basic-erase") unittest {
  auto ht = makeBasicDenseHashTable!int();
  ht.setEmptyKey(-1);
  ht.setDeletedKey(-2);
  ht.insert(3);
  ht.insert(5);
  ht.insert(8);

  Assert.equal(ht.erase(5), 1);  // Erase returns the number of items removed.
  Assert.equal(ht.erase(5), 0);  // Test duplicate erase.
  Assert.equal(ht.erase(7), 0);  // Test erase of non-existent element.
  Assert.equal(ht.size(), 2);    // Verify the new size of the DenseHashTable.
}

struct Person {
  string name;
  int age;
}

auto makePersonDenseHashTable() {
  return denseHashTable!(Person, string)(
      0,
      (string name) => typeid(name).getHash(&name),
      (Person p) => p.name,
      (ref Person p, string key) => p.name = key,
      (string k1, string k2) => k1 == k2);
}

@("struct-insert") unittest {
  auto ht = makePersonDenseHashTable();
  ht.setEmptyKey(Person("nobody", 0));
  ht.insert(Person("Jim", 32));
  ht.insert(Person("Ruffy", 5));
  ht.insert(Person("Ruffy", 5));  // Duplicate
  ht.insert(Person("Snowball", 7));
  Assert.equal(ht.size(), 3);
}

@("struct-find") unittest {
  auto ht = makePersonDenseHashTable();
  ht.setEmptyKey(Person("nobody", 0));
  ht.insert(Person("Jim", 32));
  ht.insert(Person("Ruffy", 5));
  ht.insert(Person("Snowball", 7));

  auto it = ht.find("Jim");
  Assert.equal(*it, Person("Jim", 32));
  it = ht.find("Bob");
  Assert.equal(it, ht.end());
}

@("struct-erase") unittest {
  auto ht = makePersonDenseHashTable();
  ht.setEmptyKey(Person("nobody", 0));
  ht.setDeletedKey("johndoe");
  ht.insert(Person("Jim", 32));
  ht.insert(Person("Ruffy", 5));
  ht.insert(Person("Snowball", 7));

  Assert.equal(ht.erase("Ruffy"), 1);  // Erase returns the number of items removed.
  Assert.equal(ht.erase("Ruffy"), 0);  // Test duplicate erase.
  Assert.equal(ht.erase("Fishy"), 0);  // Test erase of non-existent element.
  Assert.equal(ht.size(), 2);    // Verify the new size of the DenseHashTable.
}

class Fruit {
  string name;
  int quantity;

  this(string n, int q) {
    name = n;
    quantity = q;
  }

  override
  bool opEquals(Object o) {
    Fruit f = cast(Fruit) o;
    return f !is null && name == f.name && quantity == f.quantity;
  }
}

auto makeFruitDenseHashTable() {
  return denseHashTable!(Fruit, string)(
      0,
      (string name) => typeid(name).getHash(&name),
      (in Fruit f) => f.name,
      (ref Fruit p, string key) => p.name = key,
      (string k1, string k2) => k1 == k2);
}

@("class-insert") unittest {
  auto ht = makeFruitDenseHashTable();
  ht.setEmptyKey(new Fruit("nothing", 0));
  ht.insert(new Fruit("Apple", 32));
  ht.insert(new Fruit("Pear", 5));
  ht.insert(new Fruit("Pear", 5));  // Duplicate
  ht.insert(new Fruit("Orange", 7));
  Assert.equal(ht.size(), 3);
}

@("class-find") unittest {
  auto ht = makeFruitDenseHashTable();
  ht.setEmptyKey(new Fruit("nothing", 0));
  ht.insert(new Fruit("Apple", 32));
  ht.insert(new Fruit("Pear", 5));
  ht.insert(new Fruit("Orange", 7));

  auto it = ht.find("Apple");
  Assert.equal((*it).name, "Apple");
  Assert.equal((*it).quantity, 32);
  it = ht.find("Lychee");
  Assert.equal(it, ht.end());
}

@("class-erase") unittest {
  auto ht = makeFruitDenseHashTable();
  ht.setEmptyKey(new Fruit("nothing", 0));
  ht.setDeletedKey("spoiled");
  ht.insert(new Fruit("Apple", 32));
  ht.insert(new Fruit("Pear", 5));
  ht.insert(new Fruit("Orange", 7));

  Assert.equal(ht.erase("Pear"), 1);  // Erase returns the number of items removed.
  Assert.equal(ht.erase("Pear"), 0);  // Test duplicate erase.
  Assert.equal(ht.erase("Lychee"), 0);  // Test erase of non-existent element.
  Assert.equal(ht.size(), 2);    // Verify the new size of the DenseHashTable.
}
