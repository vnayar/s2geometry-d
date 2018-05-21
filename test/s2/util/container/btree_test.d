module s2.util.container.btree_test;

import fluent.asserts;
import s2.util.container.btree;

import std.stdio;


void printTree(T)(T btree) {
  btree.Node[] nodes = [btree.root];
  while (nodes.length > 0) {
    btree.Node[] newNodes = [];
    foreach (node; nodes) {
      foreach (i; 0 .. node.numValues()) {
        write(node.getValue(i));
      }
      write("   ");
      if (!node.isLeaf()) {
        foreach (i; 0 .. node.numChildren()) {
          newNodes ~= node.getChild(i);
        }
      }
    }
    nodes = newNodes[];
    writeln("");
  }
}


/**
 * Creates a BTRee representing a standard setup used for testing.
 * MIN_DEGREE should be 3, meaning min keys per node is 2 and max is 5.
 *
 * The Tree structure is as follows:
 *                    [    O     ]
 *                   /            \
 *      [  C     G     J ]         [     T      W   ]
 *       /   /      \    \           /        \     \
 *  [A B] [D E F] [H I] [K L M N]  [P Q R S] [U V] [X Y Z]
 */
auto createTestBTree() {
  auto btree = new BTree!(char, 62)();
  string testData = "YBZIOQPKTJUMLWHVAXGSNFCRDE";
  foreach (c; testData) {
    btree.insert(c);
  }

  Assert.equal(btree.root.getValue(0), 'O');

  Assert.equal(btree.root.getChild(0).getValue(0), 'C');
  Assert.equal(btree.root.getChild(0).getValue(1), 'G');
  Assert.equal(btree.root.getChild(0).getValue(2), 'J');

  Assert.equal(btree.root.getChild(0).getChild(0).getValue(0), 'A');
  Assert.equal(btree.root.getChild(0).getChild(0).getValue(1), 'B');

  Assert.equal(btree.root.getChild(0).getChild(1).getValue(0), 'D');
  Assert.equal(btree.root.getChild(0).getChild(1).getValue(1), 'E');
  Assert.equal(btree.root.getChild(0).getChild(1).getValue(2), 'F');

  Assert.equal(btree.root.getChild(0).getChild(2).getValue(0), 'H');
  Assert.equal(btree.root.getChild(0).getChild(2).getValue(1), 'I');

  Assert.equal(btree.root.getChild(0).getChild(3).getValue(0), 'K');
  Assert.equal(btree.root.getChild(0).getChild(3).getValue(1), 'L');
  Assert.equal(btree.root.getChild(0).getChild(3).getValue(2), 'M');
  Assert.equal(btree.root.getChild(0).getChild(3).getValue(3), 'N');

  Assert.equal(btree.root.getChild(1).getValue(0), 'T');
  Assert.equal(btree.root.getChild(1).getValue(1), 'W');

  Assert.equal(btree.root.getChild(1).getChild(0).getValue(0), 'P');
  Assert.equal(btree.root.getChild(1).getChild(0).getValue(1), 'Q');
  Assert.equal(btree.root.getChild(1).getChild(0).getValue(2), 'R');
  Assert.equal(btree.root.getChild(1).getChild(0).getValue(3), 'S');

  Assert.equal(btree.root.getChild(1).getChild(1).getValue(0), 'U');
  Assert.equal(btree.root.getChild(1).getChild(1).getValue(1), 'V');

  Assert.equal(btree.root.getChild(1).getChild(2).getValue(0), 'X');
  Assert.equal(btree.root.getChild(1).getChild(2).getValue(1), 'Y');
  Assert.equal(btree.root.getChild(1).getChild(2).getValue(2), 'Z');

  return btree;
}

unittest {
  auto btree = createTestBTree();
  assert(btree.MIN_DEGREE == 3);
}

