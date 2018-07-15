module s2.util.container.btree_test;

import fluent.asserts;
import s2.util.container.btree;

import std.stdio;

////
// BTree.Node tests
////

@("Node.findFirstGEIndex<16")
unittest {
  auto btree = new BTree!(char, 200);
  static assert(btree.MIN_DEGREE < 16);

  auto node = new btree.Node();
  node._values = "BBDDEFHH";
  node._numValues = 8;

  Assert.equal(node.findFirstGEIndex('A'), 0);
  Assert.equal(node.findFirstGEIndex('B'), 0);
  Assert.equal(node.findFirstGEIndex('C'), 2);
  Assert.equal(node.findFirstGEIndex('D'), 2);
  Assert.equal(node.findFirstGEIndex('E'), 4);
  Assert.equal(node.findFirstGEIndex('F'), 5);
  Assert.equal(node.findFirstGEIndex('G'), 6);
  Assert.equal(node.findFirstGEIndex('H'), 6);
  Assert.equal(node.findFirstGEIndex('I'), 8);
}

@("Node.findFirstGTIndex<16")
unittest {
  auto btree = new BTree!(char, 200);
  static assert(btree.MIN_DEGREE < 16);

  auto node = new btree.Node();
  node._values = "BBDDEFHH";
  node._numValues = 8;

  Assert.equal(node.findFirstGTIndex('A'), 0);
  Assert.equal(node.findFirstGTIndex('B'), 2);
  Assert.equal(node.findFirstGTIndex('C'), 2);
  Assert.equal(node.findFirstGTIndex('D'), 4);
  Assert.equal(node.findFirstGTIndex('E'), 5);
  Assert.equal(node.findFirstGTIndex('F'), 6);
  Assert.equal(node.findFirstGTIndex('G'), 6);
  Assert.equal(node.findFirstGTIndex('H'), 8);
  Assert.equal(node.findFirstGTIndex('I'), 8);
}

@("Node.findFirstGEIndex>=16")
unittest {
  import std.stdio;
  auto btree = new BTree!(char, 2000);
  assert(btree.MIN_DEGREE >= 16);

  auto node = new btree.Node();
  node._values = "BBDDEFHH";
  node._numValues = 8;

  Assert.equal(node.findFirstGEIndex('A'), 0);
  Assert.equal(node.findFirstGEIndex('B'), 0);
  Assert.equal(node.findFirstGEIndex('C'), 2);
  Assert.equal(node.findFirstGEIndex('D'), 2);
  Assert.equal(node.findFirstGEIndex('E'), 4);
  Assert.equal(node.findFirstGEIndex('F'), 5);
  Assert.equal(node.findFirstGEIndex('G'), 6);
  Assert.equal(node.findFirstGEIndex('H'), 6);
  Assert.equal(node.findFirstGEIndex('I'), 8);
}

@("Node.findFirstGTIndex>=16")
unittest {
  import std.stdio;
  auto btree = new BTree!(char, 2000);
  assert(btree.MIN_DEGREE >= 16);

  auto node = new btree.Node();
  node._values = "BBDDEFHH";
  node._numValues = 8;

  Assert.equal(node.findFirstGTIndex('A'), 0);
  Assert.equal(node.findFirstGTIndex('B'), 2);
  Assert.equal(node.findFirstGTIndex('C'), 2);
  Assert.equal(node.findFirstGTIndex('D'), 4);
  Assert.equal(node.findFirstGTIndex('E'), 5);
  Assert.equal(node.findFirstGTIndex('F'), 6);
  Assert.equal(node.findFirstGTIndex('G'), 6);
  Assert.equal(node.findFirstGTIndex('H'), 8);
  Assert.equal(node.findFirstGTIndex('I'), 8);
}

@("Node.splitChild")
unittest {
  alias BTreeType = BTree!(char, 1);
  static assert(BTreeType.MIN_DEGREE == 1);

  alias NodeType = BTreeType.Node;
  NodeType root = new NodeType();
  NodeType child = new NodeType();
  NodeType grandChild1 = new NodeType();
  NodeType grandChild2 = new NodeType();

  root._isLeaf = false;
  root._numValues = 0;
  root._children[0] = child;

  child._isLeaf = false;
  child._parent = root;
  child._values[0] = 'b';
  child._numValues = 1;
  child._children = [grandChild1, grandChild2];

  grandChild1._isLeaf = true;
  grandChild1._parent = child;
  grandChild1._numValues = 1;
  grandChild1._values[0] = 'a';

  grandChild1._isLeaf = true;
  grandChild1._parent = child;
  grandChild1._numValues = 1;
  grandChild1._values[0] = 'c';

  root.splitChild(0);

  Assert.equal(root._numValues, 1);
  Assert.equal(root._values[0], 'b');

  Assert.equal(root._children[0]._numValues, 0);
  Assert.equal(root._children[0]._parent, root);
  Assert.equal(root._children[0]._children[0], grandChild1);

  Assert.equal(root._children[1]._numValues, 0);
  Assert.equal(root._children[1]._parent, root);
  Assert.equal(root._children[1]._children[0], grandChild2);
}

@("Node.insertNonFull.leaf")
unittest {
  alias BTreeType = BTree!(char, 75);
  static assert(BTreeType.MIN_DEGREE == 2);

  BTreeType.Node node = new BTreeType.Node();
  node._isLeaf = true;
  node._numValues = 1;
  node._values = "b";

  node.insertNonFull('a');
  Assert.equal(node.getValues(), "ab");
  Assert.equal(node._numValues, 2);

  node.insertNonFull('c');
  Assert.equal(node.getValues(), "abc");
  Assert.equal(node._numValues, 3);
}

@("Node.insertNonFull.nonLeaf")
unittest {
  alias BTreeType = BTree!(char, 75);
  static assert(BTreeType.MIN_DEGREE == 2);

  BTreeType.Node node = new BTreeType.Node();
  node._isLeaf = false;

  // TODO: Resume here.
}

////
// BTree tests
////

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
  auto btree = new BTree!(char, 90)();

  static assert(btree.MIN_DEGREE == 3);

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

@("remove.1")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(0).getChild(1);
  node.remove('E');
  Assert.equal(node.getValues(), "DF");
}

@("remove.2a")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(0);
  node.remove('G');
  Assert.equal(node.getValues(), "CFJ");
  Assert.equal(node.getChild(1).getValues(), "DE");
}

@("remove.2b.2c")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(0);
  node.remove('C');
  Assert.equal(node.getValues(), "DGJ");
  Assert.equal(node.getChild(1).getValues(), "EF");

  node.remove('G');
  Assert.equal(node.getValues(), "DJ");
  Assert.equal(node.getChild(1).getValues(), "EFHI");
  Assert.equal(node.getChild(2).getValues(), "KLMN");
}

@("remove.3a")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(1);
  // Test having to grab an extra key from the right sibling.
  node.remove('U');
  Assert.equal(node.getValues(), "TX");
  Assert.equal(node.getChild(1).getValues(), "VW");
  Assert.equal(node.getChild(2).getValues(), "YZ");

  // Test having to grab an extra key from the left sibling.
  node.remove('V');
  Assert.equal(node.getValues(), "SX");
  Assert.equal(node.getChild(0).getValues(), "PQR");
  Assert.equal(node.getChild(1).getValues(), "TW");
}

@("remove.3b1")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(0);
  node.getChild(1).remove('D');
  // Node now looks like this:
  //   [  C      G      J   ]
  //    /     \     \      \
  // [ A B ] [ E F ] [H I] [K L M N]
  node.remove('F');
  Assert.equal(node.getValues(), "CJ");
  Assert.equal(node.getChild(1).getValues(), "EGHI");
  Assert.equal(node.getChild(2).getValues(), "KLMN");
}

@("remove.3b2")
unittest {
  auto btree = createTestBTree();
  btree.Node node = btree.root.getChild(0);
  node.getChild(1).remove('D');
  node.getChild(3).remove('K');
  node.getChild(3).remove('M');
  // Node now looks like this:
  //   [  C     G     J]
  //    /    \    \    \
  // [A B] [E F] [H I] [L N]
  node.remove('N');
  Assert.equal(node.getValues(), "CG");
  Assert.equal(node.getChild(2).getValues(), "HIJL");
}

@("remove.all")
unittest {
  auto btree = createTestBTree();
  foreach (c; 'A' .. cast(char)('Z' + 1)) {
    btree.remove(c);
  }
  Assert.equal(btree.root.isLeaf(), true);
  Assert.equal(btree.root.numValues(), 0);
}

@("Iterator.begin,end")
unittest {
  auto btree = createTestBTree();
  auto first = btree.begin();

  Assert.equal(btree.begin().getValue(), 'A');

  // Note that "end" is one past the last value.
  auto last = btree.end();
  last--;
  Assert.equal(last.getValue(), 'Z');
}

@("Iterator.increment")
unittest {
  auto btree = createTestBTree();
  char expected = 'A';
  for (auto iterator = btree.begin(); iterator != btree.end(); iterator++) {
    Assert.equal(iterator.getValue(), expected);
    expected++;
  }
}

@("Iterator.decrement")
unittest {
  auto btree = createTestBTree();
  char expected = 'Z';
  auto iterator = btree.end();
  do {
    iterator--;
    Assert.equal(iterator.getValue(), expected);
    expected--;
  } while (iterator != btree.begin());
}
