// Copyright 2016 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

// Original Author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.sequence_lexicon_test;

import s2.sequence_lexicon;
import fluent.asserts;

void expectSequence(T)(in T[] actual, in T[] expected, ) {
  Assert.equal(actual.length, expected.length);
  Assert.equal(actual, expected);
}

alias Seq = long[];

@("SequenceLexicon.long") unittest {
  auto lex = new SequenceLexicon!long();
  Assert.equal(lex.add(new long[0]), 0);
  Assert.equal(lex.add([5L]), 1);
  Assert.equal(lex.add(new long[0]), 0);
  Assert.equal(lex.add([5L, 5L]), 2);
  Assert.equal(lex.add([5L, 0L, -3L]), 3);
  Assert.equal(lex.add([5L]), 1);
  Assert.equal(lex.add([0x7fffffffffffffffL]), 4);
  Assert.equal(lex.add([5L, 0L, -3L]), 3);
  Assert.equal(lex.add(new long[0]), 0);
  Assert.equal(lex.size(), 5);
  expectSequence(lex.sequence(0), []);
  expectSequence(lex.sequence(1), [5L]);
  expectSequence(lex.sequence(2), [5L, 5L]);
  expectSequence(lex.sequence(3), [5L, 0L, -3L]);
  expectSequence(lex.sequence(4), [0x7fffffffffffffffL]);
}

@("SequenceLexicon.Clear") unittest {
  auto lex = new SequenceLexicon!long();
  Assert.equal(lex.add([1L]), 0);
  Assert.equal(lex.add([2L]), 1);
  lex.clear();
  Assert.equal(lex.add([2L]), 0);
  Assert.equal(lex.add([1L]), 1);
}

@("SequenceLexicon.CopyConstructor") unittest {
  auto original = new SequenceLexicon!long();
  Assert.equal(original.add([1L, 2L]), 0);
  auto lex = new SequenceLexicon!long(original);
  Assert.equal(lex.add([3L, 4L]), 1);
  expectSequence(lex.sequence(0), [1L, 2L]);
  expectSequence(lex.sequence(1), [3L, 4L]);
}
