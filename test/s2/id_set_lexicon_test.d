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

module s2.id_set_lexicon_test;

import s2.id_set_lexicon;

import fluent.asserts;

void expectIdSet(in int[] actual, in int[] expected) {
  Assert.equal(actual.length, expected.length);
  Assert.equal(actual, expected);
}

//using IdSet = IdSetLexicon::IdSet;
//using Seq = std::vector<int32>;

@("IdSetLexicon.EmptySet") unittest {
  auto lexicon = new IdSetLexicon();
  expectIdSet(lexicon.idSet(lexicon.add(new int[0])), []);
}

@("IdSetLexicon.SingletonSets") unittest {
  auto lexicon = new IdSetLexicon();
  Assert.equal(lexicon.add([5]), 5);
  Assert.equal(lexicon.add([0]), 0);
  Assert.equal(lexicon.addSingleton(1), 1);
  int m = int.max;
  Assert.equal(lexicon.add([m]), m);

  expectIdSet(lexicon.idSet(0), [0]);
  expectIdSet(lexicon.idSet(1), [1]);
  expectIdSet(lexicon.idSet(5), [5]);
  expectIdSet(lexicon.idSet(m), [m]);
}

@("IdSetLexicon.SetsAreSorted") unittest {
  auto lexicon = new IdSetLexicon();
  Assert.equal(lexicon.add([2, 5]), ~0);
  Assert.equal(lexicon.add([3, 2, 5]), ~1);
  Assert.equal(lexicon.add([5, 2]), ~0);
  Assert.equal(lexicon.add([5, 3, 2, 5]), ~1);

  expectIdSet(lexicon.idSet(~0), [2, 5]);
  expectIdSet(lexicon.idSet(~1), [2, 3, 5]);
}

@("IdSetLexicon.Clear") unittest {
  auto lexicon = new IdSetLexicon();
  Assert.equal(lexicon.add([1, 2]), ~0);
  Assert.equal(lexicon.add([3, 4]), ~1);
  lexicon.clear();
  Assert.equal(lexicon.add([3, 4]), ~0);
  Assert.equal(lexicon.add([1, 2]), ~1);
}
