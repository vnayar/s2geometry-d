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

module s2.value_lexicon_test;

import s2.value_lexicon;

import s2.s1angle;
import s2.s2point;

import core.stdc.string;
import fluent.asserts;

@("ValueLexicon.DuplicateValues") unittest {
  auto lex = new ValueLexicon!long();
  Assert.equal(lex.add(5), 0);
  Assert.equal(lex.add(0), 1);
  Assert.equal(lex.add(0), 1);
  Assert.equal(lex.add(-3), 2);
  Assert.equal(lex.add(5), 0);
  Assert.equal(lex.add(0), 1);
  Assert.equal(lex.add(0x7fffffffffffffff), 3);
  Assert.equal(lex.add(-0x8000000000000000), 4);
  Assert.equal(lex.add(0x7fffffffffffffff), 3);
  Assert.equal(lex.add(-0x8000000000000000), 4);
  Assert.equal(lex.size(), 5);
  Assert.equal(lex.value(0), 5);
  Assert.equal(lex.value(1), 0);
  Assert.equal(lex.value(2), -3);
  Assert.equal(lex.value(3), 0x7fffffffffffffff);
  Assert.equal(lex.value(4), -0x8000000000000000);
}

@("ValueLexicon.Clear") unittest {
  auto lex = new ValueLexicon!long();
  Assert.equal(lex.add(1), 0);
  Assert.equal(lex.add(2), 1);
  Assert.equal(lex.add(1), 0);
  lex.clear();
  Assert.equal(lex.add(2), 0);
  Assert.equal(lex.add(1), 1);
  Assert.equal(lex.add(2), 0);
}

@("ValueLexicon.FloatEquality") unittest {
  auto lex = new ValueLexicon!(S2Point)();
  auto a = S2Point(1, 0.0, 0.0);
  auto b = S2Point(1, -0.0, 0.0);
  auto c = S2Point(1, 0.0, -0.0);
  Assert.notEqual(memcmp(&a, &b, a.sizeof), 0);
  Assert.notEqual(memcmp(&a, &c, a.sizeof), 0);
  Assert.notEqual(memcmp(&b, &c, a.sizeof), 0);
  Assert.equal(lex.add(a), 0);
  Assert.equal(lex.add(b), 0);
  Assert.equal(lex.add(c), 0);
  Assert.equal(lex.size(), 1);
  auto result = lex.value(0);
  Assert.equal(memcmp(&result, &a, a.sizeof), 0);
}

@("ValueLexicon.CopyConstructor") unittest {
  auto original = new ValueLexicon!long();
  Assert.equal(original.add(5), 0);
  auto lex = new ValueLexicon!long(original);
  original.clear();
  Assert.equal(lex.add(10), 1);
  Assert.equal(lex.value(0), 5);
  Assert.equal(lex.value(1), 10);
}
