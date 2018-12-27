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

module s2.sequence_lexicon;

import s2.util.container.dense_hash_set;

// SequenceLexicon is a class for compactly representing sequences of values
// (e.g., tuples).  It automatically eliminates duplicates, and maps the
// remaining sequences to sequentially increasing integer ids.  See also
// ValueLexicon and IdSetLexicon.
//
// Each distinct sequence is mapped to a 32-bit integer.  The space used for
// each sequence is approximately 11 bytes plus the memory needed to represent
// the sequence elements.  For example, a sequence of three "double"s would
// need about 11 + 3*8 = 35 bytes.  Note also that sequences are referred to
// using 32-bit ids rather than 64-bit pointers.
//
// This class has the same thread-safety properties as "string": const methods
// are thread safe, and non-const methods are not thread safe.
//
// Example usage:
//
//   SequenceLexicon<string> lexicon;
//   vector<string> pets {"cat", "dog", "parrot"};
//   uint32 pets_id = lexicon.Add(pets);
//   CHECK_EQ(pets_id, lexicon.Add(pets));
//   string values;
//   for (const auto& pet : lexicon.sequence(pets_id)) {
//     values += pet;
//   }
//   CHECK_EQ("catdogparrot", values);
//
class SequenceLexicon(T) {
public:
  this() {
    _idSet = new IdSet();
    _idSet.setEmptyKey(EMPTY_KEY);
    _begins ~= 0;
  }

  this(in SequenceLexicon x) {
    _values = x._values;
    _begins = x._begins;
    _idSet = new IdSet(x._idSet.begin(), x._idSet.end(), EMPTY_KEY, 0);
  }

  // Clears all data from the lexicon.
  void clear() {
    _values.clear();
    _begins.clear();
    _idSet.clear();
    _begins ~= 0;
  }

  // Add the given sequence of values to the lexicon if it is not already
  // present, and return its integer id.  Ids are assigned sequentially
  // starting from zero.  "begin" and "end" are forward iterators over a
  // sequence of values of type T.
  //uint add(FwdIterator)(FwdIterator begin, FwdIterator end) {
  uint add(Range)(Range r)
  if (isInputRange!Range && is(r.front : uint)) {
    foreach (elem; r) {
      _values ~= elem;
    }
    _begins ~= _values.length;
    uint id = _begins.length - 2;
    auto result = _idSet.insert(id);
    if (result.second) {
      return id;
    } else {
      _begins.popBack();
      _values.length = _begins.back();
      return *result.first;
    }
  }

  // Return the number of value sequences in the lexicon.
  uint size() const {
    return _begins.length - 1;
  }

  // A representation of a sequence of values.
  alias Sequence = T[];

  // Return the value sequence with the given id.  This method can be used
  // with range-based for loops as follows:
  //   for (const auto& value : lexicon.sequence(id)) { ... }
  Sequence sequence(uint id) const {
    return _values[_begins[id] .. _begins[id + 1]];
  }

private:
  // Choose kEmptyKey to be the last key that will ever be generated.
  enum uint EMPTY_KEY = uint.max;

  alias IdSet = DenseHashSet!uint;

  T[] _values;
  uint[] _begins;
  IdSet _idSet;
}
