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

import std.range;

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
    _idSet = new IdSet(0, new IdHashFunctor(), new IdKeyEqualFunctor());
    _idSet.setEmptyKey(EMPTY_KEY);
    _begins ~= 0;
  }

  this(SequenceLexicon x) {
    _values = x._values.dup;
    _begins = x._begins.dup;
    _idSet = new IdSet(
        x._idSet.begin(), x._idSet.end(), EMPTY_KEY,
        0, new IdHashFunctor(), new IdKeyEqualFunctor());
  }

  // Clears all data from the lexicon.
  void clear() {
    _values.length = 0;
    _begins.length = 0;
    _idSet.clear();
    _begins ~= 0;
  }

  // Add the given sequence of values to the lexicon if it is not already
  // present, and return its integer id.  Ids are assigned sequentially
  // starting from zero.  "begin" and "end" are forward iterators over a
  // sequence of values of type T.
  //uint add(FwdIterator)(FwdIterator begin, FwdIterator end) {
  uint add(Range)(Range r)
  if (isInputRange!Range && is(typeof(r.front) : T)) {
    // Add all new T values to _values.
    foreach (elem; r) {
      _values ~= elem;
    }
    // Add a new "begin" which is the end of the current sequence and start of the next.
    _begins ~= cast(uint) _values.length;
    // ids start from zero, and by now, there's the initial begin of 0 and the one just added.
    // Thus to start with 0, we subtract two from the length.
    uint id = cast(uint) _begins.length - 2;
    // Equality is redefined so be if the id maps to the same sequence, they are the same.
    // This also means that when no new values are added, the id maps to an empty sequence,
    // and all these are considered to be the same.
    auto result = _idSet.insert(id);
    if (result.second) {
      // If the insert worked, keep it.
      return id;
    } else {
      // If the insert did not work, take back the id and remove the new values.
      _begins.popBack();
      _values.length = _begins.back();
      return *result.first;
    }
  }

  // Return the number of value sequences in the lexicon.
  size_t size() const {
    return _begins.length - 1;
  }

  // A representation of a sequence of values.
  alias Sequence = T[];

  // Return the value sequence with the given id.  This method can be used
  // with range-based for loops as follows:
  //   for (const auto& value : lexicon.sequence(id)) { ... }
  const(Sequence) sequence(uint id) const {
    return _values[_begins[id] .. _begins[id + 1]];
  }

private:
  // Choose kEmptyKey to be the last key that will ever be generated.
  enum uint EMPTY_KEY = uint.max;

  class IdHashFunctor {
    size_t opCall(in uint id) const {
      import s2.util.hash.mix;

      HashMix mix;
      foreach (value; this.outer.sequence(id)) {
        mix.mix(typeid(uint).getHash(&value));
      }
      return mix.get();
    }
  }

  class IdKeyEqualFunctor {
    bool opCall(in uint id1, in uint id2) const {
      import std.algorithm : equal;

      if (id1 == id2) return true;
      if (id1 == EMPTY_KEY || id2 == EMPTY_KEY) {
        return false;
      }
      auto seq1 = this.outer.sequence(id1);
      auto seq2 = this.outer.sequence(id2);
      return seq1.length == seq2.length && equal(seq1, seq2);
    }
  }

  alias IdSet = DenseHashSet!(uint, IdHashFunctor, IdKeyEqualFunctor);

  T[] _values;
  uint[] _begins;
  IdSet _idSet;
}
