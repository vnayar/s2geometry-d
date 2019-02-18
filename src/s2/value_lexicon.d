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

module s2.value_lexicon;

import s2.util.container.dense_hash_set;

import std.range;

/**
 * ValueLexicon is a class that maps distinct values to sequentially numbered
 * integer identifiers.  It automatically eliminates duplicates and uses a
 * compact representation.  See also SequenceLexicon.
 *
 * Each distinct value is mapped to a 32-bit integer.  The space used for each
 * value is approximately 7 bytes plus the space needed for the value itself.
 * For example, int64 values would need approximately 15 bytes each.  Note
 * also that values are referred to using 32-bit ids rather than 64-bit
 * pointers.
 *
 * This class has the same thread-safety properties as "string": const methods
 * are thread safe, and non-const methods are not thread safe.
 *
 * Example usage:
 *
 *   ValueLexicon<string> lexicon;
 *   uint32 cat_id = lexicon.Add("cat");
 *   EXPECT_EQ(cat_id, lexicon.Add("cat"));
 *   EXPECT_EQ("cat", lexicon.value(cat_id));
 *
 */
class ValueLexicon(T) {
public:
  this() {
    // Simply use the defaults used in the DenseHashSet.
    _hasher = &hash!T;
    _keyEqual = &equalTo!T;
    _idSet = new IdSet(0, new IdHasher(), new IdKeyEqual());
    _idSet.setEmptyKey(EMPTY_KEY);
  }

  this(ValueLexicon!T x) {
    _hasher = &hash!T;
    _keyEqual = &equalTo!T;
    _values = x._values.dup;
    _idSet = new IdSet(
        x._idSet.begin(), x._idSet.end(), EMPTY_KEY, 0, new IdHasher(), new IdKeyEqual());
  }

  /// Clears all data from the lexicon.
  void clear() {
    _values.length = 0;
    _idSet.clear();
  }

  /// Add the given value to the lexicon if it is not already present, and
  /// return its integer id.  Ids are assigned sequentially starting from zero.
  uint add(in T value) {
    if (!_values.empty() && _keyEqual(value, _values.back())) {
      return cast(uint) _values.length - 1;
    }
    _values ~= value;
    auto result = _idSet.insert(cast(uint) _values.length - 1);
    if (result.second) {
      return cast(uint) _values.length - 1;
    } else {
      _values.popBack();
      return *(result.first);
    }
  }

  /// Return the number of values in the lexicon.
  uint size() const {
    return cast(uint) _values.length;
  }

  /// Return the value with the given id.
  const(T) value(uint id) const {
    return _values[id];
  }

 private:
  // Choose kEmptyKey to be the last key that will ever be generated.
  enum int EMPTY_KEY = int.max;

  class IdHasher {
  public:
    this() { }
    size_t opCall()(uint id) const {
      return _hasher(value(id));
    }
  }

  class IdKeyEqual {
  public:
    this() { }
    bool opCall(in uint id1, in uint id2) const {
      if (id1 == id2) return true;
      if (id1 == EMPTY_KEY || id2 == EMPTY_KEY) {
        return false;
      }
      return _keyEqual(value(id1), value(id2));
    }
  }

  alias IdSet = DenseHashSet!(uint, IdHasher, IdKeyEqual);

  alias Hasher = size_t function(in T);
  alias KeyEqual = bool function(in T, in T);

  Hasher _hasher;
  KeyEqual _keyEqual;
  T[] _values;
  IdSet _idSet;
}
