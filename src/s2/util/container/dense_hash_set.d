// Copyright (c) 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// ---
//
// This is just a very thin wrapper over dense_hash_table.d, just
// like sgi stl's stl_hash_set is a very thin wrapper over
// stl_hashtable.  The major thing we define is operator[], because
// we have a concept of a data_type which stl_hashtable doesn't
// (it only has a key and a value).
//
// This is more different from dense_hash_map than you might think,
// because all iterators for sets are const (you obviously can't
// change the key, and for sets there is no value).
//
// NOTE: this is exactly like sparse_hash_set.h, with the word
// "sparse" replaced by "dense", except for the addition of
// set_empty_key().
//
//   YOU MUST CALL SET_EMPTY_KEY() IMMEDIATELY AFTER CONSTRUCTION.
//
// Otherwise your program will die in mysterious ways.  (Note if you
// use the constructor that takes an InputIterator range, you pass in
// the empty key in the constructor, rather than after.  As a result,
// this constructor differs from the standard STL version.)
//
// In other respects, we adhere mostly to the STL semantics for
// hash-map.  One important exception is that insert() may invalidate
// iterators entirely -- STL semantics are that insert() may reorder
// iterators, but they all still refer to something valid in the
// hashtable.  Not so for us.  Likewise, insert() may invalidate
// pointers into the hashtable.  (Whether insert invalidates iterators
// and pointers depends on whether it results in a hashtable resize).
// On the plus side, delete() doesn't invalidate iterators or pointers
// at all, or even change the ordering of elements.
//
// Here are a few "power user" tips:
//
//    1) set_deleted_key():
//         If you want to use erase() you must call set_deleted_key(),
//         in addition to set_empty_key(), after construction.
//         The deleted and empty keys must differ.
//
//    2) resize(0):
//         When an item is deleted, its memory isn't freed right
//         away.  This allows you to iterate over a hashtable,
//         and call erase(), without invalidating the iterator.
//         To force the memory to be freed, call resize(0).
//         For tr1 compatibility, this can also be called as rehash(0).
//
//    3) min_load_factor(0.0)
//         Setting the minimum load factor to 0.0 guarantees that
//         the hash table will never shrink.
//
// Roughly speaking:
//   (1) dense_hash_set: fastest, uses the most memory unless entries are small
//   (2) sparse_hash_set: slowest, uses the least memory
//   (3) hash_set / unordered_set (STL): in the middle
//
// Typically I use sparse_hash_set when I care about space and/or when
// I need to save the hashtable on disk.  I use hash_set otherwise.  I
// don't personally use dense_hash_set ever; some people use it for
// small sets with lots of lookups.
//
// - dense_hash_set has, typically, about 78% memory overhead (if your
//   data takes up X bytes, the hash_set uses .78X more bytes in overhead).
// - sparse_hash_set has about 4 bits overhead per entry.
// - sparse_hash_set can be 3-7 times slower than the others for lookup and,
//   especially, inserts.  See time_hash_map.cc for details.
//
// See /usr/(local/)?doc/sparsehash-*/dense_hash_set.html
// for information about how to use this class.

// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.util.container.dense_hash_set;

import s2.util.container.dense_hash_table;

// Default functions for initializing the contained DenseHashTable.

size_t hash(ValueT)(ValueT value)
if (is(typeof(value.toHash()) : size_t)) {
  return value.toHash();
}

size_t hash(ValueT)(ValueT value) {
  return typeid(value).getHash(&value);
}

class DenseHashSet(Value, alias HashFcn = hash!Value) {
public:
  // Apparently identity is not stl-standard, so we define our own
  static Value identity(Value v) {
    return v;
  }

  // Key and Value are the same type in a DenseHashSet.
  static void setKey(ref Value value, Value key) {
    value = key;
  }

  // The actual data
  alias HashTable = DenseHashTable!(Value, Value, HashFcn, identity, setKey);
  HashTable rep;
  alias rep this;

  alias KeyType = HashTable.KeyType;
  alias ValueType = HashTable.ValueType;
  alias Iterator = HashTable.Iterator;

  // Constructors
  this(size_t expected_max_items_in_table = 0) {
    rep = new HashTable(expected_max_items_in_table);
  }

  /**
  A constructor based on interators which provide the values to add to the set.
  With a DenseHashSet, the key and value types used for the DenseHashTable are the same.

  Params:
    InputIterator = Compile-time type parameter of the iterators that support the "*" operator.
    f = Iterator for the first value.
    l = Iterator for the last value.
    empty_key_val = An unused value that can be used to represent an "empty" hash slot.
    expected_max_items_in_table = Sets an initial size to help avoid additional memory allocations.
  */
  this(InputIterator)(
      InputIterator f, InputIterator l,
      Value empty_key_val, size_t expected_max_items_in_table = 0)
  if (is(typeof(*(InputIterator.init)) : Value)) {
    rep = new HashTable(expected_max_items_in_table);
    set_empty_key(empty_key_val);
    rep.insert(f, l);
  }

  // We use the default copy constructor
  // We use the default operator=()
  // We use the default destructor

  // These are tr1 methods.  bucket() is the bucket the key is or would be in.
  float loadFactor() const {
    return size() * 1.0f / bucketCount();
  }

  float max_load_factor() const {
    float shrink, grow;
    rep.getResizingParameters(shrink, grow);
    return grow;
  }

  void maxLoadFactor(float new_grow) {
    float shrink, grow;
    rep.getResizingParameters(shrink, grow);
    rep.setResizingParameters(shrink, new_grow);
  }

  // These aren't tr1 methods but perhaps ought to be.
  float minLoadFactor() const {
    float shrink, grow;
    rep.getResizingParameters(shrink, grow);
    return shrink;
  }

  void minLoadFactor(float new_shrink) {
    float shrink, grow;
    rep.getResizingParameters(shrink, grow);
    rep.setResizingParameters(new_shrink, grow);
  }

  // Comparison
  override
  bool opEquals(this ThisT)(Object o) {
    ThisT hs = cast(ThisT) o;
    if (hs is null) return false;
    return rep == hs.rep;
  }

}

void swap(Value, HashFcn)(DenseHashSet!(Value, HashFcn) hs1, DenseHashSet!(Value, HashFcn) hs2) {
  hs1.swap(hs2);
}
