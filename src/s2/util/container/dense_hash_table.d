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
//
// Converted to D:  madric@gmail.com (Vijay Nayar)

// A dense hashtable is a particular implementation of
// a hashtable: one that is meant to minimize memory allocation.
// It does this by using an array to store all the data.  We
// steal a value from the key space to indicate "empty" array
// elements (ie indices where no item lives) and another to indicate
// "deleted" elements.
//
// (Note it is possible to change the value of the delete key
// on the fly; you can even remove it, though after that point
// the hashtable is insert_only until you set it again.  The empty
// value however can't be changed.)
//
// To minimize allocation and pointer overhead, we use internal
// probing, in which the hashtable is a single table, and collisions
// are resolved by trying to insert again in another bucket.  The
// most cache-efficient internal probing schemes are linear probing
// (which suffers, alas, from clumping) and quadratic probing, which
// is what we implement by default.
//
// Type requirements: value_type is required to be Copy Constructible
// and Default Constructible. It is not required to be (and commonly
// isn't) Assignable.
//
// You probably shouldn't use this code directly.  Use dense_hash_map<>
// or dense_hash_set<> instead.

// You can change the following below:
// HT_OCCUPANCY_PCT      -- how full before we double size
// HT_EMPTY_PCT          -- how empty before we halve size
// HT_MIN_BUCKETS        -- default smallest bucket size
//
// You can also change enlarge_factor (which defaults to
// HT_OCCUPANCY_PCT), and shrink_factor (which defaults to
// HT_EMPTY_PCT) with set_resizing_parameters().
//
// How to decide what values to use?
// shrink_factor's default of .4 * OCCUPANCY_PCT, is probably good.
// HT_MIN_BUCKETS is probably unnecessary since you can specify
// (indirectly) the starting number of buckets at construct-time.
// For enlarge_factor, you can use this chart to try to trade-off
// expected lookup time to the space taken up.  By default, this
// code uses quadratic probing, though you can change it to linear
// via JUMP_ below if you really want to.
//
// From http://www.augustana.ca/~mohrj/courses/1999.fall/csc210/lecture_notes/hashing.html
// NUMBER OF PROBES / LOOKUP       Successful            Unsuccessful
// Quadratic collision resolution   1 - ln(1-L) - L/2    1/(1-L) - L - ln(1-L)
// Linear collision resolution     [1+1/(1-L)]/2         [1+1/(1-L)2]/2
//
// -- enlarge_factor --           0.10  0.50  0.60  0.75  0.80  0.90  0.99
// QUADRATIC COLLISION RES.
//    probes/successful lookup    1.05  1.44  1.62  2.01  2.21  2.85  5.11
//    probes/unsuccessful lookup  1.11  2.19  2.82  4.64  5.81  11.4  103.6
// LINEAR COLLISION RES.
//    probes/successful lookup    1.06  1.5   1.75  2.5   3.0   5.5   50.5
//    probes/unsuccessful lookup  1.12  2.5   3.6   8.5   13.0  50.0  5000.0

module s2.util.container.dense_hash_table;

import std.algorithm;
import std.array;
import std.range;
import std.typecons;

/**
Hashtable class, used to implement the hashed associative containers
hash_set and hash_map.
Params:
  Value = what is stored in the table (each bucket is a Value).
  Key = something in a 1-to-1 correspondence to a Value, that can be used
      to search for a Value in the table (find() takes a Key).
  HashFcn = Takes a Key and returns an integer, the more unique the better.
  ExtractKey = given a Value, returns the unique Key associated with it.
      Must inherit from unary_function, or at least have a
      result_type enum indicating the return type of operator().
  SetKey = given a ref Value and a Key, modifies the value such that
      ExtractKey(value) == key.  We guarantee this is only called
      with key == deleted_key or key == empty_key.
*/
class DenseHashTable(Value, Key, alias HashFcn, alias ExtractKey, alias SetKey) {
public:
  // ACCESSOR FUNCTIONS for the things we templatize on, basically.
  // The functions validate at compile time that the accessor functions have the right signature.
  static size_t hash(in Key key) {
    return HashFcn(key);
  }

  static Key getKey(in Value v) {
    return ExtractKey(v);
  }

  static void setKey(ref Value v, Key k) {
    SetKey(v, k);
  }

  alias KeyType = Key;
  alias ValueType = Value;

  alias ThisT = DenseHashTable!(Value, Key, HashFcn, ExtractKey, SetKey);
  alias Iterator = DenseHashTableIterator!(Value, Key, HashFcn, ExtractKey, SetKey);

  /**
  How full we let the table get before we resize.  Knuth says .8 is
  good -- higher causes us to probe too much, though saves memory.
  However, we go with .5, getting better performance at the cost of
  more space (a trade-off densehashtable explicitly chooses to make).
  Feel free to play around with different values, though, via
  max_load_factor() and/or set_resizing_parameters().
  */
  enum size_t HT_OCCUPANCY_PCT = 50;

  /**
  How empty we let the table get before we resize lower, by default.
  (0.0 means never resize lower.)
  It should be less than OCCUPANCY_PCT / 2 or we thrash resizing.
  */
  enum size_t HT_EMPTY_PCT = cast(size_t)(0.4 * HT_OCCUPANCY_PCT);

  /**
  Minimum size we're willing to let hashtables be.
  Must be a power of two, and at least 4.
  Note, however, that for a given hashtable, the initial size is a
  function of the first constructor arg, and may be >HT_MIN_BUCKETS.
  */
  enum size_t HT_MIN_BUCKETS = 4;

  /**
  By default, if you don't specify a hashtable size at
  construction-time, we use this size.  Must be a power of two, and
  at least HT_MIN_BUCKETS.
  */
  enum size_t HT_DEFAULT_STARTING_BUCKETS = 32;

  // ITERATOR FUNCTIONS
  Iterator begin() {
    return Iterator(this, _table[0 .. _numBuckets], true);
  }

  Iterator end() {
    return Iterator(this, _table[_numBuckets .. _numBuckets], true);
  }

  // Accessor function for statistics gathering.
  int numTableCopies() const {
    return _settings.numHtCopies;
  }

private:
  void destroyBuckets(size_t first, size_t last) {
    for ( ; first != last; ++first)
      _table[first] = Value.init;
  }

  // DELETE HELPER FUNCTIONS
  // This lets the user describe a key that will indicate deleted
  // table entries.  This key should be an "impossible" entry --
  // if you try to insert it for real, you won't be able to retrieve it!
  // (NB: while you pass in an entire value, only the key part is looked
  // at.  This is just because I don't know how to assign just a key.)
  void squashDeleted() {
    if (_numDeleted) {
      auto tmp = new DenseHashTable(this);  // copying will get rid of deleted
      swap(tmp);                   // now we are tmp
    }
    assert(_numDeleted == 0);
  }

  // Test if the given key is the deleted indicator.  Requires
  // num_deleted > 0, for correctness of read(), and because that
  // guarantees that key_info.delkey is valid.
  bool testDeletedKey(in Key key) const {
    assert(_numDeleted > 0);
    return _keyValInfo.delKey == key;
  }

public:
  void setDeletedKey(in Key key) {
    // the empty indicator (if specified) and the deleted indicator
    // must be different
    assert(!_settings.useEmpty || key != getKey(_keyValInfo.emptyValue),
        "The deleted-key must be distinct from the the empty-key.");
    // It's only safe to change what "deleted" means if we purge deleted guys
    squashDeleted();
    _settings.useDeleted = true;
    _keyValInfo.delKey = key;
  }

  void clearDeletedKey() {
    squashDeleted();
    _settings.useDeleted = false;
  }
  Key deletedKey() const {
    assert(_settings.useDeleted, "Must set deleted key before calling deletedKey().");
    return _keyValInfo.delKey;
  }

  // These are public so the iterators can use them
  // True if the item at position bucknum is "deleted" marker
  bool testDeleted(size_t bucknum) const {
    // Invariant: !use_deleted() implies num_deleted is 0.
    assert(_settings.useDeleted || _numDeleted == 0);
    return _numDeleted > 0 && testDeletedKey(getKey(_table[bucknum]));
  }
  bool testDeleted(in Iterator it) const {
    // Invariant: !use_deleted() implies num_deleted is 0.
    assert(_settings.useDeleted || _numDeleted == 0);
    return _numDeleted > 0 && testDeletedKey(getKey(*it));
  }

private:
  void checkUseDeleted(string caller) {
    assert(_settings.useDeleted, caller);
  }

  // Set it so test_deleted is true.  true if object didn't used to be deleted.
  bool setDeleted(Iterator it) {
    checkUseDeleted("set_deleted()");
    bool retval = !testDeleted(it);
    // &* converts from iterator to value-type.
    setKey(*it, _keyValInfo.delKey);
    return retval;
  }

  // Set it so test_deleted is false.  true if object used to be deleted.
  bool clearDeleted(Iterator it) {
    checkUseDeleted("clear_deleted()");
    // Happens automatically when we assign something else in its place.
    return testDeleted(it);
  }

  // EMPTY HELPER FUNCTIONS
  // This lets the user describe a key that will indicate empty (unused)
  // table entries.  This key should be an "impossible" entry --
  // if you try to insert it for real, you won't be able to retrieve it!
  // (NB: while you pass in an entire value, only the key part is looked
  // at.  This is just because I don't know how to assign just a key.)
public:
  // These are public so the iterators can use them
  // True if the item at position bucknum is "empty" marker
  bool testEmpty(size_t bucknum) const {
    assert(_settings.useEmpty, "HashTable must call setEmpty() before use!");
    return getKey(_keyValInfo.emptyValue) == getKey(_table[bucknum]);
  }

  bool testEmpty(in Iterator it) const {
    assert(_settings.useEmpty, "HashTable must call setEmpty() before use!");
    return getKey(_keyValInfo.emptyValue) == getKey(*it);
  }

private:
  void fillRangeWithEmpty(size_t table_start, size_t table_end) {
    _table[table_start..table_end] = _keyValInfo.emptyValue;
  }

public:
  // TODO(csilvers): change all callers of this to pass in a key instead,
  //                 and take a const key_type instead of const value_type.
  void setEmptyKey(Value val) {
    // Once you set the empty key, you can't change it
    assert(!_settings.useEmpty, "Calling setEmptyKey() multiple times, which is invalid!");
    // The deleted indicator (if specified) and the empty indicator
    // must be different.
    assert(!_settings.useDeleted || getKey(val) != _keyValInfo.delKey,
        "setEmptyKey() must be called with a key distinct from the deleted-key!");
    _settings.useEmpty = true;
    _keyValInfo = new KeyValInfo();
    _keyValInfo.emptyValue = val;

    assert(_table is null);            // must set before first use
    // num_buckets was set in constructor even though table was not initialized.
    _table = new Value[_numBuckets];
    assert(_table);
    fillRangeWithEmpty(0, _numBuckets);
  }

  Key emptyKey() const {
    assert(_settings.useEmpty);
    return getKey(_keyValInfo.emptyValue);
  }

  // FUNCTIONS CONCERNING SIZE
public:
  size_t size() const {
    return _numElements - _numDeleted;
  }

  size_t maxSize() const {
    return size_t.max;
  }

  bool empty() const {
    return size() == 0;
  }

  size_t bucketCount() const {
    return _numBuckets;
  }

  size_t maxBucketCount() const {
    return maxSize();
  }

  size_t nonemptyBucketCount() const {
    return _numElements;
  }

  // These are tr1 methods.  Their idea of 'bucket' doesn't map well to
  // what we do.  We just say every bucket has 0 or 1 items in it.
  size_t bucketSize(size_t i) /*const*/ {
    return begin() == end() ? 0 : 1;
  }

private:
  // Because of the above, size_t(-1) is never legal; use it for errors
  enum size_t ILLEGAL_BUCKET = cast(size_t) -1;

  // Used after a string of deletes.  Returns true if we actually shrunk.
  // TODO(csilvers): take a delta so we can take into account inserts
  // done after shrinking.  Maybe make part of the Settings class?
  bool maybeShrink()
  in {
    assert(_numElements >= _numDeleted);
    assert((bucketCount() & (bucketCount() - 1)) == 0); // is a power of two
    assert(bucketCount() >= HT_MIN_BUCKETS);
  } body {
    bool retval = false;

    // If you construct a hashtable with < HT_DEFAULT_STARTING_BUCKETS,
    // we'll never shrink until you get relatively big, and we'll never
    // shrink below HT_DEFAULT_STARTING_BUCKETS.  Otherwise, something
    // like "dense_hash_set<int> x; x.insert(4); x.erase(4);" will
    // shrink us down to HT_MIN_BUCKETS buckets, which is too small.
    const size_t num_remain = _numElements - _numDeleted;
    const size_t shrink_threshold = _settings.shrinkThreshold;
    if (shrink_threshold > 0 && num_remain < shrink_threshold
        && bucketCount() > HT_DEFAULT_STARTING_BUCKETS) {
      const float shrink_factor = _settings.shrinkFactor;
      size_t sz = bucketCount() / 2;    // find how much we should shrink
      while (sz > HT_DEFAULT_STARTING_BUCKETS && num_remain < sz * shrink_factor) {
        sz /= 2;                            // stay a power of 2
      }
      auto tmp = new DenseHashTable(this, sz);   // Do the actual resizing
      swap(tmp);                                 // now we are tmp
      retval = true;
    }
    _settings.considerShrink = false;    // because we just considered it
    return retval;
  }

  // We'll let you resize a hashtable -- though this makes us copy all!
  // When you resize, you say, "make it big enough for this many more elements"
  // Returns true if we actually resized, false if size was already ok.
  bool resizeDelta(size_t delta) {
    bool did_resize = false;
    if (_settings.considerShrink) {  // see if lots of deletes happened
      if (maybeShrink())
        did_resize = true;
    }
    if (_numElements >= size_t.max - delta) {
      throw new Exception("resize overflow");
    }
    if ( bucketCount() >= HT_MIN_BUCKETS
        && (_numElements + delta) <= _settings.enlargeThreshold )
      return did_resize;                          // we're ok as we are

    // Sometimes, we need to resize just to get rid of all the
    // "deleted" buckets that are clogging up the hashtable.  So when
    // deciding whether to resize, count the deleted buckets (which
    // are currently taking up room).  But later, when we decide what
    // size to resize to, *don't* count deleted buckets, since they
    // get discarded during the resize.
    size_t needed_size = _settings.minBuckets(_numElements + delta, 0);
    if ( needed_size <= bucketCount() )      // we have enough buckets
      return did_resize;

    size_t resize_to =
      _settings.minBuckets(_numElements - _numDeleted + delta, bucketCount());

    // When num_deleted is large, we may still grow but we do not want to
    // over expand.  So we reduce needed_size by a portion of num_deleted
    // (the exact portion does not matter).  This is especially helpful
    // when min_load_factor is zero (no shrink at all) to avoid doubling
    // the bucket count to infinity.  See also test ResizeWithoutShrink.
    needed_size = _settings.minBuckets(_numElements - _numDeleted / 4 + delta, 0);
    if (resize_to < needed_size &&    // may double resize_to
        resize_to < size_t.max / 2) {
      // This situation means that we have enough deleted elements,
      // that once we purge them, we won't actually have needed to
      // grow.  But we may want to grow anyway: if we just purge one
      // element, say, we'll have to grow anyway next time we
      // insert.  Might as well grow now, since we're already going
      // through the trouble of copying (in order to purge the
      // deleted elements).
      const size_t target = _settings.shrinkSize(resize_to * 2);
      if (_numElements - _numDeleted + delta >= target) {
        // Good, we won't be below the shrink threshhold even if we double.
        resize_to *= 2;
      }
    }
    auto tmp = new DenseHashTable(this, resize_to);
    swap(tmp);                             // now we are tmp
    return true;
  }

  // We require table be not-NULL and empty before calling this.
  void resizeTable(size_t new_size) {
    _table.length = new_size;
  }

  void resizeTable(size_t old_size, size_t new_size) {
    _table = new Value[new_size];
  }

  // Used to actually do the rehashing when we grow/shrink a hashtable
  void copyFrom(/*in*/ DenseHashTable ht, size_t min_buckets_wanted) {
    clearToSize(_settings.minBuckets(ht.size(), min_buckets_wanted));

    // We use a normal iterator to get non-deleted bcks from ht
    // We could use insert() here, but since we know there are
    // no duplicates and no deleted items, we can be more efficient
    assert((bucketCount() & (bucketCount() - 1)) == 0);      // a power of two
    for ( Iterator it = ht.begin(); it != ht.end(); ++it ) {
      size_t num_probes = 0;              // how many times we've probed
      size_t bucknum;
      const size_t bucket_count_minus_one = bucketCount() - 1;
      for (bucknum = hash(getKey(*it)) & bucket_count_minus_one;
           !testEmpty(bucknum);                               // not empty
           bucknum = (bucknum + num_probes) & bucket_count_minus_one) {
        ++num_probes;
        assert(num_probes < bucketCount(),
               "Hashtable is full: an error in key_equal<> or hash<>");
      }
      _table[bucknum] = *it;       // copies the value to here
      _numElements++;
    }
    _settings.numHtCopies++;
  }

  // Required by the spec for hashed associative container
public:
  // Though the docs say this should be num_buckets, I think it's much
  // more useful as num_elements.  As a special feature, calling with
  // req_elements==0 will cause us to shrink if we can, saving space.
  void resize(size_t req_elements) {       // resize to this or larger
    if ( _settings.considerShrink || req_elements == 0 )
      maybeShrink();
    if ( req_elements > _numElements )
      resizeDelta(req_elements - _numElements);
  }

  // Get and change the value of shrink_factor and enlarge_factor.  The
  // description at the beginning of this file explains how to choose
  // the values.  Setting the shrink parameter to 0.0 ensures that the
  // table never shrinks.
  void getResizingParameters(out float shrink, out float grow) const {
    shrink = _settings.shrinkFactor;
    grow = _settings.enlargeFactor;
  }
  void setResizingParameters(float shrink, float grow) {
    _settings.setResizingParameters(shrink, grow);
    _settings.resetThresholds(bucketCount());
  }

  // CONSTRUCTORS -- as required by the specs, we take a size,
  // but also let you specify a hashfunction, key comparator,
  // and key extractor.  We also define a copy constructor and =.
  // DESTRUCTOR -- needs to free the table
  this(size_t expected_max_items_in_table = 0) {
    _settings = new Settings();
    _numDeleted = 0;
    _numElements = 0;
    _numBuckets = expected_max_items_in_table == 0
        ? HT_DEFAULT_STARTING_BUCKETS
        : _settings.minBuckets(expected_max_items_in_table, 0);
    _table = null;
    // table is NULL until emptyval is set.  However, we set num_buckets
    // here so we know how much space to allocate once emptyval is set
    _settings.resetThresholds(bucketCount());
  }

  // As a convenience for resize(), we allow an optional second argument
  // which lets you make this new hashtable a different size than ht
  this(this DenseHashTableT)(DenseHashTableT ht, size_t min_buckets_wanted = HT_DEFAULT_STARTING_BUCKETS) {
    _settings = ht._settings;
    _keyValInfo = ht._keyValInfo;
    _numDeleted = 0;
    _numElements = 0;
    _numBuckets = 0;
    _table = null;

    if (!ht._settings.useEmpty) {
      // If use_empty isn't set, copy_from will crash, so we do our own copying.
      assert(ht.empty());
      _numBuckets = _settings.minBuckets(ht.size(), min_buckets_wanted);
      _settings.resetThresholds(bucketCount());
      return;
    }
    _settings.resetThresholds(bucketCount());
    copyFrom(ht, min_buckets_wanted);   // copy_from() ignores deleted entries
  }

  ~this() {
    if (_table) {
      destroyBuckets(0, _numBuckets);
      //_keyValInfo.deallocate(table, _numBuckets);
    }
  }

  // Many STL algorithms use swap instead of copy constructors
  void swap(DenseHashTable ht) {
    .swap(_settings, ht._settings);
    .swap(_keyValInfo, ht._keyValInfo);
    .swap(_numDeleted, ht._numDeleted);
    .swap(_numElements, ht._numElements);
    .swap(_numBuckets, ht._numBuckets);
    .swap(_keyValInfo, ht._keyValInfo);
    .swap(_table, ht._table);
    _settings.resetThresholds(bucketCount());  // also resets consider_shrink
    ht._settings.resetThresholds(ht.bucketCount());
  }

private:
  void clearToSize(size_t new_num_buckets) {
    if (!_table) {
      // TODO: Use a custom allocator here.
      _table = new Value[new_num_buckets];
    } else {
      destroyBuckets(0, _numBuckets);
      if (new_num_buckets != _numBuckets) {   // resize, if necessary
        resizeTable(_numBuckets, new_num_buckets);
      }
    }
    assert(_table);
    fillRangeWithEmpty(0, new_num_buckets);
    _numElements = 0;
    _numDeleted = 0;
    _numBuckets = new_num_buckets;          // our new size
    _settings.resetThresholds(bucketCount());
  }

public:
  // It's always nice to be able to clear a table without deallocating it
  void clear() {
    // If the table is already empty, and the number of buckets is
    // already as we desire, there's nothing to do.
    const size_t new_num_buckets = _settings.minBuckets(0, 0);
    if (_numElements == 0 && new_num_buckets == _numBuckets) {
      return;
    }
    clearToSize(new_num_buckets);
  }

  // Clear the table without resizing it.
  // Mimicks the stl_hashtable's behaviour when clear()-ing in that it
  // does not modify the bucket count
  void clearNoResize() {
    if (_numElements > 0) {
      assert(_table);
      destroyBuckets(0, _numBuckets);
      fillRangeWithEmpty(0, _numBuckets);
    }
    // don't consider to shrink before another erase()
    _settings.resetThresholds(bucketCount());
    _numElements = 0;
    _numDeleted = 0;
  }

  // LOOKUP ROUTINES
private:

  struct Pair(T1, T2) {
    T1 first;
    T2 second;
  }

  // Returns a pair of positions: 1st where the object is, 2nd where
  // it would go if you wanted to insert it.  1st is ILLEGAL_BUCKET
  // if object is not found; 2nd is ILLEGAL_BUCKET if it is.
  // Note: because of deletions where-to-insert is not trivial: it's the
  // first deleted bucket we see, as long as we don't find the key later
  Pair!(size_t, size_t) findPosition(in Key key) const {
    size_t num_probes = 0;              // how many times we've probed
    const size_t bucket_count_minus_one = bucketCount() - 1;
    size_t bucknum = hash(key) & bucket_count_minus_one;
    size_t insert_pos = ILLEGAL_BUCKET; // where we would insert
    while ( 1 ) {                          // probe until something happens
      if ( testEmpty(bucknum) ) {         // bucket is empty
        if ( insert_pos == ILLEGAL_BUCKET )   // found no prior place to insert
          return Pair!(size_t, size_t)(ILLEGAL_BUCKET, bucknum);
        else
          return Pair!(size_t, size_t)(ILLEGAL_BUCKET, insert_pos);

      } else if ( testDeleted(bucknum) ) { // keep searching, but mark to insert
        if ( insert_pos == ILLEGAL_BUCKET )
          insert_pos = bucknum;

      } else if ( key == getKey(_table[bucknum]) ) {
        return Pair!(size_t, size_t)(bucknum, ILLEGAL_BUCKET);
      }
      ++num_probes;                        // we're doing another probe
      bucknum = (bucknum + num_probes) & bucket_count_minus_one;
      assert(num_probes < bucketCount(), "Hashtable is full: an error in key_equal<> or hash<>");
    }
  }

public:

  Iterator find(in Key key) {
    if (size() == 0) return end();
    auto pos = findPosition(key);
    if ( pos.first == ILLEGAL_BUCKET )     // alas, not there
      return end();
    else
      return Iterator(this, _table[pos.first .. _numBuckets], false);
  }

  // This is a tr1 method: the bucket a given key is in, or what bucket
  // it would be put in, if it were to be inserted.  Shrug.
  size_t bucket(in Key key) const {
    auto pos = findPosition(key);
    return pos.first == ILLEGAL_BUCKET ? pos.second : pos.first;
  }

  // Counts how many elements have key key.  For maps, it's either 0 or 1.
  size_t count(in Key key) const {
    auto pos = findPosition(key);
    return pos.first == ILLEGAL_BUCKET ? 0 : 1;
  }

  // Likewise, equal_range doesn't really make sense for us.  Oh well.
  Pair!(Iterator, Iterator) equalRange(in Key key) {
    Iterator pos = find(key);      // either an iterator or end
    if (pos == end()) {
      return Pair!(Iterator, Iterator)(pos, pos);
    } else {
      Iterator startpos = pos++;
      return Pair!(Iterator, Iterator)(startpos, pos);
    }
  }

  // INSERTION ROUTINES
private:
  // Private method used by insert_noresize and find_or_insert.
  Iterator insertAt(Value obj, size_t pos) {
    if (size() >= maxSize()) {
      throw new Exception("insert overflow");
    }
    if ( testDeleted(pos) ) {      // just replace if it's been del.
      // shrug: shouldn't need to be const.
      auto delpos = Iterator(this, _table[pos .. _numBuckets], false);
      clearDeleted(delpos);
      assert( _numDeleted > 0);
      --_numDeleted;                // used to be, now it isn't
    } else {
      ++_numElements;               // replacing an empty bucket
    }
    _table[pos] = obj;
    return Iterator(this, _table[pos .. _numBuckets], false);
  }

  // If you know *this is big enough to hold obj, use this routine
  Pair!(Iterator, bool) insertNoResize(Value obj) {
    // First, double-check we're not inserting delkey or emptyval
    assert(!_settings.useEmpty || getKey(obj) != getKey(_keyValInfo.emptyValue),
        "Inserting the empty key");
    assert(!_settings.useDeleted || getKey(obj) != _keyValInfo.delKey,
        "Inserting the deleted key");
    const Pair!(size_t, size_t) pos = findPosition(getKey(obj));
    if ( pos.first != ILLEGAL_BUCKET) {      // object was already there
      return Pair!(Iterator, bool)(
          Iterator(this, _table[pos.first .. _numBuckets], false),
          false);          // false: we didn't insert
    } else {               // pos.second says where to put it
      return Pair!(Iterator, bool)(insertAt(obj, pos.second), true);
    }
  }

public:
  // This is the normal insert routine, used by the outside world
  Pair!(Iterator, bool) insert(Value obj) {
    resizeDelta(1);                      // adding an object, grow if need be
    return insertNoResize(obj);
  }

  // Specializations of insert(it, it) depending on the power of the iterator:
  void insert(Iterator f, Iterator l) {
    for (; f != l; f++)
      insert(*f);
  }

  // DefaultValue is a functor that takes a key and returns a value_type
  // representing the default value to be inserted if none is found.
  Value findOrInsert(Value function(Key) defaultValue)(in Key key) {
    // First, double-check we're not inserting emptykey or delkey
    assert(!_settings.useEmpty || key != getKey(_keyValInfo.emptyValue),
        "Inserting the empty key");
    assert(!_settings.useDeleted || key != _keyValInfo.delKey,
        "Inserting the deleted key");
    Pair!(size_t,size_t) pos = findPosition(key);
    if (pos.first != ILLEGAL_BUCKET) {  // object was already there
      return _table[pos.first];
    } else if (resizeDelta(1)) {        // needed to rehash to make room
      // Since we resized, we can't use pos, so recalculate where to insert.
      return insertNoResize(defaultValue(key)).first;
    } else {                             // no need to rehash, insert right here
      return insertAt(defaultValue(key), pos.second);
    }
  }

  // DELETION ROUTINES
  size_t erase(in Key key) {
    // First, double-check we're not trying to erase delkey or emptyval.
    assert(!_settings.useEmpty || key != getKey(_keyValInfo.emptyValue),
        "Erasing the empty key");
    assert(!_settings.useDeleted || key != _keyValInfo.delKey,
        "Erasing the deleted key");
    Iterator pos = find(key);   // shrug: shouldn't need to be const
    if ( pos != end() ) {
      assert(!testDeleted(pos));  // or find() shouldn't have returned it
      setDeleted(pos);
      ++_numDeleted;
      _settings.considerShrink = true; // will think about shrink after next insert
      return 1;                    // because we deleted one thing
    } else {
      return 0;                    // because we deleted nothing
    }
  }

  // We return the iterator past the deleted item.
  void erase(Iterator pos) {
    if ( pos == end() ) return;    // sanity check
    if ( setDeleted(pos) ) {      // true if object has been newly deleted
      ++_numDeleted;
      _settings.considerShrink = true; // will think about shrink after next insert
    }
  }

  void erase(Iterator f, Iterator l) {
    for ( ; f != l; ++f) {
      if (setDeleted(f))       // should always be true
        ++_numDeleted;
    }
    _settings.considerShrink = true; // will think about shrink after next insert
  }

  // COMPARISON
  override
  bool opEquals(in Object o) {
    ThisT ht = cast(ThisT) o;
    if (size() != ht.size()) {
      return false;
    } else if (this is ht) {
      return true;
    } else {
      // Iterate through the elements in "this" and see if the
      // corresponding element is in ht
      for ( Iterator it = begin(); it != end(); ++it ) {
        Iterator it2 = ht.find(getKey(*it));
        if (it2 == ht.end() || *it != *it2) {
          return false;
        }
      }
      return true;
    }
  }

  // I/O
  // We support reading and writing hashtables to disk.  Alas, since
  // I don't know how to write a hasher or key_equal, you have to make
  // sure everything but the table is the same.  We compact before writing.
private:
  // Every time the disk format changes, this should probably change too
  alias MagicNumberType = ulong;
  enum MagicNumberType MAGIC_NUMBER = 0x13578642;

public:
  // I/O -- this is an add-on for writing hash table to disk
  //
  // INPUT and OUTPUT must be either a FILE, *or* a C++ stream
  //    (istream, ostream, etc) *or* a class providing
  //    Read(void*, size_t) and Write(const void*, size_t)
  //    (respectively), which writes a buffer into a stream
  //    (which the INPUT/OUTPUT instance presumably owns).

  // typedef sparsehash_internal::pod_serializer<value_type> NopointerSerializer;

  // TODO: Convert if serialization logic is needed.
  // ValueSerializer: a functor.  operator()(OUTPUT*, const value_type&)
  // template <typename ValueSerializer, typename OUTPUT>
  // bool serialize(ValueSerializer serializer, OUTPUT *fp) {
  //   squash_deleted();           // so we don't have to worry about delkey
  //   if ( !sparsehash_internal::write_bigendian_number(fp, MAGIC_NUMBER, 4) )
  //     return false;
  //   if ( !sparsehash_internal::write_bigendian_number(fp, num_buckets, 8) )
  //     return false;
  //   if ( !sparsehash_internal::write_bigendian_number(fp, num_elements, 8) )
  //     return false;
  //   // Now write a bitmap of non-empty buckets.
  //   for ( size_t i = 0; i < num_buckets; i += 8 ) {
  //     unsigned char bits = 0;
  //     for ( int bit = 0; bit < 8; ++bit ) {
  //       if ( i + bit < num_buckets && !test_empty(i + bit) )
  //         bits |= (1 << bit);
  //     }
  //     if ( !sparsehash_internal::write_data(fp, &bits, sizeof(bits)) )
  //       return false;
  //     for ( int bit = 0; bit < 8; ++bit ) {
  //       if ( bits & (1 << bit) ) {
  //         if ( !serializer(fp, table[i + bit]) ) return false;
  //       }
  //     }
  //   }
  //   return true;
  // }

  // TODO: Convert if serialization logic is needed.
  // INPUT: anything we've written an overload of read_data() for.
  // ValueSerializer: a functor.  operator()(INPUT*, value_type*)
  // template <typename ValueSerializer, typename INPUT>
  // bool unserialize(ValueSerializer serializer, INPUT *fp) {
  //   assert(settings.use_empty() && "empty_key not set for read");

  //   clear();                        // just to be consistent
  //   MagicNumberType magic_read;
  //   if ( !sparsehash_internal::read_bigendian_number(fp, &magic_read, 4) )
  //     return false;
  //   if ( magic_read != MAGIC_NUMBER ) {
  //     return false;
  //   }
  //   size_t new_num_buckets;
  //   if ( !sparsehash_internal::read_bigendian_number(fp, &new_num_buckets, 8) )
  //     return false;
  //   clear_to_size(new_num_buckets);
  //   if ( !sparsehash_internal::read_bigendian_number(fp, &num_elements, 8) )
  //     return false;

  //   // Read the bitmap of non-empty buckets.
  //   for (size_t i = 0; i < num_buckets; i += 8) {
  //     unsigned char bits;
  //     if ( !sparsehash_internal::read_data(fp, &bits, sizeof(bits)) )
  //       return false;
  //     for ( int bit = 0; bit < 8; ++bit ) {
  //       if ( i + bit < num_buckets && (bits & (1 << bit)) ) {  // not empty
  //         if ( !serializer(fp, &table[i + bit]) ) return false;
  //       }
  //     }
  //   }
  //   return true;
  // }

private:
  // Package functors with another class to eliminate memory needed for zero-size functors.
  static class Settings {
    size_t enlargeThreshold = 0;  // table.size() * enlarge_factor
    size_t shrinkThreshold = 0;   // table.size() * shrink_factor
    float enlargeFactor = HT_OCCUPANCY_PCT / 100.0f;  // how full before resize
    float shrinkFactor = HT_EMPTY_PCT / 100.0f;  // how empty before resize
    // consider_shrink=true if we should try to shrink before next insert
    bool considerShrink = false;
    bool useEmpty = false;    // used only by densehashtable, not sparsehashtable
    bool useDeleted = false;  // false until delkey has been set
    // num_ht_copies is a counter incremented every Copy/Move
    uint numHtCopies = 0;

    size_t enlargeSize(size_t x) const {
      return cast(size_t)(x * enlargeFactor);
    }

    size_t shrinkSize(size_t x) const {
      return cast(size_t)(x * shrinkFactor);
    }

    // Reset the enlarge and shrink thresholds
    void resetThresholds(size_t num_buckets) {
      enlargeThreshold = enlargeSize(num_buckets);
      shrinkThreshold = shrinkSize(num_buckets);
      // whatever caused us to reset already considered
      considerShrink = false;
    }

    // Caller is resposible for calling reset_threshold right after
    // set_resizing_parameters.
    void setResizingParameters(float shrink, float grow) {
      assert(shrink >= 0.0);
      assert(grow <= 1.0);
      if (shrink > grow / 2.0f)
        shrink = grow / 2.0f;     // otherwise we thrash hashtable size
      shrinkFactor = shrink;
      enlargeFactor = grow;
    }

    // This is the smallest size a hashtable can be without being too crowded
    // If you like, you can give a min #buckets as well as a min #elts
    size_t minBuckets(size_t num_elts, size_t min_buckets_wanted) {
      float enlarge = enlargeFactor;
      size_t sz = HT_MIN_BUCKETS;             // min buckets allowed
      while ( sz < min_buckets_wanted || num_elts >= cast(size_t)(sz * enlarge) ) {
        // This just prevents overflowing size_t, since sz can exceed
        // max_size() here.
        if (cast(size_t)(sz * 2) < sz) {
          throw new Exception("resize overflow");  // protect against overflow
        }
        sz *= 2;
      }
      return sz;
    }

  }

  // A class is passed by reference, meaning large numbers of hash-tables can share the same data.
  static class KeyValInfo {
    Value emptyValue;
    Key delKey;
  }

private:
  // Actual data
  Settings _settings;
  KeyValInfo _keyValInfo;

  size_t _numDeleted;  // how many occupied buckets are marked deleted
  size_t _numElements;
  size_t _numBuckets;
  Value[] _table;
}

/**
 * A basic iterator type for finding entries and iterating.
 * We're just an array, but we need to skip over empty and deleted elements.
 *
 * TODO(vnayar): Coonvert DenseHashTable to be range based after S2Builder is modified to use it.
 */
struct DenseHashTableIterator(Value, Key, alias HashFcn, alias ExtractKey, alias SetKey) {
public:
  alias Iterator = DenseHashTableIterator!(Value, Key, HashFcn, ExtractKey, SetKey);

  // "Real" constructor and default constructor
  this(in DenseHashTable!(Value, Key, HashFcn, ExtractKey, SetKey) h, Value[] data, bool advance) {
    _ht = h;
    _data = data;
    if (advance) advancePastEmptyAndDeleted();
  }

  // Happy dereferencer
  ref inout(Value) opUnary(string op)() inout
    if (op == "*")
  in {
    assert(!_data.empty(), "Iterator is already at end!");
  } body {
    // return _data.front(); -- Bug: See https://issues.dlang.org/show_bug.cgi?id=19518
    return _data[0];
  }

  // Arithmetic.  The only hard part is making sure that
  // we're not on an empty or marked-deleted array element
  void advancePastEmptyAndDeleted() {
    while ( !_data.empty() && (_ht.testEmpty(this) || _ht.testDeleted(this)) )
      _data.popFront();
  }

  DenseHashTableIterator opUnary(string op)() if (op == "++")
  in {
    assert(!_data.empty());
  } body {
    _data.popFront();
    advancePastEmptyAndDeleted();
    return this;
  }

  // Comparison.
  bool opEquals(this ThisT)(ThisT o) const {
    return _data == o._data;
  }

  // The actual data
  Rebindable!(const DenseHashTable!(Value, Key, HashFcn, ExtractKey, SetKey)) _ht;
  Value[] _data;
}
