// Copyright 2013 Google Inc. All Rights Reserved.
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

// Original author: jyrki@google.com (Jyrki Alakuijala)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.util.hash.mix;

// Fast mixing of hash values -- not strong enough for fingerprinting.
// May change from time to time.
//
// Values given are expected to be hashes from good hash functions.
// What constitutes a good hash may depend on your application. As a rule of
// thumb, if std::hash<int> is strong enough for your hashing need if
// your data were just ints, it will most likely be the correct choice
// for a mixed hash of data members. HashMix does one round of multiply and
// rotate mixing, so you get some additional collision avoidance guarantees
// compared to just using std::hash<int> directly.
//
// Possible use:
//
// struct Xyzzy {
//   int x;
//   int y;
//
//   size_t toHash() const {
//     auto mix = HashMix(x);
//     mix.mix(y);
//     return mix.get()
//   }
// }
//

struct HashMix {
private:
  static size_t MUL = 0xdc3eb94af8ab4c93;

  size_t _hash = 1;

public:
  this(size_t val) nothrow @safe {
    _hash = val + 83;
  }

  HashMix mix(size_t val) nothrow @safe {
    _hash *= MUL;
    _hash = ((_hash << 19) |
        (_hash >> size_t.sizeof * 8 - 19)) + val;
    return this;
  }

  size_t get() const nothrow @safe {
    return _hash;
  }
}
