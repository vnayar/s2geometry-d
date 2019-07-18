/**
   A mutex-free lock based on noop loops.

   Copyright: 2018 Google Inc. All Rights Reserved.

   License:
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   $(LINK http://www.apache.org/licenses/LICENSE-2.0)

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS-IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Authors: madric@gmail.com (Vijay Nayar)
*/
module s2.base.spinlock;

import core.atomic;

/// A mutex-free lock based on noop loops.
class SpinLock {
public:

  /// Blocks until `unlock()` is called.
  void lock() shared @safe nothrow @nogc {
    while (!cas(&_flag, false, true)) {
      // Spin.
      continue;
    }
  }

  /// Allows code blocking on `lock()` to continue to execute.
  void unlock() shared @safe nothrow @nogc
  in {
    assert(_flag == true);
  } do {
    while (!cas(&_flag, true, false)) {
      // Spin.
      continue;
    }
  }

private:
  shared bool _flag;
}
