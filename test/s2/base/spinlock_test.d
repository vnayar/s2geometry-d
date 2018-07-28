module s2.base.spinlock_test;

import fluent.asserts;
import s2.base.spinlock;

import std.concurrency : spawnLinked, receiveOnly, LinkTerminated;
import core.thread : Thread;
import core.time : msecs;

class SharedValue {
  int counter;
}

void testLock(shared SpinLock spinLock, shared SharedValue sharedValue) {
  foreach (i; 0 .. 200) {
    spinLock.lock();
    Thread.sleep(msecs(5));
    // An unsafe operation for shared variables.
    sharedValue.counter = sharedValue.counter + 1;
    spinLock.unlock();
  }
}

@("lock")
unittest {
  shared SpinLock spinLock = new SpinLock();
  shared SharedValue sharedValue = new SharedValue();

  auto tid1 = spawnLinked(&testLock, spinLock, sharedValue);
  auto tid2 = spawnLinked(&testLock, spinLock, sharedValue);
  receiveOnly!LinkTerminated();
  receiveOnly!LinkTerminated();

  Assert.equal(sharedValue.counter, 400);
}
