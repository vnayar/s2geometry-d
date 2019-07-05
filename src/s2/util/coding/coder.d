// Copyright 2000 Google Inc. All Rights Reserved.
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

//
//
// This holds the encoding/decoding routines that used to live in netutil

module s2.util.coding.coder;

import std.array;
import std.bitmanip : append, read;
import std.range;
import std.traits;

auto encoder() {
  return new Encoder!(Appender!(ubyte[]))(appender(new ubyte[0]));
}

auto encoder(RangeT)(RangeT r) {
  return new Encoder!RangeT(r);
}

/**
 * Class for encoding data into a memory buffer
 */
class Encoder(RangeT)
if (isOutputRange!(RangeT, ubyte) && isOutputRange!(RangeT, ubyte[])) {
public:
  // Initialize encoder to encode into "buf"
  this(RangeT buf, size_t maxn = size_t.max) {
    _buf = buf;
    _limit = maxn;
  }

  final void reset(RangeT buf, size_t maxn) {
    _buf = buf;
    _limit = maxn;
  }

  // Encoding routines.  Note that these do not check bounds
  final void put8(ubyte v)
  in (avail() >= v.sizeof) {
    _pos += v.sizeof;
    _buf.append!ubyte(v);
  }

  final void put16(ushort v)
  in (avail() >= v.sizeof) {
    _pos += v.sizeof;
    _buf.append!ushort(v);
  }

  final void put32(uint v)
  in (avail() >= v.sizeof) {
    _pos += v.sizeof;
    _buf.append!uint(v);
  }

  final void put64(ulong v)
  in (avail() >= v.sizeof) {
    _pos += v.sizeof;
    _buf.append!ulong(v);
  }

  final void putdouble(double v)
  in (avail() >= v.sizeof) {
    _pos += v.sizeof;
    _buf.append!double(v);
  }

  /**
   * Encodes into the OutputRange a raw object of a type with no indirections,
   * e.g. no arrays, no classes, no pointers, etc.
   *
   * Note: The encoded byte order is machine dependent.
   */
  final void putRaw(T)(T item) @trusted
  if (!hasIndirections!T) {
    ubyte[] itemBytes = *(cast(ubyte[typeof(item).sizeof]*) &item);
    put(_buf, itemBytes);
    _pos += item.sizeof;
  }

  /**
   * Encodes into the OutputRange a range of objects in raw form.
   *
   * Note: The encoded byte order is machine dependent.
   */
  final void putRaw(R)(R items)
  if (isInputRange!R) {
    foreach (ref item; items) {
      ubyte[] itemBytes = *(cast(ubyte[typeof(item).sizeof]*) &item);
      put(_buf, itemBytes);
      _pos += item.sizeof;
    }
  }

  // Return number of bytes encoded so far
  final size_t length() const {
    return _pos;
  }

  // Return number of bytes of space remaining in buffer
  final size_t avail() const {
    return _limit - _pos;
  }

  // REQUIRES: Encoder was created with the 0-argument constructor interface.
  //
  // This interface ensures that at least "N" more bytes are available
  // in the underlying buffer by resizing the buffer (if necessary).
  //
  // Note that no bounds checking is done on any of the put routines,
  // so it is the client's responsibility to call Ensure() at
  // appropriate intervals to ensure that enough space is available
  // for the data being added.
  static if (hasMember!(RangeT, "reserve")) {
    void ensure(size_t n) {
      _buf.reserve(n);
    }
  }

  // Return ptr to start of encoded data.  This pointer remains valid
  // until reset or Ensure is called.
  final RangeT buffer() {
    return _buf;
  }

private:
  // buf_ points into the orig_ buffer, just past the last encoded byte.
  RangeT _buf;

  // limits_ points just past the last allocated byte in the orig_ buffer.
  size_t _limit;
  size_t _pos;
}

@("Encoder.put")
unittest {
  auto enc = encoder();
  // Use hex just to make the bytes easier to identify.
  enc.put8(0x0a);
  assert(enc.buffer().data == [0x0a]);
  assert(enc.length() == 1);

  enc.put16(0x0b0c);
  assert(enc.buffer().data == [0x0a, 0x0b, 0x0c]);
  assert(enc.length() == 3);
}

@("Encoder.putRaw")
unittest {
  auto enc = encoder();
  struct Thing {
    int a;
    double b;
  }

  enc.putRaw(Thing(1, 4.1));
  assert(enc.length() == Thing.sizeof);

  auto things = [Thing(3, 2.1), Thing(2, 3.1)];
  enc.putRaw(things);

  assert(enc.length() == (things.length + 1) * Thing.sizeof);
}

auto decoder(RangeT)(RangeT r) {
  static if (hasLength!RangeT) {
    return new Decoder!RangeT(r, r.length);
  } else {
    return new Decoder!RangeT(r, size_t.max);
  }
}

auto decoder(RangeT)(RangeT r, size_t maxn) {
  return new Decoder!RangeT(r, maxn);
}

/* Class for decoding data from a memory buffer */
class Decoder(RangeT)
if (isInputRange!RangeT && is(ElementType!RangeT == ubyte)) {
public:
  // Initialize decoder to decode from "buf"
  this(RangeT buf, size_t maxn = size_t.max) {
    reset(buf, maxn);
  }

  void reset(RangeT buf, size_t maxn) {
    _buf = buf;
    _limit = maxn;
  }

  // Decoding routines.  Note that these do not check bounds
  ubyte get8()
  in (avail() >= ubyte.sizeof) {
    _pos += ubyte.sizeof;
    return _buf.read!ubyte();
  }

  ushort get16()
  in (avail() >= ushort.sizeof) {
    _pos += ushort.sizeof;
    return _buf.read!ushort();
  }

  uint get32()
  in (avail() >= uint.sizeof) {
    _pos += uint.sizeof;
    return _buf.read!uint();
  }

  ulong get64()
  in (avail() >= ulong.sizeof) {
    _pos += ulong.sizeof;
    return _buf.read!ulong();
  }

  double getdouble()
  in (avail() >= double.sizeof) {
    _pos += double.sizeof;
    return _buf.read!double();
  }

  /**
   * Decodes an object without indirections from raw bytes in the InputRange.
   *
   * Note: The decoding byte order is machine dependent.
   */
  T getRaw(T)() @trusted
  if (!hasIndirections!T) {
    ubyte[T.sizeof] itemBytes = _buf.takeExactly(T.sizeof);
    _buf = _buf.dropExactly(T.sizeof);
    _pos += T.sizeof;
    return *(cast(T*) &itemBytes);
  }

  T[] getRaw(T)(size_t n) {
    T[] output;
    output.reserve(n);
    foreach (i; 0 .. n) {
      output ~= getRaw!T();
    }
    return output;
  }

  void skip(size_t n)
  in (avail() >= n) {
    _buf = _buf.dropExactly(n);
    _pos += n;
  }

  /// Returns number of bytes decoded so far.
  size_t pos() const {
    return _pos;
  }

  /// Returns number of available bytes to read.
  size_t avail() const {
    return _limit - _pos;
  }

private:
  RangeT _buf;

  size_t _pos;
  size_t _limit;
}

version(unittest) {
  auto toNativeBytes(T)(T t) {
    import std.system : endian, Endian;
    import std.bitmanip : nativeToLittleEndian, nativeToBigEndian;

    if (endian == Endian.littleEndian) {
      return nativeToLittleEndian(t);
    } else {
      return nativeToBigEndian(t);
    }
  }
}

@("Decoder.get")
unittest {
  auto dec = decoder(cast(ubyte[]) [0x0a, 0x0b, 0x0c]);
  assert(dec.avail() == 3);
  assert(dec.get8() == 0x0a);
  assert(dec.avail() == 2);
  assert(dec.get16() == 0x0b0c);
  assert(dec.avail() == 0);
  assert(dec.pos() == 3);
}

@("Decoder.getRaw")
unittest {
  struct Thing {
    short a;
    byte b;
    short c;
  }

  // Pick some test data with clearly written bytes in big-endian order (most significant first).
  short testA = 0x0102;
  byte testB = 0x03;
  short testC = 0x0405;

  auto dec = decoder(chain(
          toNativeBytes(testA)[],
          toNativeBytes(testB)[], [cast(ubyte) 0x00],
          toNativeBytes(testC)[])
      .array);

  assert(dec.avail() == Thing.sizeof);
  assert(dec.getRaw!Thing() == Thing(testA, testB, testC));
  assert(dec.avail() == 0);
  assert(dec.pos() == Thing.sizeof);
}

@("Decoder.getRaw")
unittest {
  short[] testData = [0x0102, 0x0304, 0x0506];
  ubyte[] buffer;
  foreach (data; testData) {
    buffer ~= toNativeBytes(data);
  }
  auto dec = decoder(buffer);
  short[] result = dec.getRaw!short(3);
  assert(result == testData);
}
