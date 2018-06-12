// Copyright 2018 Google Inc. All Rights Reserved.
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

// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.util.math.vector;

// Simple classes to handle vectors in 2D, 3D, and 4D.
//
// Maintainers: Please be mindful of extreme degradations in unoptimized
// performance here.

import algorithm = std.algorithm;
import conv = std.conv;
import format = std.format;
import std.math;
import range = std.range;
import traits = std.traits;
import s2.util.hash.mix;

// CRTP base class for all Vector templates.
struct Vector(ElemT, size_t SizeV)
if (SizeV >= 1) {
protected:

  alias ThisT = Vector!(ElemT, SizeV);
  // FloatType is the type returned by Norm() and Angle().  These methods are
  // special because they return floating-point values even when VType is an
  // integer.
  static if (traits.isIntegral!ElemT) {
    alias FloatT = double;
  } else {
    alias FloatT = ElemT;
  }

  static if (traits.isNumeric!ElemT) {
    ElemT[SizeV] _data = 0;
  } else {
    ElemT[SizeV] _data;
  }

public:

  this(ThisT v) {
    _data = v.data;
  }

  this(ElemT v) {
    _data[] = v;
  }

  this(ElemT[SizeV] v ...) {
    _data = v;
  }

  this(RangeT)(RangeT data)
  if (range.isInputRange!RangeT) {
    _data[] = conv.to!(ElemT[])(range.take(range.array(data), SizeV));
  }

  @property
  ref ElemT[SizeV] data() {
    return _data;
  }

  static if (SizeV >= 1) {
    @property
    void x(ElemT v) {
      _data[0] = v;
    }

    @property
    ElemT x() const {
      return _data[0];
    }
  }

  static if (SizeV >= 2) {
    @property
    void y(ElemT v) {
      _data[1] = v;
    }

    @property
    ElemT y() const {
      return _data[1];
    }
  }

  static if (SizeV >= 3) {
    @property
    void z(ElemT v) {
      _data[2] = v;
    }

    @property
    ElemT z() const {
      return _data[2];
    }
  }

  static if (SizeV >= 4) {
    @property
    void w(ElemT v) {
      _data[3] = v;
    }

    @property
    ElemT w() const {
      return _data[3];
    }
  }

  static size_t size() { return SizeV; }

  // Convert from another vector type
  static ThisT from(Elem2T)(in Vector!(Elem2T, SizeV) b) {
    return ThisT(algorithm.map!(elem => cast(ElemT) elem)(b._data[]));
  }

  static if (is(typeof(ElemT.nan) : ElemT)) {
    // A Vector populated with all NaN values.
    static ThisT nan() {
      ThisT retvalue = ThisT();
      retvalue._data[] = ElemT.nan;
      return retvalue;
    }
  }

  void clear() { _data = _data.init; }

  // Implement index assignment of the form: v[i1] = value
  void opIndexAssign(ElemT value, size_t i1)
  in {
    assert(i1 >= 0);
    assert(i1 < SizeV);
  } body {
    _data[i1] = value;
  }

  // Implement indexing of the form: v[i1]
  ref inout(ElemT) opIndex(size_t i1) inout
  in {
    assert(i1 >= 0);
    assert(i1 < SizeV);
  } body {
    return _data[i1];
  }

  // Support the == and != operators.
  bool opEquals(in ThisT v) const {
    return _data == v._data;
  }

  // Support the <, <=, >, and >= operators.
  int opCmp(in ThisT v) const {
    foreach (size_t i; 0 .. SizeV) {
      if (_data[i] == v._data[i]) {
        continue;
      }
      return _data[i] > v._data[i] ? 1 : -1;
    }
    return 0;
  }

  // Support the +=, -=, *=, and /= operators.
  ThisT opOpAssign(string op)(ThisT v) {
    static if (op == "+") {
      _data[] += v._data[];
    } else static if (op == "-") {
      _data[] -= v._data[];
    } else static if (op == "*") {
      _data[] *= v._data[];
    } else static if (op == "/") {
      _data[] /= v._data[];
    } else {
      static assert(false, "Operator " ~ op ~ "= not implemented!");
    }
    return cast(ThisT) this;
  }

  ThisT opOpAssign(string op, ScalarT)(ScalarT v) {
    mixin("_data[] " ~ op ~ "= v;");
    return this;
  }

  // Support the +, -, *, and / operators.
  ThisT opBinary(string op)(ThisT v) const {
    ThisT retval = ThisT(this);
    static if (op == "+") {
      retval += v;
    } else static if (op == "-") {
      retval -= v;
    } else static if (op == "*") {
      retval *= v;
    } else static if (op == "/") {
      retval /= v;
    } else {
      static assert(0, "Operator " ~ op ~ " not implemented.");
    }
    return retval;
  }

  // Support the +, -, *, and / operators.
  ThisT opBinary(string op)(in ElemT k) const
  if (op == "+" || op == "-" || op == "*" || op == "/") {
    ThisT retval = ThisT();
    static foreach(i; 0..SizeV) {
      mixin("retval._data[i] = _data[i] " ~ op ~ " k;");
    }
    return retval;
  }

  // Support the +, -, *, and / operators.
  ThisT opBinaryRight(string op)(in ElemT k) const
  if (op == "+" || op == "-" || op == "*" || op == "/") {
    ThisT retval = ThisT();
    static foreach(i; 0..SizeV) {
      mixin("retval._data[i] = k " ~ op ~ " _data[i];");
    }
    return retval;
  }

  // Support negation.
  ThisT opUnary(string op)() const {
    ThisT retval = ThisT(this);
    static if (op == "-") {
      retval._data[] = -retval._data[];
    } else {
      static assert(0, "Unary operator " ~ op ~ " not implemented.");
    }
    return retval;
  }

  // Element-wise max.  {max(a[0],b[0]), max(a[1],b[1]), ...}
  static ThisT max(in ThisT a, in ThisT b) {
    return ThisT(algorithm.map!(vals => algorithm.max(vals[0], vals[1]))(
        range.zip(a._data[], b._data[])));
  }

  // Element-wise min.  {min(a[0],b[0]), min(a[1],b[1]), ...}
  static ThisT min(in ThisT a, in ThisT b) {
    return ThisT(algorithm.map!(vals => algorithm.min(vals[0], vals[1]))(
        range.zip(a._data[], b._data[])));
  }

  ElemT dotProd(in ThisT b) const {
    return cast(ElemT) algorithm.sum(algorithm.map!(vals => vals[0] * vals[1])(
        range.zip(_data[], b._data[])));
  }

  // Squared Euclidean norm (the dot product with itself).
  ElemT norm2() const {
    return cast(ElemT) dotProd(cast(ThisT) this);
  }

  // Normalized vector if the norm is nonzero. Not for integer types.
  static if (traits.isFloatingPoint!ElemT) {
    // Euclidean norm. For integer T, correct only if Norm2 does not overflow.
    ElemT norm() const {
      return norm2().sqrt();
    }

    ThisT normalize() const {
      ElemT n = norm();
      if (n != cast(ElemT) 0.0) {
        n = cast(ElemT) 1.0 / n;
      }
      return ThisT(algorithm.map!(a => a * n)(_data[]));
    }

    // Compose a vector from the sqrt of each component.
    ThisT sqrt() const {
      return ThisT(algorithm.map!(a => a.sqrt())(_data[]));
    }

    // Take the floor of each component.
    ThisT floor() const {
      return ThisT(algorithm.map!(a => a.floor())(_data[]));
    }

    // Take the ceil of each component.
    ThisT ceil() const {
      return ThisT(algorithm.map!(a => a.ceil())(_data[]));
    }

    // Round of each component. Formerly FRound().
    ThisT round() const {
      return ThisT(algorithm.map!(a => a.rint())(_data[]));
    }

    // Round of each component and return an integer vector. Formerly IRound().
    Vector!(int, SizeV) toIntVector() const {
      return Vector!(int, SizeV)(algorithm.map!(a => cast(int) a.lrint())(_data[]));
    }

    // True if any of the components is not a number.
    bool isNaN() const {
      return algorithm.any!(a => a.isNaN())(_data[]);
    }

    ThisT fabs() const {
      return ThisT(algorithm.map!(a => a.fabs())(_data[]));
    }
  }

  // Compute the absolute value of each component.
  ThisT abs() const {
    return ThisT(algorithm.map!(a => a.abs())(_data[]));
  }

  // Approximately equal.
  bool aequal(in ThisT vb, FloatT margin) const {
    return algorithm.all!(vals => (vals[0] - vals[1]).abs() < margin)(
        range.zip(_data[], vb._data[]));
  }

  string toString() const {
    string val = "[";
    string sep = "";
    foreach (int i; 0 .. SizeV) {
      val ~= sep ~ conv.to!string(_data[i]);
      sep = ", ";
    }
    return val ~ "]";
  }

  static if (traits.isNumeric!ElemT) {
    size_t toHash() const nothrow @safe {
      HashMix h = HashMix(cast(size_t) _data[0]);
      foreach (d; _data[1 .. $]) {
        h.mix(cast(size_t) d);
      }
      return h.get();
    }
  } else {
    size_t toHash() const nothrow @safe {
      HashMix h = HashMix(_data[0].toHash());
      foreach (d; _data[1 .. $]) {
        h.mix(d.toHash());
      }
      return h.get();
    }
  }

  // Function implementations specific to 2-dimentional vectors.
  static if (SizeV == 2 && traits.isSigned!ElemT) {
    // Cross product.  Be aware that if T is an integer type, the high bits
    // of the result are silently discarded.
    ElemT crossProd(in ThisT vb) const {
      return x * vb.y - y * vb.x;
    }

    // return the angle between "this" and v in radians
    FloatT angle(in ThisT v) const {
      return atan2(cast(FloatT) crossProd(v), cast(FloatT) dotProd(v));
    }

    // return a vector orthogonal to the current one
    // with the same norm and counterclockwise to it
    ThisT ortho() const {
      return ThisT([-y, x]);
    }
  }

  // Function implementations speficif to 3-dimentional vectors.
  static if (SizeV == 3) {
    // Cross product.  Be aware that if VType is an integer type, the high bits
    // of the result are silently discarded.
    ThisT crossProd(in ThisT vb) const {
      return ThisT([
        y * vb.z - z * vb.y,
        z * vb.x - x * vb.z,
        x * vb.y - y * vb.x]);
    }

    static if (traits.isFloatingPoint!ElemT) {
      // Returns a unit vector orthogonal to this one.
      ThisT ortho() const {
        int k = largestAbsComponent() - 1;
        if (k < 0) {
          k = 2;
        }
        ThisT temp = ThisT(0);
        temp[k] = cast(ElemT) 1;
        return crossProd(temp).normalize();
      }

      // return the angle between two vectors in radians
      FloatT angle(in ThisT va) const {
        return cast(FloatT) crossProd(va).norm().atan2(dotProd(va));
      }
    }

    // return the index of the largest component (fabs)
    int largestAbsComponent() const {
      ThisT temp = abs();
      return temp[0] > temp[1]
          ? temp[0] > temp[2] ? 0 : 2
          : temp[1] > temp[2] ? 1 : 2;
    }

    // return the index of the smallest, median ,largest component of the vector
    Vector!(int, 3) componentOrder() const {
      auto temp = Vector!(int, 3)([0, 1, 2]);
      if (_data[temp[0]] > _data[temp[1]]) {
        algorithm.swap(temp.data[0], temp.data[1]);
      }
      if (_data[temp[1]] > _data[temp[2]]) {
        algorithm.swap(temp.data[1], temp.data[2]);
      }
      if (_data[temp[0]] > _data[temp[1]]) {
        algorithm.swap(temp.data[0], temp.data[1]);
      }
      return temp;
    }

  }
}

template Vector2(ElemT) {
  alias Vector2 = Vector!(ElemT, 2);
}
alias Vector2_b = Vector!(ubyte, 2);
alias Vector2_i = Vector!(int, 2);
alias Vector2_f = Vector!(float, 2);
alias Vector2_d = Vector!(double, 2);
alias Vector2_r = Vector!(real, 2);

template Vector3(ElemT) {
  alias Vector3 = Vector!(ElemT, 3);
}
alias Vector3_b = Vector!(ubyte, 3);
alias Vector3_i = Vector!(int, 3);
alias Vector3_f = Vector!(float, 3);
alias Vector3_d = Vector!(double, 3);
alias Vector3_r = Vector!(real, 3);

template Vector4(ElemT) {
  alias Vector4 = Vector!(ElemT, 4);
}
alias Vector4_b = Vector!(ubyte, 4);
alias Vector4_i = Vector!(int, 4);
alias Vector4_f = Vector!(float, 4);
alias Vector4_d = Vector!(double, 4);
alias Vector4_r = Vector!(real, 4);
