module s2.util.math.vector;

// Simple classes to handle vectors in 2D, 3D, and 4D.
//
// Maintainers: Please be mindful of extreme degradations in unoptimized
// performance here.

import algorithm = std.algorithm;
import conv = std.conv;
import format = std.format;
import math = std.math;
import range = std.range;
import traits = std.traits;

// CRTP base class for all Vector templates.
class BasicVector(ThisT : VectorT!ElemT, size_t SizeV, alias VectorT, ElemT)
if (traits.isNumeric!ElemT) {
protected:

  ElemT[SIZE] _data;

  alias SIZE = SizeV;
  // FloatType is the type returned by Norm() and Angle().  These methods are
  // special because they return floating-point values even when VType is an
  // integer.
  static if (traits.isIntegral!ElemT) {
    alias FloatT = double;
  } else {
    alias FloatT = ElemT;
  }

public:
  this() {
    _data = _data.init;
  }

  this(ThisT v) {
    _data = v.data;
  }

  @property
  ElemT[SIZE] data() {
    return _data;
  }

  @property
  ThisT data(RangeT)(RangeT data)
  if (range.isInputRange!RangeT) {
    _data[] = range.take(range.array(data), SIZE);
    return cast(ThisT) this;
  }

  static size_t size() { return SIZE; }

  // Convert from another vector type
  static ThisT from(Elem2T)(in VectorT!Elem2T b) {
    return new ThisT(algorithm.map!(elem => cast(ElemT) elem)(b._data));
  }

  static if (is(ElemT.nan : ElemT)) {
    // A Vector populated with all NaN values.
    static ThisT NaN() {
      ThisT retvalue = new ThisT();
      retvalue._data[] = ElemT.nan;
    }
  }

  void clear() { _data = _data.init; }

  // Implement index assignment of the form: v[i1] = value
  void opIndexAssign(ElemT value, size_t i1)
  in {
    assert(i1 >= 0);
    assert(i1 < SIZE);
  } body {
    _data[i1] = value;
  }

  // Implement indexing of the form: v[i1]
  ref ElemT opIndex(size_t i1)
  in {
    assert(i1 >= 0);
    assert(i1 < SIZE);
  } body {
    return _data[i1];
  }

  // Support the == and != operators.
  // TODO(user): Relationals should be nonmembers.
  override
  bool opEquals(Object o) {
    if (auto v = cast(typeof(this)) o) {
      return _data == v._data;
    }
    return false;
  }

  // Support the <, <=, >, and >= operators.
  override
  int opCmp(Object o) {
    if (auto v = cast(typeof(this)) o) {
      foreach (size_t i; 0 .. SIZE) {
        if (_data[i] == v._data[i]) {
          continue;
        }
        return _data[i] > v._data[i] ? 1 : -1;
      }
      return 0;
    }
    return false;
  }

  // Support the +=, -=, *=, and /= operators.
  ThisT opOpAssign(string op)(ThisT v) {
    static if (op == "+=") {
      _data += v._data;
    } else static if (op == "-=") {
      _data -= v._data;
    } else static if (op == "*=") {
      _data *= v._data;
    } else static if (op == "/=") {
      _data /= v._data;
    } else {
      static assert(false, "Operator " ~ op ~ " not implemented!");
    }
    return cast(ThisT) this;
  }

  // Support the +, -, *, and / operators.
  ThisT opBinary(string op)(ThisT v) const {
    ThisT retval = new ThisT(this);
    static if (op == "+") {
      retval += v;
    } else static if (op == "-") {
      retval -= v;
    } else static if (op == "*") {
      retval *= v;
    } else static if (op == "/") {
      retval /= v;
    } else {
      static assert(0, "Operator " ~ op ~ " not implemented");
    }
    return retval;
  }

  // Support negation.
  ThisT opUnary(string op)() const {
    ThisT retval = new ThisT(this);
    static if (op == "-") {
      retval._data = -retval._data;
    } else {
      static assert(0, "Unary operator " ~ op ~ " not implemented.");
    }
    return retval;
  }

  // Element-wise max.  {max(a[0],b[0]), max(a[1],b[1]), ...}
  ThisT max(in ThisT a, in ThisT b) {
    return new ThisT().data(
        algorithm.map!(vals => algorithm.max(vals[0], vals[1]))(
                range.zip(a._data[], b._data[])));
  }

  // Element-wise min.  {min(a[0],b[0]), min(a[1],b[1]), ...}
  ThisT min(in ThisT a, in ThisT b) {
    return new ThisT().data(
        algorithm.map!(vals => algorithm.min(vals[0], vals[1]))(
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

  // Euclidean norm. For integer T, correct only if Norm2 does not overflow.
  FloatT norm() const {
    return math.sqrt(cast(FloatT) norm2());
  }

  // Normalized vector if the norm is nonzero. Not for integer types.
  static if (traits.isFloatingPoint!ElemT) {
    ThisT normalize() {
      ElemT n = norm();
      if (n != cast(ElemT) 0.0) {
        n = cast(ElemT) 1.0 / n;
      }
      _data[] = _data[] * n;
      return cast(ThisT) this;
    }

    // Compose a vector from the sqrt of each component.
    ThisT sqrt() const {
      return new ThisT().data(algorithm.map!(math.sqrt)(_data[]));
    }

    // Take the floor of each component.
    ThisT floor() const {
      return new ThisT().data(algorithm.map!(math.floor)(_data[]));
    }

    // Take the ceil of each component.
    ThisT ceil() const {
      return new ThisT().data(algorithm.map!(math.ceil)(_data[]));
    }

    // Round of each component. Formerly FRound().
    ThisT round() const {
      return new ThisT().data(algorithm.map!(math.rint)(_data[]));
    }

    // Round of each component and return an integer vector. Formerly IRound().
    VectorT!int toIntVector() const {
      return new VectorT!int().data(algorithm.map!(val => cast(int) math.lrint(val))(_data[]));
    }

    // True if any of the components is not a number.
    bool isNaN() const {
      return algorithm.any!(math.isNaN)(_data[]);
    }
  }

  // Compute the absolute value of each component.
  ThisT abs() const {
    return new ThisT().data(algorithm.map!(math.abs)(_data[]));
  }

  // Approximately equal.
  bool aequal(in Vector2!ElemT vb, FloatT margin) const {
    return algorithm.all!(vals => math.abs(vals[0] - vals[1]) < margin)(
        range.zip(_data[], vb._data[]));
  }

  override
  string toString() const {
    string val = "[";
    string sep = "";
    foreach (int i; 0 .. SIZE) {
      val ~= sep ~ conv.to!string(_data[i]);
      sep = ", ";
    }
    return val ~ "]";
  }

}

// ======================================================================
class Vector2(ElemT) : BasicVector!(Vector2!ElemT, 2) {
 private:
  alias BaseT = BasicVector!(Vector2!ElemT, 2);

 public:
  alias FloatT = BaseT.FloatT;
  alias SIZE = BaseT.SIZE;

  this() {
    super();
  }

  this(ElemT x, ElemT y) {
    _data[0] = x;
    _data[1] = y;
  }

  this(in Vector3!ElemT b) {
    this(b.x(), b.y());
  }

  this(in Vector4!ElemT b) {
    this(b.x(), b.y());
  }

  @property
  void x(ElemT v) {
    _data[0] = v;
  }

  @property
  void y(ElemT v) {
    _data[1] = v;
  }

  @property
  ElemT x() const {
    return _data[0];
  }

  @property
  ElemT y() const {
    return _data[1];
  }

  void set(ElemT x, ElemT y) {
    this.x = x;
    this.y = y;
  }

  static if (traits.isSigned!ElemT) {
    // Cross product.  Be aware that if T is an integer type, the high bits
    // of the result are silently discarded.
    ElemT crossProd(in Vector2!ElemT vb) const {
      return x * vb.y - y * vb.x;
    }

    // return the angle between "this" and v in radians
    FloatT angle(in Vector2!ElemT v) const {
      return math.atan2(cast(FloatT) crossProd(v), cast(FloatT) dotProd(v));
    }

    // return a vector orthogonal to the current one
    // with the same norm and counterclockwise to it
    Vector2!ElemT ortho() const {
      return new Vector2!ElemT(-y, x);
    }
  }

}


class Vector3(ElemT) : BasicVector!(Vector3!ElemT, 3) {
 private:
  alias BaseT = BasicVector!(Vector3!ElemT, 3);

 public:
  alias FloatT = BaseT.FloatT;

  this() {
    super();
  }

  this(ElemT x, ElemT y, ElemT z) {
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
  }

  this(in Vector2!ElemT b, ElemT z) {
    this(b.x, b.y, z);
  }

  this(in Vector4!ElemT b) {
    this(b.x, b.y, b.z);
  }

  @property
  void x(ElemT v) {
    _data[0] = v;
  }

  @property
  void y(ElemT v) {
    _data[1] = v;
  }

  @property
  void z(ElemT v) {
    _data[2] = v;
  }

  @property
  ElemT x() const {
    return _data[0];
  }

  @property
  ElemT y() const {
    return _data[1];
  }

  @property
  ElemT z() const {
    return _data[2];
  }

  void set(ElemT x, ElemT y, ElemT z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  // Cross product.  Be aware that if VType is an integer type, the high bits
  // of the result are silently discarded.
  Vector3!ElemT crossProd(in Vector3!ElemT vb) const {
    return new Vector3(
        cast(ElemT)(y * vb.z - z * vb.y),
        cast(ElemT)(z * vb.x - x * vb.z),
        cast(ElemT)(x * vb.y - y * vb.x));
  }

  static if (traits.isFloatingPoint!ElemT) {
    // Returns a unit vector orthogonal to this one.
    Vector3!ElemT ortho() const {
      int k = largestAbsComponent() - 1;
      if (k < 0) {
        k = 2;
      }
      Vector3!ElemT temp = new Vector3!ElemT();
      temp[k] = cast(ElemT) 1;
      return crossProd(temp).normalize();
    }
  }

  // return the angle between two vectors in radians
  FloatT angle(in Vector3!ElemT va) const {
    return math.atan2(crossProd(va).norm(), cast(FloatT) dotProd(va));
  }

  // return the index of the largest component (fabs)
  int largestAbsComponent() const {
    Vector3!ElemT temp = abs();
    return temp[0] > temp[1]
        ? temp[0] > temp[2] ? 0 : 2
        : temp[1] > temp[2] ? 1 : 2;
  }

  // return the index of the smallest, median ,largest component of the vector
  Vector3!int componentOrder() const {
    auto temp = new Vector3!int(0, 1, 2);
    if (_data[temp[0]] > _data[temp[1]]) {
      algorithm.swap(temp[0], temp[1]);
    }
    if (_data[temp[1]] > _data[temp[2]]) {
      algorithm.swap(temp[1], temp[2]);
    }
    if (_data[temp[0]] > _data[temp[1]]) {
      algorithm.swap(temp[0], temp[1]);
    }
    return temp;
  }

}

class Vector4(ElemT) : BasicVector!(Vector4!ElemT, 4) {
 private:
  alias BaseT = BasicVector!(Vector4!ElemT, 4);

 public:
  alias FloatT = BaseT.FloatT;

  this() {
    super();
  }

  this(ElemT x, ElemT y, ElemT z, ElemT w) {
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
    _data[3] = w;
  }

  this(in Vector2!ElemT b, ElemT z, ElemT w) {
    this(b.x, b.y, z, w);
  }

  this(in Vector2!ElemT a, in Vector2!ElemT b) {
    this(a.x, a.y, b.x, b.y);
  }

  this(in Vector3!ElemT b, ElemT w) {
    this(b.x, b.y, b.z, w);
  }

  @property
  void x(ElemT v) {
    _data[0] = v;
  }
  @property
  void y(ElemT v) {
    _data[1] = v;
  }

  @property
  void z(ElemT v) {
    _data[2] = v;
  }

  @property
  void w(ElemT v) {
    _data[3] = v;
  }

  @property
  ElemT x() const {
    return _data[0];
  }

  @property
  ElemT y() const {
    return _data[1];
  }

  @property
  ElemT z() const {
    return _data[2];
  }

  @property
  ElemT w() const {
    return _data[3];
  }

  void set(ElemT x, ElemT y, ElemT z, ElemT w) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
  }

};

alias Vector2_b = Vector2!ubyte;
alias Vector2_i = Vector2!int;
alias Vector2_f = Vector2!float;
alias Vector2_d = Vector2!double;

alias Vector3_b = Vector3!ubyte;
alias Vector3_i = Vector3!int;
alias Vector3_f = Vector3!float;
alias Vector3_d = Vector3!double;

alias Vector4_b = Vector4!ubyte;
alias Vector4_i = Vector4!int;
alias Vector4_f = Vector4!float;
alias Vector4_d = Vector4!double;
