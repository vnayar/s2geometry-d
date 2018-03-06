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

version(unittest) {
  import fluent.asserts;
  import std.stdio;
}

// CRTP base class for all Vector templates.
struct Vector(ElemT, size_t SizeV)
if (traits.isNumeric!ElemT && SizeV >= 1) {
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

  ElemT[SizeV] _data;

public:

  this(ThisT v) {
    _data = v.data;
  }

  this(ElemT v) {
    _data[] = v;
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
  ElemT opIndex(size_t i1) const
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
      return ThisT(algorithm.map!(math.sqrt)(_data[]));
    }

    // Take the floor of each component.
    ThisT floor() const {
      return ThisT(algorithm.map!(math.floor)(_data[]));
    }

    // Take the ceil of each component.
    ThisT ceil() const {
      return ThisT(algorithm.map!(math.ceil)(_data[]));
    }

    // Round of each component. Formerly FRound().
    ThisT round() const {
      return ThisT(algorithm.map!(math.rint)(_data[]));
    }

    // Round of each component and return an integer vector. Formerly IRound().
    Vector!(int, SizeV) toIntVector() const {
      return Vector!(int, SizeV)(algorithm.map!(val => cast(int) math.lrint(val))(_data[]));
    }

    // True if any of the components is not a number.
    bool isNaN() const {
      return algorithm.any!(math.isNaN)(_data[]);
    }
  }

  // Compute the absolute value of each component.
  ThisT abs() const {
    return ThisT(algorithm.map!(math.abs)(_data[]));
  }

  // Approximately equal.
  bool aequal(in ThisT vb, FloatT margin) const {
    return algorithm.all!(vals => math.abs(vals[0] - vals[1]) < margin)(
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

  // Function implementations specific to 2-dimentional vectors.
  static if (SizeV == 2 && traits.isSigned!ElemT) {
    // Cross product.  Be aware that if T is an integer type, the high bits
    // of the result are silently discarded.
    ElemT crossProd(in ThisT vb) const {
      return x * vb.y - y * vb.x;
    }

    // return the angle between "this" and v in radians
    FloatT angle(in ThisT v) const {
      return math.atan2(cast(FloatT) crossProd(v), cast(FloatT) dotProd(v));
    }

    // return a vector orthogonal to the current one
    // with the same norm and counterclockwise to it
    ThisT ortho() const {
      return ThisT([-y, x]);
    }
  }

  // Function implementations speficif to 3-dimentional vectors.
  static if (SizeV == 3 && traits.isSigned!ElemT) {
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
    }

    // return the angle between two vectors in radians
    FloatT angle(in ThisT va) const {
      return math.atan2(crossProd(va).norm(), cast(FloatT) dotProd(va));
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

alias Vector2_b = Vector!(ubyte, 2);
alias Vector2_i = Vector!(int, 2);
alias Vector2_f = Vector!(float, 2);
alias Vector2_d = Vector!(double, 2);

alias Vector3_b = Vector!(ubyte, 3);
alias Vector3_i = Vector!(int, 3);
alias Vector3_f = Vector!(float, 3);
alias Vector3_d = Vector!(double, 3);

alias Vector4_b = Vector!(ubyte, 4);
alias Vector4_i = Vector!(int, 4);
alias Vector4_f = Vector!(float, 4);
alias Vector4_d = Vector!(double, 4);

/// Construction
unittest {
  // Test a constructor with ranges.
  auto v1 = Vector3_i([1, 2, 3]);
  assert(v1._data[0] == 1);
  assert(v1._data[1] == 2);
  assert(v1._data[2] == 3);

  // Make sure copy constructor works.
  auto v2 = Vector3_i(v1);
  assert(v1._data[0] == v2._data[0]);
  assert(v1._data[1] == v2._data[1]);
  assert(v1._data[2] == v2._data[2]);

  // Make sure the new object has unique references.
  v2._data[0] = 4;
  assert(v1._data[0] != v2._data[0]);
}

/// Member access
unittest {
  auto v1 = Vector2_f([2.3f, 1.5f]);

  static assert(traits.hasMember!(Vector2_f, "x"));
  static assert(traits.hasMember!(Vector2_f, "y"));
  static assert(!traits.hasMember!(Vector2_f, "z"));
  static assert(!traits.hasMember!(Vector2_f, "w"));

  static assert(traits.hasMember!(Vector3_d, "x"));
  static assert(traits.hasMember!(Vector3_d, "y"));
  static assert(traits.hasMember!(Vector3_d, "z"));
  static assert(!traits.hasMember!(Vector3_d, "w"));

  static assert(traits.hasMember!(Vector4_i, "x"));
  static assert(traits.hasMember!(Vector4_i, "y"));
  static assert(traits.hasMember!(Vector4_i, "z"));
  static assert(traits.hasMember!(Vector4_i, "w"));

  // Access as an array through data.
  assert(v1.size() == 2);
  assert(v1.data.length == 2);
  assert(v1.data == [2.3f, 1.5f]);

  // Access through named properties.
  assert(v1.x == 2.3f);
  assert(v1.y == 1.5f);

  // Access through the index operator.
  assert(v1[0] == 2.3f);
  assert(v1[1] == 1.5f);
}

/// Mutation
unittest {
  auto v1 = Vector2_i([2, 3]);
  v1.x = 12;
  v1[1] = 13;

  assert(v1.data == [12, 13]);
}

/// Comparison operators
unittest {
  auto v1 = Vector2_b([3u, 1u]);
  auto v2 = Vector2_b([3u, 6u]);
  auto v3 = Vector2_b([4u, 2u]);
  auto v4 = Vector2_b([4u, 2u]);

  assert(v1 < v2);
  assert(v2 < v3);
  assert(v3 == v4);
  assert(v3 > v2);
  assert(v2 > v1);
}

/// Arithmetic operators
unittest {
  assert(Vector2_i([3, 4]) == -Vector2_i([-3, -4]));
  assert(Vector2_i([3, 4]) + Vector2_i([-2, 4]) == Vector2_i([1, 8]));
  assert(Vector2_i([3, 4]) - Vector2_i([-2, 4]) == Vector2_i([5, 0]));
  assert(Vector2_i([3, 4]) * Vector2_i([-2, 4]) == Vector2_i([-6, 16]));
  assert(Vector2_i([3, 4]) / Vector2_i([-2, 4]) == Vector2_i([-1, 1]));

  auto v1 = Vector2_i([-5, 3]);
  v1 += Vector2_i([6, -2]);
  assert(v1 == Vector2_i([1, 1]));

  v1 = Vector2_i([-5, 3]);
  v1 -= Vector2_i([-6, 2]);
  assert(v1 == Vector2_i([1, 1]));

  v1 = Vector2_i([-5, 3]);
  v1 *= Vector2_i([-2, -3]);
  assert(v1 == Vector2_i([10, -9]));

  v1 = Vector2_i([-5, 3]);
  v1 /= Vector2_i([-2, 2]);
  assert(v1 == Vector2_i([2, 1]));

}

/// Other helper functions
unittest {
  auto v = Vector2_i.from(Vector2_f([3.3f, 4.7f]));
  assert(v[0] == 3);
  assert(v[1] == 4);

  v.clear();
  assert(v.data == [0, 0]);

  auto v2 = Vector2_f.nan();
  assert(math.isNaN(v2.x));
  assert(math.isNaN(v2.y));
}

/// Min and Max
unittest {
  assert(Vector2_i.max(Vector2_i([4, 7]), Vector2_i([7, 1])) == Vector2_i([7, 7]));
  assert(Vector2_i.min(Vector2_i([4, 7]), Vector2_i([7, 1])) == Vector2_i([4, 1]));
}

/// dotProduct(), norm2() and norm()
unittest {
  Assert.equal(Vector3_i([2, -3, 5]).dotProd(Vector3_i([-1, 2, -3])), -23);
  Assert.equal(Vector3_i([2, -3, 5]).norm2(), 38);
  Assert.equal(Vector2_i([3, 4]).norm(), 5.0);
}

/// Floating point operations: normalize(), sqrt(), floor(), ceil(), round()
unittest {
  float delta = 0.0001;
  auto vn = Vector3_f([6.7, -2.3, 8.9]).normalize();
  Assert.approximately(vn.x, 0.589012411578, delta);
  Assert.approximately(vn.y, -0.202198290542, delta);
  Assert.approximately(vn.z, 0.782419472096, delta);

  auto vs = Vector3_f([6.7, 2.3, 8.9]).sqrt();
  Assert.approximately(vs.x, 2.58843582111, delta);
  Assert.approximately(vs.y, 1.51657508881, delta);
  Assert.approximately(vs.z, 2.98328677804, delta);

  auto vf = Vector3_f([6.7, -2.3, 8.9]).floor();
  Assert.approximately(vf.x, 6.0, delta);
  Assert.approximately(vf.y, -3.0, delta);
  Assert.approximately(vf.z, 8.0, delta);

  auto vc = Vector3_f([6.7, -2.3, 8.9]).ceil();
  Assert.approximately(vc.x, 7.0, delta);
  Assert.approximately(vc.y, -2.0, delta);
  Assert.approximately(vc.z, 9.0, delta);

  auto vr = Vector3_f([6.2, -2.8, 8.9]).round();
  Assert.approximately(vr.x, 6.0, delta);
  Assert.approximately(vr.y, -3.0, delta);
  Assert.approximately(vr.z, 9.0, delta);
  Assert.approximately(vr.z, 9.0, delta);
}

/// Floating point operations: toIntVector(), isNaN()
unittest {
  auto vi = Vector3_d([2.6, 9.1, 5.5]).toIntVector();
  Assert.equal(is(typeof(vi) == Vector3_i), true);
  Assert.equal(Vector3_d([2.6, 9.1, 5.5]).isNaN(), false);
  Assert.equal(Vector3_d([2.6, double.nan, 5.5]).isNaN(), true);
}

/// Test abs() and aequal()
unittest {
  Assert.equal(Vector3_f([-4.7, 3.4, -8.2]).abs() == Vector3_f([4.7, 3.4, 8.2]), true);

  Assert.equal(
      Vector3_f([-4.7, 3.4, -8.2]).aequal(Vector3_f([-4.693, 3.402, -8.207]), 0.01), true);
}

// Special 2-dimentional functions: crossProd(), angle(), ortho()
unittest {
  double delta = 0.0001;

  // Cross product is the norm^2 for perpendicular vectors.
  Assert.equal(Vector2_d([3.0, -4.0]).crossProd(Vector2_d([4.0, 3.0])), 25);
  // Cross product is 0 for colinear vectors.
  Assert.equal(Vector2_d([3.0, -4.0]).crossProd(Vector2_d([-3.0, 4.0])), 0);

  Assert.approximately(Vector2_d([3.0, -4.0]).angle(Vector2_d([4.0, 3.0])), math.PI/2.0, delta);
  Assert.approximately(Vector2_d([3.0, -4.0]).angle(Vector2_d([-3.0, 4.0])), math.PI, delta);

  // Counterclockwise orthoginal vector.
  Assert.equal(Vector2_d([-3.0, 4.0]).ortho() == Vector2_d([-4.0, -3.0]), true);
}

// Special 3-dimentional functions: crossProd(), angle(), ortho()
unittest {
  Assert.equal(
      Vector3_d([3.0, -4.0, 2.0]).crossProd(Vector3_d([-2.0, -3.0, 4.0]))
          == Vector3_d([-10, -16, -17]),
      true);

  Assert.approximately(Vector3_d([-2, 2, 2]).angle(Vector3_d([0, 2, 2])), 0.61548, 0.00001);

  // The orthogonal vector's axis 1 less than the largest magnitude should be zero.
  // E.g. y is largest, so x should be zero.
  Assert.equal(Vector3_d([-2, 3, 2]).ortho().x, 0);
}







