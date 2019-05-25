// Copyright 2003 Google Inc. All Rights Reserved.
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

module s2.util.math.matrix3x3;

import s2.util.math.mathutil;
import s2.util.math.s2const;
import s2.util.math.vector;
import mathutil = s2.util.math.mathutil;
import conv = std.conv;
import math = std.math;
import traits = std.traits;

/**
 * A simple class to handle 3x3 matrices
 * The aim of this class is to be able to manipulate 3x3 matrices
 * and 3D vectors as naturally as possible and make calculations
 * readable.
 * For that reason, the operators +, -, * are overloaded.
 * (Reading a = a + b*2 - c is much easier to read than
 * a = Sub(Add(a, Mul(b,2)),c)   )
 *
 * Please be careful about overflows when using those matrices wth integer types
 * The calculations are carried with ElemT. eg : if you are using ubyte as the
 * base type, all values will be modulo 256.
 * This feature is necessary to use the class in a more general framework with
 * ElemT != plain old data type.
 */
struct Matrix3x3(ElemT)
if (traits.isNumeric!ElemT) {
private:
  ElemT[3][3] _m;

public:
  alias ThisT = Matrix3x3!ElemT;
  alias VectorT = Vector!(ElemT, 3);

  static if (traits.isIntegral!ElemT) {
    alias FloatT = double;
  } else {
    alias FloatT = ElemT;
  }

  this(ThisT v) {
    _m = v._m;
  }

  /// Constructor explicitly setting the values of all the coefficient of the matrix.
  this(
      in ElemT m00, in ElemT m01, in ElemT m02,
      in ElemT m10, in ElemT m11, in ElemT m12,
      in ElemT m20, in ElemT m21, in ElemT m22) {
    _m[0][0] = m00;
    _m[0][1] = m01;
    _m[0][2] = m02;

    _m[1][0] = m10;
    _m[1][1] = m11;
    _m[1][2] = m12;

    _m[2][0] = m20;
    _m[2][1] = m21;
    _m[2][2] = m22;
  }

  /// Casting constructor
  MatrixT opCast(MatrixT : Matrix3x3!ElemT2, ElemT2)() const {
    return MatrixT(
        conv.to!ElemT2(_m[0][0]),
        conv.to!ElemT2(_m[0][1]),
        conv.to!ElemT2(_m[0][2]),
        conv.to!ElemT2(_m[1][0]),
        conv.to!ElemT2(_m[1][1]),
        conv.to!ElemT2(_m[1][2]),
        conv.to!ElemT2(_m[2][0]),
        conv.to!ElemT2(_m[2][1]),
        conv.to!ElemT2(_m[2][2]));
  }

  /// Change the value of all the coefficients of the matrix.
  ThisT set(
      in ElemT m00, in ElemT m01, in ElemT m02,
      in ElemT m10, in ElemT m11, in ElemT m12,
      in ElemT m20, in ElemT m21, in ElemT m22) {
    _m[0][0] = m00;
    _m[0][1] = m01;
    _m[0][2] = m02;

    _m[1][0] = m10;
    _m[1][1] = m11;
    _m[1][2] = m12;

    _m[2][0] = m20;
    _m[2][1] = m21;
    _m[2][2] = m22;

    return this;
  }

  /// Support Matrix operators +=, -=.
  ThisT opOpAssign(string op)(in ThisT mb)
  if (op == "+" || op == "-") {
    static foreach (x; 0..3) {
      static foreach (y; 0..3) {
        mixin("_m[x][y] " ~ op ~ "= mb._m[x][y];");
      }
    }
    return this;
  }

  // Matrix multiplication
  ThisT opOpAssign(string op)(in ThisT mb)
  if (op == "*") {
    this = this * mb;
    return this;
  }

  // Element-wise assignment operators.
  ThisT opOpAssign(string op)(in ElemT k) {
    static foreach (x; 0..3) {
      static foreach (y; 0..3) {
         mixin("_m[x][y] " ~ op ~ "= k;");
      }
    }
    return this;
  }

  // Support Matrix operators +, -
  ThisT opBinary(string op)(in ThisT mb) const
  if (op == "+" || op == "-") {
    ThisT v = ThisT(this);
    return mixin("v " ~ op ~ "= mb");
  }

  // Matrix multiplication.
  ThisT opBinary(string op)(in ThisT mb) const
  if (op == "*") {
    ThisT v = ThisT(
        _m[0][0] * mb._m[0][0] + _m[0][1] * mb._m[1][0] + _m[0][2] * mb._m[2][0],
        _m[0][0] * mb._m[0][1] + _m[0][1] * mb._m[1][1] + _m[0][2] * mb._m[2][1],
        _m[0][0] * mb._m[0][2] + _m[0][1] * mb._m[1][2] + _m[0][2] * mb._m[2][2],

        _m[1][0] * mb._m[0][0] + _m[1][1] * mb._m[1][0] + _m[1][2] * mb._m[2][0],
        _m[1][0] * mb._m[0][1] + _m[1][1] * mb._m[1][1] + _m[1][2] * mb._m[2][1],
        _m[1][0] * mb._m[0][2] + _m[1][1] * mb._m[1][2] + _m[1][2] * mb._m[2][2],

        _m[2][0] * mb._m[0][0] + _m[2][1] * mb._m[1][0] + _m[2][2] * mb._m[2][0],
        _m[2][0] * mb._m[0][1] + _m[2][1] * mb._m[1][1] + _m[2][2] * mb._m[2][1],
        _m[2][0] * mb._m[0][2] + _m[2][1] * mb._m[1][2] + _m[2][2] * mb._m[2][2]);

    return v;
  }

  // Change the sign of all the coefficients in the matrix
  ThisT opUnary(string op)() const
  if (op == "-") {
    return ThisT(
        -_m[0][0], -_m[0][1], -_m[0][2],
        -_m[1][0], -_m[1][1], -_m[1][2],
        -_m[2][0], -_m[2][1], -_m[2][2]);
  }

  // Matrix scalar binary operators.
  ThisT opBinary(string op)(in ElemT k) const {
    ThisT v = ThisT(this);
    return mixin("v " ~ op ~ "= k");
  }

  ThisT opBinaryRight(string op)(in ElemT k) const {
    ThisT v = ThisT();
    static foreach (x; 0..3) {
      static foreach (y; 0..3) {
        mixin("v._m[x][y] = k " ~ op ~ " _m[x][y];");
      }
    }
    return v;
  }

  // Multiplication of a matrix by a vector
  VectorT opBinary(string op)(in VectorT v) const
  if (op == "*") {
    return VectorT(
      _m[0][0] * v[0] + _m[0][1] * v[1] + _m[0][2] * v[2],
      _m[1][0] * v[0] + _m[1][1] * v[1] + _m[1][2] * v[2],
      _m[2][0] * v[0] + _m[2][1] * v[1] + _m[2][2] * v[2]);
  }

  // Return the determinant of the matrix
  ElemT det() const {
    return _m[0][0] * _m[1][1] * _m[2][2]
         + _m[0][1] * _m[1][2] * _m[2][0]
         + _m[0][2] * _m[1][0] * _m[2][1]
         - _m[2][0] * _m[1][1] * _m[0][2]
         - _m[2][1] * _m[1][2] * _m[0][0]
         - _m[2][2] * _m[1][0] * _m[0][1];
  }

  // Return the trace of the matrix
  ElemT trace() const {
    return _m[0][0] + _m[1][1] + _m[2][2];
  }

  // Return a pointer to the data array for interface with other libraries
  // like opencv
  @property
  ref ElemT[3][3] data() {
    return _m;
  }

  // Return matrix element (i,j) with 0<=i<=2 0<=j<=2
  ref inout(ElemT) opIndex(in size_t i, in size_t j) inout
  in {
    assert(i >= 0 && i < 3);
    assert(j >= 0 && j < 3);
  } do {
    return _m[i][j];
  }

  void opIndexAssign(ElemT value, in size_t i, in size_t j)
  in {
    assert(i >= 0 && i < 3);
    assert(j >= 0 && j < 3);
  } do {
    _m[i][j] = value;
  }

  // Return matrix element (i/3,i%3) with 0<=i<=8 (access concatenated rows).
  ElemT opIndex(in size_t i)
  in {
    assert(i >= 0 && i < 9);
  } do {
    return _m[i/3][i%3];
  }

  ElemT opIndexAssign(ElemT value, in size_t i)
  in {
    assert(i >= 0 && i < 9);
  } do {
    return _m[i/3][i%3] = value;
  }

  // Return the transposed matrix
  ThisT transpose() const {
    return ThisT(
        _m[0][0], _m[1][0], _m[2][0],
        _m[0][1], _m[1][1], _m[2][1],
        _m[0][2], _m[1][2], _m[2][2]);
  }

  // Return the transposed of the matrix of the cofactors
  // (Useful for inversion for example)
  ThisT cofactorMatrixTransposed() const {
    return ThisT(
      _m[1][1] * _m[2][2] - _m[2][1] * _m[1][2],
      _m[2][1] * _m[0][2] - _m[0][1] * _m[2][2],
      _m[0][1] * _m[1][2] - _m[1][1] * _m[0][2],

      _m[1][2] * _m[2][0] - _m[2][2] * _m[1][0],
      _m[2][2] * _m[0][0] - _m[0][2] * _m[2][0],
      _m[0][2] * _m[1][0] - _m[1][2] * _m[0][0],

      _m[1][0] * _m[2][1] - _m[2][0] * _m[1][1],
      _m[2][0] * _m[0][1] - _m[0][0] * _m[2][1],
      _m[0][0] * _m[1][1] - _m[1][0] * _m[0][1]);
  }

  // Matrix inversion
  ThisT inverse() const {
    ElemT det = det();
    assert(det != cast(ElemT) 0, "Can't inverse. Determinant = 0.");
    return (cast(ElemT) 1 / det) * cofactorMatrixTransposed();
  }

  // Return the vector 3D at row i
  VectorT row(in int i) const
  in {
    assert(i >= 0 && i < 3);
  } do {
    return VectorT(_m[i][0], _m[i][1], _m[i][2]);
  }

  // Return the vector 3D at col i
  VectorT col(in int i) const
  in {
    assert(i >= 0 && i < 3);
  } do {
    return VectorT(_m[0][i], _m[1][i], _m[2][i]);
  }

  // Create a matrix from 3 row vectors
  static ThisT fromRows(in VectorT v1, in VectorT v2, in VectorT v3) {
    ThisT temp;
    temp.set(
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2]);
    return temp;
  }

  // Create a matrix from 3 column vectors
  static ThisT fromCols(in VectorT v1, in VectorT v2, in VectorT v3) {
    ThisT temp;
    temp.set(
        v1[0], v2[0], v3[0],
        v1[1], v2[1], v3[1],
        v1[2], v2[2], v3[2]);
    return temp;
  }

  // Set the vector in row i to be v1
  void setRow(int i, in VectorT v1)
  in {
    assert(i >= 0 && i < 3);
  } do {
    _m[i][0] = v1[0];
    _m[i][1] = v1[1];
    _m[i][2] = v1[2];
  }

  // Set the vector in column i to be v1
  void setCol(int i, in VectorT v1)
  in {
    assert(i >= 0 && i < 3);
  } do {
    _m[0][i] = v1[0];
    _m[1][i] = v1[1];
    _m[2][i] = v1[2];
  }

  // Return a matrix M close to the original but verifying MtM = I
  // (useful to compensate for errors in a rotation matrix)
  static if (traits.isFloatingPoint!ElemT) {
    ThisT orthogonalize() const {
      VectorT r1, r2, r3;
      r1 = row(0).normalize();
      r2 = (row(2).crossProd(r1)).normalize();
      r3 = (r1.crossProd(r2)).normalize();
      return fromRows(r1, r2, r3);
    }
  }

  // Return the identity matrix
  static ThisT identity() {
    ThisT temp;
    temp.set(conv.to!ElemT(1), conv.to!ElemT(0), conv.to!ElemT(0),
             conv.to!ElemT(0), conv.to!ElemT(1), conv.to!ElemT(0),
             conv.to!ElemT(0), conv.to!ElemT(0), conv.to!ElemT(1));
    return temp;
  }

  // Return a matrix full of zeros
  static ThisT zero() {
    return ThisT();
  }

  // Return a diagonal matrix with the coefficients in v
  static ThisT diagonal(in VectorT v) {
    return ThisT(
        v[0], conv.to!ElemT(0), conv.to!ElemT(0),
        conv.to!ElemT(0), v[1], conv.to!ElemT(0),
        conv.to!ElemT(0), conv.to!ElemT(0), v[2]);
  }

  // Return the matrix vvT
  static ThisT sym3(in VectorT v) {
    return ThisT(
      v[0]*v[0], v[0]*v[1], v[0]*v[2],
      v[1]*v[0], v[1]*v[1], v[1]*v[2],
      v[2]*v[0], v[2]*v[1], v[2]*v[2]);
  }

  // Return a matrix M such that:
  // for each u,  M * u = v.CrossProd(u)
  static ThisT antiSym3(in VectorT v) {
    return ThisT(
        conv.to!ElemT(0),    -v[2],     v[1],
        v[2],     conv.to!ElemT(0),    -v[0],
        -v[1],       v[0],   conv.to!ElemT(0));
  }

  // Returns matrix that rotates |rot| radians around axis rot.
  static if (traits.isFloatingPoint!ElemT) {
    static ThisT rodrigues(in VectorT rot) {
      ThisT R;
      ElemT theta = rot.norm();
      VectorT w = VectorT(rot).normalize();
      ThisT Wv = ThisT.antiSym3(w);
      ThisT I = ThisT.identity();
      ThisT A = ThisT.sym3(w);
      R = (1 - math.cos(theta)) * A + math.sin(theta) * Wv + math.cos(theta) * I;
      return R;
    }
  }

  // Returns v.Transpose() * (*this) * u
  ElemT mulBothSides(in VectorT v, in VectorT u) const {
    return (this * u).dotProd(v);
  }

  // Use the 3x3 matrix as a projective transform for 2d points
  Vector2!ElemT project(in Vector!(ElemT, 2) v) const {
    VectorT temp = Matrix3x3!ElemT(this) * VectorT(v[0], v[1], 1);
    return Vector2!ElemT(temp[0] / temp[2], temp[1] / temp[2]);
  }

  // Return the Frobenius norm of the matrix: sqrt(sum(aij^2))
  ElemT frobeniusNorm() const {
    ElemT sum = conv.to!ElemT(0);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        sum += _m[i][j] * _m[i][j];
      }
    }
    return conv.to!ElemT(math.sqrt(conv.to!FloatT(sum)));
  }

  // Finds the eigen values of the matrix. Return the number of real eigenvalues
  // found
  int eigenValues(out VectorT eig_val) const {
    real r1, r2, r3;  // NOLINT
    // characteristic polynomial is
    // x^3 + (a11*a22+a22*a33+a33*a11)*x^2 - trace(A)*x - det(A)
    ElemT a = -trace();
    ElemT b = _m[0][0]*_m[1][1] + _m[1][1]*_m[2][2] + _m[2][2]*_m[0][0]
            - _m[1][0]*_m[0][1] - _m[2][1]*_m[1][2] - _m[0][2]*_m[2][0];
    ElemT c = -det();
    bool res = mathutil.realRootsForCubic(a, b, c, r1, r2, r3);
    eig_val[0] = conv.to!ElemT(r1);
    if (res) {
      eig_val[1] = conv.to!ElemT(r2);
      eig_val[2] = conv.to!ElemT(r3);
      return 3;
    }
    return 1;
  }

  // Finds the eigen values and associated eigen vectors of a
  // symmetric positive definite 3x3 matrix,eigen values are
  // sorted in decreasing order, eig_val corresponds to the
  // columns of the eig_vec matrix.
  // Note: The routine will only use the lower diagonal part
  // of the matrix, i.e.
  // |  a00,          |
  // |  a10, a11,     |
  // |  a20, a21, a22 |
  static if (traits.isFloatingPoint!ElemT) {
    void symmetricEigenSolver(out VectorT eig_val, out Matrix3x3!ElemT eig_vec) const {
      // Compute characteristic polynomial coefficients
      double c2 = -(_m[0][0] + _m[1][1] + _m[2][2]);
      double c1 = -(_m[1][0] * _m[1][0] - _m[0][0] * _m[1][1]
          - _m[0][0] * _m[2][2] - _m[1][1] * _m[2][2]
          + _m[2][0] * _m[2][0] + _m[2][1] * _m[2][1]);
      double c0 = -(_m[0][0] * _m[1][1] * _m[2][2] - _m[2][0]
          * _m[2][0] * _m[1][1] - _m[1][0] * _m[1][0]
          * _m[2][2] - _m[0][0] * _m[2][1] * _m[2][1]
          + 2 * _m[1][0] * _m[2][0] * _m[2][1]);

      // Root finding
      double q = (c2*c2-3*c1)/9.0;
      double r = (2*c2*c2*c2-9*c2*c1+27*c0)/54.0;
      // Assume R^3 <Q^3 so there are three real roots
      if (q < 0) q = 0;
      double sqrt_q = -2.0 * math.sqrt(q);
      double theta = math.acos(r / math.sqrt(q * q * q));
      double c2_3 = c2 / 3;
      eig_val[0] = conv.to!ElemT(sqrt_q * math.cos(theta / 3.0) - c2_3);
      eig_val[1] = conv.to!ElemT(sqrt_q * math.cos((theta + 2.0 * M_PI)/3.0) - c2_3);
      eig_val[2] = conv.to!ElemT(sqrt_q * math.cos((theta - 2.0 * M_PI)/3.0) - c2_3);

      // Sort eigen value in decreasing order
      Vector3!int d_order = eig_val.componentOrder();
      eig_val = VectorT(eig_val[d_order[2]], eig_val[d_order[1]], eig_val[d_order[0]]);
      // Compute eigenvectors
      for (int i = 0; i < 3; ++i) {
        VectorT r1 , r2 , r3 , e1 , e2 , e3;
        r1[0] = _m[0][0] - eig_val[i];
        r2[0] = _m[1][0];
        r3[0] = _m[2][0];
        r1[1] = _m[1][0];
        r2[1] = _m[1][1] - eig_val[i];
        r3[1] = _m[2][1];
        r1[2] = _m[2][0];
        r2[2] = _m[2][1];
        r3[2] = _m[2][2] - eig_val[i];

        e1 = r1.crossProd(r2);
        e2 = r2.crossProd(r3);
        e3 = r3.crossProd(r1);

        // Make e2 and e3 point in the same direction as e1
        if (e2.dotProd(e1) < 0) e2 = -e2;
        if (e3.dotProd(e1) < 0) e3 = -e3;
        VectorT e = (e1 + e2 + e3).normalize();
        eig_vec.setCol(i, e);
      }
    }
  }

  // Return true is one of the elements of the matrix is NaN
  static if (traits.isFloatingPoint!ElemT) {
    bool isNaN() const {
      for (int i = 0; i < 3; ++i ) {
        for (int j = 0; j < 3; ++j ) {
          if (math.isNaN(_m[i][j]) ) {
            return true;
          }
        }
      }
      return false;
    }
  }

  bool opEquals(in Matrix3x3!ElemT v) const {
    return _m[0][0] == v._m[0][0]
        && _m[0][1] == v._m[0][1]
        && _m[0][2] == v._m[0][2]
        && _m[1][0] == v._m[1][0]
        && _m[1][1] == v._m[1][1]
        && _m[1][2] == v._m[1][2]
        && _m[2][0] == v._m[2][0]
        && _m[2][1] == v._m[2][1]
        && _m[2][2] == v._m[2][2];
  }

  string toString() const {
    import std.format;
    string val = "";
    int i, j;
    for (i = 0; i < 3; i++) {
      if (i ==0) {
        val ~= "[";
      } else {
        val ~= " ";
      }
      for (j = 0; j < 3; j++) {
        val ~= conv.to!string(_m[i][j]) ~ " ";
      }
      if (i == 2) {
        val ~= "]";
      } else {
        val ~= "\n";
      }
    }
    return val;
  }
}

alias Matrix3x3_i = Matrix3x3!int;
alias Matrix3x3_f = Matrix3x3!float;
alias Matrix3x3_d = Matrix3x3!double;
