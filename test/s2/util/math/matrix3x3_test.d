module s2.util.math.matrix3x3_test;

import s2.util.math.matrix3x3;
import s2.util.math.vector;
import fluent.asserts;

@("construction")
unittest {
  Matrix3x3_i mi = Matrix3x3_i(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Assert.equal(mi[1, 1], 5);
  Matrix3x3_i mi2 = Matrix3x3_i(mi);
  Assert.equal(mi2[1, 1], 5);
  Matrix3x3_f mf = Matrix3x3_f(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9);
  Assert.equal(mf[1, 1], 0.5f);
}

@("opCast")
unittest {
  auto mi = cast(Matrix3x3_i) Matrix3x3_f(
      1.2, 1.6, 2.0,
      2.4, 2.8, 3.2,
      3.6, 4.0, 4.4);

  Assert.equal(mi[1, 1], 2);
  Assert.equal(mi[1, 2], 3);
}

@("set")
unittest {
  Matrix3x3_i mi = Matrix3x3_i(-1, 2, -3, 4, -5, 6, -7, 8, -9);
  Matrix3x3_i mi2 = Matrix3x3_i(mi);
  mi.set(1, -2, 3, -4, 5, -6, 7, -8, 9);
  Assert.equal(mi[2, 2], 9);
  Assert.equal(mi2[2, 2], -9);
}

@("opEquals")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(-1, 2, -3, 4, -5, 6, -7, 8, -9);
  Matrix3x3_i mi2 = Matrix3x3_i(-1, 2, -3, 4, -5, 6, -7, 8, -9);
  Matrix3x3_i mi3 = Matrix3x3_i(-1, 2, -3, -4, -5, 6, -7, 8, -9);

  Assert.equal(mi1, mi2);
  Assert.notEqual(mi1, mi3);
}

@("opOpAssign(ThisT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);
  Matrix3x3_i mi2 = Matrix3x3_i(-1, 0, 1, 0, 1, -1, 1, -1, 0);

  // Check return types of the operators.
  Assert.equal(Matrix3x3_i(mi1) += mi2, Matrix3x3_i(0, 2, 4, 2, 4, 0, 4, 0, 2));
  Assert.equal(Matrix3x3_i(mi1) -= mi2, Matrix3x3_i(2, 2, 2, 2, 2, 2, 2, 2, 2));

  // Check actual mutation.
  mi1 += mi2;
  Assert.equal(mi1, Matrix3x3_i(0, 2, 4, 2, 4, 0, 4, 0, 2));
}

@("opOpAssign(\"*\")(ThisT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, -2, 3, 4, -1, -1, 1, -3);
  Matrix3x3_i mi2 = Matrix3x3_i(2, 1, 3, -1, 2, 4, -2, 1, -1);

  Assert.equal(Matrix3x3_i(mi1) *= mi2, Matrix3x3_i(4, 3, 13, 4, 10, 26, 3, -2, 4));
  mi1 *= mi2;
  Assert.equal(mi1, Matrix3x3_i(4, 3, 13, 4, 10, 26, 3, -2, 4));
}

@("opOpAssign(ElemT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);

  // Check return types of the operators.
  Assert.equal(Matrix3x3_i(mi1) += 2, Matrix3x3_i(3, 4, 5, 4, 5, 3, 5, 3, 4));
  Assert.equal(Matrix3x3_i(mi1) -= 2, Matrix3x3_i(-1, 0, 1, 0, 1, -1, 1, -1, 0));

  // Check actual mutation.
  mi1 += 2;
  Assert.equal(mi1, Matrix3x3_i(3, 4, 5, 4, 5, 3, 5, 3, 4));
}

@("opBinary(ThisT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);
  Matrix3x3_i mi2 = Matrix3x3_i(-1, 0, 1, 0, 1, -1, 1, -1, 0);

  Assert.equal(mi1 + mi2, Matrix3x3_i(0, 2, 4, 2, 4, 0, 4, 0, 2));
  Assert.equal(mi1 - mi2, Matrix3x3_i(2, 2, 2, 2, 2, 2, 2, 2, 2));
}

@("opBinary(\"*\")(ThisT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, -2, 3, 4, -1, -1, 1, -3);
  Matrix3x3_i mi2 = Matrix3x3_i(2, 1, 3, -1, 2, 4, -2, 1, -1);

  Assert.equal(mi1 * mi2, Matrix3x3_i(4, 3, 13, 4, 10, 26, 3, -2, 4));
}

@("opBinary(ElemT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);

  Assert.equal(mi1 + 2, Matrix3x3_i(3, 4, 5, 4, 5, 3, 5, 3, 4));
  Assert.equal(mi1 - 2, Matrix3x3_i(-1, 0, 1, 0, 1, -1, 1, -1, 0));

  Assert.equal(2 + mi1, Matrix3x3_i(3, 4, 5, 4, 5, 3, 5, 3, 4));
  Assert.equal(2 - mi1, Matrix3x3_i(1, 0, -1, 0, -1, 1, -1, 1, 0));

  Assert.equal(mi1 * 2, Matrix3x3_i(2, 4, 6, 4, 6, 2, 6, 2, 4));
  Assert.equal(mi1 / 2, Matrix3x3_i(0, 1, 1, 1, 1, 0, 1, 0, 1));

  Assert.equal(2 * mi1, Matrix3x3_i(2, 4, 6, 4, 6, 2, 6, 2, 4));
  Assert.equal(2 / mi1, Matrix3x3_i(2, 1, 0, 1, 0, 2, 0, 2, 1));
}

@("opBinary(VectorT)")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);

  Assert.equal(mi1 * Vector3_i(1, 2, 3), Vector3_i(14, 11, 11));
}

@("opUnary")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 2, 3, 1, 3, 1, 2);
  Assert.equal(-mi1, Matrix3x3_i(-1, -2, -3, -2, -3, -1, -3, -1, -2));
}

@("det")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(6, 1, 1, 4, -2, 5, 2, 8, 7);
  Assert.equal(mi1.det(), -306);
}

@("trace")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(6, 1, 1, 4, -2, 5, 2, 8, 7);
  Assert.equal(mi1.trace(), 11);
}

@("data")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(6, 1, 1, 4, -2, 5, 2, 8, 7);
  Assert.equal(mi1.data[0][2], 1);
  Assert.equal(mi1.data[1][1], -2);
  Assert.equal(mi1.data[2][0], 2);
}

@("opIndex")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(6, 1, 1, 4, -2, 5, 2, 8, 7);
  Assert.equal(mi1[0, 2], 1);
  Assert.equal(mi1[1, 1], -2);
  Assert.equal(mi1[2, 0], 2);

  Assert.equal(mi1[2], 1);
  Assert.equal(mi1[4], -2);
  Assert.equal(mi1[6], 2);
}

@("opIndexAssign")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(6, 1, 1, 4, -2, 5, 2, 8, 7);
  mi1[0, 2] = 3;
  mi1[1, 1] = 2;
  mi1[2, 0] = 1;
  Assert.equal(mi1, Matrix3x3_i(6, 1, 3, 4, 2, 5, 1, 8, 7));

  mi1[2] = -1;
  mi1[4] = -2;
  mi1[6] = -3;
  Assert.equal(mi1, Matrix3x3_i(6, 1, -1, 4, -2, 5, -3, 8, 7));
}

@("transpose")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(
      6,  1, 1,
      4, -2, 5,
      2,  8, 7);
  Assert.equal(mi1.transpose(), Matrix3x3_i(
      6,  4, 2,
      1, -2, 8,
      1,  5, 7));
}

@("cofactorMatrixTransposed")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(
      1, 2, 3,
      0, 4, 5,
      1, 0, 6);
  // The cofactor matrix is:
  //   24,  5, -4
  //  -12,  3,  2
  //   -2, -5,  4
  Assert.equal(mi1.cofactorMatrixTransposed(), Matrix3x3_i(
          24, -12, -2,
           5,   3, -5,
          -4,   2,  4));
}

@("inverse")
unittest {
  Matrix3x3_f mf1 = Matrix3x3_f(
      3.0, 0.0, 2.0,
      2.0, 0.0, -2.0,
      0.0, 1.0, 1.0);
  Assert.equal(mf1.inverse(), Matrix3x3_f(
      0.2, 0.2, 0.0,
      -0.2, 0.3, 1.0,
      0.2, -0.3, 0));
}

@("row")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6);
  Assert.equal(mi1.row(0), Vector3_i(1, 2, 3));
  Assert.equal(mi1.row(1), Vector3_i(0, 4, 5));
  Assert.equal(mi1.row(2), Vector3_i(1, 0, 6));
}

@("col")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6);
  Assert.equal(mi1.col(0), Vector3_i(1, 0, 1));
  Assert.equal(mi1.col(1), Vector3_i(2, 4, 0));
  Assert.equal(mi1.col(2), Vector3_i(3, 5, 6));
}

@("fromRows")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i.fromRows(
      Vector3_i(1, 2, 3),
      Vector3_i(0, 4, 5),
      Vector3_i(1, 0, 6));
  Assert.equal(mi1, Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6));
}

@("fromCols")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i.fromCols(
      Vector3_i(1, 0, 1),
      Vector3_i(2, 4, 0),
      Vector3_i(3, 5, 6));
  Assert.equal(mi1, Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6));
}

@("setRow")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6);
  mi1.setRow(2, Vector3_i(5, 4, 3));
  Assert.equal(mi1, Matrix3x3_i(1, 2, 3, 0, 4, 5, 5, 4, 3));
}

@("setCol")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i(1, 2, 3, 0, 4, 5, 1, 0, 6);
  mi1.setCol(2, Vector3_i(5, 4, 3));
  Assert.equal(mi1, Matrix3x3_i(1, 2, 5, 0, 4, 4, 1, 0, 3));
}

@("identity")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i.identity();
  Assert.equal(mi1, Matrix3x3_i(1, 0, 0, 0, 1, 0, 0, 0, 1));
}

@("zero")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i.zero();
  Assert.equal(mi1, Matrix3x3_i(0, 0, 0, 0, 0, 0, 0, 0, 0));
}

@("diagonal")
unittest {
  Matrix3x3_i mi1 = Matrix3x3_i.diagonal(Vector3_i(3, 4, 2));
  Assert.equal(mi1, Matrix3x3_i(3, 0, 0, 0, 4, 0, 0, 0, 2));
}

@("sym3")
unittest {
  Vector3_i v = Vector3_i(-1, 2, 3);
  Assert.equal(Matrix3x3_i.sym3(v), Matrix3x3_i(1, -2, -3, -2, 4, 6, -3, 6, 9));
}
