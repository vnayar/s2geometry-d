module s2.util.math.vector_test;

import fluent.asserts;
import s2.util.math.vector;
import traits = std.traits;
import math = std.math;

/// Construction
unittest {
  // Test a constructor with fixed number of parameters.
  auto v1 = Vector3_i(1, 2, 3);
  v1.data[0].should.equal(1);
  v1.data[1].should.equal(2);
  v1.data[2].should.equal(3);

  // Test a constructor with ranges.
  auto v2 = Vector3_i([1, 2, 3]);
  v2.data[0].should.equal(1);
  v2.data[1].should.equal(2);
  v2.data[2].should.equal(3);

  // Make sure copy constructor works.
  auto v3 = Vector3_i(v1);
  v3.data[0].should.equal(v3[0]);
  v3.data[1].should.equal(v3[1]);
  v3.data[2].should.equal(v3[2]);

  // Make sure the new object has unique references.
  v3.data[0] = 4;
  v2.data[0].should.not.equal(v3.data[0]);
}

/// Member access
unittest {
  auto v1 = Vector2_f([2.3f, 1.5f]);

  traits.hasMember!(Vector2_f, "x").should.equal(true);
  traits.hasMember!(Vector2_f, "y").should.equal(true);
  traits.hasMember!(Vector2_f, "z").should.equal(false);
  traits.hasMember!(Vector2_f, "w").should.equal(false);

  traits.hasMember!(Vector3_d, "x").should.equal(true);
  traits.hasMember!(Vector3_d, "y").should.equal(true);
  traits.hasMember!(Vector3_d, "z").should.equal(true);
  traits.hasMember!(Vector3_d, "w").should.equal(false);

  traits.hasMember!(Vector4_i, "x").should.equal(true);
  traits.hasMember!(Vector4_i, "y").should.equal(true);
  traits.hasMember!(Vector4_i, "z").should.equal(true);
  traits.hasMember!(Vector4_i, "w").should.equal(true);

  // Access as an array through data.
  v1.size().should.equal(2);
  v1.data.length.should.equal(2);
  v1.data.should.equal([2.3f, 1.5f]);

  // Access through named properties.
  v1.x.should.equal(2.3f);
  v1.y.should.equal(1.5f);

  // Access through the index operator.
  v1[0].should.equal(2.3f);
  v1[1].should.equal(1.5f);
}

/// Mutation
unittest {
  auto v1 = Vector2_i(2, 3);
  v1.x = 12;
  v1[1] = 13;

  v1.data.should.equal([12, 13]);
}

/// Comparison operators
unittest {
  auto v1 = Vector2_b(3u, 1u);
  auto v2 = Vector2_b(3u, 6u);
  auto v3 = Vector2_b(4u, 2u);
  auto v4 = Vector2_b(4u, 2u);

  v1.should.be.lessThan(v2);
  v2.should.be.lessThan(v3);
  v3.should.equal(v4);
  v3.should.be.greaterThan(v2);
  v2.should.be.greaterThan(v1);
}

/// Arithmetic operators
unittest {
  Vector2_i(3, 4).should.equal(-Vector2_i(-3, -4));
  (Vector2_i(3, 4) + Vector2_i(-2, 4)).should.equal(Vector2_i(1, 8));
  (Vector2_i(3, 4) - Vector2_i(-2, 4)).should.equal(Vector2_i(5, 0));
  (Vector2_i(3, 4) * Vector2_i(-2, 4)).should.equal(Vector2_i(-6, 16));
  (Vector2_i(3, 4) / Vector2_i(-2, 4)).should.equal(Vector2_i(-1, 1));

  auto v1 = Vector2_i(-5, 3);
  v1 += Vector2_i(6, -2);
  v1.should.equal(Vector2_i(1, 1));

  v1 = Vector2_i(-5, 3);
  v1 -= Vector2_i(-6, 2);
  v1.should.equal(Vector2_i(1, 1));

  v1 = Vector2_i(-5, 3);
  v1 *= Vector2_i(-2, -3);
  v1.should.equal(Vector2_i(10, -9));

  v1 = Vector2_i(-5, 3);
  v1 /= Vector2_i(-2, 2);
  v1.should.equal(Vector2_i(2, 1));

}

/// Other helper functions
unittest {
  auto v = Vector2_i.from(Vector2_f(3.3f, 4.7f));
  v[0].should.equal(3);
  v[1].should.equal(4);

  v.clear();
  v.data.should.equal([0, 0]);

  auto v2 = Vector2_f.nan();
  math.isNaN(v2.x).should.equal(true);
  math.isNaN(v2.y).should.equal(true);
}

/// Min and Max
unittest {
  Vector2_i.max(Vector2_i(4, 7), Vector2_i(7, 1)).should.equal(Vector2_i(7, 7));
  Vector2_i.min(Vector2_i(4, 7), Vector2_i(7, 1)).should.equal(Vector2_i(4, 1));
}

/// dotProduct(), norm2() and norm()
unittest {
  Assert.equal(Vector3_i(2, -3, 5).dotProd(Vector3_i(-1, 2, -3)), -23);
  Assert.equal(Vector3_i(2, -3, 5).norm2(), 38);
  Assert.equal(Vector2_i(3, 4).norm(), 5.0);
}

/// Floating point operations: normalize(), sqrt(), floor(), ceil(), round()
unittest {
  float delta = 0.0001;
  auto vn = Vector3_f(6.7, -2.3, 8.9).normalize();
  Assert.approximately(vn.x, 0.589012411578, delta);
  Assert.approximately(vn.y, -0.202198290542, delta);
  Assert.approximately(vn.z, 0.782419472096, delta);

  auto vs = Vector3_f(6.7, 2.3, 8.9).sqrt();
  Assert.approximately(vs.x, 2.58843582111, delta);
  Assert.approximately(vs.y, 1.51657508881, delta);
  Assert.approximately(vs.z, 2.98328677804, delta);

  auto vf = Vector3_f(6.7, -2.3, 8.9).floor();
  Assert.approximately(vf.x, 6.0, delta);
  Assert.approximately(vf.y, -3.0, delta);
  Assert.approximately(vf.z, 8.0, delta);

  auto vc = Vector3_f(6.7, -2.3, 8.9).ceil();
  Assert.approximately(vc.x, 7.0, delta);
  Assert.approximately(vc.y, -2.0, delta);
  Assert.approximately(vc.z, 9.0, delta);

  auto vr = Vector3_f(6.2, -2.8, 8.9).round();
  Assert.approximately(vr.x, 6.0, delta);
  Assert.approximately(vr.y, -3.0, delta);
  Assert.approximately(vr.z, 9.0, delta);
  Assert.approximately(vr.z, 9.0, delta);
}

/// Floating point operations: toIntVector(), isNaN()
unittest {
  auto vi = Vector3_d(2.6, 9.1, 5.5).toIntVector();
  Assert.equal(is(typeof(vi) == Vector3_i), true);
  Assert.equal(Vector3_d(2.6, 9.1, 5.5).isNaN(), false);
  Assert.equal(Vector3_d(2.6, double.nan, 5.5).isNaN(), true);
}

/// Test abs() and aequal()
unittest {
  Assert.equal(Vector3_f(-4.7, 3.4, -8.2).abs() == Vector3_f(4.7, 3.4, 8.2), true);

  Assert.equal(
      Vector3_f(-4.7, 3.4, -8.2).aequal(Vector3_f(-4.693, 3.402, -8.207), 0.01), true);
}

// Special 2-dimentional functions: crossProd(), angle(), ortho()
unittest {
  double delta = 0.0001;

  // Cross product is the norm^2 for perpendicular vectors.
  Assert.equal(Vector2_d(3.0, -4.0).crossProd(Vector2_d(4.0, 3.0)), 25);
  // Cross product is 0 for colinear vectors.
  Assert.equal(Vector2_d(3.0, -4.0).crossProd(Vector2_d(-3.0, 4.0)), 0);

  Assert.approximately(Vector2_d(3.0, -4.0).angle(Vector2_d(4.0, 3.0)), math.PI/2.0, delta);
  Assert.approximately(Vector2_d(3.0, -4.0).angle(Vector2_d(-3.0, 4.0)), math.PI, delta);

  // Counterclockwise orthoginal vector.
  Assert.equal(Vector2_d(-3.0, 4.0).ortho(), Vector2_d(-4.0, -3.0));
}

// Special 3-dimentional functions: crossProd(), angle(), ortho()
unittest {
  Assert.equal(
      Vector3_d(3.0, -4.0, 2.0).crossProd(Vector3_d(-2.0, -3.0, 4.0)),
      Vector3_d(-10, -16, -17));

  Assert.approximately(Vector3_d(-2, 2, 2).angle(Vector3_d(0, 2, 2)), 0.61548, 0.00001);

  // The orthogonal vector's axis 1 less than the largest magnitude should be zero.
  // E.g. y is largest, so x should be zero.
  Assert.equal(Vector3_d(-2, 3, 2).ortho().x, 0);
}
