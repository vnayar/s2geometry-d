module s2.s1chordangle_test;

import fluent.asserts;
import math = std.math;
import s2.s1angle;
import s2.s1chordangle;
import s2.s2edgedistances;
import s2.s2point;
import s2.s2predicates;
import s2.s2testing;
import std.stdio;


enum double EPSILON = 0.0001;

@("DefaultConstructor")
unittest {
  // Check that the default constructor returns an angle of 0.
  S1ChordAngle a;
  Assert.equal(S1ChordAngle.zero(), a);
}

@("TwoPointConstructor")
unittest {
  foreach (int iter; 0 .. 100) {
    S2Point x, y, z;
    S2Testing.getRandomFrame(x, y, z);
    Assert.equal(S1Angle.zero(), S1ChordAngle(z, z).toS1Angle());
    Assert.approximately(S1ChordAngle(-z, z).radians(), math.PI, 1e-7);
    Assert.approximately(S1ChordAngle(x, z).radians(), math.PI_2, EPSILON);
    S2Point w = (y + z).normalize();
    Assert.approximately(S1ChordAngle(w, z).radians(), math.PI_4, EPSILON);
  }
}

@("FromLength2")
unittest {
  Assert.approximately(S1ChordAngle.fromLength2(0).degrees(), 0, EPSILON);
  Assert.approximately(S1ChordAngle.fromLength2(1).degrees(), 60, EPSILON);
  Assert.approximately(S1ChordAngle.fromLength2(2).degrees(), 90, EPSILON);
  Assert.approximately(S1ChordAngle.fromLength2(4).degrees(), 180, EPSILON);
  Assert.approximately(S1ChordAngle.fromLength2(5).degrees(), 180, EPSILON);
}

@("Zero")
unittest {
  Assert.equal(S1Angle.zero(), S1ChordAngle.zero().toS1Angle());
}

@("Right")
unittest {
  Assert.approximately(S1ChordAngle.right().degrees(), 90, EPSILON);
}

@("Straight")
unittest {
  Assert.equal(S1Angle.fromDegrees(180), S1ChordAngle.straight().toS1Angle());
}

@("Infinity")
unittest {
  Assert.lessThan(S1ChordAngle.straight(), S1ChordAngle.infinity());
  Assert.equal(S1ChordAngle.infinity(), S1ChordAngle.infinity());
  Assert.equal(S1Angle.infinity(), S1ChordAngle.infinity().toS1Angle());
}

@("Negative")
unittest {
  Assert.lessThan(S1ChordAngle.negative(), S1ChordAngle.zero());
  Assert.equal(S1ChordAngle.negative(), S1ChordAngle.negative());
  Assert.lessThan(S1ChordAngle.negative().toS1Angle(), S1Angle.zero());
}

@("Predicates")
unittest {
  Assert.equal(S1ChordAngle.zero().isZero(), true);
  Assert.equal(S1ChordAngle.zero().isNegative(), false);
  Assert.equal(S1ChordAngle.zero().isSpecial(), false);
  Assert.equal(S1ChordAngle.straight().isSpecial(), false);
  Assert.equal(S1ChordAngle.negative().isNegative(), true);
  Assert.equal(S1ChordAngle.negative().isSpecial(), true);
  Assert.equal(S1ChordAngle.infinity().isInfinity(), true);
  Assert.equal(S1ChordAngle.infinity().isSpecial(), true);
}

@("ToFromS1Angle")
unittest {
  Assert.equal(S1ChordAngle(S1Angle.zero()).radians(), 0);
  Assert.equal(S1ChordAngle(S1Angle.fromRadians(math.PI)).length2(), 4);
  Assert.equal(S1ChordAngle(S1Angle.fromRadians(math.PI)).radians(), math.PI);
  Assert.equal(S1Angle.infinity(), S1ChordAngle(S1Angle.infinity()).toS1Angle());
  Assert.lessThan(S1ChordAngle(S1Angle.fromRadians(-1)).radians(), 0);
  Assert.equal(1.0, S1ChordAngle(S1Angle.fromRadians(1.0)).radians());
}

@("Successor")
unittest {
  Assert.equal(S1ChordAngle.zero(), S1ChordAngle.negative().successor());
  Assert.equal(S1ChordAngle.infinity(), S1ChordAngle.straight().successor());
  Assert.equal(S1ChordAngle.infinity(), S1ChordAngle.infinity().successor());
  S1ChordAngle x = S1ChordAngle.negative();
  for (int i = 0; i < 10; ++i) {
    Assert.lessThan(x, x.successor());
    x = x.successor();
  }
}

@("Predecessor")
unittest {
  Assert.equal(S1ChordAngle.straight(), S1ChordAngle.infinity().predecessor());
  Assert.equal(S1ChordAngle.negative(), S1ChordAngle.zero().predecessor());
  Assert.equal(S1ChordAngle.negative(), S1ChordAngle.negative().predecessor());
  S1ChordAngle x = S1ChordAngle.infinity();
  for (int i = 0; i < 10; ++i) {
    Assert.greaterThan(x, x.predecessor());
    x = x.predecessor();
  }
}

@("Arithmetic")
unittest {
  S1ChordAngle zero = S1ChordAngle.zero();
  S1ChordAngle degree30 = S1ChordAngle.fromDegrees(30);
  S1ChordAngle degree60 = S1ChordAngle.fromDegrees(60);
  S1ChordAngle degree90 = S1ChordAngle.fromDegrees(90);
  S1ChordAngle degree120 = S1ChordAngle.fromDegrees(120);
  S1ChordAngle degree180 = S1ChordAngle.straight();
  Assert.equal((zero + zero).degrees(), 0);
  Assert.equal((zero - zero).degrees(), 0);
  Assert.equal((degree60 - degree60).degrees(), 0);
  Assert.equal((degree180 - degree180).degrees(), 0);
  Assert.equal((zero - degree60).degrees(), 0);
  Assert.equal((degree30 - degree90).degrees(), 0);
  Assert.approximately((degree60 + zero).degrees(), 60, EPSILON);
  Assert.approximately((degree60 - zero).degrees(), 60, EPSILON);
  Assert.approximately((zero + degree60).degrees(), 60, EPSILON);
  Assert.approximately((degree30 + degree60).degrees(), 90, EPSILON);
  Assert.approximately((degree60 + degree30).degrees(), 90, EPSILON);
  Assert.approximately((degree90 - degree30).degrees(), 60, EPSILON);
  Assert.approximately((degree90 - degree60).degrees(), 30, EPSILON);
  Assert.approximately((degree180 + zero).degrees(), 180, EPSILON);
  Assert.approximately((degree180 - zero).degrees(), 180, EPSILON);
  Assert.approximately((degree90 + degree90).degrees(), 180, EPSILON);
  Assert.approximately((degree120 + degree90).degrees(), 180, EPSILON);
  Assert.approximately((degree120 + degree120).degrees(), 180, EPSILON);
  Assert.approximately((degree30 + degree180).degrees(), 180, EPSILON);
  Assert.approximately((degree180 + degree180).degrees(), 180, EPSILON);
}

@("Trigonometry")
unittest {
  immutable int kIters = 20;
  foreach (int iter; 0 .. kIters) {
    double rads = math.PI * iter / kIters;
    S1ChordAngle angle = S1ChordAngle(S1Angle.fromRadians(rads));
    Assert.approximately(math.sin(rads), angle.sin(), 1e-15);
    Assert.approximately(math.cos(rads), angle.cos(), 1e-15);
    // Since the tan(x) is unbounded near Pi/4, we map the result back to an
    // angle before comparing.  (The assertion is that the result is equal to
    // the tangent of a nearby angle.)
    Assert.approximately(math.atan(math.tan(rads)), math.atan(angle.tan()), 1e-15);
  }

  // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
  S1ChordAngle angle90 = S1ChordAngle.fromLength2(2);
  S1ChordAngle angle180 = S1ChordAngle.fromLength2(4);
  Assert.equal(angle90.sin(), 1);
  Assert.equal(angle90.cos(), 0);
  Assert.equal(angle90.tan(), double.infinity);
  Assert.equal(angle180.sin(), 0);
  Assert.equal(angle180.cos(), -1);
  Assert.equal(angle180.tan(), 0);
}

@("PlusError")
unittest {
  Assert.equal(S1ChordAngle.negative(), S1ChordAngle.negative().plusError(5));
  Assert.equal(S1ChordAngle.infinity(), S1ChordAngle.infinity().plusError(-5));
  Assert.equal(S1ChordAngle.straight(), S1ChordAngle.straight().plusError(5));
  Assert.equal(S1ChordAngle.zero(), S1ChordAngle.zero().plusError(-5));
  Assert.equal(S1ChordAngle.fromLength2(1.25), S1ChordAngle.fromLength2(1).plusError(0.25));
  Assert.equal(S1ChordAngle.fromLength2(0.75), S1ChordAngle.fromLength2(1).plusError(-0.25));
}

@("GetS2PointConstructorMaxError")
unittest {
  // Check that the error bound returned by GetS2PointConstructorMaxError() is
  // large enough.
  auto rnd = S2Testing.rnd;
  foreach (iter; 0 .. 100_000) {
    rnd.reset(iter);  // Easier to reproduce a specific case.
    S2Point x = S2Testing.randomPoint();
    S2Point y = S2Testing.randomPoint();
    if (rnd.oneIn(10)) {
      // Occasionally test a point pair that is nearly identical or antipodal.
      S1Angle r = S1Angle.fromRadians(1e-15 * rnd.randDouble());
      y = interpolateAtDistance(r, x, y);
      if (rnd.oneIn(2)) y = -y;
    }
    S1ChordAngle dist = S1ChordAngle(x, y);
    double error = dist.getS2PointConstructorMaxError();
    Assert.notGreaterThan(compareDistance(x, y, dist.plusError(error)), 0);
    Assert.notLessThan(compareDistance(x, y, dist.plusError(-error)), 0);
  }
}
