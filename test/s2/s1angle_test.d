module s2.s1angle_test;

// Original author: ericv@google.com (Eric Veach)

import fluent.asserts;
import math = std.math;
import s2.s1angle;
import s2.s2latlng;
import s2.s2point;
import s2.s2testing;

enum double DOUBLE_ERR = 0.0001;
// Run timing tests with this many iterations.
enum int iters = 100_000_000;

@("S1Angle.DefaultConstructor")
unittest {
  // Check that the default constructor returns an angle of 0.
  S1Angle a;
  Assert.equal(a.radians(), 0);
}

@("S1Angle.Infinity")
unittest {
  Assert.lessThan(S1Angle.fromRadians(1e30), S1Angle.infinity());
  Assert.lessThan(-S1Angle.infinity(), S1Angle.zero());
  Assert.equal(S1Angle.infinity(), S1Angle.infinity());
}

@("S1Angle.Zero")
unittest {
  Assert.equal(S1Angle.fromRadians(0), S1Angle.zero());
}

@("S1Angle.PiRadiansExactly180Degrees")
unittest {
  // Check that the conversion between Pi radians and 180 degrees is exact.
  Assert.equal(S1Angle.fromRadians(math.PI).radians(), math.PI);
  Assert.equal(S1Angle.fromRadians(math.PI).degrees(), 180.0);
  Assert.equal(S1Angle.fromDegrees(180).radians(), math.PI);
  Assert.equal(S1Angle.fromDegrees(180).degrees(), 180.0);

  Assert.equal(S1Angle.fromRadians(math.PI_2).degrees(), 90.0);

  // Check negative angles.
  Assert.equal(S1Angle.fromRadians(-math.PI_2).degrees(), -90.0);
  Assert.equal(S1Angle.fromDegrees(-45).radians(), -math.PI_4);
}

@("S1Angle.E5E6E7Representations")
unittest {
  // Check that E5/E6/E7 representations work as expected.
  Assert.approximately(
      S1Angle.fromDegrees(-45).radians(),
      S1Angle.fromE5(-4500000).radians(), DOUBLE_ERR);
  Assert.approximately(
      S1Angle.fromDegrees(-60).radians(),
      S1Angle.fromE6(-60000000).radians(), DOUBLE_ERR);
  Assert.approximately(
      S1Angle.fromDegrees(75).radians(),
      S1Angle.fromE7(750000000).radians(), DOUBLE_ERR);
  Assert.equal(S1Angle.fromDegrees(-172.56123).e5(), -17256123);
  Assert.equal(S1Angle.fromDegrees(12.345678).e6(), 12345678);
  Assert.equal(S1Angle.fromDegrees(-12.3456789).e7(), -123456789);
}

@("S1Angle.E6E7RepresentationsUnsigned")
unittest {
  // Check that unsigned E6/E7 representations work as expected.
  Assert.approximately(
      S1Angle.fromDegrees(60).radians(),
      S1Angle.fromUnsignedE6(60000000u).radians(), DOUBLE_ERR);
  Assert.approximately(
      S1Angle.fromDegrees(-60).radians(),
      S1Angle.fromUnsignedE6(-60000000u).radians(), DOUBLE_ERR);
  Assert.approximately(
      S1Angle.fromDegrees(75).radians(),
      S1Angle.fromUnsignedE7(750000000u).radians(), DOUBLE_ERR);
  Assert.approximately(
      S1Angle.fromDegrees(-75).radians(),
      S1Angle.fromUnsignedE7(-750000000u).radians(), DOUBLE_ERR);
}

@("S1Angle.NormalizeCorrectlyCanonicalizesAngles")
unittest {
  Assert.approximately(S1Angle.fromDegrees(360.0).normalized().degrees(), 0.0, DOUBLE_ERR);
  Assert.approximately(S1Angle.fromDegrees(-180.0).normalized().degrees(), 180.0, DOUBLE_ERR);
  Assert.approximately(S1Angle.fromDegrees(180.0).normalized().degrees(), 180.0, DOUBLE_ERR);
  Assert.approximately(S1Angle.fromDegrees(540.0).normalized().degrees(), 180.0, DOUBLE_ERR);
  Assert.approximately(S1Angle.fromDegrees(-270.0).normalized().degrees(), 90.0, DOUBLE_ERR);
}

@("S1Angle.ArithmeticOperationsOnAngles")
unittest {
  Assert.approximately(S1Angle.fromRadians(-0.3).abs().radians(), 0.3, DOUBLE_ERR);
  Assert.approximately((-S1Angle.fromRadians(0.1)).radians(), -0.1, DOUBLE_ERR);
  Assert.approximately(
      (S1Angle.fromRadians(0.1) + S1Angle.fromRadians(0.3)).radians(), 0.4, DOUBLE_ERR);
  Assert.approximately(
      (S1Angle.fromRadians(0.1) - S1Angle.fromRadians(0.3)).radians(), -0.2, DOUBLE_ERR);
  Assert.approximately((2 * S1Angle.fromRadians(0.3)).radians(), 0.6, DOUBLE_ERR);
  Assert.approximately((S1Angle.fromRadians(0.3) * 2).radians(), 0.6, DOUBLE_ERR);
  Assert.approximately((S1Angle.fromRadians(0.3) / 2).radians(), 0.15, DOUBLE_ERR);
  Assert.approximately(
      (S1Angle.fromRadians(0.3) / S1Angle.fromRadians(0.6)).radians(), 0.5, DOUBLE_ERR);

  S1Angle tmp = S1Angle.fromRadians(1.0);
  tmp += S1Angle.fromRadians(0.5);
  Assert.approximately(tmp.radians(), 1.5, DOUBLE_ERR);
  tmp -= S1Angle.fromRadians(1.0);
  Assert.approximately(tmp.radians(), 0.5, DOUBLE_ERR);
  tmp *= 5;
  Assert.approximately(tmp.radians(), 2.5, DOUBLE_ERR);
  tmp /= 2;
  Assert.approximately(tmp.radians(), 1.25, DOUBLE_ERR);
}

@("S1Angle.Trigonometry")
unittest {
  // Spot check a few angles to ensure that the correct function is called.
  Assert.approximately(cos(S1Angle.fromDegrees(0)), 1, DOUBLE_ERR);
  Assert.approximately(sin(S1Angle.fromDegrees(90)), 1, DOUBLE_ERR);
  Assert.approximately(tan(S1Angle.fromDegrees(45)), 1, DOUBLE_ERR);
}

@("S1Angle.ConstructorsThatMeasureAngles")
unittest {
  Assert.approximately(
      S1Angle(S2Point(1, 0, 0), S2Point(0, 0, 2)).radians(), math.PI_2, DOUBLE_ERR);
  Assert.approximately(S1Angle(S2Point(1, 0, 0), S2Point(1, 0, 0)).radians(), 0.0, DOUBLE_ERR);
  Assert.approximately(
      S1Angle(S2LatLng.fromDegrees(20, 20), S2LatLng.fromDegrees(70, 20)).degrees(),
      50.0, 1e-13);
}

@("S1Angle.TestFormatting")
unittest {
  Assert.equal(S1Angle.fromDegrees(180.0).toString(), "180.0000000");
}

/+
// TODO: After S2Testing.getCpuTime() is implemented.
@("S1Angle.TestPerformance")
unittest {
  // Verify that the conversion to E5/E6/E7 is not much slower than the
  // conversion from E5/E6/E7.  (Float-to-integer conversions can be quite
  // slow on some platforms.)  We only check the times for E6; the times for
  // E5/E7 should be similar.

  // To reduce the impact of loop overhead, we do kOpsPerLoop ops per loop.
  static const int kOpsPerLoop = 8;

  // Time conversion from E6 to radians.
  double rad_sum = 0;
  const double from_e6_start = S2Testing::GetCpuTime();
  for (int i = FLAGS_iters; i > 0; i -= kOpsPerLoop) {
    // We structure both loops so that all the conversions can be done in
    // parallel.  Otherwise on some platforms the optimizer happens to do a
    // much better job of parallelizing one loop than the other.
    double r0 = S1Angle::E6(i-0).radians();
    double r1 = S1Angle::E6(i-1).radians();
    double r2 = S1Angle::E6(i-2).radians();
    double r3 = S1Angle::E6(i-3).radians();
    double r4 = S1Angle::E6(i-4).radians();
    double r5 = S1Angle::E6(i-5).radians();
    double r6 = S1Angle::E6(i-6).radians();
    double r7 = S1Angle::E6(i-7).radians();
    rad_sum += ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7));
  }
  const double from_e6_time = S2Testing::GetCpuTime() - from_e6_start;
  EXPECT_NE(rad_sum, 0);  // Don't let the sum get optimized away.
  LOG(INFO) << "From E6: "
            << (FLAGS_iters / from_e6_time)
            << " values per second";

  // Time conversion from radians to E6.
  const double delta = (2 * M_PI) / (FLAGS_iters - 1);
  double angle = -M_PI;
  int64 e6_sum = 0;
  const double to_e6_start = S2Testing::GetCpuTime();
  for (int i = FLAGS_iters; i > 0; i -= kOpsPerLoop) {
    int64 r0 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r1 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r2 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r3 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r4 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r5 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r6 = S1Angle::Radians(angle).e6(); angle += delta;
    int64 r7 = S1Angle::Radians(angle).e6(); angle += delta;
    e6_sum += ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7));
  }
  const double to_e6_time = S2Testing::GetCpuTime() - to_e6_start;
  EXPECT_NE(e6_sum + angle, 0);  // Don't let them get optimized away.
  LOG(INFO) << "  To E6: "
            << (FLAGS_iters / to_e6_time)
            << " values per second";

  // Make sure that the To/From E6 times are not much different.
  // The difference factor slightly less than 2 on an x86_64.
  EXPECT_LE(from_e6_time / to_e6_time, 3);
  EXPECT_LE(to_e6_time / from_e6_time, 3);
}
+/

// The current implementation guarantees exact conversions between
// Degrees() and E6() when the Degrees() argument is an integer.
@("S1Angle.DegreesVsE6")
unittest {
  foreach (i; 0 .. 180) {
    Assert.equal(S1Angle.fromDegrees(i), S1Angle.fromE6(1000000 * i));
  }
}

// The current implementation guarantees exact conversions between
// Degrees() and E7() when the Degrees() argument is an integer.
@("S1Angle.DegreesVsE7")
unittest {
  foreach (i; 0 .. 180) {
    Assert.equal(S1Angle.fromDegrees(i), S1Angle.fromE7(10000000 * i));
  }
}

// The current implementation guarantees exact conversions between
// E6() and E7() when the E6() argument is an integer.
@("S1Angle.E6VsE7")
unittest {
  S2Testing.rnd.reset(s2RandomSeed);
  foreach (iter; 0 .. 1000) {
    int i = S2Testing.rnd.uniform(180000000);
    Assert.equal(S1Angle.fromE6(i), S1Angle.fromE7(10 * i));
  }
}

// The current implementation guarantees certain exact conversions between
// degrees and radians (see the header file for details).
@("S1Angle.DegreesVsRadians")
unittest {
  for (int k = -8; k <= 8; ++k) {
    Assert.equal(S1Angle.fromDegrees(45 * k), S1Angle.fromRadians(k * S1Angle.PI / 4));
    Assert.equal(S1Angle.fromDegrees(45 * k).degrees(), 45 * k);
  }
  for (int k = 0; k <= 30; ++k) {
    int n = 1 << k;
    Assert.equal(S1Angle.fromDegrees(180.0 / n), S1Angle.fromRadians(S1Angle.PI / n));
    Assert.equal(S1Angle.fromDegrees(60.0 / n), S1Angle.fromRadians(S1Angle.PI / (3.0 * n)));
    Assert.equal(S1Angle.fromDegrees(36.0 / n), S1Angle.fromRadians(S1Angle.PI / (5.0 * n)));
    Assert.equal(S1Angle.fromDegrees(20.0 / n), S1Angle.fromRadians(S1Angle.PI / (9.0 * n)));
    Assert.equal(S1Angle.fromDegrees(4.0 / n), S1Angle.fromRadians(S1Angle.PI / (45.0 * n)));
  }
}
