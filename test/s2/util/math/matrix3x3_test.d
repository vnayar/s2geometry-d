module s2.util.math.matrix3x3_test;

import s2.util.math.matrix3x3;
import fluent.asserts;

@("construction")
unittest {
  Matrix3x3_i mi = Matrix3x3_i(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Matrix3x3_i mi2 = Matrix3x3_i(mi);
  Matrix3x3_f mf = Matrix3x3_f(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9);
  // No assertions, just expect to compile and not throw.
}
