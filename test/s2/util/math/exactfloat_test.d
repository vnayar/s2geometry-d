module s2.util.math.exactfloat_test;

import fluent.asserts;
import s2.util.math.exactfloat;
import std.bigint;
import std.stdio;

@("createFromInt")
unittest {
  ExactFloat xf = -3;
  writeln("xf has value: ", xf);
  Assert.equal(xf.toDouble(), -3.0);
}

@("createFromFloat")
unittest {
  ExactFloat xf = -3.0;
  writeln("xf has value: ", xf);
  Assert.equal(xf.toDouble(), -3.0);
}
