// Copyright 2001 Google Inc. All Rights Reserved.
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

module s2.util.math.mathutil;

import math = std.math;

bool realRootsForCubic(
    in real a,
    in real b,
    in real c,
    out real r1,
    out real r2,
    out real r3) {
  // According to Numerical Recipes (pp. 184-5), what
  // follows is an arrangement of computations to
  // compute the roots of a cubic that minimizes
  // roundoff error (as pointed out by A.J. Glassman).

  const real a_squared = a * a;
  const real a_third = a / 3.0;
  const real b_tripled = 3.0 * b;
  const real Q = (a_squared - b_tripled) / 9.0;
  const real R = (2.0 * a_squared * a - 3.0 * a * b_tripled + 27.0 * c) / 54.0;

  const real R_squared = R * R;
  const real Q_cubed = Q * Q * Q;

  if (R_squared < Q_cubed) {
    const real root_Q = math.sqrt(Q);
    const real two_pi_third = 2.0 * math.PI / 3.0;
    const real theta_third = math.acos(R / math.sqrt(Q_cubed)) / 3.0;
    const real minus_two_root_Q = -2.0 * root_Q;

    r1 = minus_two_root_Q * math.cos(theta_third) - a_third;
    r2 = minus_two_root_Q * math.cos(theta_third + two_pi_third) - a_third;
    r3 = minus_two_root_Q * math.cos(theta_third - two_pi_third) - a_third;

    return true;
  }

  const real A = -math.signbit(R)
      * math.pow(math.abs(R) + math.sqrt(R_squared - Q_cubed), 1.0 / 3.0L);

  if (A != 0.0) {  // in which case, B from NR is zero
    r1 = A + Q / A - a_third;
    return false;
  }

  r1 = r2 = r3 = -a_third;
  return true;
}
