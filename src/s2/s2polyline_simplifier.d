// Copyright 2016 Google Inc. All Rights Reserved.
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
//

// Original Author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2polyline_simplifier;

import s2.s1chord_angle;
import s2.s1interval;
import s2.s2point;

import std.math;

/**
 * This is a helper class for simplifying polylines.  It allows you to compute
 * a maximal edge that intersects a sequence of discs, and that optionally
 * avoids a different sequence of discs.  The results are conservative in that
 * the edge is guaranteed to intersect or avoid the specified discs using
 * exact arithmetic (see s2predicates.h).
 *
 * Note that S2Builder can also simplify polylines and supports more features
 * (e.g., snapping to S2CellId centers), so it is only recommended to use this
 * class if S2Builder does not meet your needs.
 *
 * Here is a simple example showing how to simplify a polyline into a sequence
 * of edges that stay within "max_error" of the original edges:
 *
 *   vector<S2Point> v = { ... };
 *   S2PolylineSimplifier simplifier;
 *   simplifier.Init(v[0]);
 *   for (int i = 1; i < v.size(); ++i) {
 *     if (!simplifier.Extend(v[i])) {
 *       OutputEdge(simplifier.src(), v[i-1]);
 *       simplifier.Init(v[i-1]);
 *     }
 *     simplifier.TargetDisc(v[i], max_error);
 *   }
 *   OutputEdge(simplifer.src(), v.back());
 *
 * Note that the points targeted by TargetDisc do not need to be the same as
 * the candidate endpoints passed to Extend.  So for example, you could target
 * the original vertices of a polyline, but only consider endpoints that are
 * snapped to E7 coordinates or S2CellId centers.
 *
 * Please be aware that this class works by maintaining a range of acceptable
 * angles (bearings) from the start vertex to the hypothetical destination
 * vertex.  It does not keep track of distances to any of the discs to be
 * targeted or avoided.  Therefore to use this class correctly, constraints
 * should be added in increasing order of distance.  (The actual requirement
 * is slightly weaker than this, which is why it is not enforced, but
 * basically you should only call TargetDisc() and AvoidDisc() with arguments
 * that you want to constrain the immediately following call to Extend().)
 */
class S2PolylineSimplifier {
public:
  this() {}

  /// Starts a new simplified edge at "src".
  void initialize(in S2Point src) {
    _src = src;
    _window = S1Interval.full();

    // Precompute basis vectors for the tangent space at "src".  This is similar
    // to GetFrame() except that we don't normalize the vectors.  As it turns
    // out, the two basis vectors below have the same magnitude (up to the
    // length error in S2Point::Normalize).

    // Find the index of the component whose magnitude is smallest.
    S2Point tmp = src.abs();
    int i = (tmp[0] < tmp[1] ?
        (tmp[0] < tmp[2] ? 0 : 2) : (tmp[1] < tmp[2] ? 1 : 2));

    // We define the "y" basis vector as the cross product of "src" and the
    // basis vector for axis "i".  Let "j" and "k" be the indices of the other
    // two components in cyclic order.
    int j = (i == 2 ? 0 : i + 1), k = (i == 0 ? 2 : i - 1);
    _yDir[i] = 0;
    _yDir[j] = src[k];
    _yDir[k] = -src[j];

    // Compute the cross product of "y_dir" and "src".  We write out the cross
    // product here mainly for documentation purposes; it also happens to save a
    // few multiplies because unfortunately the optimizer does *not* get rid of
    // multiplies by zero (since these multiplies propagate NaN, for example).
    _xDir[i] = src[j] * src[j] + src[k] * src[k];
    _xDir[j] = -src[j] * src[i];
    _xDir[k] = -src[k] * src[i];
  }

  // Returns the source vertex of the output edge.
  //S2Point src() const;

  /**
   * Returns true if the edge (src, dst) satisfies all of the targeting
   * requirements so far.  Returns false if the edge would be longer than
   * 90 degrees (such edges are not supported).
   */
  bool extend(in S2Point dst) const {
    // We limit the maximum edge length to 90 degrees in order to simplify the
    // error bounds.  (The error gets arbitrarily large as the edge length
    // approaches 180 degrees.)
    if (S1ChordAngle(_src, dst) > S1ChordAngle.right()) return false;

    // Otherwise check whether this vertex is in the acceptable angle range.
    return _window.contains(getAngle(dst));
  }

  /// Requires that the output edge must pass through the given disc.
  bool targetDisc(in S2Point p, S1ChordAngle r) {
    // Shrink the target interval by the maximum error from all sources.  This
    // guarantees that the output edge will intersect the given disc.
    double semiwidth = getSemiwidth(p, r, -1 /*round down*/);
    if (semiwidth >= M_PI) {
      // The target disc contains "src", so there is nothing to do.
      return true;
    }
    if (semiwidth < 0) {
      _window = S1Interval.empty();
      return false;
    }
    // Otherwise compute the angle interval corresponding to the target disc and
    // intersect it with the current window.
    double center = getAngle(p);
    S1Interval target = S1Interval.fromPoint(center).expanded(semiwidth);
    _window = _window.intersection(target);
    return !_window.isEmpty();
  }

  /**
   * Requires that the output edge must avoid the given disc.  "disc_on_left"
   * specifies whether the disc must be to the left or right of the edge.
   * (This feature allows the simplified edge to preserve the topology of the
   * original polyline with respect to other nearby points.)
   *
   * If your input is a polyline, you can compute "disc_on_left" as follows.
   * Let the polyline be ABCDE and assume that it already avoids a set of
   * points X_i.  Suppose that you have aleady added ABC to the simplifer, and
   * now want to extend the edge chain to D.  First find the X_i that are near
   * the edge CD, then discard the ones such that AX_i <= AC or AX_i >= AD
   * (since these points have either already been considered or aren't
   * relevant yet).  Now X_i is to the left of the polyline if and only if
   * s2pred::OrderedCCW(A, D, X, C) (in other words, if X_i is to the left of
   * the angle wedge ACD).
   */
  bool avoidDisc(in S2Point p, S1ChordAngle r, bool disc_on_left) {
    // Expand the interval by the maximum error from all sources.  This
    // guarantees that the final output edge will avoid the given disc.
    double semiwidth = getSemiwidth(p, r, 1 /*round up*/);
    if (semiwidth >= M_PI) {
      // The avoidance disc contains "src", so it is impossible to avoid.
      _window = S1Interval.empty();
      return false;
    }
    double center = getAngle(p);
    double opposite = (center > 0) ? center - M_PI : center + M_PI;
    S1Interval target =
        (disc_on_left ? S1Interval(opposite, center) : S1Interval(center, opposite));
    _window = _window.intersection(target.expanded(-semiwidth));
    return !_window.isEmpty();
  }

private:
  double getAngle(in S2Point p) const {
    return atan2(p.dotProd(_yDir), p.dotProd(_xDir));
  }
  double getSemiwidth(in S2Point p, S1ChordAngle r, int round_direction) const {
    enum double DBL_ERR = 0.5 * double.epsilon;

    // Using spherical trigonometry,
    //
    //   sin(semiwidth) = sin(r) / sin(a)
    //
    // where "a" is the angle between "src" and "p".  Rather than measuring
    // these angles, instead we measure the squared chord lengths through the
    // interior of the sphere (i.e., Cartersian distance).  Letting "r2" be the
    // squared chord distance corresponding to "r", and "a2" be the squared
    // chord distance corresponding to "a", we use the relationships
    //
    //    sin^2(r) = r2 (1 - r2 / 4)
    //    sin^2(a) = d2 (1 - d2 / 4)
    //
    // which follow from the fact that r2 = (2 * sin(r / 2)) ^ 2, etc.

    // "a2" has a relative error up to 5 * DBL_ERR, plus an absolute error of up
    // to 64 * DBL_ERR^2 (because "src" and "p" may differ from unit length by
    // up to 4 * DBL_ERR).  We can correct for the relative error later, but for
    // the absolute error we use "round_direction" to account for it now.
    double r2 = r.length2();
    double a2 = S1ChordAngle(_src, p).length2();
    a2 -= 64 * DBL_ERR * DBL_ERR * round_direction;
    if (a2 <= r2) return M_PI;  // The given disc contains "src".

    double sin2_r = r2 * (1 - 0.25 * r2);
    double sin2_a = a2 * (1 - 0.25 * a2);
    double semiwidth = asin(sqrt(sin2_r / sin2_a));

    // We compute bounds on the errors from all sources:
    //
    //   - The call to GetSemiwidth (this call).
    //   - The call to GetAngle that computes the center of the interval.
    //   - The call to GetAngle in Extend that tests whether a given point
    //     is an acceptable destination vertex.
    //
    // Summary of the errors in GetAngle:
    //
    // y_dir_ has no error.
    //
    // x_dir_ has a relative error of DBL_ERR in two components, a relative
    // error of 2 * DBL_ERR in the other component, plus an overall relative
    // length error of 4 * DBL_ERR (compared to y_dir_) because "src" is assumed
    // to be normalized only to within the tolerances of S2Point::Normalize().
    //
    // p.DotProd(y_dir_) has a relative error of 1.5 * DBL_ERR and an
    // absolute error of 1.5 * DBL_ERR * y_dir_.Norm().
    //
    // p.DotProd(x_dir_) has a relative error of 5.5 * DBL_ERR and an absolute
    // error of 3.5 * DBL_ERR * y_dir_.Norm() (noting that x_dir_ and y_dir_
    // have the same length to within a relative error of 4 * DBL_ERR).
    //
    // It's possible to show by taking derivatives that these errors can affect
    // the angle atan2(y, x) by up 7.093 * DBL_ERR radians.  Rounding up and
    // including the call to atan2 gives a final error bound of 10 * DBL_ERR.
    //
    // Summary of the errors in GetSemiwidth:
    //
    // The distance a2 has a relative error of 5 * DBL_ERR plus an absolute
    // error of 64 * DBL_ERR^2 because the points "src" and "p" may differ from
    // unit length (by up to 4 * DBL_ERR).  We have already accounted for the
    // absolute error above, leaving only the relative error.
    //
    // sin2_r has a relative error of 2 * DBL_ERR.
    //
    // sin2_a has a relative error of 12 * DBL_ERR assuming that a2 <= 2,
    // i.e. distance(src, p) <= 90 degrees.  (The relative error gets
    // arbitrarily larger as this distance approaches 180 degrees.)
    //
    // semiwidth has a relative error of 17 * DBL_ERR.
    //
    // Finally, (center +/- semiwidth) has a rounding error of up to 4 * DBL_ERR
    // because in theory, the result magnitude may be as large as 1.5 * M_PI
    // which is larger than 4.0.  This gives a total error of:
    double error = (2 * 10 + 4) * DBL_ERR + 17 * DBL_ERR * semiwidth;
    return semiwidth + round_direction * error;
  }

  S2Point _src;
  S2Point _xDir, _yDir;
  S1Interval _window;
}
