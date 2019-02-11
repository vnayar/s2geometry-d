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

module s2.builder.util.snap_functions;

import s2.s1angle;
import s2.s2builder;
import s2.s2builder;
import s2.s2cell_id;
import s2.s2latlng;
import s2.s2metrics;
import s2.s2point;

import std.algorithm;
import std.math;

/**
 * A SnapFunction that snaps every vertex to itself.  It should be used when
 * vertices do not need to be snapped to a discrete set of locations (such as
 * E7 lat/lngs), or when maximum accuracy is desired.
 *
 * If the given "snap_radius" is zero, then all input vertices are preserved
 * exactly.  Otherwise, S2Builder merges nearby vertices to ensure that no
 * vertex pair is closer than "snap_radius".  Furthermore, vertices are
 * separated from non-incident edges by at least "min_edge_vertex_separation",
 * equal to (0.5 * snap_radius).  For example, if the snap_radius is 1km, then
 * vertices will be separated from non-incident edges by at least 500m.
 */
class IdentitySnapFunction : S2Builder.SnapFunction {
public:
  /// The default constructor uses a snap_radius of zero (i.e., no snapping).
  this() {
    _snapRadius = S1Angle.zero();
  }

  /// Convenience constructor that calls set_snap_radius().
  this(S1Angle snap_radius) {
    setSnapRadius(snap_radius);
  }

  this(in IdentitySnapFunction other) {
    _snapRadius = other._snapRadius;
  }

  /// REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
  void setSnapRadius(S1Angle snap_radius)
  in {
    assert(snap_radius <= kMaxSnapRadius());
  } do {
    _snapRadius = snap_radius;
  }

  override
  S1Angle snapRadius() const {
    return _snapRadius;
  }

  /// For the identity snap function, all vertex pairs are separated by at
  /// least snap_radius().
  override
  S1Angle minVertexSeparation() const {
    // Since SnapFunction does not move the input point, output vertices are
    // separated by the full snap_radius().
    return _snapRadius;
  }

  // For the identity snap function, edges are separated from all non-incident
  // vertices by at least 0.5 * snap_radius().
  override
  S1Angle minEdgeVertexSeparation() const {
    // In the worst case configuration, the edge separation is half of the
    // vertex separation.
    return 0.5 * _snapRadius;
  }

  override
  S2Point snapPoint(in S2Point point) const {
    return point;
  }

  override
  S2Builder.SnapFunction clone() const {
    return new IdentitySnapFunction(this);
  }

private:
  // Copying and assignment are allowed.
  S1Angle _snapRadius;
}

/**
 * A SnapFunction that snaps vertices to S2CellId centers.  This can be useful
 * if you want to encode your geometry compactly using S2Polygon::Encode(),
 * for example.  You can snap to the centers of cells at any level.
 *
 * Every snap level has a corresponding minimum snap radius, which is simply
 * the maximum distance that a vertex can move when snapped.  It is
 * approximately equal to half of the maximum diagonal length for cells at the
 * chosen level.  You can also set the snap radius to a larger value; for
 * example, you could snap to the centers of leaf cells (1cm resolution) but
 * set the snap_radius() to 10m.  This would result in significant extra
 * simplification, without moving vertices unnecessarily (i.e., vertices that
 * are at least 10m away from all other vertices will move by less than 1cm).
 */
class S2CellIdSnapFunction : S2Builder.SnapFunction {
public:
  /// The default constructor snaps to S2CellId::kMaxLevel (i.e., the centers
  /// of leaf cells), and uses the minimum allowable snap radius at that level.
  this() {
    setLevel(S2CellId.MAX_LEVEL);
  }

  /// Convenience constructor equivalent to calling set_level(level).
  this(int level) {
    setLevel(level);
  }

  this(in S2CellIdSnapFunction other) {
    _level = other._level;
    _snapRadius = other._snapRadius;
  }

  /**
   * Snaps vertices to S2Cell centers at the given level.  As a side effect,
   * this method also resets "snap_radius" to the minimum value allowed at
   * this level:
   *
   *   set_snap_radius(MinSnapRadiusForLevel(level))
   *
   * This means that if you want to use a larger snap radius than the minimum,
   * you must call set_snap_radius() *after* calling set_level().
   */
  void setLevel(int level)
  in {
    assert(level >= 0);
    assert(level <= S2CellId.MAX_LEVEL);
  } do {
    _level = level;
    setSnapRadius(minSnapRadiusForLevel(level));
  }

  int level() const {
    return _level;
  }

  /**
   * Defines the snap radius to be used (see s2builder.h).  The snap radius
   * must be at least the minimum value for the current level(), but larger
   * values can also be used (e.g., to simplify the geometry).
   *
   * REQUIRES: snap_radius >= MinSnapRadiusForLevel(level())
   * REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
   */
  void setSnapRadius(S1Angle snap_radius)
  in {
    assert(snap_radius >= minSnapRadiusForLevel(level()));
    assert(snap_radius <= kMaxSnapRadius());
  } do {
    _snapRadius = snap_radius;
  }

  override
  S1Angle snapRadius() const {
      return _snapRadius;
  }

  /// Returns the minimum allowable snap radius for the given S2Cell level
  /// (approximately equal to half of the maximum cell diagonal length).
  static S1Angle minSnapRadiusForLevel(int level) {
    // snap_radius() needs to be an upper bound on the true distance that a
    // point can move when snapped, taking into account numerical errors.
    //
    // The maximum error when converting from an S2Point to an S2CellId is
    // S2::kMaxDiag.deriv() * DBL_EPSILON.  The maximum error when converting an
    // S2CellId center back to an S2Point is 1.5 * DBL_EPSILON.  These add up to
    // just slightly less than 4 * DBL_EPSILON.
    return S1Angle.fromRadians(0.5 * MAX_DIAG.getValue(level) + 4 * double.epsilon);
  }

  /**
   * Returns the minimum S2Cell level (i.e., largest S2Cells) such that
   * vertices will not move by more than "snap_radius".  This can be useful
   * when choosing an appropriate level to snap to.  The return value is
   * always a valid level (out of range values are silently clamped).
   *
   * If you want to choose the snap level based on a distance, and then use
   * the minimum possible snap radius for the chosen level, do this:
   *
   *   S2CellIdSnapFunction f(
   *       S2CellIdSnapFunction::LevelForMaxSnapRadius(distance));
   */
  static int levelForMaxSnapRadius(S1Angle snap_radius) {
    // When choosing a level, we need to acount for the error bound of
    // 4 * DBL_EPSILON that is added by MinSnapRadiusForLevel().
    return MAX_DIAG.getLevelForMaxValue(2 * (snap_radius.radians() - 4 * double.epsilon));
  }

  // TODO: Resume here.

  /**
   * For S2CellId snapping, the minimum separation between vertices depends on
   * level() and snap_radius().  It can vary between 0.5 * snap_radius()
   * and snap_radius().
   */
  override
  S1Angle minVertexSeparation() const {
    // We have three different bounds for the minimum vertex separation: one is
    // a constant bound, one is proportional to snap_radius, and one is equal to
    // snap_radius minus a constant.  These bounds give the best results for
    // small, medium, and large snap radii respectively.  We return the maximum
    // of the three bounds.
    //
    // 1. Constant bound: Vertices are always separated by at least
    //    kMinEdge(level), the minimum edge length for the chosen snap level.
    //
    // 2. Proportional bound: It can be shown that in the plane, the worst-case
    //    configuration has a vertex separation of 2 / sqrt(13) * snap_radius.
    //    This is verified in the unit test, except that on the sphere the ratio
    //    is slightly smaller at cell level 2 (0.54849 vs. 0.55470).  We reduce
    //    that value a bit more below to be conservative.
    //
    // 3. Best asymptotic bound: This bound bound is derived by observing we
    //    only select a new site when it is at least snap_radius() away from all
    //    existing sites, and the site can move by at most 0.5 * kMaxDiag(level)
    //    when snapped.
    S1Angle min_edge = S1Angle.fromRadians(MIN_EDGE.getValue(_level));
    S1Angle max_diag = S1Angle.fromRadians(MAX_DIAG.getValue(_level));
    return max(min_edge,
        max(0.548 * _snapRadius,  // 2 / sqrt(13) in the plane
            _snapRadius - 0.5 * max_diag));
  }

  /**
   * For S2CellId snapping, the minimum separation between edges and
   * non-incident vertices depends on level() and snap_radius().  It can
   * be as low as 0.219 * snap_radius(), but is typically 0.5 * snap_radius()
   * or more.
   */
  override
  S1Angle minEdgeVertexSeparation() const {
    // Similar to min_vertex_separation(), in this case we have four bounds: a
    // constant bound that holds only at the minimum snap radius, a constant
    // bound that holds for any snap radius, a bound that is proportional to
    // snap_radius, and a bound that is equal to snap_radius minus a constant.
    //
    // 1. Constant bounds:
    //
    //    (a) At the minimum snap radius for a given level, it can be shown that
    //    vertices are separated from edges by at least 0.5 * kMinDiag(level) in
    //    the plane.  The unit test verifies this, except that on the sphere the
    //    worst case is slightly better: 0.5652980068 * kMinDiag(level).
    //
    //    (b) Otherwise, for arbitrary snap radii the worst-case configuration
    //    in the plane has an edge-vertex separation of sqrt(3/19) *
    //    kMinDiag(level), where sqrt(3/19) is about 0.3973597071.  The unit
    //    test verifies that the bound is slighty better on the sphere:
    //    0.3973595687 * kMinDiag(level).
    //
    // 2. Proportional bound: In the plane, the worst-case configuration has an
    //    edge-vertex separation of 2 * sqrt(3/247) * snap_radius, which is
    //    about 0.2204155075.  The unit test verifies this, except that on the
    //    sphere the bound is slightly worse for certain large S2Cells: the
    //    minimum ratio occurs at cell level 6, and is about 0.2196666953.
    //
    // 3. Best asymptotic bound: If snap_radius() is large compared to the
    //    minimum snap radius, then the best bound is achieved by 3 sites on a
    //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
    //    apart.  An input edge passing just to one side of the center of the
    //    circle intersects the Voronoi regions of the two end sites but not the
    //    Voronoi region of the center site, and gives an edge separation of
    //    (min_vertex_separation ** 2) / (2 * snap_radius).  This bound
    //    approaches 0.5 * snap_radius for large snap radii, i.e.  the minimum
    //    edge-vertex separation approaches half of the minimum vertex
    //    separation as the snap radius becomes large compared to the cell size.

    S1Angle min_diag = S1Angle.fromRadians(MIN_DIAG.getValue(_level));
    if (snapRadius() == minSnapRadiusForLevel(_level)) {
      // This bound only holds when the minimum snap radius is being used.
      return 0.565 * min_diag;            // 0.500 in the plane
    }
    // Otherwise, these bounds hold for any snap_radius().
    S1Angle vertex_sep = minVertexSeparation();
    return max(0.397 * min_diag,          // sqrt(3 / 19) in the plane
        max(0.219 * _snapRadius,  // 2 * sqrt(3 / 247) in the plane
            0.5 * (vertex_sep / _snapRadius) * vertex_sep));
  }

  override
  S2Point snapPoint(in S2Point point) const {
    return S2CellId(point).parent(_level).toS2Point();
  }

  override
  S2Builder.SnapFunction clone() const {
    return new S2CellIdSnapFunction(this);
  }

private:
  // Copying and assignment are allowed.
  int _level;
  S1Angle _snapRadius;
}

/**
 * A SnapFunction that snaps vertices to S2LatLng E5, E6, or E7 coordinates.
 * These coordinates are expressed in degrees multiplied by a power of 10 and
 * then rounded to the nearest integer.  For example, in E6 coordinates the
 * point (23.12345651, -45.65432149) would become (23123457, -45654321).
 *
 * The main argument of the SnapFunction is the exponent for the power of 10
 * that coordinates should be multipled by before rounding.  For example,
 * IntLatLngSnapFunction(7) is a function that snaps to E7 coordinates.  The
 * exponent can range from 0 to 10.
 *
 * Each exponent has a corresponding minimum snap radius, which is simply the
 * maximum distance that a vertex can move when snapped.  It is approximately
 * equal to 1/sqrt(2) times the nominal point spacing; for example, for
 * snapping to E7 the minimum snap radius is (1e-7 / sqrt(2)) degrees.
 * You can also set the snap radius to any value larger than this; this can
 * result in significant extra simplification (similar to using a larger
 * exponent) but does not move vertices unnecessarily.
 */
class IntLatLngSnapFunction : S2Builder.SnapFunction {
public:
  /// The default constructor yields an invalid snap function.  You must set
  /// the exponent explicitly before using it.
  this() {
    _exponent = -1;
    _snapRadius = S1Angle();
    _fromDegrees = 0;
    _toDegrees = 0;
  }

  /// Convenience constructor equivalent to calling set_exponent(exponent).
  this(int exponent) {
      setExponent(exponent);
  }

  this(in IntLatLngSnapFunction other) {
    _exponent = other._exponent;
    _snapRadius = other._snapRadius;
    _fromDegrees = other._fromDegrees;
    _toDegrees = other._toDegrees;
  }

  /**
   * Snaps vertices to points whose (lat, lng) coordinates are integers after
   * converting to degrees and multiplying by 10 raised to the given exponent.
   * For example, (exponent == 7) yields E7 coordinates.  As a side effect,
   * this method also resets "snap_radius" to the minimum value allowed for
   * this exponent:
   *
   *   set_snap_radius(MinSnapRadiusForExponent(exponent))
   *
   * This means that if you want to use a larger snap radius than the minimum,
   * you must call set_snap_radius() *after* calling set_exponent().
   *
   * REQUIRES: kMinExponent <= exponent <= kMaxExponent
   */
  void setExponent(int exponent)
  in {
    assert(exponent >= MIN_EXPONENT);
    assert(exponent <= MAX_EXPONENT);
  } do {
    _exponent = exponent;
    setSnapRadius(minSnapRadiusForExponent(exponent));

    // Precompute the scale factors needed for snapping.  Note that these
    // calculations need to exactly match the ones in s1angle.h to ensure
    // that the same S2Points are generated.
    double power = 1;
    for (int i = 0; i < exponent; ++i) power *= 10;
    _fromDegrees = power;
    _toDegrees = 1 / power;
  }

  int exponent() const {
    return _exponent;
  }

  /// The minum exponent supported for snapping.
  enum int MIN_EXPONENT = 0;

  /// The maximum exponent supported for snapping.
  enum int MAX_EXPONENT = 10;

  /**
   * Defines the snap radius to be used (see s2builder.h).  The snap radius
   * must be at least the minimum value for the current exponent(), but larger
   * values can also be used (e.g., to simplify the geometry).
   *
   * REQUIRES: snap_radius >= MinSnapRadiusForExponent(exponent())
   * REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
   */
  void setSnapRadius(S1Angle snap_radius)
  in {
    assert(snap_radius >= minSnapRadiusForExponent(exponent()));
    assert(snap_radius <= kMaxSnapRadius());
  } do {
    _snapRadius = snap_radius;
  }

  override
  S1Angle snapRadius() const {
    return _snapRadius;
  }

  // Returns the minimum allowable snap radius for the given exponent
  // (approximately equal to (pow(10, -exponent) / sqrt(2)) degrees).
  static S1Angle minSnapRadiusForExponent(int exponent) {
    // snap_radius() needs to be an upper bound on the true distance that a
    // point can move when snapped, taking into account numerical errors.
    //
    // The maximum errors in latitude and longitude can be bounded as
    // follows (as absolute errors in terms of DBL_EPSILON):
    //
    //                                      Latitude      Longitude
    // Convert to S2LatLng:                    1.000          1.000
    // Convert to degrees:                     1.032          2.063
    // Scale by 10**exp:                       0.786          1.571
    // Round to integer: 0.5 * S1Angle::Degrees(to_degrees_)
    // Scale by 10**(-exp):                    1.375          2.749
    // Convert to radians:                     1.252          1.503
    // ------------------------------------------------------------
    // Total (except for rounding)             5.445          8.886
    //
    // The maximum error when converting the S2LatLng back to an S2Point is
    //
    //   sqrt(2) * (maximum error in latitude or longitude) + 1.5 * DBL_EPSILON
    //
    // which works out to (9 * sqrt(2) + 1.5) * DBL_EPSILON radians.  Finally
    // we need to consider the effect of rounding to integer coordinates
    // (much larger than the errors above), which can change the position by
    // up to (sqrt(2) * 0.5 * to_degrees_) radians.
    double power = 1;
    for (int i = 0; i < exponent; ++i) power *= 10;
    return (S1Angle.fromDegrees(SQRT1_2 / power) +
        S1Angle.fromRadians((9 * SQRT2 + 1.5) * double.epsilon));
  }

  /**
   * Returns the minimum exponent such that vertices will not move by more
   * than "snap_radius".  This can be useful when choosing an appropriate
   * exponent for snapping.  The return value is always a valid exponent
   * (out of range values are silently clamped).
   *
   * If you want to choose the exponent based on a distance, and then use
   * the minimum possible snap radius for that exponent, do this:
   *
   *   IntLatLngSnapFunction f(
   *       IntLatLngSnapFunction::ExponentForMaxSnapRadius(distance));
   */
  static int exponentForMaxSnapRadius(S1Angle snap_radius) {
    // When choosing an exponent, we need to acount for the error bound of
    // (9 * sqrt(2) + 1.5) * DBL_EPSILON added by MinSnapRadiusForExponent().
    snap_radius -= S1Angle.fromRadians((9 * SQRT2 + 1.5) * double.epsilon);
    snap_radius = max(snap_radius, S1Angle.fromRadians(1e-30));
    double exponent = log10(SQRT1_2 / snap_radius.degrees());

    // There can be small errors in the calculation above, so to ensure that
    // this function is the inverse of MinSnapRadiusForExponent() we subtract a
    // small error tolerance.
    return max(MIN_EXPONENT,
        min(MAX_EXPONENT, cast(int) ceil(exponent - 2 * double.epsilon)));
  }

  /**
   * For IntLatLng snapping, the minimum separation between vertices depends on
   * exponent() and snap_radius().  It can vary between snap_radius()
   * and snap_radius().
   */
  override
  S1Angle minVertexSeparation() const {
    // We have two bounds for the minimum vertex separation: one is proportional
    // to snap_radius, and one is equal to snap_radius minus a constant.  These
    // bounds give the best results for small and large snap radii respectively.
    // We return the maximum of the two bounds.
    //
    // 1. Proportional bound: It can be shown that in the plane, the worst-case
    //    configuration has a vertex separation of (sqrt(2) / 3) * snap_radius.
    //    This is verified in the unit test, except that on the sphere the ratio
    //    is slightly smaller (0.471337 vs. 0.471404).  We reduce that value a
    //    bit more below to be conservative.
    //
    // 2. Best asymptotic bound: This bound bound is derived by observing we
    //    only select a new site when it is at least snap_radius() away from all
    //    existing sites, and snapping a vertex can move it by up to
    //    ((1 / sqrt(2)) * to_degrees_) degrees.
    return max(0.471 * _snapRadius,        // sqrt(2) / 3 in the plane
        _snapRadius - S1Angle.fromDegrees(SQRT1_2 * _toDegrees));
  }

  /**
   * For IntLatLng snapping, the minimum separation between edges and
   * non-incident vertices depends on level() and snap_radius().  It can
   * be as low as 0.222 * snap_radius(), but is typically 0.39 * snap_radius()
   * or more.
   */
  override
  S1Angle minEdgeVertexSeparation() const {
    // Similar to min_vertex_separation(), in this case we have three bounds:
    // one is a constant bound, one is proportional to snap_radius, and one is
    // equal to snap_radius minus a constant.
    //
    // 1. Constant bound: In the plane, the worst-case configuration has an
    //    edge-vertex separation of ((1 / sqrt(13)) * to_degrees_) degrees.
    //    The unit test verifies this, except that on the sphere the ratio is
    //    slightly lower when small exponents such as E1 are used
    //    (0.2772589 vs 0.2773501).
    //
    // 2. Proportional bound: In the plane, the worst-case configuration has an
    //    edge-vertex separation of (2 / 9) * snap_radius (0.222222222222).  The
    //    unit test verifies this, except that on the sphere the bound can be
    //    slightly worse with large exponents (e.g., E9) due to small numerical
    //    errors (0.222222126756717).
    //
    // 3. Best asymptotic bound: If snap_radius() is large compared to the
    //    minimum snap radius, then the best bound is achieved by 3 sites on a
    //    circular arc of radius "snap_radius", spaced "min_vertex_separation"
    //    apart (see S2CellIdSnapFunction::min_edge_vertex_separation).  This
    //    bound approaches 0.5 * snap_radius as the snap radius becomes large
    //    relative to the grid spacing.

    S1Angle vertex_sep = minVertexSeparation();
    return max(0.277 * S1Angle.fromDegrees(_toDegrees),  // 1/sqrt(13) in the plane
        max(0.222 * _snapRadius,               // 2/9 in the plane
            0.5 * (vertex_sep / _snapRadius) * vertex_sep));
  }

  override
  S2Point snapPoint(in S2Point point) const
  in {
    assert(_exponent >= 0);  // Make sure the snap function was initialized.
  } do {
    auto input = S2LatLng(point);
    long lat = lround(input.lat().degrees() * _fromDegrees);
    long lng = lround(input.lng().degrees() * _fromDegrees);
    return S2LatLng.fromDegrees(lat * _toDegrees, lng * _toDegrees).toS2Point();
  }

  override
  S2Builder.SnapFunction clone() const {
    return new IntLatLngSnapFunction(this);
  }

private:
  // Copying and assignment are allowed.
  int _exponent;
  S1Angle _snapRadius;
  double _fromDegrees;
  double _toDegrees;
}
