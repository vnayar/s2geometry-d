// Copyright 2005 Google Inc. All Rights Reserved.
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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2testing;

import fluent.asserts;
import math = std.math;
import random = std.random;
import s2.r2point;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2loop;
import s2.s2point;
import s2.s2pointutil;
import s2.s2polygon;
import s2.s2region;
import s2.util.math.vector;
import s2.util.math.matrix3x3;

import std.algorithm;
import std.math;
import std.range;
import std.stdio;

// You can optionally call S2Testing::rnd.Reset(FLAGS_s2_random_seed) at the
// start of a test or benchmark to ensure that its results do not depend on
// which other tests of benchmarks have run previously.  This can help with
// debugging.
//
// This flag currently does *not* affect the initial seed value for
// S2Testing::rnd.  TODO(user): Fix this.
package int s2RandomSeed = 1;

enum double DOUBLE_ERR = 0.0001;

// This class defines various static functions that are useful for writing
// unit tests.
struct S2Testing {
public:

  // The Earth's mean radius in kilometers (according to NASA).
  static immutable double EARTH_RADIUS_KM = 6371.01;

  /**
   * Returns a vector of points shaped as a regular polygon with
   * num_vertices vertices, all on a circle of the specified angular
   * radius around the center.  The radius is the actual distance from
   * the center to the circle along the sphere.
   *
   * If you want to construct a regular polygon, try this:
   *   S2Polygon polygon(S2Loop::MakeRegularLoop(center, radius, num_vertices));
   */
  static S2Point[] makeRegularPoints(in S2Point center, in S1Angle radius, int num_vertices) {
    auto loop = S2Loop.makeRegularLoop(center, radius, num_vertices);
    S2Point[] points;
    for (int i = 0; i < loop.numVertices(); i++) {
      points ~= loop.vertex(i);
    }
    return points;
  }

  /** Append the vertices of "loop" to "vertices". */
  //static void AppendLoopVertices(const S2Loop& loop, std::vector<S2Point>* vertices);

  /**
   * A simple class that generates "Koch snowflake" fractals (see Wikipedia
   * for an introduction).  There is an option to control the fractal
   * dimension (between 1.0 and 2.0); values between 1.02 and 1.50 are
   * reasonable simulations of various coastlines.  The default dimension
   * (about 1.26) corresponds to the standard Koch snowflake.  (The west coast
   * of Britain has a fractal dimension of approximately 1.25.)
   *
   * The fractal is obtained by starting with an equilateral triangle and
   * recursively subdividing each edge into four segments of equal length.
   * Therefore the shape at level "n" consists of 3*(4**n) edges.  Multi-level
   * fractals are also supported: if you set min_level() to a non-negative
   * value, then the recursive subdivision has an equal probability of
   * stopping at any of the levels between the given min and max (inclusive).
   * This yields a fractal where the perimeter of the original triangle is
   * approximately equally divided between fractals at the various possible
   * levels.  If there are k distinct levels {min,..,max}, the expected number
   * of edges at each level "i" is approximately 3*(4**i)/k.
   */
  static class Fractal {
  public:
    // You must call set_max_level() or SetLevelForApproxMaxEdges() before
    // calling MakeLoop().
    this() {
      _maxLevel = -1;
      _minLevelArg = -1;
      _minLevel = -1;
      _dimension = log(4.0) / log(3.0);  /* standard Koch curve */
      _edgeFraction = 0;
      _offsetFraction = 0;
      computeOffsets();
    }

    // Set the maximum subdivision level for the fractal (see above).
    // REQUIRES: max_level >= 0
    void setMaxLevel(int max_level)
    in {
      assert(max_level >= 0);
    } body {
      _maxLevel = max_level;
      computeMinLevel();
    }

    int maxLevel() const {
      return _maxLevel;
    }

    // Set the minimum subdivision level for the fractal (see above).  The
    // default value of -1 causes the min and max levels to be the same.  A
    // min_level of 0 should be avoided since this creates a significant
    // chance that none of the three original edges will be subdivided at all.
    //
    // DEFAULT: max_level()
    void setMinLevel(int min_level_arg)
    in {
      assert(min_level_arg >= -1);
    } body {
      _minLevelArg = min_level_arg;
      computeMinLevel();
    }

    int minLevel() const {
      return _minLevelArg;
    }

    // Set the min and/or max level to produce approximately the given number
    // of edges.  (The values are rounded to a nearby value of 3*(4**n).)
    void setLevelForApproxMinEdges(int min_edges) {
      // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
      setMinLevel(cast(int) round(0.5 * log2(min_edges / 3)));
    }

    void setLevelForApproxMaxEdges(int max_edges) {
      // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
      setMaxLevel(cast(int) round(0.5 * log2(max_edges / 3)));
    }

    // Set the fractal dimension.  The default value of approximately 1.26
    // corresponds to the stardard Koch curve.  The value must lie in the
    // range [1.0, 2.0).
    //
    // DEFAULT: log(4) / log(3) ~= 1.26
    void setFractalDimension(double dimension)
    in {
      assert(dimension >= 1.0);
      assert(dimension < 2.0);
    } body {
      _dimension = dimension;
      computeOffsets();
    }

    double fractalDimension() const {
      return _dimension;
    }

    // Return a lower bound on ratio (Rmin / R), where "R" is the radius
    // passed to MakeLoop() and "Rmin" is the minimum distance from the
    // fractal boundary to its center, where all distances are measured in the
    // tangent plane at the fractal's center.  This can be used to inscribe
    // another geometric figure within the fractal without intersection.
    double minRadiusFactor() const {
      // The minimum radius is attained at one of the vertices created by the
      // first subdivision step as long as the dimension is not too small (at
      // least kMinDimensionForMinRadiusAtLevel1, see below).  Otherwise we fall
      // back on the incircle radius of the original triangle, which is always a
      // lower bound (and is attained when dimension = 1).
      //
      // The value below was obtained by letting AE be an original triangle edge,
      // letting ABCDE be the corresponding polyline after one subdivision step,
      // and then letting BC be tangent to the inscribed circle at the center of
      // the fractal O.  This gives rise to a pair of similar triangles whose edge
      // length ratios can be used to solve for the corresponding "edge fraction".
      // This method is slightly conservative because it is computed using planar
      // rather than spherical geometry.  The value below is equal to
      // -log(4)/log((2 + cbrt(2) - cbrt(4))/6).
      const double kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407;
      if (_dimension >= kMinDimensionForMinRadiusAtLevel1) {
        return sqrt(1 + 3 * _edgeFraction * (_edgeFraction - 1));
      }
      return 0.5;
    }

    // Return the ratio (Rmax / R), where "R" is the radius passed to
    // MakeLoop() and "Rmax" is the maximum distance from the fractal boundary
    // to its center, where all distances are measured in the tangent plane at
    // the fractal's center.  This can be used to inscribe the fractal within
    // some other geometric figure without intersection.
    double maxRadiusFactor() const {
      // The maximum radius is always attained at either an original triangle
      // vertex or at a middle vertex from the first subdivision step.
      return max(1.0, _offsetFraction * sqrt(3.0) + 0.5);
    }

    // Return a fractal loop centered around the z-axis of the given
    // coordinate frame, with the first vertex in the direction of the
    // positive x-axis.  In order to avoid self-intersections, the fractal is
    // generated by first drawing it in a 2D tangent plane to the unit sphere
    // (touching at the fractal's center point) and then projecting the edges
    // onto the sphere.  This has the side effect of shrinking the fractal
    // slightly compared to its nominal radius.
    S2Loop makeLoop(in Matrix3x3_d frame, in S1Angle nominal_radius) const {
      R2Point[] r2vertices;
      getR2Vertices(r2vertices);
      S2Point[] vertices;
      double r = nominal_radius.radians();
      foreach (const R2Point v; r2vertices) {
        auto p = S2Point(v[0] * r, v[1] * r, 1);
        vertices ~= fromFrame(frame, p).normalize();
      }
      return new S2Loop(vertices);
    }

  private:
    void computeMinLevel() {
      if (_minLevelArg >= 0 && _minLevelArg <= _maxLevel) {
        _minLevel = _minLevelArg;
      } else {
        _minLevel = _maxLevel;
      }
    }

    void computeOffsets() {
      _edgeFraction = pow(4.0, -1.0 / _dimension);
      _offsetFraction = sqrt(_edgeFraction - 0.25);
    }

    void getR2Vertices(ref R2Point[] vertices) const {
      // The Koch "snowflake" consists of three Koch curves whose initial edges
      // form an equilateral triangle.
      auto v0 = R2Point(1.0, 0.0);
      auto v1 = R2Point(-0.5, sqrt(3.0)/2);
      auto v2 = R2Point(-0.5, -sqrt(3.0)/2);
      getR2VerticesHelper(v0, v1, 0, vertices);
      getR2VerticesHelper(v1, v2, 0, vertices);
      getR2VerticesHelper(v2, v0, 0, vertices);
    }

    // Given the two endpoints (v0,v4) of an edge, recursively subdivide the edge
    // to the desired level, and insert all vertices of the resulting curve up to
    // but not including the endpoint "v4".
    void getR2VerticesHelper(
        in R2Point v0, in R2Point v4, int level, ref R2Point[] vertices) const {
      if (level >= _minLevel && S2Testing.rnd.oneIn(_maxLevel - level + 1)) {
        // Stop subdivision at this level.
        vertices ~= v0;
        return;
      }
      // Otherwise compute the intermediate vertices v1, v2, and v3.
      Vector2_d dir = v4 - v0;
      R2Point v1 = v0 + _edgeFraction * dir;
      R2Point v2 = 0.5 * (v0 + v4) - _offsetFraction * dir.ortho();
      R2Point v3 = v4 - _edgeFraction * dir;

      // And recurse on the four sub-edges.
      getR2VerticesHelper(v0, v1, level + 1, vertices);
      getR2VerticesHelper(v1, v2, level + 1, vertices);
      getR2VerticesHelper(v2, v3, level + 1, vertices);
      getR2VerticesHelper(v3, v4, level + 1, vertices);
    }

    int _maxLevel;
    int _minLevelArg;  // Value set by user
    int _minLevel;      // Actual min level (depends on max_level_)
    double _dimension;

    // The ratio of the sub-edge length to the original edge length at each
    // subdivision step.
    double _edgeFraction;

    // The distance from the original edge to the middle vertex at each
    // subdivision step, as a fraction of the original edge length.
    double _offsetFraction;
  }

  // Convert a distance on the Earth's surface to an angle.
  // Do not use these methods in non-testing code; use s2earth.h instead.
  static S1Angle metersToAngle(double meters) {
    return kmToAngle(0.001 * meters);
  }

  static S1Angle kmToAngle(double km) {
    return S1Angle.fromRadians(km / EARTH_RADIUS_KM);
  }

  // Convert an area in steradians (as returned by the S2 area methods) to
  // square meters or square kilometers.
  //static double AreaToMeters2(double steradians);
  //static double AreaToKm2(double steradians);

  // The Earth's mean radius in kilometers (according to NASA).
  //static const double kEarthRadiusKm;

  // A deterministically-seeded random number generator.
  static Random rnd = Random();

  // Return a random unit-length vector.
  static S2Point randomPoint() {
    // The order of evaluation of function arguments is unspecified,
    // so we may not just call S2Point with three RandDouble-based args.
    // Use temporaries to induce sequence points between calls.
    double x = rnd.uniformDouble(-1.0, 1.0);
    double y = rnd.uniformDouble(-1.0, 1.0);
    double z = rnd.uniformDouble(-1.0, 1.0);
    return S2Point(x, y, z).normalize();
  }

  // Return a right-handed coordinate frame (three orthonormal vectors).
  static void getRandomFrame(out S2Point x, out S2Point y, out S2Point z) {
    z = randomPoint();
    getRandomFrameAt(z, x, y);
  }

  static Matrix3x3_d getRandomFrame() {
    return getRandomFrameAt(randomPoint());
  }

  // Given a unit-length z-axis, compute x- and y-axes such that (x,y,z) is a
  // right-handed coordinate frame (three orthonormal vectors).
  static void getRandomFrameAt(in S2Point z, out S2Point x, out S2Point y) {
    x = z.crossProd(randomPoint()).normalize();
    y = z.crossProd(x).normalize();
  }

  static Matrix3x3_d getRandomFrameAt(in S2Point z) {
    S2Point x, y;
    getRandomFrameAt(z, x, y);
    return Matrix3x3_d.fromCols(x, y, z);
  }

  // Return a cap with a random axis such that the log of its area is
  // uniformly distributed between the logs of the two given values.
  // (The log of the cap angle is also approximately uniformly distributed.)
  static S2Cap getRandomCap(double min_area, double max_area) {
    double cap_area = max_area * math.pow(min_area / max_area, rnd.randDouble());
    Assert.notLessThan(cap_area, min_area);
    Assert.notGreaterThan(cap_area, max_area);

    // The surface area of a cap is 2*Pi times its height.
    return S2Cap.fromCenterArea(randomPoint(), cap_area);
  }

  // Return a point chosen uniformly at random (with respect to area)
  // from the given cap.
  static S2Point samplePoint(in S2Cap cap) {
    // We consider the cap axis to be the "z" axis.  We choose two other axes to
    // complete the coordinate frame.

    Matrix3x3_d m;
    getFrame(cap.center(), m);

    // The surface area of a spherical cap is directly proportional to its
    // height.  First we choose a random height, and then we choose a random
    // point along the circle at that height.

    double h = rnd.randDouble() * cap.height();
    double theta = 2 * math.PI * rnd.randDouble();
    double r = math.sqrt(h * (2 - h));  // Radius of circle.

    // The result should already be very close to unit-length, but we might as
    // well make it accurate as possible.
    return fromFrame(m, S2Point(math.cos(theta) * r, math.sin(theta) * r, 1 - h)).normalize();
  }

  // Return a random cell id at the given level or at a randomly chosen
  // level.  The distribution is uniform over the space of cell ids,
  // but only approximately uniform over the surface of the sphere.
  static S2CellId getRandomCellId(int level) {
    int face = rnd.uniform(S2CellId.NUM_FACES);
    ulong pos = rnd.rand64() & ((1uL << S2CellId.POS_BITS) - 1);
    return S2CellId.fromFacePosLevel(face, pos, level);
  }

  static S2CellId getRandomCellId() {
    return getRandomCellId(rnd.uniform(S2CellId.MAX_LEVEL + 1));
  }

  /// Return a polygon with the specified center, number of concentric loops
  /// and vertices per loop.
  static void concentricLoopsPolygon(
      in S2Point center, int num_loops, int num_vertices_per_loop, S2Polygon polygon) {
    writeln("s2testing.concentricLoopsPolygon >");
    scope(exit) writeln("s2testing.concentricLoopsPolygon <");
    writeln("s2testing.concentricLoopsPolygon 1: num_loops=", num_loops,
        ", num_vertices_per_loop=", num_vertices_per_loop);
    Matrix3x3_d m;
    getFrame(center, m);
    /**/m = Matrix3x3_d(
        0.418314995042454, 0.701198425261565, 0.577350269189626,
        0.816413146891386, -0.0116721998599441, -0.577350269189626,
        -0.398098151848932, 0.712870625121509, -0.577350269189626
    );
    writeln("s2testing.concentricLoopsPolygon 2: m=", m);
    S2Loop[] loops;
    for (int li = 0; li < num_loops; ++li) {
      writeln("s2testing.concentricLoopsPolygon 3:");
      S2Point[] vertices;
      double radius = 0.005 * (li + 1) / num_loops;
      double radian_step = 2 * M_PI / num_vertices_per_loop;
      writeln("s2testing.concentricLoopsPolygon 4: radius=", radius, ", radian_step=", radian_step);
      for (int vi = 0; vi < num_vertices_per_loop; ++vi) {
        double angle = vi * radian_step;
        auto p = S2Point(radius * math.cos(angle), radius * math.sin(angle), 1);
        vertices ~= fromFrame(m, p.normalize());
      }
      writeln("s2testing.concentricLoopsPolygon 5:");
      loops ~= new S2Loop(vertices);
    }
    polygon.initializeNested(loops);
  }

  // Checks that "covering" completely covers the given region.  If
  // "check_tight" is true, also checks that it does not contain any cells
  // that do not intersect the given region.  ("id" is only used internally.)
  static void checkCovering(
      S2Region region, in S2CellUnion covering, bool check_tight, S2CellId id = S2CellId()) {
    if (!id.isValid()) {
      for (int face = 0; face < 6; ++face) {
        checkCovering(region, covering, check_tight, S2CellId.fromFace(face));
      }
      return;
    }

    if (!region.mayIntersect(new S2Cell(id))) {
      // If region does not intersect id, then neither should the covering.
      if (check_tight) Assert.equal(covering.intersects(id), false);

    } else if (!covering.contains(id)) {
      // The region may intersect id, but we can't assert that the covering
      // intersects id because we may discover that the region does not actually
      // intersect upon further subdivision.  (MayIntersect is not exact.)
      Assert.equal(region.contains(new S2Cell(id)), false);
      Assert.equal(id.isLeaf(), false);
      S2CellId end = id.childEnd();
      S2CellId child;
      for (child = id.childBegin(); child != end; child = child.next()) {
        checkCovering(region, covering, check_tight, child);
      }
    }
  }

  // // Returns the user time consumed by this process, in seconds.
  // static double GetCpuTime();
}

// Functions in this class return random numbers that are as good as random()
// is.  The results are reproducible since the seed is deterministic.  This
// class is *NOT* thread-safe; it is only intended for testing purposes.
struct Random {
private:
  random.Random _rnd = random.Random(1);

public:
  // Initialize using a deterministic seed.

  // Reset the generator state using the given seed.
  void reset(int seed) {
    _rnd.seed(seed);
  }

  // Return a uniformly distributed 64-bit unsigned integer.
  ulong rand64() {
    return random.uniform(0, ulong.max, _rnd);
  }

  // Return a uniformly distributed 32-bit unsigned integer.
  uint rand32() {
    return random.uniform(0, uint.max, _rnd);
  }

  private ulong getBits(int num_bits)
  in {
    assert(num_bits >= 0);
    assert(num_bits <= 64);
  } body {
    return rand64() & ((1uL << num_bits) - 1);
  }

  // Return a uniformly distributed "double" in the range [0,1).  Note that
  // the values returned are all multiples of 2**-53, which means that not all
  // possible values in this range are returned.
  double randDouble() {
    const int NUM_BITS = 53;
    return math.ldexp(cast(double) getBits(NUM_BITS), -NUM_BITS);
  }

  // Return a uniformly distributed integer in the range [0,n).
  int uniform(int n) {
    if (n == 0) {
      return 0;
    }
    return random.uniform(0, n, _rnd);
  }

  // Return a uniformly distributed "double" in the range [min, limit).
  double uniformDouble(double min, double limit) {
    return min + random.uniform01(_rnd) * (limit - min);
  }

  // A functor-style version of Uniform, so that this class can be used with
  // STL functions that require a RandomNumberGenerator concept.
  //int32 operator() (int32 n) {
  //  return Uniform(n);
  //}

  // Return true with probability 1 in n.
  bool oneIn(int n) {
    return uniform(n) == 0;
  }

  // Skewed: pick "base" uniformly from range [0,max_log] and then
  // return "base" random bits.  The effect is to pick a number in the
  // range [0,2^max_log-1] with bias towards smaller numbers.
  int skewed(int max_log)
  in {
    assert(max_log >= 0);
    assert(max_log <= 31);
  } body {
    int base = uniform(max_log + 1);
    return cast(int) getBits(31) & ((1u << base) - 1);
  }
}

struct Pair(T1, T2) {
  T1 first;
  T2 second;
}

// Compare two sets of "closest" items, where "expected" is computed via brute
// force (i.e., considering every possible candidate) and "actual" is computed
// using a spatial data structure.  Here "max_size" is a bound on the maximum
// number of items, "max_distance" is a limit on the distance to any item, and
// "max_error" is the maximum error allowed when selecting which items are
// closest (see S2ClosestEdgeQuery::Options::max_error).
bool checkDistanceResults(Id)(
    in Pair!(S1Angle, Id)[] expected,
    in Pair!(S1Angle, Id)[] actual,
    int max_size, S1Angle max_distance, S1Angle max_error) {
  // This is a conservative bound on the error in computing the distance from
  // the target geometry to an S2Cell.  Such errors can cause candidates to be
  // pruned from the result set even though they may be slightly closer.
  enum S1Angle kMaxPruningError = S1Angle.fromRadians(1e-15);
  return (checkResultSet(
          actual, expected, max_size, max_distance, max_error,
          kMaxPruningError, "Missing") & /*not &&*/
      checkResultSet(
          expected, actual, max_size, max_distance, max_error,
          S1Angle.zero(), "Extra"));
}


// Check that result set "x" contains all the expected results from "y", and
// does not include any duplicate results.
bool checkResultSet(Id)(
    in Pair!(S1Angle, Id)[] x,
    in Pair!(S1Angle, Id)[] y,
    int max_size, S1Angle max_distance, S1Angle max_error,
    S1Angle max_pruning_error, string label) {
  alias Result = Pair!(S1Angle, Id);
  // Results should be sorted by distance, but not necessarily then by Id.
  Assert.equal(isSorted!"a.first < b.first"(x), true);

  // Result set X should contain all the items from U whose distance is less
  // than "limit" computed below.
  S1Angle limit = S1Angle.zero();
  if (x.length < max_size) {
    // Result set X was not limited by "max_size", so it should contain all
    // the items up to "max_distance", except that a few items right near the
    // distance limit may be missed because the distance measurements used for
    // pruning S2Cells are not conservative.
    limit = max_distance - max_pruning_error;
  } else if (x.length != 0) {
    // Result set X contains only the closest "max_size" items, to within a
    // tolerance of "max_error + max_pruning_error".
    limit = x.back().first - max_error - max_pruning_error;
  }
  bool result = true;
  foreach (const yp; y) {
    // Note that this test also catches duplicate values.
    size_t count = count!"a.second == b.second"(x, yp);
    if (yp.first < limit && count != 1) {
      result = false;
      writeln((count > 1 ? "Duplicate" : label), " distance = ",
          yp.first, ", id = ", yp.second);
    }
  }
  return result;
}
