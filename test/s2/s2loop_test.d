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
//

// Original Author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2loop_test;

import s2.s2loop;

import fluent.asserts;
import s2.logger;
import s2.r1interval;
import s2.s1angle;
import s2.s1interval;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2debug;
import s2.s2edge_crossings;
import s2.s2edge_distances;
import s2.s2error;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2measures;
//import s2.s2point_compression;
import s2.s2point;
import s2.s2pointutil;
import s2.s2predicates;
import s2.s2testing;
import s2.s2text_format;
import s2.util.math.matrix3x3;
import s2.util.math.vector;

import std.algorithm;
import std.math;
import std.conv;

import std.stdio;
import std.exception : enforce;
import std.algorithm : find;
import std.range;

class S2LoopTestBase {
public:
  // The set of all loops declared below.
  S2Loop[] allLoops;

  // Some standard loops to use in the tests (see descriptions below).
  S2Loop empty;
  S2Loop full;
  S2Loop northHemi;
  S2Loop northHemi3;
  S2Loop southHemi;
  S2Loop westHemi;
  S2Loop eastHemi;
  S2Loop nearHemi;
  S2Loop farHemi;
  S2Loop candyCane;
  S2Loop smallNeCw;
  S2Loop arctic80;
  S2Loop antarctic80;
  S2Loop lineTriangle;
  S2Loop skinnyChevron;
  S2Loop loopA;
  S2Loop loopB;
  S2Loop aIntersectB;
  S2Loop aUnionB;
  S2Loop aMinusB;
  S2Loop bMinusA;
  S2Loop loopC;
  S2Loop loopD;
  S2Loop loopE;
  S2Loop loopF;
  S2Loop loopG;
  S2Loop loopH;
  S2Loop loopI;
  S2Loop snappedLoopA;

public:
  this() {
    // The set of all loops declared below.
    allLoops = [];

    // The empty loop.
    empty = addLoop(new S2Loop(S2Loop.empty()));

    // The full loop.
    full = addLoop(new S2Loop(S2Loop.full()));

    // The northern hemisphere, defined using two pairs of antipodal points.
    northHemi = addLoop("0:-180, 0:-90, 0:0, 0:90");

    // The northern hemisphere, defined using three points 120 degrees apart.
    northHemi3 = addLoop("0:-180, 0:-60, 0:60");

    // The southern hemisphere, defined using two pairs of antipodal points.
    southHemi = addLoop("0:90, 0:0, 0:-90, 0:-180");

    // The western hemisphere, defined using two pairs of antipodal points.
    westHemi = addLoop("0:-180, -90:0, 0:0, 90:0");

    // The eastern hemisphere, defined using two pairs of antipodal points.
    eastHemi = addLoop("90:0, 0:0, -90:0, 0:-180");

    // The "near" hemisphere, defined using two pairs of antipodal points.
    nearHemi = addLoop("0:-90, -90:0, 0:90, 90:0");

    // The "far" hemisphere, defined using two pairs of antipodal points.
    farHemi = addLoop("90:0, 0:90, -90:0, 0:-90");

    // A spiral stripe that slightly over-wraps the equator.
    candyCane = addLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70");

    // A small clockwise loop in the northern & eastern hemisperes.
    smallNeCw = addLoop("35:20, 45:20, 40:25");

    // Loop around the north pole at 80 degrees.
    arctic80 = addLoop("80:-150, 80:-30, 80:90");

    // Loop around the south pole at 80 degrees.
    antarctic80 = addLoop("-80:120, -80:0, -80:-120");

    // A completely degenerate triangle along the equator that Sign()
    // considers to be CCW.
    lineTriangle = addLoop("0:1, 0:2, 0:3");

    // A nearly-degenerate CCW chevron near the equator with very long sides
    // (about 80 degrees).  Its area is less than 1e-640, which is too small
    // to represent in double precision.
    skinnyChevron = addLoop("0:0, -1e-320:80, 0:1e-320, 1e-320:80");

    // A diamond-shaped loop around the point 0:180.
    loopA = addLoop("0:178, -1:180, 0:-179, 1:-180");

    // Another diamond-shaped loop around the point 0:180.
    loopB = addLoop("0:179, -1:180, 0:-178, 1:-180");

    // The intersection of A and B.
    aIntersectB = addLoop("0:179, -1:180, 0:-179, 1:-180");

    // The union of A and B.
    aUnionB = addLoop("0:178, -1:180, 0:-178, 1:-180");

    // A minus B (concave).
    aMinusB = addLoop("0:178, -1:180, 0:179, 1:-180");

    // B minus A (concave).
    bMinusA = addLoop("0:-179, -1:180, 0:-178, 1:-180");

    // A shape gotten from A by adding a triangle to one edge, and
    // subtracting a triangle from the opposite edge.
    loopC = addLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180");

    // A shape gotten from A by adding a triangle to one edge, and
    // adding another triangle to the opposite edge.
    loopD = addLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180");

    //   3------------2
    //   |            |               ^
    //   |  7-8  b-c  |               |
    //   |  | |  | |  |      Latitude |
    //   0--6-9--a-d--1               |
    //   |  | |       |               |
    //   |  f-e       |               +----------->
    //   |            |                 Longitude
    //   4------------5
    //
    // Important: It is not okay to skip over collinear vertices when
    // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
    // uses symbolic perturbations to ensure that no three vertices are
    // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
    // collinear).  In other words, it is unpredictable (modulo knowing the
    // details of the symbolic perturbations) whether 0123 contains 06123,
    // for example.
    //
    // Loop E:  0,6,9,a,d,1,2,3
    // Loop F:  0,4,5,1,d,a,9,6
    // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
    // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
    // Loop I:  7,6,f,e,9,8
    loopE = addLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30");
    loopF = addLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34");
    loopG = addLoop(
        "0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30");
    loopH = addLoop(
        "0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30");
    loopI = addLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36");

    // Like loop_a, but the vertices are at leaf cell centers.
    S2Point[] snapped_loop_a_vertices = [
        S2CellId(makePointOrDie("0:178")).toS2Point(),
        S2CellId(makePointOrDie("-1:180")).toS2Point(),
        S2CellId(makePointOrDie("0:-179")).toS2Point(),
        S2CellId(makePointOrDie("1:-180")).toS2Point()];
    snappedLoopA = new S2Loop(snapped_loop_a_vertices);
  }

  S2Loop addLoop(string str) {
    return addLoop(makeLoopOrDie(str));
  }

  S2Loop addLoop(S2Loop loop) {
    allLoops ~= loop;
    return loop;
  }

  // Wrapper function that encodes "loop" into "encoder" using the private
  // EncodeCompressed() method.
  // void TestEncodeCompressed(const S2Loop& loop, int level, Encoder* encoder) {
  //   absl::FixedArray<S2XYZFaceSiTi> points(loop.num_vertices());
  //   loop.GetXYZFaceSiTiVertices(points.data());
  //   loop.EncodeCompressed(encoder, points.data(), level);
  // }

  // Wrapper function that decodes the contents of "encoder" into "loop" using
  // the private DecodeCompressed() method.
  // void TestDecodeCompressed(const Encoder& encoder, int level, S2Loop* loop) {
  //   Decoder decoder(encoder.base(), encoder.length());
  //   ASSERT_TRUE(loop->DecodeCompressed(&decoder, level));
  // }
}

private const S2LatLng kRectError = S2LatLngRectBounder.maxErrorForTests();

@("S2LoopTestBase.GetRectBound") unittest {
  auto t = new S2LoopTestBase();
  Assert.equal(t.empty.getRectBound().isEmpty(), true);
  Assert.equal(t.full.getRectBound().isFull(), true);
  Assert.equal(t.candyCane.getRectBound().lng().isFull(), true);
  Assert.lessThan(t.candyCane.getRectBound().latLo().degrees(), -20);
  Assert.greaterThan(t.candyCane.getRectBound().latHi().degrees(), 10);
  Assert.equal(t.smallNeCw.getRectBound().isFull(), true);
  Assert.equal(t.arctic80.getRectBound().approxEquals(
      new S2LatLngRect(S2LatLng.fromDegrees(80, -180), S2LatLng.fromDegrees(90, 180)), kRectError),
      true);
  Assert.equal(t.antarctic80.getRectBound().approxEquals(
      new S2LatLngRect(S2LatLng.fromDegrees(-90, -180), S2LatLng.fromDegrees(-80, 180)), kRectError),
      true);

  // Create a loop that contains the complement of the "arctic_80" loop.
  auto arctic_80_inv = t.arctic80.clone();
  arctic_80_inv.invert();
  // The highest latitude of each edge is attained at its midpoint.
  S2Point mid = 0.5 * (arctic_80_inv.vertex(0) + arctic_80_inv.vertex(1));
  Assert.approximately(
      arctic_80_inv.getRectBound().latHi().radians(),
      S2LatLng(mid).lat().radians(),
      kRectError.lat().radians());

  Assert.equal(t.southHemi.getRectBound().lng().isFull(), true);
  Assert.equal(t.southHemi.getRectBound().lat().approxEquals(
      R1Interval(-PI_2, 0), kRectError.lat().radians()), true);
}

private void rotate(ref S2Loop ptr) {
  S2Loop loop = ptr;
  S2Point[] vertices;
  for (int i = 1; i <= loop.numVertices(); ++i) {
    vertices ~= loop.vertex(i);
  }
  ptr = new S2Loop(vertices);
}

@("S2LoopTestBase.AreaConsistentWithTurningAngle") unittest {
  // Check that the area computed using GetArea() is consistent with the
  // turning angle of the loop computed using GetTurnAngle().  According to
  // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
  // minus its turning angle.
  auto t = new S2LoopTestBase();
  foreach (S2Loop loop; t.allLoops) {
    double area = loop.getArea();
    double gauss_area = 2 * M_PI - loop.getTurningAngle();
    // TODO(ericv): The error bound below is much larger than it should be.
    // Need to improve the error minimization analysis in S2::Area().
    Assert.notGreaterThan(fabs(area - gauss_area), 1e-9,
        "Failed loop: " ~ loop.toString()
        ~ "\nArea = " ~ area.to!string ~ ", Gauss Area = " ~ gauss_area.to!string);
  }
}

// TODO: Resume here, test crashes.
@("S2LoopTestBase.GetAreaConsistentWithSign") unittest {
  // Test that GetArea() returns an area near 0 for degenerate loops that
  // contain almost no points, and an area near 4*Pi for degenerate loops that
  // contain almost all points.
  Random rnd = S2Testing.rnd;
  enum int kMaxVertices = 6;
  for (int i = 0; i < 50; ++i) {
    int num_vertices = 3 + rnd.uniform(kMaxVertices - 3 + 1);
    // Repeatedly choose N vertices that are exactly on the equator until we
    // find some that form a valid loop.
    auto loop = new S2Loop();
    loop.s2DebugOverride(S2Debug.DISABLE);
    do {
      S2Point[] vertices;
      for (int j = 0; j < num_vertices; ++j) {
        // We limit longitude to the range [0, 90] to ensure that the loop is
        // degenerate (as opposed to following the entire equator).
        vertices ~= S2LatLng.fromRadians(0, rnd.randDouble() * M_PI_2).toS2Point();
      }
      loop.initialize(vertices);
    } while (!loop.isValid());
    bool ccw = loop.isNormalized();
    // TODO(ericv): The error bound below is much larger than it should be.
    // Need to improve the error minimization analysis in S2::Area().
    Assert.approximately(loop.getArea(), ccw ? 0 : 4 * M_PI, 2e-8,
        "Failed loop " ~ i.to!string ~ ": " ~ toString(loop));
    Assert.equal(!ccw, loop.contains(S2Point(0, 0, 1)));
  }
}

@("S2LoopTestBase.GetAreaAccuracy") unittest {
  // TODO(ericv): Test that GetArea() has an accuracy significantly better
  // than 1e-15 on loops whose area is small.
}

@("S2LoopTestBase.GetAreaAndCentroid") unittest {
  auto t = new S2LoopTestBase();
  Assert.equal(0.0, t.empty.getArea());
  Assert.equal(4 * M_PI, t.full.getArea());
  Assert.equal(S2Point(0, 0, 0), t.empty.getCentroid());
  Assert.equal(S2Point(0, 0, 0), t.full.getCentroid());

  Assert.approximately(t.northHemi.getArea(), 2 * M_PI, DOUBLE_ERR);
  Assert.notGreaterThan(t.eastHemi.getArea(), 2 * M_PI + 1e-12);
  Assert.notLessThan(t.eastHemi.getArea(), 2 * M_PI - 1e-12);

  // Construct spherical caps of random height, and approximate their boundary
  // with closely spaces vertices.  Then check that the area and centroid are
  // correct.

  for (int i = 0; i < 50; ++i) {
    // Choose a coordinate frame for the spherical cap.
    Vector3_d x, y, z;
    S2Testing.getRandomFrame(x, y, z);

    // Given two points at latitude phi and whose longitudes differ by dtheta,
    // the geodesic between the two points has a maximum latitude of
    // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
    // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
    //
    // We want to position the vertices close enough together so that their
    // maximum distance from the boundary of the spherical cap is kMaxDist.
    // Thus we want fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
    static const double kMaxDist = 1e-6;
    double height = 2 * S2Testing.rnd.randDouble();
    double phi = asin(1 - height);
    double max_dtheta = 2 * acos(tan(fabs(phi)) / tan(fabs(phi) + kMaxDist));
    max_dtheta = min(M_PI, max_dtheta);  // At least 3 vertices.

    S2Point[] vertices;
    for (double theta = 0; theta < 2 * M_PI; theta += S2Testing.rnd.randDouble() * max_dtheta) {
      vertices ~= cos(theta) * cos(phi) * x
          + sin(theta) * cos(phi) * y
          + sin(phi) * z;
    }
    auto loop = new S2Loop(vertices);
    double area = loop.getArea();
    S2Point centroid = loop.getCentroid();
    double expected_area = 2 * M_PI * height;
    Assert.notGreaterThan(fabs(area - expected_area), 2 * M_PI * kMaxDist);
    S2Point expected_centroid = expected_area * (1 - 0.5 * height) * z;
    Assert.notGreaterThan((centroid - expected_centroid).norm(), 2 * kMaxDist);
  }
}

// Check that the turning angle is *identical* when the vertex order is
// rotated, and that the sign is inverted when the vertices are reversed.
private void checkTurningAngleInvariants(in S2Loop loop) {
  double expected = loop.getTurningAngle();
  S2Loop loop_copy = loop.clone();
  for (int i = 0; i < loop.numVertices(); ++i) {
    loop_copy.invert();
    Assert.equal(loop_copy.getTurningAngle(), -expected);
    loop_copy.invert();
    rotate(loop_copy);
    Assert.equal(loop_copy.getTurningAngle(), expected);
  }
}

@("S2LoopTestBase.GetTurningAngle") unittest {
  auto t = new S2LoopTestBase();

  Assert.equal(t.empty.getTurningAngle(), 2 * M_PI);
  Assert.equal(t.full.getTurningAngle(), -2 * M_PI);

  Assert.approximately(t.northHemi3.getTurningAngle(), 0, 1e-15);
  checkTurningAngleInvariants(t.northHemi3);

  Assert.approximately(t.westHemi.getTurningAngle(), 0, 1e-15);
  checkTurningAngleInvariants(t.westHemi);

  // We don't have an easy way to estimate the turning angle of this loop, but
  // we can still check that the expected invariants hold.
  checkTurningAngleInvariants(t.candyCane);

  Assert.approximately(t.lineTriangle.getTurningAngle(), 2 * M_PI, DOUBLE_ERR);
  checkTurningAngleInvariants(t.lineTriangle);

  Assert.approximately(t.skinnyChevron.getTurningAngle(), 2 * M_PI, DOUBLE_ERR);
  checkTurningAngleInvariants(t.skinnyChevron);

  // Build a narrow spiral loop starting at the north pole.  This is designed
  // to test that the error in GetTurningAngle is linear in the number of
  // vertices even when the partial sum of the turning angles gets very large.
  // The spiral consists of two "arms" defining opposite sides of the loop.
  const int kArmPoints = 10000;    // Number of vertices in each "arm"
  const double kArmRadius = 0.01;  // Radius of spiral.
  auto vertices = new S2Point[](2 * kArmPoints);
  vertices[kArmPoints] = S2Point(0, 0, 1);
  for (int i = 0; i < kArmPoints; ++i) {
    double angle = (2 * M_PI / 3) * i;
    double x = cos(angle);
    double y = sin(angle);
    double r1 = i * kArmRadius / kArmPoints;
    double r2 = (i + 1.5) * kArmRadius / kArmPoints;
    vertices[kArmPoints - i - 1] = S2Point(r1 * x, r1 * y, 1).normalize();
    vertices[kArmPoints + i] = S2Point(r2 * x, r2 * y, 1).normalize();
  }
  // This is a pathological loop that contains many long parallel edges, and
  // takes tens of seconds to validate in debug mode.
  auto spiral = new S2Loop(vertices, S2Debug.DISABLE);

  // Check that GetTurningAngle() is consistent with GetArea() to within the
  // error bound of the former.  We actually use a tiny fraction of the
  // worst-case error bound, since the worst case only happens when all the
  // roundoff errors happen in the same direction and this test is not
  // designed to achieve that.  The error in GetArea() can be ignored for the
  // purposes of this test since it is generally much smaller.
  Assert.approximately(
      spiral.getTurningAngle(),
      2 * M_PI - spiral.getArea(),
      0.01 * spiral.getTurningAngleMaxError());
}

// Checks that if a loop is normalized, it doesn't contain a
// point outside of it, and vice versa.
static void checkNormalizeAndContains(S2Loop loop) {
  S2Point p = makePointOrDie("40:40");

  auto flip = loop.clone();
  flip.invert();
  Assert.equal(loop.isNormalized() ^ loop.contains(p), true);
  Assert.equal(flip.isNormalized() ^ flip.contains(p), true);

  Assert.equal(loop.isNormalized() ^ flip.isNormalized(), true);

  flip.normalize();
  Assert.equal(flip.contains(p), false);
}

@("S2LoopTestBase.NormalizedCompatibleWithContains") unittest {
  auto t = new S2LoopTestBase();

  checkNormalizeAndContains(t.lineTriangle);
  checkNormalizeAndContains(t.skinnyChevron);
}

@("S2LoopTestBase.Contains") unittest {
  auto t = new S2LoopTestBase();

  // Check the full and empty loops have the correct containment relationship
  // with the special "vertex" that defines them.
  Assert.equal(t.empty.contains(S2Loop.empty()[0]), false);
  Assert.equal(t.full.contains(S2Loop.full()[0]), true);

  Assert.equal(t.candyCane.contains(S2LatLng.fromDegrees(5, 71).toS2Point()), true);

  // Create copies of these loops so that we can change the vertex order.
  S2Loop north_copy = t.northHemi.clone();
  S2Loop south_copy = t.southHemi.clone();
  S2Loop west_copy = t.westHemi.clone();
  S2Loop east_copy = t.eastHemi.clone();
  for (int i = 0; i < 4; ++i) {
    Assert.equal(north_copy.contains(S2Point(0, 0, 1)), true);
    Assert.equal(north_copy.contains(S2Point(0, 0, -1)), false);
    Assert.equal(south_copy.contains(S2Point(0, 0, 1)), false);
    Assert.equal(south_copy.contains(S2Point(0, 0, -1)), true);
    Assert.equal(west_copy.contains(S2Point(0, 1, 0)), false);
    Assert.equal(west_copy.contains(S2Point(0, -1, 0)), true);
    Assert.equal(east_copy.contains(S2Point(0, 1, 0)), true);
    Assert.equal(east_copy.contains(S2Point(0, -1, 0)), false);
    rotate(north_copy);
    rotate(south_copy);
    rotate(east_copy);
    rotate(west_copy);
  }

  // This code checks each cell vertex is contained by exactly one of
  // the adjacent cells.
  for (int level = 0; level < 3; ++level) {
    S2Loop[] loops;
    S2Point[] loop_vertices;
    bool[S2Point] points;
    for (S2CellId id = S2CellId.begin(level);
         id != S2CellId.end(level); id = id.next()) {
      auto cell = new S2Cell(id);
      points[cell.getCenter()] = true;
      for (int k = 0; k < 4; ++k) {
        loop_vertices ~= cell.getVertex(k);
        points[cell.getVertex(k)] = true;
      }
      loops ~= new S2Loop(loop_vertices);
      loop_vertices.length = 0;
    }
    foreach (point; points.keys()) {
      int count = 0;
      foreach (loop; loops) {
        if (loop.contains(point)) ++count;
      }
      Assert.equal(count, 1);
    }
  }
}

@("S2Loop.ContainsMatchesCrossingSign") unittest {
  // This test demonstrates a former incompatibility between CrossingSign()
  // and Contains(const S2Point&).  It constructs an S2Cell-based loop L and
  // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
  // previously, Contains() returned false for both endpoints of E.
  //
  // The reason for the bug was that the loop bound was sometimes too tight.
  // The Contains() code for a0 bailed out early because a0 was found not to
  // be inside the bound of L.

  // Start with a cell that ends up producing the problem.
  const S2CellId cell_id = S2CellId(S2Point(1, 1, 1)).parent(21);

  S2Cell[4] children;
  (new S2Cell(cell_id)).subdivide(children);

  auto points = new S2Point[](4);
  for (int i = 0; i < 4; ++i) {
    // Note extra normalization. GetCenter() is already normalized.
    // The test results will no longer be inconsistent if the extra
    // Normalize() is removed.
    points[i] = children[i].getCenter().normalize();
  }

  auto loop = new S2Loop(points);

  // Get a vertex from a grandchild cell.
  // +---------------+---------------+
  // |               |               |
  // |    points[3]  |   points[2]   |
  // |       v       |       v       |
  // |       +-------+------ +       |
  // |       |       |       |       |
  // |       |       |       |       |
  // |       |       |       |       |
  // +-------+-------+-------+-------+
  // |       |       |       |       |
  // |       |    <----------------------- grandchild_cell
  // |       |       |       |       |
  // |       +-------+------ +       |
  // |       ^       |       ^       | <-- cell
  // | points[0]/a0  |     points[1] |
  // |               |               |
  // +---------------+---------------+
  auto grandchild_cell = new S2Cell(cell_id.child(0).child(2));
  const S2Point a0 = grandchild_cell.getVertex(0);

  // If this doesn't hold, the rest of the test is pointless.
  Assert.notEqual(points[0], a0,
      "This test depends on rounding errors that should make "
      ~ "a0 slightly different from points[0]"
      ~ "\npoints[0]:" ~ points[0].to!string
      ~ "\n       a0:" ~ a0.to!string);

  // The edge from a0 to the origin crosses one boundary.
  Assert.equal(crossingSign(a0, origin(), loop.vertex(0), loop.vertex(1)), -1);
  Assert.equal(crossingSign(a0, origin(), loop.vertex(1), loop.vertex(2)), 1);
  Assert.equal(crossingSign(a0, origin(), loop.vertex(2), loop.vertex(3)), -1);
  Assert.equal(crossingSign(a0, origin(), loop.vertex(3), loop.vertex(4)), -1);

  // Contains should return false for the origin, and true for a0.
  Assert.equal(loop.contains(origin()), false);
  Assert.equal(loop.contains(a0), true);

  // Since a0 is inside the loop, it should be inside the bound.
  S2LatLngRect bound = loop.getRectBound();
  Assert.equal(bound.contains(a0), true);
}

// Given a pair of loops where A contains B, check various identities.
private void checkOneNestedPair(S2Loop a, S2Loop b) {
  Assert.equal(a.contains(b), true);
  Assert.equal(a.boundaryEquals(b), b.contains(a));
  Assert.equal(a.intersects(b), !b.isEmpty());
  Assert.equal(b.intersects(a), !b.isEmpty());
}

// Given a pair of disjoint loops A and B, check various identities.
private void checkOneDisjointPair(S2Loop a, S2Loop b) {
  Assert.equal(a.intersects(b), false);
  Assert.equal(b.intersects(a), false);
  Assert.equal(b.isEmpty(), a.contains(b));
  Assert.equal(a.isEmpty(), b.contains(a));
}

// Given loops A and B whose union covers the sphere, check various identities.
private void checkOneCoveringPair(S2Loop a, S2Loop b) {
  Assert.equal(a.isFull(), a.contains(b));
  Assert.equal(b.isFull(), b.contains(a));
  auto a1 = a.clone();
  a1.invert();
  bool complementary = a1.boundaryEquals(b);
  Assert.equal(a.intersects(b), !complementary);
  Assert.equal(b.intersects(a), !complementary);
}

// Given loops A and B such that both A and its complement intersect both B
// and its complement, check various identities.
private void checkOneOverlappingPair(S2Loop a, S2Loop b) {
  Assert.equal(a.contains(b), false);
  Assert.equal(b.contains(a), false);
  Assert.equal(a.intersects(b), true);
  Assert.equal(b.intersects(a), true);
}

// Given a pair of loops where A contains B, test various identities
// involving A, B, and their complements.
private void checkNestedPair(S2Loop a, S2Loop b) {
  auto a1 = a.clone();
  auto b1 = b.clone();
  a1.invert();
  b1.invert();
  checkOneNestedPair(a, b);
  checkOneNestedPair(b1, a1);
  checkOneDisjointPair(a1, b);
  checkOneCoveringPair(a, b1);
}

// Given a pair of disjoint loops A and B, test various identities
// involving A, B, and their complements.
private void checkDisjointPair(S2Loop a, S2Loop b) {
  auto a1 = a.clone();
  a1.invert();
  checkNestedPair(a1, b);
}

// Given loops A and B whose union covers the sphere, test various identities
// involving A, B, and their complements.
private void checkCoveringPair(S2Loop a, S2Loop b) {
  auto b1 = b.clone();
  b1.invert();
  checkNestedPair(a, b1);
}

// Given loops A and B such that both A and its complement intersect both B
// and its complement, test various identities involving these four loops.
private void checkOverlappingPair(S2Loop a, S2Loop b) {
  auto a1 = a.clone();
  auto b1 = b.clone();
  a1.invert();
  b1.invert();
  checkOneOverlappingPair(a, b);
  checkOneOverlappingPair(a1, b1);
  checkOneOverlappingPair(a1, b);
  checkOneOverlappingPair(a, b1);
}

enum RelationFlags {
  CONTAINS =  0x01,  // A contains B
  CONTAINED = 0x02,  // B contains A
  DISJOINT =  0x04,  // A and B are disjoint (intersection is empty)
  COVERS =    0x08,  // (A union B) covers the entire sphere
}

// Verify the relationship between two loops A and B.  "flags" is the set of
// RelationFlags that apply.  "shared_edge" means that the loops share at
// least one edge (possibly reversed).
private void checkRelationWithDesc(
    S2Loop a, S2Loop b, int flags, bool shared_edge, string test_description) {
  logger.logTrace(test_description);
  if (flags & RelationFlags.CONTAINS) {
    checkNestedPair(a, b);
  }
  if (flags & RelationFlags.CONTAINED) {
    checkNestedPair(b, a);
  }
  if (flags & RelationFlags.COVERS) {
    checkCoveringPair(a, b);
  }
  if (flags & RelationFlags.DISJOINT) {
    checkDisjointPair(a, b);
  } else if (!(flags & (RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.COVERS))) {
    checkOverlappingPair(a, b);
  }
  if (!shared_edge
      && (flags & (RelationFlags.CONTAINS | RelationFlags.CONTAINED | RelationFlags.DISJOINT))) {
    Assert.equal(a.contains(b), a.containsNested(b));
  }
  // A contains the boundary of B if either A contains B, or the two loops
  // contain each other's boundaries and there are no shared edges (since at
  // least one such edge must be reversed, and therefore is not considered to
  // be contained according to the rules of CompareBoundary).
  int comparison = 0;
  if ((flags & RelationFlags.CONTAINS) || ((flags & RelationFlags.COVERS) && !shared_edge)) {
    comparison = 1;
  }
  // Similarly, A excludes the boundary of B if either A and B are disjoint,
  // or B contains A and there are no shared edges (since A is considered to
  // contain such edges according to the rules of CompareBoundary).
  if ((flags & RelationFlags.DISJOINT) || ((flags & RelationFlags.CONTAINED) && !shared_edge)) {
    comparison = -1;
  }
  // CompareBoundary requires that neither loop is empty.
  if (!a.isEmpty() && !b.isEmpty()) {
    Assert.equal(comparison, a.compareBoundary(b));
  }
}

private void checkRelation(S2Loop a, S2Loop b, int flags, bool shared_edge) {
  checkRelationWithDesc(a, b, flags, shared_edge, "args " ~ toString(a) ~ ", " ~ toString(b));
}

@("S2LoopTestBase.LoopRelations") unittest {
  auto t = new S2LoopTestBase();

  // Check full and empty relationships with normal loops and each other.
  checkRelation(
      t.full, t.full, RelationFlags.CONTAINS|RelationFlags.CONTAINED|RelationFlags.COVERS, true);
  checkRelation(t.full, t.northHemi, RelationFlags.CONTAINS|RelationFlags.COVERS, false);
  checkRelation(
      t.full, t.empty, RelationFlags.CONTAINS|RelationFlags.DISJOINT|RelationFlags.COVERS, false);
  checkRelation(t.northHemi, t.full, RelationFlags.CONTAINED|RelationFlags.COVERS, false);
  checkRelation(t.northHemi, t.empty, RelationFlags.CONTAINS|RelationFlags.DISJOINT, false);
  checkRelation(
      t.empty, t.full, RelationFlags.CONTAINED|RelationFlags.DISJOINT|RelationFlags.COVERS, false);
  checkRelation(t.empty, t.northHemi, RelationFlags.CONTAINED|RelationFlags.DISJOINT, false);
  checkRelation(t.empty, t.empty,
      RelationFlags.CONTAINS|RelationFlags.CONTAINED|RelationFlags.DISJOINT, false);

  checkRelation(t.northHemi, t.northHemi, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.northHemi, t.southHemi, RelationFlags.DISJOINT|RelationFlags.COVERS, true);
  checkRelation(t.northHemi, t.eastHemi, 0, false);
  checkRelation(t.northHemi, t.arctic80, RelationFlags.CONTAINS, false);
  checkRelation(t.northHemi, t.antarctic80, RelationFlags.DISJOINT, false);
  checkRelation(t.northHemi, t.candyCane, 0, false);

  // We can't compare north_hemi3 vs. north_hemi or south_hemi because the
  // result depends on the "simulation of simplicity" implementation details.
  checkRelation(t.northHemi3, t.northHemi3, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.northHemi3, t.eastHemi, 0, false);
  checkRelation(t.northHemi3, t.arctic80, RelationFlags.CONTAINS, false);
  checkRelation(t.northHemi3, t.antarctic80, RelationFlags.DISJOINT, false);
  checkRelation(t.northHemi3, t.candyCane, 0, false);

  checkRelation(t.southHemi, t.northHemi, RelationFlags.DISJOINT|RelationFlags.COVERS, true);
  checkRelation(t.southHemi, t.southHemi, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.southHemi, t.farHemi, 0, false);
  checkRelation(t.southHemi, t.arctic80, RelationFlags.DISJOINT, false);
  checkRelation(t.southHemi, t.antarctic80, RelationFlags.CONTAINS, false);
  checkRelation(t.southHemi, t.candyCane, 0, false);

  checkRelation(t.candyCane, t.northHemi, 0, false);
  checkRelation(t.candyCane, t.southHemi, 0, false);
  checkRelation(t.candyCane, t.arctic80, RelationFlags.DISJOINT, false);
  checkRelation(t.candyCane, t.antarctic80, RelationFlags.DISJOINT, false);
  checkRelation(t.candyCane, t.candyCane, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);

  checkRelation(t.nearHemi, t.westHemi, 0, false);

  checkRelation(t.smallNeCw, t.southHemi, RelationFlags.CONTAINS, false);
  checkRelation(t.smallNeCw, t.westHemi, RelationFlags.CONTAINS, false);

  checkRelation(t.smallNeCw, t.northHemi, RelationFlags.COVERS, false);
  checkRelation(t.smallNeCw, t.eastHemi, RelationFlags.COVERS, false);

  checkRelation(t.loopA, t.loopA, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.loopA, t.loopB, 0, false);
  checkRelation(t.loopA, t.aIntersectB, RelationFlags.CONTAINS, true);
  checkRelation(t.loopA, t.aUnionB, RelationFlags.CONTAINED, true);
  checkRelation(t.loopA, t.aMinusB, RelationFlags.CONTAINS, true);
  checkRelation(t.loopA, t.bMinusA, RelationFlags.DISJOINT, true);

  checkRelation(t.loopB, t.loopA, 0, false);
  checkRelation(t.loopB, t.loopB, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.loopB, t.aIntersectB, RelationFlags.CONTAINS, true);
  checkRelation(t.loopB, t.aUnionB, RelationFlags.CONTAINED, true);
  checkRelation(t.loopB, t.aMinusB, RelationFlags.DISJOINT, true);
  checkRelation(t.loopB, t.bMinusA, RelationFlags.CONTAINS, true);

  checkRelation(t.aIntersectB, t.loopA, RelationFlags.CONTAINED, true);
  checkRelation(t.aIntersectB, t.loopB, RelationFlags.CONTAINED, true);
  checkRelation(t.aIntersectB, t.aIntersectB, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.aIntersectB, t.aUnionB, RelationFlags.CONTAINED, false);
  checkRelation(t.aIntersectB, t.aMinusB, RelationFlags.DISJOINT, true);
  checkRelation(t.aIntersectB, t.bMinusA, RelationFlags.DISJOINT, true);

  checkRelation(t.aUnionB, t.loopA, RelationFlags.CONTAINS, true);
  checkRelation(t.aUnionB, t.loopB, RelationFlags.CONTAINS, true);
  checkRelation(t.aUnionB, t.aIntersectB, RelationFlags.CONTAINS, false);
  checkRelation(t.aUnionB, t.aUnionB, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.aUnionB, t.aMinusB, RelationFlags.CONTAINS, true);
  checkRelation(t.aUnionB, t.bMinusA, RelationFlags.CONTAINS, true);

  checkRelation(t.aMinusB, t.loopA, RelationFlags.CONTAINED, true);
  checkRelation(t.aMinusB, t.loopB, RelationFlags.DISJOINT, true);
  checkRelation(t.aMinusB, t.aIntersectB, RelationFlags.DISJOINT, true);
  checkRelation(t.aMinusB, t.aUnionB, RelationFlags.CONTAINED, true);
  checkRelation(t.aMinusB, t.aMinusB, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
  checkRelation(t.aMinusB, t.bMinusA, RelationFlags.DISJOINT, false);

  checkRelation(t.bMinusA, t.loopA, RelationFlags.DISJOINT, true);
  checkRelation(t.bMinusA, t.loopB, RelationFlags.CONTAINED, true);
  checkRelation(t.bMinusA, t.aIntersectB, RelationFlags.DISJOINT, true);
  checkRelation(t.bMinusA, t.aUnionB, RelationFlags.CONTAINED, true);
  checkRelation(t.bMinusA, t.aMinusB, RelationFlags.DISJOINT, false);
  checkRelation(t.bMinusA, t.bMinusA, RelationFlags.CONTAINS|RelationFlags.CONTAINED, true);
}

// Make sure the relations are correct if the loop crossing happens on
// two ends of a shared boundary segment.
@("S2LoopTestBase.LoopRelationsWhenSameExceptPiecesStickingOutAndIn") unittest {
  auto t = new S2LoopTestBase();

  checkRelation(t.loopA, t.loopC, 0, true);
  checkRelation(t.loopC, t.loopA, 0, true);
  checkRelation(t.loopA, t.loopD, RelationFlags.CONTAINED, true);
  checkRelation(t.loopD, t.loopA, RelationFlags.CONTAINS, true);
  checkRelation(t.loopE, t.loopF, RelationFlags.DISJOINT, true);
  checkRelation(t.loopE, t.loopG, RelationFlags.CONTAINS, true);
  checkRelation(t.loopE, t.loopH, 0, true);
  checkRelation(t.loopE, t.loopI, 0, false);
  checkRelation(t.loopF, t.loopG, RelationFlags.DISJOINT, true);
  checkRelation(t.loopF, t.loopH, 0, true);
  checkRelation(t.loopF, t.loopI, 0, false);
  checkRelation(t.loopG, t.loopH, RelationFlags.CONTAINED, true);
  checkRelation(t.loopH, t.loopG, RelationFlags.CONTAINS, true);
  checkRelation(t.loopG, t.loopI, RelationFlags.DISJOINT, true);
  checkRelation(t.loopH, t.loopI, RelationFlags.CONTAINS, true);
}

private S2Loop makeCellLoop(S2CellId begin, S2CellId end) {
  // Construct a CCW polygon whose boundary is the union of the cell ids
  // in the range [begin, end).  We add the edges one by one, removing
  // any edges that are already present in the opposite direction.

  bool[S2Point][S2Point] edges;
  for (S2CellId id = begin; id != end; id = id.next()) {
    auto cell = new S2Cell(id);
    for (int k = 0; k < 4; ++k) {
      S2Point a = cell.getVertex(k);
      S2Point b = cell.getVertex(k + 1);
      if (b in edges && a in edges[b]) {
        edges[b].remove(a);
        if (edges[b].length == 0)
          edges.remove(b);
      } else {
        edges[a][b] = true;
      }
    }
  }

  // The remaining edges form a single loop.  We simply follow it starting
  // at an arbitrary vertex and build up a list of vertices.

  S2Point[] vertices;
  S2Point p = edges.keys[0];
  while (edges.length > 0) {
    enforce(edges[p].length == 1,
        "Error, found " ~ edges[p].length.to!string ~ " edges for point " ~ p.to!string);
    S2Point next = edges[p].keys[0];
    vertices ~= p;
    edges.remove(p);
    p = next;
  }

  return new S2Loop(vertices);
}

@("S2Loop.LoopRelations2") unittest {
  // Construct polygons consisting of a sequence of adjacent cell ids
  // at some fixed level.  Comparing two polygons at the same level
  // ensures that there are no T-vertices.
  auto rnd = S2Testing.rnd;
  foreach (int iter; 0..1000) {
    auto begin = S2CellId(rnd.rand64() | 1);
    if (!begin.isValid()) continue;
    begin = begin.parent(rnd.uniform(S2CellId.MAX_LEVEL));
    S2CellId a_begin = begin.advance(rnd.skewed(6));
    S2CellId a_end = a_begin.advance(rnd.skewed(6) + 1);
    S2CellId b_begin = begin.advance(rnd.skewed(6));
    S2CellId b_end = b_begin.advance(rnd.skewed(6) + 1);
    if (!a_end.isValid() || !b_end.isValid()) continue;

    S2Loop a = makeCellLoop(a_begin, a_end);
    S2Loop b = makeCellLoop(b_begin, b_end);
    if (a && b) {
      bool contained = (a_begin <= b_begin && b_end <= a_end);
      bool intersects = (a_begin < b_end && b_begin < a_end);
      logger.logTrace("Checking ", a.numVertices(), " vs. ", b.numVertices(), ", contained = ",
          contained, ", intersects = ", intersects);
      Assert.equal(a.contains(b), contained);
      Assert.equal(a.intersects(b), intersects);
    } else {
      logger.logTrace("MakeCellLoop failed to create a loop.");
    }
  }
}

@("S2Loop.BoundsForLoopContainment") unittest {
  import s2.s2edge_distances : interpolate;
  import s2.s2predicates : sign;
  // To reliably test whether one loop contains another, the bounds of the
  // outer loop are expanded slightly.  This test constructs examples where
  // this expansion is necessary and verifies that it is sufficient.
  auto rnd = S2Testing.rnd;
  foreach (int iter; 0..1000) {
    // We construct a triangle ABC such that A,B,C are nearly colinear, B is
    // the point of maximum latitude, and the edge AC passes very slightly
    // below B (i.e., ABC is CCW).
    S2Point b = (S2Testing.randomPoint() + S2Point(0, 0, 1)).normalize();
    S2Point v = b.crossProd(S2Point(0, 0, 1)).normalize();
    S2Point a = interpolate(rnd.randDouble(), -v, b);
    S2Point c = interpolate(rnd.randDouble(), b, v);
    if (sign(a, b, c) < 0) {
      --iter; continue;
    }
    // Now construct another point D directly below B, and create two loops
    // ABCD and ACD.
    S2Point d = S2Point(b.x(), b.y(), 0).normalize();
    S2Point[] vertices = [ c, d, a, b ];  // Reordered for convenience
    auto outer = new S2Loop(vertices[0..4]);
    auto inner = new S2Loop(vertices[0..3]);
    // Now because the bounds calculation is less accurate when the maximum is
    // attained along an edge (rather than at a vertex), sometimes the inner
    // loop will have a *larger* bounding box than the outer loop.  We look
    // only for those cases.
    if (outer.getRectBound().contains(inner.getRectBound())) {
      --iter; continue;
    }
    Assert.equal(outer.contains(inner), true);
  }
}

void debugDumpCrossings(S2Loop loop) {
  import s2.s2pointutil : ortho, origin;
  import s2.s2predicates : sign, orderedCCW;
  import s2.s2edge_crossings : edgeOrVertexCrossing;

  // This function is useful for debugging.
  logger.logTrace("Ortho(v1): ", ortho(loop.vertex(1)));
  writefln("Contains(kOrigin): %d", loop.contains(origin()));
  for (int i = 1; i <= loop.numVertices(); ++i) {
    S2Point a = ortho(loop.vertex(i));
    S2Point b = loop.vertex(i-1);
    S2Point c = loop.vertex(i+1);
    S2Point o = loop.vertex(i);
    writefln("Vertex %d: [%.17g, %.17g, %.17g], %d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %d",
           i, loop.vertex(i).x(), loop.vertex(i).y(), loop.vertex(i).z(),
           i - 1, i, sign(b, o, a),
           i + 1, i, i - 1, sign(c, o, b),
           i, i + 1, sign(a, o, c),
           orderedCCW(a, b, c, o));
  }
  for (int i = 0; i < loop.numVertices() + 2; ++i) {
    S2Point orig = origin();
    S2Point dest;
    if (i < loop.numVertices()) {
      dest = loop.vertex(i);
      writef("Origin->%d crosses:", i);
    } else {
      dest = S2Point(0, 0, 1);
      if (i == loop.numVertices() + 1) orig = loop.vertex(1);
      writef("Case %d:", i);
    }
    for (int j = 0; j < loop.numVertices(); ++j) {
      writef(" %d", edgeOrVertexCrossing(orig, dest, loop.vertex(j), loop.vertex(j+1)));
    }
    writeln();
  }
  for (int i = 0; i <= 2; i += 2) {
    writef("Origin->v1 crossing v%d->v1: ", i);
    S2Point a = ortho(loop.vertex(1));
    S2Point b = loop.vertex(i);
    S2Point c = origin();
    S2Point o = loop.vertex(1);
    writefln("%d1R=%d, M1%d=%d, R1M=%d, crosses: %d",
           i, sign(b, o, a),
           i, sign(c, o, b),
           sign(a, o, c),
           edgeOrVertexCrossing(c, o, b, a));
  }
}

private void checkNear(string a_str, string b_str, S1Angle max_error, bool expected) {
  auto a = makeLoopOrDie(a_str);
  auto b = makeLoopOrDie(b_str);
  Assert.equal(a.boundaryNear(b, max_error), expected);
  Assert.equal(b.boundaryNear(a, max_error), expected);
}

@("S2Loop.BoundaryNear") unittest {
  S1Angle degree = S1Angle.fromDegrees(1);

  checkNear(
      "0:0, 0:10, 5:5",
      "0:0.1, -0.1:9.9, 5:5.2",
      0.5 * degree, true);
  checkNear(
      "0:0, 0:3, 0:7, 0:10, 3:7, 5:5",
      "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1",
      S1Angle.fromRadians(1e-3), true);

  // All vertices close to some edge, but not equivalent.
  checkNear(
      "0:0, 0:2, 2:2, 2:0",
      "0:0, 1.9999:1, 0:2, 2:2, 2:0",
      0.5 * degree, false);

  // Two triangles that backtrack a bit on different edges.  A simple
  // greedy matching algorithm would fail on this example.
  string t1 = "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, "
      ~ "2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2";
  string t2 = "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, "
      ~ "0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4";
  checkNear(t1, t2, 1.5 * degree, true);
  checkNear(t1, t2, 0.5 * degree, false);
}

private void checkIdentical(S2Loop loop, S2Loop loop2) {
  import s2.s2pointutil : origin;

  Assert.equal(loop.depth(), loop2.depth());
  Assert.equal(loop.numVertices(), loop2.numVertices());
  for (int i = 0; i < loop.numVertices(); ++i) {
    Assert.equal(loop.vertex(i), loop2.vertex(i));
  }
  Assert.equal(loop.isEmpty(), loop2.isEmpty());
  Assert.equal(loop.isFull(), loop2.isFull());
  Assert.equal(loop.depth(), loop2.depth());
  Assert.equal(loop.isNormalized(), loop2.isNormalized());
  Assert.equal(loop.contains(origin()), loop2.contains(origin()));
  Assert.equal(loop.getRectBound(), loop2.getRectBound());
}

/+ TODO: Convert when encode/decode are implemented.
private void TestEncodeDecode(const S2Loop& loop) {
  Encoder encoder;
  loop.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2Loop loop2;
  loop2.set_s2debug_override(loop.s2debug_override());
  ASSERT_TRUE(loop2.Decode(&decoder));
  CheckIdentical(loop, loop2);
}

TEST(S2Loop, EncodeDecode) {
  unique_ptr<S2Loop> l(s2textformat::MakeLoop("30:20, 40:20, 39:43, 33:35"));
  l->set_depth(3);
  TestEncodeDecode(*l);

  S2Loop empty(S2Loop::kEmpty());
  TestEncodeDecode(empty);
  S2Loop full(S2Loop::kFull());
  TestEncodeDecode(full);

  S2Loop uninitialized;
  TestEncodeDecode(uninitialized);
}
+/

private void checkEmptyFullSnapped(S2Loop loop, int level) {
  Assert.equal(loop.isEmptyOrFull(), true);
  auto cellid = S2CellId(loop.vertex(0)).parent(level);
  S2Point[] vertices = [cellid.toS2Point()];
  auto loop2 = new S2Loop(vertices);
  Assert.equal(loop.boundaryEquals(loop2), true);
  Assert.equal(loop.boundaryApproxEquals(loop2), true);
  Assert.equal(loop.boundaryNear(loop2), true);
}

// Test converting the empty/full loops to S2LatLng representations.  (We
// don't bother testing E5/E6/E7 because that test is less demanding.)
private void checkEmptyFullLatLng(S2Loop loop) {
  Assert.equal(loop.isEmptyOrFull(), true);
  S2Point[] vertices = [S2LatLng(loop.vertex(0)).toS2Point()];
  auto loop2 = new S2Loop(vertices);
  Assert.equal(loop.boundaryEquals(loop2), true);
  Assert.equal(loop.boundaryApproxEquals(loop2), true);
  Assert.equal(loop.boundaryNear(loop2), true);
}

private void checkEmptyFullConversions(S2Loop loop) {
  checkEmptyFullSnapped(loop, S2CellId.MAX_LEVEL);
  checkEmptyFullSnapped(loop, 1);  // Worst case for approximation
  checkEmptyFullSnapped(loop, 0);
  checkEmptyFullLatLng(loop);
}

@("S2Loop.EmptyFullLossyConversions") unittest {
  // Verify that the empty and full loops can be encoded lossily.
  auto empty = new S2Loop(S2Loop.empty());
  checkEmptyFullConversions(empty);

  auto full = new S2Loop(S2Loop.full());
  checkEmptyFullConversions(full);
}

/+ TODO: Add when encode/decode are supported.

TEST(S2Loop, EncodeDecodeWithinScope) {
  unique_ptr<S2Loop> l(s2textformat::MakeLoop("30:20, 40:20, 39:43, 33:35"));
  l->set_depth(3);
  Encoder encoder;
  l->Encode(&encoder);
  Decoder decoder1(encoder.base(), encoder.length());

  // Initialize the loop using DecodeWithinScope and check that it is the
  // same as the original loop.
  S2Loop loop1;
  ASSERT_TRUE(loop1.DecodeWithinScope(&decoder1));
  EXPECT_TRUE(l->BoundaryEquals(&loop1));
  EXPECT_EQ(l->depth(), loop1.depth());
  EXPECT_EQ(l->GetRectBound(), loop1.GetRectBound());

  // Initialize the same loop using Init with a vector of vertices, and
  // check that it doesn't deallocate the original memory.
  vector<S2Point> vertices = {loop1.vertex(0), loop1.vertex(2),
                              loop1.vertex(3)};
  loop1.Init(vertices);
  Decoder decoder2(encoder.base(), encoder.length());
  S2Loop loop2;
  ASSERT_TRUE(loop2.DecodeWithinScope(&decoder2));
  EXPECT_TRUE(l->BoundaryEquals(&loop2));
  EXPECT_EQ(l->vertex(1), loop2.vertex(1));
  EXPECT_NE(loop1.vertex(1), loop2.vertex(1));

  // Initialize loop2 using Decode with a decoder on different data.
  // Check that the original memory is not deallocated or overwritten.
  unique_ptr<S2Loop> l2(s2textformat::MakeLoop("30:40, 40:75, 39:43, 80:35"));
  l2->set_depth(2);
  Encoder encoder2;
  l2->Encode(&encoder2);
  Decoder decoder3(encoder2.base(), encoder2.length());
  ASSERT_TRUE(loop2.Decode(&decoder3));
  Decoder decoder4(encoder.base(), encoder.length());
  ASSERT_TRUE(loop1.DecodeWithinScope(&decoder4));
  EXPECT_TRUE(l->BoundaryEquals(&loop1));
  EXPECT_EQ(l->vertex(1), loop1.vertex(1));
  EXPECT_NE(loop1.vertex(1), loop2.vertex(1));
}

TEST_F(S2LoopTestBase, FourVertexCompressedLoopRequires36Bytes) {
  Encoder encoder;
  TestEncodeCompressed(*snapped_loop_a_, S2CellId::kMaxLevel, &encoder);

  // 1 byte for num_vertices
  // 1 byte for origin_inside and boolean indicating we did not
  //   encode the bound
  // 1 byte for depth
  // Vertices:
  // 1 byte for faces
  // 8 bytes for each vertex.
  // 1 byte indicating that there is no unsnapped vertex.
  EXPECT_EQ(37, encoder.length());
}

TEST_F(S2LoopTestBase, CompressedEncodedLoopDecodesApproxEqual) {
  unique_ptr<S2Loop> loop(snapped_loop_a_->Clone());
  loop->set_depth(3);

  Encoder encoder;
  TestEncodeCompressed(*loop, S2CellId::kMaxLevel, &encoder);
  S2Loop decoded_loop;
  TestDecodeCompressed(encoder, S2CellId::kMaxLevel, &decoded_loop);
  CheckIdentical(*loop, decoded_loop);
}

+/

// This test checks that S2Loops created directly from S2Cells behave
// identically to S2Loops created from the vertices of those cells; this
// previously was not the case, because S2Cells calculate their bounding
// rectangles slightly differently, and S2Loops created from them just copied
// the S2Cell bounds.
@("S2Loop.S2CellConstructorAndContains") unittest {
  auto cell = new S2Cell(S2CellId(S2LatLng.fromE6(40565459, -74645276)));
  auto cell_as_loop = new S2Loop(cell);

  S2Point[] vertices;
  for (int i = 0; i < cell_as_loop.numVertices(); ++i) {
    vertices ~= cell_as_loop.vertex(i);
  }
  auto loop_copy = new S2Loop(vertices);
  Assert.equal(loop_copy.contains(cell_as_loop), true);
  Assert.equal(cell_as_loop.contains(loop_copy), true);

  // Demonstrates the reason for this test; the cell bounds are more
  // conservative than the resulting loop bounds.
  Assert.equal(loop_copy.getRectBound().contains(cell.getRectBound()), false);
}

// Construct a loop using s2textformat::MakeLoop(str) and check that it
// produces a validation error that includes "snippet".
private void checkLoopIsInvalid(string str, string snippet) {
  auto loop = makeLoopOrDie(str);
  S2Error error;
  Assert.equal(loop.findValidationError(error), true);
  Assert.equal(error.text().find(snippet).empty(), false);
}

private void checkLoopIsInvalid(in S2Point[] points, string snippet) {
  auto l = new S2Loop(points);
  S2Error error;
  Assert.equal(l.findValidationError(error), true);
  Assert.equal(error.text().find(snippet).empty(), false);
}

/+ TODO: Determine proper debug behavior.
@("S2Loop.IsValidDetectsInvalidLoops") unittest {
  //google::FlagSaver flag_saver;
  //FLAGS_s2debug = false;

  // Not enough vertices.  Note that all single-vertex loops are valid; they
  // are interpreted as being either "empty" or "full".
  writeln("Test 0:");
  checkLoopIsInvalid("", "at least 3 vertices");
  writeln("Test 1:");
  checkLoopIsInvalid("20:20, 21:21", "at least 3 vertices");

  // There is a degenerate edge
  writeln("Test 2:");
  checkLoopIsInvalid("20:20, 20:20, 20:21", "degenerate");
  writeln("Test 3:");
  checkLoopIsInvalid("20:20, 20:21, 20:20", "degenerate");

  // There is a duplicate vertex
  writeln("Test 4:");
  checkLoopIsInvalid("20:20, 21:21, 21:20, 20:20, 20:21", "duplicate vertex");

  // Some edges cross
  writeln("Test 5:");
  checkLoopIsInvalid("20:20, 21:21, 21:20.5, 21:20, 20:21", "crosses");

  // Points with non-unit length (triggers DCHECK failure in debug)
  //EXPECT_DEBUG_DEATH(
  writeln("Test 6:");
  checkLoopIsInvalid([S2Point(2, 0, 0), S2Point(0, 1, 0), S2Point(0, 0, 1)], "unit length");
  //"IsUnitLength");

  // Adjacent antipodal vertices
  writeln("Test 7:");
  checkLoopIsInvalid([S2Point(1, 0, 0), S2Point(-1, 0, 0), S2Point(0, 0, 1)], "antipodal");
}
+/

// Helper function for testing the distance methods.  "boundary_x" is the
// expected result of projecting "x" onto the loop boundary.  For convenience
// it can be set to S2Point() to indicate that (boundary_x == x).
private void checkDistanceMethods(S2Loop loop, in S2Point x, S2Point boundary_x) {
  // This error is not guaranteed by the implementation but is okay for tests.
  const S1Angle kMaxError = S1Angle.fromRadians(1e-15);

  if (boundary_x == S2Point()) boundary_x = x;
  Assert.notGreaterThan(S1Angle(boundary_x, loop.projectToBoundary(x)), kMaxError);

  if (loop.isEmptyOrFull()) {
    Assert.equal(S1Angle.infinity(), loop.getDistanceToBoundary(x));
  } else {
    // EXPECT_NEAR only works with doubles.
    Assert.approximately(S1Angle(x, boundary_x).degrees(),
        loop.getDistanceToBoundary(x).degrees(), kMaxError.degrees());
  }
  if (loop.contains(x)) {
    Assert.equal(S1Angle.zero(), loop.getDistance(x));
    Assert.equal(x, loop.project(x));
  } else {
    Assert.equal(loop.getDistanceToBoundary(x), loop.getDistance(x));
    Assert.equal(loop.projectToBoundary(x), loop.project(x));
  }
}

@("S2LoopTestBase.DistanceMethods") unittest {
  auto t = new S2LoopTestBase();
  // S2ClosestEdgeQuery is already tested, so just do a bit of sanity checking.

  // The empty and full loops don't have boundaries.
  checkDistanceMethods(t.empty, S2Point(0, 1, 0), S2Point());
  checkDistanceMethods(t.full, S2Point(0, 1, 0), S2Point());

  // A CCW square around the S2LatLng point (0,0).  Note that because lines of
  // latitude are curved on the sphere, it is not straightforward to project
  // points onto any edge except along the equator.  (The equator is the only
  // line of latitude that is also a geodesic.)
  auto square = makeLoopOrDie("-1:-1, -1:1, 1:1, 1:-1");
  Assert.equal(square.isNormalized(), true);

  // A vertex.
  checkDistanceMethods(square, S2LatLng.fromDegrees(1, -1).toS2Point(), S2Point());
  // A point on one of the edges.
  checkDistanceMethods(square, S2LatLng.fromDegrees(0.5, 1).toS2Point(), S2Point());
  // A point inside the square.
  checkDistanceMethods(square, S2LatLng.fromDegrees(0, 0.5).toS2Point(),
      S2LatLng.fromDegrees(0, 1).toS2Point());
  // A point outside the square that projects onto an edge.
  checkDistanceMethods(square, S2LatLng.fromDegrees(0, -2).toS2Point(),
      S2LatLng.fromDegrees(0, -1).toS2Point());
  // A point outside the square that projects onto a vertex.
  checkDistanceMethods(square, S2LatLng.fromDegrees(3, 4).toS2Point(),
      S2LatLng.fromDegrees(1, 1).toS2Point());
}

@("S2LoopTestBase.MakeRegularLoop") unittest {
  S2Point center = S2LatLng.fromDegrees(80, 135).toS2Point();
  S1Angle radius = S1Angle.fromDegrees(20);
  auto loop = S2Loop.makeRegularLoop(center, radius, 4);

  Assert.equal(loop.numVertices(), 4);
  S2Point p0 = loop.vertex(0);
  S2Point p1 = loop.vertex(1);
  S2Point p2 = loop.vertex(2);
  S2Point p3 = loop.vertex(3);
  // Make sure that the radius is correct.
  Assert.approximately(20.0, S2LatLng(center).getDistance(S2LatLng(p0)).degrees(), DOUBLE_ERR);
  Assert.approximately(20.0, S2LatLng(center).getDistance(S2LatLng(p1)).degrees(), DOUBLE_ERR);
  Assert.approximately(20.0, S2LatLng(center).getDistance(S2LatLng(p2)).degrees(), DOUBLE_ERR);
  Assert.approximately(20.0, S2LatLng(center).getDistance(S2LatLng(p3)).degrees(), DOUBLE_ERR);
  // Make sure that all angles of the polygon are the same.
  Assert.approximately((p1 - p0).angle(p3 - p0), M_PI_2, DOUBLE_ERR);
  Assert.approximately((p2 - p1).angle(p0 - p1), M_PI_2, DOUBLE_ERR);
  Assert.approximately((p3 - p2).angle(p1 - p2), M_PI_2, DOUBLE_ERR);
  Assert.approximately((p0 - p3).angle(p2 - p3), M_PI_2, DOUBLE_ERR);
  // Make sure that all edges of the polygon have the same length.
  Assert.approximately(
      S2LatLng(p0).getDistance(S2LatLng(p1)).degrees(), 27.990890717782829, DOUBLE_ERR);
  Assert.approximately(
      S2LatLng(p1).getDistance(S2LatLng(p2)).degrees(), 27.990890717782829, DOUBLE_ERR);
  Assert.approximately(
      S2LatLng(p2).getDistance(S2LatLng(p3)).degrees(), 27.990890717782829, DOUBLE_ERR);
  Assert.approximately(
      S2LatLng(p3).getDistance(S2LatLng(p0)).degrees(), 27.990890717782829, DOUBLE_ERR);

  // Check actual coordinates. This may change if we switch the algorithm
  // intentionally.
  Assert.approximately(S2LatLng(p0).lat().degrees(), 62.162880741097204, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p0).lng().degrees(), 103.11051028343407, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p1).lat().degrees(), 61.955157772928345, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p1).lng().degrees(), 165.25681963683536, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p2).lat().degrees(), 75.139812547718478, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p2).lng().degrees(), -119.13042521187423, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p3).lat().degrees(), 75.524190079054392, DOUBLE_ERR);
  Assert.approximately(S2LatLng(p3).lng().degrees(), 26.392175948257943, DOUBLE_ERR);
}

@("S2LoopShape.Basic") unittest {
  S2Loop loop = makeLoopOrDie("0:0, 0:1, 1:0");
  auto shape = new S2Loop.Shape(loop);
  Assert.equal(shape.loop(), loop);
  Assert.equal(shape.numEdges(), 3);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, 3);
  auto edge2 = shape.edge(2);
  Assert.equal(.toString(edge2.v0), "1:0");
  Assert.equal(.toString(edge2.v1), "0:0");
  Assert.equal(shape.dimension(), 2);
  Assert.equal(shape.hasInterior(), true);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2LoopShape.EmptyLoop") unittest {
  auto loop = new S2Loop(S2Loop.empty());
  auto shape = new S2Loop.Shape(loop);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
}
