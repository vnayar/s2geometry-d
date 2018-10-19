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

module s2.s2polyline_test;

import s2.s2polyline;

import s2.logger;
import s2.s1angle;
import s2.s2cell;
import s2.s2debug;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2pointutil : approxEquals;
import s2.s2pointutil;
import s2.s2testing;
import s2.s2text_format;
import s2.util.math.vector;

import fluent.asserts;

import std.math;
import std.range;
import std.conv;

// TODO: Add then Encoder/Decoder is completed.
/+
// Wraps s2textformat::MakePolyline in order to test Encode/Decode.
unique_ptr<S2Polyline> MakePolyline(const string& str) {
  unique_ptr<S2Polyline> polyline(s2textformat::MakePolyline(str));
  Encoder encoder;
  polyline->Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  auto decoded_polyline = make_unique<S2Polyline>();
  decoded_polyline->Decode(&decoder);
  return decoded_polyline;
}
+/

@("S2Polyline.Basic") unittest {
  S2Point[] vertices;
  auto emptyLine = new S2Polyline(vertices);
  Assert.equal(S2LatLngRect.empty(), emptyLine.getRectBound());
  emptyLine.reverse();
  Assert.equal(emptyLine.numVertices(), 0);

  S2LatLng[] latlngs = [
      S2LatLng.fromDegrees(0.0, 0.0),
      S2LatLng.fromDegrees(0.0, 90.0),
      S2LatLng.fromDegrees(0.0, 180.0)];
  auto semi_equator = new S2Polyline(latlngs);
  Assert.equal(approxEquals(semi_equator.interpolate(0.5), S2Point(0, 1, 0)), true);
  semi_equator.reverse();
  Assert.equal(S2Point(1, 0, 0), semi_equator.vertex(2));
}

@("S2Polyline.GetLengthAndCentroid") unittest {
  // Construct random great circles and divide them randomly into segments.
  // Then make sure that the length and centroid are correct.  Note that
  // because of the way the centroid is computed, it does not matter how
  // we split the great circle into segments.

  for (int i = 0; i < 100; ++i) {
    // Choose a coordinate frame for the great circle.
    Vector3_d x, y, z;
    S2Testing.getRandomFrame(x, y, z);

    S2Point[] vertices;
    for (double theta = 0; theta < 2 * PI;
         theta += pow(S2Testing.rnd.randDouble(), 10.0)) {
      S2Point p = cos(theta) * x + sin(theta) * y;
      if (vertices.empty() || p != vertices.back())
        vertices ~= p;
    }
    // Close the circle.
    vertices ~= vertices[0];
    auto line = new S2Polyline(vertices);
    S1Angle length = line.getLength();
    Assert.notGreaterThan(fabs(length.radians() - 2 * PI), 2e-14);
    S2Point centroid = line.getCentroid();
    Assert.notGreaterThan(centroid.norm(), 2e-14);
  }
}

@("S2Polyline.MayIntersect") unittest {
  S2Point[] vertices = [S2Point(1, -1.1, 0.8).normalize(), S2Point(1, -0.8, 1.1).normalize()];
  auto line = new S2Polyline(vertices);
  for (int face = 0; face < 6; ++face) {
    S2Cell cell = S2Cell.fromFace(face);
    Assert.equal((face & 1) == 0, line.mayIntersect(cell));
  }
}

@("S2Polyline.Interpolate") unittest {
  S2Point[] vertices = [
      S2Point(1, 0, 0),
      S2Point(0, 1, 0),
      S2Point(0, 1, 1).normalize(),
      S2Point(0, 0, 1)];
  auto line = new S2Polyline(vertices);
  Assert.equal(vertices[0], line.interpolate(-0.1));
  Assert.equal(
      approxEquals(line.interpolate(0.1), S2Point(1, tan(0.2 * PI / 2), 0).normalize()), true);
  Assert.equal(approxEquals(line.interpolate(0.25), S2Point(1, 1, 0).normalize()), true);
  Assert.equal(vertices[1], line.interpolate(0.5));
  Assert.equal(approxEquals(vertices[2], line.interpolate(0.75)), true);
  int next_vertex;
  Assert.equal(vertices[0], line.getSuffix(-0.1, next_vertex));
  Assert.equal(next_vertex, 1);
  Assert.equal(approxEquals(vertices[2], line.getSuffix(0.75, next_vertex)), true);
  Assert.equal(next_vertex, 3);
  Assert.equal(vertices[3], line.getSuffix(1.1, next_vertex));
  Assert.equal(next_vertex, 4);

  // Check the case where the interpolation fraction is so close to 1 that
  // the interpolated point is identical to the last vertex.
  vertices = [
      S2Point(1, 1, 1).normalize(),
      S2Point(1, 1, 1 + 1e-15).normalize(),
      S2Point(1, 1, 1 + 2e-15).normalize()];
  auto short_line = new S2Polyline(vertices);
  Assert.equal(vertices[2], short_line.getSuffix(1.0 - 2e-16, next_vertex));
  Assert.equal(next_vertex, 3);
}

@("S2Polyline.UnInterpolate") unittest {
  S2Point[] vertices = [S2Point(1, 0, 0)];
  auto point_line = new S2Polyline(vertices);
  Assert.approximately(point_line.unInterpolate(S2Point(0, 1, 0), 1), 0.0, DOUBLE_ERR);

  vertices ~= S2Point(0, 1, 0);
  vertices ~= S2Point(0, 1, 1).normalize();
  vertices ~= S2Point(0, 0, 1);
  auto line = new S2Polyline(vertices);

  S2Point interpolated;
  int next_vertex;
  interpolated = line.getSuffix(-0.1, next_vertex);
  Assert.approximately(line.unInterpolate(interpolated, next_vertex), 0.0, DOUBLE_ERR);
  interpolated = line.getSuffix(0.0, next_vertex);
  Assert.approximately(line.unInterpolate(interpolated, next_vertex), 0.0, DOUBLE_ERR);
  interpolated = line.getSuffix(0.5, next_vertex);
  Assert.approximately(line.unInterpolate(interpolated, next_vertex), 0.5, DOUBLE_ERR);
  interpolated = line.getSuffix(0.75, next_vertex);
  Assert.approximately(line.unInterpolate(interpolated, next_vertex), 0.75, DOUBLE_ERR);
  interpolated = line.getSuffix(1.1, next_vertex);
  Assert.approximately(line.unInterpolate(interpolated, next_vertex), 1.0, DOUBLE_ERR);

  // Check that the return value is clamped to 1.0.
  Assert.approximately(
      line.unInterpolate(S2Point(0, 1, 0), cast(int) vertices.length), 1.0, DOUBLE_ERR);
}

@("S2Polyline.Project") unittest {
  S2LatLng[] latlngs = [
      S2LatLng.fromDegrees(0.0, 0.0), S2LatLng.fromDegrees(0.0, 1.0),
      S2LatLng.fromDegrees(0.0, 2.0), S2LatLng.fromDegrees(1.0, 2.0)];
  auto line = new S2Polyline(latlngs);

  int next_vertex;
  Assert.equal(approxEquals(
          line.project(S2LatLng.fromDegrees(0.5, -0.5).toS2Point(), next_vertex),
          S2LatLng.fromDegrees(0, 0).toS2Point()),
      true);
  Assert.equal(next_vertex, 1);
  Assert.equal(approxEquals(
          line.project(S2LatLng.fromDegrees(0.5, 0.5).toS2Point(), next_vertex),
          S2LatLng.fromDegrees(0, 0.5).toS2Point()),
      true);
  Assert.equal(next_vertex, 1);
  Assert.equal(approxEquals(
          line.project(S2LatLng.fromDegrees(0.5, 1).toS2Point(), next_vertex),
          S2LatLng.fromDegrees(0, 1).toS2Point()),
      true);
  Assert.equal(next_vertex, 2);
  Assert.equal(approxEquals(
          line.project(S2LatLng.fromDegrees(-0.5, 2.5).toS2Point(), next_vertex),
          S2LatLng.fromDegrees(0, 2).toS2Point()),
      true);
  Assert.equal(next_vertex, 3);
  Assert.equal(approxEquals(
          line.project(S2LatLng.fromDegrees(2, 2).toS2Point(), next_vertex),
          S2LatLng.fromDegrees(1, 2).toS2Point()),
      true);
  Assert.equal(next_vertex, 4);
}

@("S2Polyline.IsOnRight") unittest {
  S2LatLng[] latlngs = [
      S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
      S2LatLng.fromDegrees(0, 2), S2LatLng.fromDegrees(1, 2)];
  auto line = new S2Polyline(latlngs);

  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(-0.5, 0.5).toS2Point()), true);
  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(0.5, -0.5).toS2Point()), false);
  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(0.5, 0.5).toS2Point()), false);
  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(0.5, 1).toS2Point()), false);
  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(-0.5, 2.5).toS2Point()), true);
  Assert.equal(line.isOnRight(S2LatLng.fromDegrees(1.5, 2.5).toS2Point()), true);

  // Explicitly test the case where the closest point is an interior vertex.
  latlngs = [S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1), S2LatLng.fromDegrees(-1, 0)];
  auto line2 = new S2Polyline(latlngs);

  // The points are chosen such that they are on different sides of the two
  // edges that the interior vertex is on.
  Assert.equal(line2.isOnRight(S2LatLng.fromDegrees(-0.5, 5).toS2Point()), false);
  Assert.equal(line2.isOnRight(S2LatLng.fromDegrees(5.5, 5).toS2Point()), false);
}

@("S2Polyline.IntersectsEmptyPolyline") unittest {
  auto line1 = makePolylineOrDie("1:1, 4:4");
  auto empty_polyline = new S2Polyline();
  Assert.equal(empty_polyline.intersects(line1), false);
}

@("S2Polyline.IntersectsOnePointPolyline") unittest {
  auto line1 = makePolylineOrDie("1:1, 4:4");
  auto line2 = makePolylineOrDie("1:1");
  Assert.equal(line1.intersects(line2), false);
}

@("S2Polyline.Intersects") unittest {
  auto line1 = makePolylineOrDie("1:1, 4:4");
  auto small_crossing = makePolylineOrDie("1:2, 2:1");
  auto small_noncrossing = makePolylineOrDie("1:2, 2:3");
  auto big_crossing = makePolylineOrDie("1:2, 2:3, 4:3");

  Assert.equal(line1.intersects(small_crossing), true);
  Assert.equal(line1.intersects(small_noncrossing), false);
  Assert.equal(line1.intersects(big_crossing), true);
}

@("S2Polyline.IntersectsAtVertex") unittest {
  auto line1 = makePolylineOrDie("1:1, 4:4, 4:6");
  auto line2 = makePolylineOrDie("1:1, 1:2");
  auto line3 = makePolylineOrDie("5:1, 4:4, 2:2");
  Assert.equal(line1.intersects(line2), true);
  Assert.equal(line1.intersects(line3), true);
}

@("S2Polyline.IntersectsVertexOnEdge") unittest {
  auto horizontal_left_to_right = makePolylineOrDie("0:1, 0:3");
  auto vertical_bottom_to_top = makePolylineOrDie("-1:2, 0:2, 1:2");
  auto horizontal_right_to_left = makePolylineOrDie("0:3, 0:1");
  auto vertical_top_to_bottom = makePolylineOrDie("1:2, 0:2, -1:2");
  Assert.equal(horizontal_left_to_right.intersects(vertical_bottom_to_top), true);
  Assert.equal(horizontal_left_to_right.intersects(vertical_top_to_bottom), true);
  Assert.equal(horizontal_right_to_left.intersects(vertical_bottom_to_top), true);
  Assert.equal(horizontal_right_to_left.intersects(vertical_top_to_bottom), true);
}

@("S2Polyline.SpaceUsedEmptyPolyline") unittest {
  auto line = makePolylineOrDie("");
  Assert.greaterThan(line.spaceUsed(), 0);
}

@("S2Polyline.SpaceUsedNonEmptyPolyline") unittest {
  auto line = makePolylineOrDie("1:1, 4:4, 4:6");
  Assert.greaterThan(line.spaceUsed(), 3 * S2Point.sizeof);
}

string joinInts(int[] ints) {
  string result;
  size_t n = ints.length;
  for (size_t i = 0; i + 1 < n; ++i) {
    result ~= to!string(ints[i]) ~ ",";
  }
  if (n > 0) {
    result ~= to!string(ints[n - 1]);
  }
  return result;
}

void checkSubsample(string polyline_str, double tolerance_degrees, string expected_str) {
  logger.logTrace("\"" ~ polyline_str ~ "\", tolerance ", tolerance_degrees);
  auto polyline = makePolylineOrDie(polyline_str);
  int[] indices;
  polyline.subsampleVertices(S1Angle.fromDegrees(tolerance_degrees), indices);
  Assert.equal(expected_str, joinInts(indices));
}

@("S2Polyline.SubsampleVerticesTrivialInputs") unittest {
  // No vertices.
  checkSubsample("", 1.0, "");
  // One vertex.
  checkSubsample("0:1", 1.0, "0");
  // Two vertices.
  checkSubsample("10:10, 11:11", 5.0, "0,1");
  // Three points on a straight line.
  // In theory, zero tolerance should work, but in practice there are floating
  // point errors.
  checkSubsample("-1:0, 0:0, 1:0", 1e-15, "0,2");
  // Zero tolerance on a non-straight line.
  checkSubsample("-1:0, 0:0, 1:1", 0.0, "0,1,2");
  // Negative tolerance should return all vertices.
  checkSubsample("-1:0, 0:0, 1:1", -1.0, "0,1,2");
  // Non-zero tolerance with a straight line.
  checkSubsample("0:1, 0:2, 0:3, 0:4, 0:5", 1.0, "0,4");

  // And finally, verify that we still do something reasonable if the client
  // passes in an invalid polyline with two or more adjacent vertices.
  checkSubsample("0:1, 0:1, 0:1, 0:2", 0.0, "0,3");
}

@("S2Polyline.SubsampleVerticesSimpleExample") unittest {
  string poly_str = "0:0, 0:1, -1:2, 0:3, 0:4, 1:4, 2:4.5, 3:4, 3.5:4, 4:4";
  checkSubsample(poly_str, 3.0, "0,9");
  checkSubsample(poly_str, 2.0, "0,6,9");
  checkSubsample(poly_str, 0.9, "0,2,6,9");
  checkSubsample(poly_str, 0.4, "0,1,2,3,4,6,9");
  checkSubsample(poly_str, 0, "0,1,2,3,4,5,6,7,8,9");
}

@("S2Polyline.SubsampleVerticesGuarantees") unittest {
  // Check that duplicate vertices are never generated.
  checkSubsample("10:10, 12:12, 10:10", 5.0, "0");
  checkSubsample("0:0, 1:1, 0:0, 0:120, 0:130", 5.0, "0,3,4");

  // Check that points are not collapsed if they would create a line segment
  // longer than 90 degrees, and also that the code handles original polyline
  // segments longer than 90 degrees.
  checkSubsample("90:0, 50:180, 20:180, -20:180, -50:180, -90:0, 30:0, 90:0", 5.0, "0,2,4,5,6,7");

  // Check that the output polyline is parametrically equivalent and not just
  // geometrically equivalent, i.e. that backtracking is preserved.  The
  // algorithm achieves this by requiring that the points must be encountered
  // in increasing order of distance along each output segment, except for
  // points that are within "tolerance" of the first vertex of each segment.
  checkSubsample("10:10, 10:20, 10:30, 10:15, 10:40", 5.0, "0,2,3,4");
  checkSubsample("10:10, 10:20, 10:30, 10:10, 10:30, 10:40", 5.0, "0,2,3,5");
  checkSubsample("10:10, 12:12, 9:9, 10:20, 10:30", 5.0, "0,4");
}

static bool testEquals(string a_str, string b_str, S1Angle max_error) {
  auto a = makePolylineOrDie(a_str);
  auto b = makePolylineOrDie(b_str);
  return a.approxEquals(b, max_error);
}

@("S2Polyline.ApproxEquals") unittest {
  S1Angle degree = S1Angle.fromDegrees(1.0);

  // Close lines, differences within max_error.
  Assert.equal(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.5 * degree), true);

  // Close lines, differences outside max_error.
  Assert.equal(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.01 * degree), false);

  // Same line, but different number of vertices.
  Assert.equal(testEquals("0:0, 0:10, 0:20", "0:0, 0:20", 0.1 * degree), false);

  // Same vertices, in different order.
  Assert.equal(testEquals("0:0, 5:5, 0:10", "5:5, 0:10, 0:0", 0.1 * degree), false);
}

/+ TODO: Add when Encode/Decode are added.
@("S2Polyline.EncodeDecode") unittest {
  auto polyline = makePolylineOrDie("0:0, 0:10, 10:20, 20:30");
  Encoder encoder;
  polyline->Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2Polyline decoded_polyline;
  EXPECT_TRUE(decoded_polyline.Decode(&decoder));
  EXPECT_TRUE(decoded_polyline.ApproxEquals(*polyline, S1Angle::Zero()));
}
+/

@("S2PolylineShape.Basic") unittest {
  auto polyline = makePolylineOrDie("0:0, 1:0, 1:1, 2:1");
  auto shape = new S2Polyline.Shape(polyline);
  Assert.equal(polyline, shape.polyline());
  Assert.equal(shape.numEdges(), 3);
  Assert.equal(shape.numChains(), 1);
  Assert.equal(shape.chain(0).start, 0);
  Assert.equal(shape.chain(0).length, 3);
  auto edge2 = shape.edge(2);
  Assert.equal(S2LatLng.fromDegrees(1, 1).toS2Point(), edge2.v0);
  Assert.equal(S2LatLng.fromDegrees(2, 1).toS2Point(), edge2.v1);
  Assert.equal(shape.dimension(), 1);
  Assert.equal(shape.hasInterior(), false);
  Assert.equal(shape.getReferencePoint().contained, false);
}

@("S2PolylineShape.EmptyPolyline") unittest {
  auto polyline = new S2Polyline();
  auto shape = new S2Polyline.Shape(polyline);
  Assert.equal(shape.numEdges(), 0);
  Assert.equal(shape.numChains(), 0);
}

void checkNearlyCovers(
    string a_str, string b_str, double max_error_degrees,
    bool expect_b_covers_a, bool expect_a_covers_b) {
  logger.logTrace("a=\"", a_str, "\", b=\"", b_str, "\", max error=", max_error_degrees);
  auto a = makePolylineOrDie(a_str);
  auto b = makePolylineOrDie(b_str);
  auto max_error = S1Angle.fromDegrees(max_error_degrees);
  Assert.equal(expect_b_covers_a, b.nearlyCovers(a, max_error));
  Assert.equal(expect_a_covers_b, a.nearlyCovers(b, max_error));
}

@("S2PolylineCoveringTest.PolylineOverlapsSelf") unittest {
  string pline = "1:1, 2:2, -1:10";
  checkNearlyCovers(pline, pline, 1e-10, true, true);
}

@("S2PolylineCoveringTest.PolylineDoesNotOverlapReverse") unittest {
  checkNearlyCovers("1:1, 2:2, -1:10", "-1:10, 2:2, 1:1", 1e-10, false, false);
}

@("S2PolylineCoveringTest.PolylineOverlapsEquivalent") unittest {
  // These two polylines trace the exact same polyline, but the second one uses
  // three points instead of two.
  checkNearlyCovers("1:1, 2:1", "1:1, 1.5:1, 2:1", 1e-10, true, true);
}

@("S2PolylineCoveringTest.ShortCoveredByLong") unittest {
  // The second polyline is always within 0.001 degrees of the first polyline,
  // but the first polyline is too long to be covered by the second.
  checkNearlyCovers("-5:1, 10:1, 10:5, 5:10", "9:1, 9.9995:1, 10.0005:5", 1e-3, false, true);
}

@("S2PolylineCoveringTest.PartialOverlapOnly") unittest {
  // These two polylines partially overlap each other, but neither fully
  // overlaps the other.
  checkNearlyCovers("-5:1, 10:1", "0:1, 20:1", 1.0, false, false);
}

@("S2PolylineCoveringTest.ShortBacktracking") unittest {
  // Two lines that backtrack a bit (less than 1.5 degrees) on different edges.
  // A simple greedy matching algorithm would fail on this example.
  string t1 = "0:0, 0:2, 0:1, 0:4, 0:5";
  string t2 = "0:0, 0:2, 0:4, 0:3, 0:5";
  checkNearlyCovers(t1, t2, 1.5, true, true);
  checkNearlyCovers(t1, t2, 0.5, false, false);
}

@("S2PolylineCoveringTest.LongBacktracking") unittest {
  // Two arcs with opposite direction do not overlap if the shorter arc is
  // longer than max_error, but do if the shorter arc is shorter than max-error.
  checkNearlyCovers("5:1, -5:1", "1:1, 3:1", 1.0, false, false);
  checkNearlyCovers("5:1, -5:1", "1:1, 3:1", 2.5, false, true);
}

@("S2PolylineCoveringTest.IsResilientToDuplicatePoints") unittest {
  // S2Polyines are not generally supposed to contain adjacent, identical
  // points, but it happens in practice.  When --s2debug=true, debug-mode
  // binaries abort on such polylines, so we also set --s2debug=false.
  checkNearlyCovers("0:1, 0:2, 0:2, 0:3", "0:1, 0:1, 0:1, 0:3", 1e-10, true, true);
}

@("S2PolylineCoveringTest.CanChooseBetweenTwoPotentialStartingPoints") unittest {
  // Can handle two possible starting points, only one of which leads to finding
  // a correct path.  In the first polyline, the edge from 0:1.1 to 0:0 and the
  // edge from 0:0.9 to 0:2 might be lucrative starting states for covering the
  // second polyline, because both edges are with the max_error of 1.5 degrees
  // from 0:10.  However, only the latter is actually effective.
  checkNearlyCovers("0:11, 0:0, 0:9, 0:20", "0:10, 0:15", 1.5, false, true);
}

@("S2PolylineCoveringTest.StraightAndWigglyPolylinesCoverEachOther") unittest {
  checkNearlyCovers(
      "40:1, 20:1",
      "39.9:0.9, 40:1.1, 30:1.15, 29:0.95, 28:1.1, 27:1.15, "
          ~ "26:1.05, 25:0.85, 24:1.1, 23:0.9, 20:0.99",
      0.2, true, true);
}

@("S2PolylineCoveringTest.MatchStartsAtLastVertex") unittest {
  // The first polyline covers the second, but the matching segment starts at
  // the last vertex of the first polyline.
  checkNearlyCovers("0:0, 0:2", "0:2, 0:3", 1.5, false, true);
}

@("S2PolylineCoveringTest.MatchStartsAtDuplicatedLastVertex") unittest {
  checkNearlyCovers(
      "0:0, 0:2, 0:2, 0:2", "0:2, 0:3", 1.5, false, true);
}

@("S2PolylineCoveringTest.EmptyPolylines") unittest {
  // We expect:
  //    anything.covers(empty) = true
  //    empty.covers(nonempty) = false
  checkNearlyCovers("0:1, 0:2", "", 0.0, false, true);
  checkNearlyCovers("", "", 0.0, true, true);
}
