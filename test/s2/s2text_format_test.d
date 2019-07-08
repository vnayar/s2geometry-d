// Copyright 2018 Google Inc. All Rights Reserved.
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

// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2text_format_test;

import s2.s2text_format;

import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s2latlng;
import s2.s2lax_polygon_shape;
import s2.s2lax_polyline_shape;
import s2.s2loop;
//import s2.s2polygon;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2polyline;
import s2.s2testing;

import std.algorithm : find;
import std.array : split;
import std.math : pow, lround;
import std.range : empty, back;

import fluent.asserts;

private enum int kIters = 10000;

// Verify that s2textformat::ToString() formats the given lat/lng with at most
// "max_digits" after the decimal point and has no trailing zeros.
void expectMaxDigits(in S2LatLng ll, int max_digits) {
  string result = toString(ll.toS2Point());
  string[] values = result.split(':');
  Assert.equal(values.length, 2, result);
  foreach (value; values) {
    size_t num_digits = 0;
    if (!value.find('.').empty) {
      num_digits = value.find('.').length - 1;
      Assert.notEqual(value.back(), '0');
    }
    Assert.notGreaterThan(num_digits, max_digits, value);
  }
}

void expectString(string expected, in S2LatLng ll) {
  Assert.equal(expected, toString(ll.toS2Point()));
}

@("ToString.SpecialCases") unittest {
  expectString("0:0", S2LatLng.fromDegrees(0, 0));
  expectString("1e-20:1e-30", S2LatLng.fromDegrees(1e-20, 1e-30));
}

@("ToString.MinimalDigitsE5") unittest {
  S2Testing.rnd.reset(1);
  for (int iter = 0; iter < kIters; ++iter) {
    auto ll = S2LatLng(S2Testing.randomPoint());
    auto ll_e5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5());
    expectMaxDigits(ll_e5, 5);
  }
}

@("ToString.MinimalDigitsE6") unittest {
  for (int iter = 0; iter < kIters; ++iter) {
    auto ll = S2LatLng(S2Testing.randomPoint());
    auto ll_e6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6());
    expectMaxDigits(ll_e6, 6);
  }
}

@("ToString.MinimalDigitsE7") unittest {
  expectMaxDigits(S2LatLng.fromDegrees(0, 0), 7);
  for (int iter = 0; iter < kIters; ++iter) {
    auto ll = S2LatLng(S2Testing.randomPoint());
    auto ll_e7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7());
    expectMaxDigits(ll_e7, 7);
  }
}

@("ToString.MinimalDigitsDoubleConstants") unittest {
  // Verify that points specified as floating-point literals in degrees using
  // up to 10 digits after the decimal point are formatted with the minimal
  // number of digits.
  for (int iter = 0; iter < kIters; ++iter) {
    int max_digits = S2Testing.rnd.uniform(11);
    long scale = lround(pow(10.0, max_digits));
    long lat = lround(S2Testing.rnd.uniformDouble(-90 * scale, 90 * scale));
    long lng = lround(S2Testing.rnd.uniformDouble(-180 * scale, 180 * scale));
    S2LatLng ll = S2LatLng.fromDegrees(lat / cast(double) scale, lng / cast(double) scale);
    expectMaxDigits(ll, max_digits);
  }
}

@("ToString.UninitializedLoop") unittest {
  auto loop = new S2Loop();
  Assert.equal(toString(loop), "");
}

@("ToString.EmptyLoop") unittest {
  auto empty = new S2Loop(S2Loop.empty());
  Assert.equal("empty", toString(empty));
}

@("ToString.FullLoop") unittest {
  auto full = new S2Loop(S2Loop.full());
  Assert.equal("full", toString(full));
}

@("ToString.EmptyPolyline") unittest {
  auto polyline = new S2Polyline();
  Assert.equal(toString(polyline), "");
}

@("ToString.EmptyPointVector") unittest {
  S2Point[] points;
  Assert.equal(toString(points), "");
}

// @("ToString.EmptyPolygon") unittest {
//   S2Polygon empty;
//   Assert.equal(toString(empty), "empty");
// }

// TEST(ToString, FullPolygon) {
//   S2Polygon full(absl::make_unique<S2Loop>(S2Loop::kFull()));
//   EXPECT_EQ("full", s2textformat::ToString(full));
// }

@("MakeLaxPolygon.Empty") unittest {
  // Verify that "" and "empty" both create empty polygons.
  auto shape = makeLaxPolygonOrDie("");
  Assert.equal(shape.numLoops(), 0);
  shape = makeLaxPolygonOrDie("empty");
  Assert.equal(0, shape.numLoops());
}

@("MakeLaxPolygon.Full") unittest {
  auto shape = makeLaxPolygonOrDie("full");
  Assert.equal(shape.numLoops(), 1);
  Assert.equal(shape.numLoopVertices(0), 0);
}

@("MakeLaxPolygon.FullWithHole") unittest {
  auto shape = makeLaxPolygonOrDie("full; 0:0");
  Assert.equal(shape.numLoops(), 2);
  Assert.equal(shape.numLoopVertices(0), 0);
  Assert.equal(shape.numLoopVertices(1), 1);
  Assert.equal(shape.numEdges(), 1);
}

void checkS2ShapeIndex(in string str) {
  Assert.equal(toString(makeIndexOrDie(str)), str);
}

@("ToString.S2ShapeIndex") unittest {
  checkS2ShapeIndex("# #");
  checkS2ShapeIndex("0:0 # #");
  checkS2ShapeIndex("0:0 | 1:0 # #");
  checkS2ShapeIndex("0:0 | 1:0 # #");
  checkS2ShapeIndex("# 0:0, 0:0 #");
  checkS2ShapeIndex("# 0:0, 0:0 | 1:0, 2:0 #");
  checkS2ShapeIndex("# # 0:0");
  checkS2ShapeIndex("# # 0:0, 0:1");
  checkS2ShapeIndex("# # 0:0, 0:1, 1:0");
  checkS2ShapeIndex("# # 0:0, 0:1, 1:0; 2:2");
}

@("MakePoint.ValidInput") unittest {
  S2Point point;
  Assert.equal(makePoint("-20:150", point), true);
  Assert.equal(S2LatLng.fromDegrees(-20, 150).toS2Point(), point);
}

@("MakePoint.InvalidInput") unittest {
  S2Point point;
  Assert.equal(makePoint("blah", point), false);
}

@("SafeParseLatLngs.ValidInput") unittest {
  S2LatLng[] latlngs;
  Assert.equal(parseLatLngs("-20:150, -20:151, -19:150", latlngs), true);
  Assert.equal(latlngs.length, 3);
  Assert.equal(S2LatLng.fromDegrees(-20, 150), latlngs[0]);
  Assert.equal(S2LatLng.fromDegrees(-20, 151), latlngs[1]);
  Assert.equal(S2LatLng.fromDegrees(-19, 150), latlngs[2]);
}

@("SafeParseLatLngs.InvalidInput") unittest {
  S2LatLng[] latlngs;
  Assert.equal(parseLatLngs("blah", latlngs), false);
}

@("SafeParsePoints.ValidInput") unittest {
  S2Point[] vertices;
  Assert.equal(parsePoints("-20:150, -20:151, -19:150", vertices), true);
  Assert.equal(vertices.length, 3);
  Assert.equal(S2LatLng.fromDegrees(-20, 150).toS2Point(), vertices[0]);
  Assert.equal(S2LatLng.fromDegrees(-20, 151).toS2Point(), vertices[1]);
  Assert.equal(S2LatLng.fromDegrees(-19, 150).toS2Point(), vertices[2]);
}

@("SafeParsePoints.InvalidInput") unittest {
  S2Point[] vertices;
  Assert.equal(parsePoints("blah", vertices), false);
}

@("SafeMakeLatLngRect.ValidInput") unittest {
  S2LatLngRect rect;
  Assert.equal(makeLatLngRect("-10:-10, 10:10", rect), true);
  Assert.equal(new S2LatLngRect(S2LatLng.fromDegrees(-10, -10), S2LatLng.fromDegrees(10, 10)), rect);
}

@("SafeMakeLatLngRect.InvalidInput") unittest {
  S2LatLngRect rect;
  Assert.equal(makeLatLngRect("blah", rect), false);
}

@("SafeMakeLatLng.ValidInput") unittest {
  S2LatLng latlng;
  Assert.equal(makeLatLng("-12.3:45.6", latlng), true);
  Assert.equal(S2LatLng.fromDegrees(-12.3, 45.6), latlng);
}

@("SafeMakeLatLng.InvalidInput") unittest {
  S2LatLng latlng;
  Assert.equal(makeLatLng("blah", latlng), false);
}

// @("SafeMakeLoop.ValidInput") unittest {
//   std::unique_ptr<S2Loop> loop;
//   EXPECT_TRUE(s2textformat::MakeLoop("-20:150, -20:151, -19:150", &loop));
//   EXPECT_TRUE(loop->BoundaryApproxEquals(
//       S2Loop({S2LatLng::FromDegrees(-20, 150).ToPoint(),
//               S2LatLng::FromDegrees(-20, 151).ToPoint(),
//               S2LatLng::FromDegrees(-19, 150).ToPoint()})));
// }

// TEST(SafeMakeLoop, InvalidInput) {
//   std::unique_ptr<S2Loop> loop;
//   EXPECT_FALSE(s2textformat::MakeLoop("blah", &loop));
// }

@("SafeMakePolyline.ValidInput") unittest {
  S2Polyline polyline;
  Assert.equal(
      makePolyline("-20:150, -20:151, -19:150", polyline),
      true);
  auto expected = new S2Polyline([
          S2LatLng.fromDegrees(-20, 150).toS2Point(),
          S2LatLng.fromDegrees(-20, 151).toS2Point(),
          S2LatLng.fromDegrees(-19, 150).toS2Point()]);
  Assert.equal(polyline, expected);
}

@("SafeMakePolyline.InvalidInput") unittest {
  S2Polyline polyline;
  Assert.equal(makePolyline("blah", polyline), false);
}

@("SafeMakeLaxPolyline.ValidInput") unittest {
  S2LaxPolylineShape lax_polyline;
  Assert.equal(makeLaxPolyline("-20:150, -20:151, -19:150", lax_polyline), true);
  // No easy equality check for LaxPolylines; check vertices instead.
  Assert.equal(lax_polyline.numVertices(), 3);
  Assert.equal(
      S2LatLng(lax_polyline.vertex(0)).approxEquals(S2LatLng.fromDegrees(-20, 150)), true);
  Assert.equal(
      S2LatLng(lax_polyline.vertex(1)).approxEquals(S2LatLng.fromDegrees(-20, 151)), true);
  Assert.equal(
      S2LatLng(lax_polyline.vertex(2)).approxEquals(S2LatLng.fromDegrees(-19, 150)), true);
}

@("SafeMakeLaxPolyline.InvalidInput") unittest {
  S2LaxPolylineShape lax_polyline;
  Assert.equal(makeLaxPolyline("blah", lax_polyline), false);
}

// @("SafeMakePolygon.ValidInput") unittest {
//   S2Polygon polygon;
//   Assert.equal(makePolygon("-20:150, -20:151, -19:150", polygon), true);
//   S2Point[] vertices = [
//       S2LatLng.fromDegrees(-20, 150).toS2Point(),
//       S2LatLng.fromDegrees(-20, 151).toS2Point(),
//       S2LatLng.fromDegrees(-19, 150).toS2Point()];
//   S2Polygon expected(absl::make_unique<S2Loop>(vertices));
//   EXPECT_TRUE(polygon->Equals(&expected));
// }

// TEST(SafeMakePolygon, InvalidInput) {
//   std::unique_ptr<S2Polygon> polygon;
//   EXPECT_FALSE(s2textformat::MakePolygon("blah", &polygon));
// }

// TEST(SafeMakePolygon, Empty) {
//   // Verify that "" and "empty" both create empty polygons.
//   std::unique_ptr<S2Polygon> polygon;
//   EXPECT_TRUE(s2textformat::MakePolygon("", &polygon));
//   EXPECT_TRUE(polygon->is_empty());
//   EXPECT_TRUE(s2textformat::MakePolygon("empty", &polygon));
//   EXPECT_TRUE(polygon->is_empty());
// }

// TEST(SafeMakePolygon, Full) {
//   // Verify that "full" creates the full polygon.
//   std::unique_ptr<S2Polygon> polygon;
//   EXPECT_TRUE(s2textformat::MakePolygon("full", &polygon));
//   EXPECT_TRUE(polygon->is_full());
// }

// TEST(SafeMakeVerbatimPolygon, ValidInput) {
//   std::unique_ptr<S2Polygon> polygon;
//   EXPECT_TRUE(
//       s2textformat::MakeVerbatimPolygon("-20:150, -20:151, -19:150", &polygon));
//   std::vector<S2Point> vertices({S2LatLng::FromDegrees(-20, 150).ToPoint(),
//                                  S2LatLng::FromDegrees(-20, 151).ToPoint(),
//                                  S2LatLng::FromDegrees(-19, 150).ToPoint()});
//   S2Polygon expected(absl::make_unique<S2Loop>(vertices));
//   EXPECT_TRUE(polygon->Equals(&expected));
// }

// TEST(SafeMakeVerbatimPolygon, InvalidInput) {
//   std::unique_ptr<S2Polygon> polygon;
//   EXPECT_FALSE(s2textformat::MakeVerbatimPolygon("blah", &polygon));
// }

@("SafeMakeLaxPolygon.ValidInput") unittest {
  S2LaxPolygonShape lax_polygon;
  Assert.equal(makeLaxPolygon("-20:150, -20:151, -19:150", lax_polygon), true);
  // No easy equality check for LaxPolygons; check vertices instead.
  Assert.equal(lax_polygon.numLoops(), 1);
  Assert.equal(lax_polygon.numVertices(), 3);
  Assert.equal(
      S2LatLng(lax_polygon.loopVertex(0, 0)).approxEquals(S2LatLng.fromDegrees(-20, 150)), true);
  Assert.equal(
      S2LatLng(lax_polygon.loopVertex(0, 1)).approxEquals(S2LatLng.fromDegrees(-20, 151)), true);
  Assert.equal(
      S2LatLng(lax_polygon.loopVertex(0, 2)).approxEquals(S2LatLng.fromDegrees(-19, 150)), true);
}

@("SafeMakeLaxPolygon.InvalidInput") unittest {
  S2LaxPolygonShape lax_polygon;
  Assert.equal(makeLaxPolygon("blah", lax_polygon), false);
}

@("SafeMakeIndex.ValidInput") unittest {
  auto index = new MutableS2ShapeIndex();
  Assert.equal(makeIndex("# 0:0, 0:0 | 1:0, 2:0 #", index), true);
  Assert.equal("# 0:0, 0:0 | 1:0, 2:0 #", toString(index));
}

@("SafeMakeIndex.InvalidInput") unittest {
  auto index = new MutableS2ShapeIndex();
  Assert.equal(makeIndex("# blah #", index), false);
}
