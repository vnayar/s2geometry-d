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

module s2.s2text_format;

// s2text_format contains a collection of functions for converting
// geometry to and from a human-readable format.  It is mainly
// intended for testing and debugging.  Be aware that the
// human-readable format is *not* designed to preserve the full
// precision of the original object, so it should not be used
// for data storage.

import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2point;
import s2.strings.serialize;
import s2.s2polyline;
import s2.s2shape;

import std.conv;
import std.exception;
import std.range;
import std.string;

// Returns an S2Point corresponding to the given a latitude-longitude
// coordinate in degrees.  Example of the input format:
//     "-20:150"
S2Point makePointOrDie(string str) {
  S2Point point;
  enforce(makePoint(str, point), ": str == \"" ~ str ~ "\"");
  return point;
}

// As above, but do not CHECK-fail on invalid input. Returns true if conversion
// is successful.
bool makePoint(string str, out S2Point point) {
  S2Point[] vertices;
  if (!parsePoints(str, vertices) || vertices.length != 1) return false;
  point = vertices[0];
  return true;
}

deprecated("Use MakePointOrDie.")
S2Point makePoint(string str) {
  return makePointOrDie(str);
}

// Parses a string of one or more latitude-longitude coordinates in degrees,
// and return the corresponding vector of S2LatLng points.
// Examples of the input format:
//     ""                            // no points
//     "-20:150"                     // one point
//     "-20:150, -20:151, -19:150"   // three points
S2LatLng[] parseLatLngsOrDie(string str) {
  S2LatLng[] latlngs;
  enforce(parseLatLngs(str, latlngs), ": str == \"" ~ str ~ "\"");
  return latlngs;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool parseLatLngs(string str, out S2LatLng[] latlngs) {
  string[2][] ps;
  if (!dictionaryParse(str, ps)) return false;
  foreach (p; ps) {
    try {
      double lat = to!double(strip(p[0]));
      double lng = to!double(strip(p[1]));
      latlngs ~= S2LatLng.fromDegrees(lat, lng);
    } catch (ConvException e) {
      return false;
    }
  }
  return true;
}

deprecated("Use ParseLatLngsOrDie.")
S2LatLng[] parseLatLngs(string str) {
  return parseLatLngsOrDie(str);
}

// Parses a string in the same format as ParseLatLngs, and return the
// corresponding vector of S2Point values.
S2Point[] parsePointsOrDie(string str) {
  S2Point[] vertices;
  enforce(parsePoints(str, vertices), ": str == \"" ~ str ~ "\"");
  return vertices;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool parsePoints(string str, out S2Point[] vertices) {
  S2LatLng[] latlngs;
  if (!parseLatLngs(str, latlngs)) return false;
  foreach (latlng; latlngs) {
    vertices ~= latlng.toS2Point();
  }
  return true;
}

deprecated("Use ParsePointsOrDie.")
S2Point[] parsePoints(string str) {
  return parsePointsOrDie(str);
}

bool makeLatLng(string str, out S2LatLng latlng) {
  S2LatLng[] latlngs;
  if (!parseLatLngs(str, latlngs) || latlngs.length != 1) return false;
  latlng = latlngs[0];
  return true;
}

// Given a string in the same format as ParseLatLngs, returns a single S2LatLng
S2LatLng makeLatLngOrDie(string str) {
  S2LatLng latlng;
  enforce(makeLatLng(str, latlng), ": str == \"" ~ str ~ "\"");
  return latlng;
}

// Given a string in the same format as ParseLatLngs, returns the minimal
// bounding S2LatLngRect that contains the coordinates.
S2LatLngRect makeLatLngRectOrDie(string str) {
  S2LatLngRect rect;
  enforce(makeLatLngRect(str, rect), ": str == \"" ~ str ~ "\"");
  return rect;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeLatLngRect(string str, out S2LatLngRect rect) {
  S2LatLng[] latlngs;
  if (!parseLatLngs(str, latlngs) || latlngs.empty()) return false;
  rect = S2LatLngRect.fromPoint(latlngs[0]);
  for (int i = 1; i < latlngs.length; ++i) {
    rect.addPoint(latlngs[i]);
  }
  return true;
}

deprecated("Use MakeLatLngRectOrDie.")
S2LatLngRect makeLatLngRect(string str) {
  return makeLatLngRectOrDie(str);
}

// Given a string of latitude-longitude coordinates in degrees,
// returns a newly allocated loop.  Example of the input format:
//     "-20:150, 10:-120, 0.123:-170.652"
// The strings "empty" or "full" create an empty or full loop respectively.
//std::unique_ptr<S2Loop> MakeLoopOrDie(absl::string_view str);

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakeLoop(absl::string_view str,
//                                   std::unique_ptr<S2Loop>* loop);

//ABSL_DEPRECATED("Use MakeLoopOrDie.")
//std::unique_ptr<S2Loop> MakeLoop(absl::string_view str);

// Similar to MakeLoop(), but returns an S2Polyline rather than an S2Loop.
S2Polyline makePolylineOrDie(string str) {
  S2Polyline polyline;
  enforce(makePolyline(str, polyline), ": str == \"" ~ str ~ "\"");
  return polyline;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makePolyline(string str, out S2Polyline polyline) {
  S2Point[] vertices;
  if (!parsePoints(str, vertices)) return false;
  polyline = new S2Polyline(vertices);
  return true;
}

deprecated("Use MakePolylineOrDie.")
S2Polyline makePolyline(string str) {
  return makePolylineOrDie(str);
}

// Like MakePolyline, but returns an S2LaxPolylineShape instead.
//std::unique_ptr<S2LaxPolylineShape> MakeLaxPolylineOrDie(absl::string_view str);

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakeLaxPolyline(
//    absl::string_view str, std::unique_ptr<S2LaxPolylineShape>* lax_polyline);

//ABSL_DEPRECATED("Use MakeLaxPolylineOrDie.")
//std::unique_ptr<S2LaxPolylineShape> MakeLaxPolyline(absl::string_view str);

// Given a sequence of loops separated by semicolons, returns a newly
// allocated polygon.  Loops are automatically normalized by inverting them
// if necessary so that they enclose at most half of the unit sphere.
// (Historically this was once a requirement of polygon loops.  It also
// hides the problem that if the user thinks of the coordinates as X:Y
// rather than LAT:LNG, it yields a loop with the opposite orientation.)
//
// Examples of the input format:
//     "10:20, 90:0, 20:30"                                  // one loop
//     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"   // two loops
//     ""       // the empty polygon (consisting of no loops)
//     "empty"  // the empty polygon (consisting of no loops)
//     "full"   // the full polygon (consisting of one full loop).
//std::unique_ptr<S2Polygon> MakePolygonOrDie(absl::string_view str);

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakePolygon(absl::string_view str, std::unique_ptr<S2Polygon>* polygon);

//ABSL_DEPRECATED("Use MakePolygonOrDie.")
//std::unique_ptr<S2Polygon> MakePolygon(absl::string_view str);

// Like MakePolygon(), except that it does not normalize loops (i.e., it
// gives you exactly what you asked for).
//std::unique_ptr<S2Polygon> MakeVerbatimPolygonOrDie(absl::string_view str);

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakeVerbatimPolygon(
//    absl::string_view str, std::unique_ptr<S2Polygon>* polygon);

//ABSL_DEPRECATED("Use MakeVerbatimPolygonOrDie.")
//std::unique_ptr<S2Polygon> MakeVerbatimPolygon(absl::string_view str);

// Parses a string in the same format as MakePolygon, except that loops must
// be oriented so that the interior of the loop is always on the left, and
// polygons with degeneracies are supported.  As with MakePolygon, "full" and
// denotes the full polygon and "" or "empty" denote the empty polygon.
//std::unique_ptr<S2LaxPolygonShape> MakeLaxPolygonOrDie(absl::string_view str);

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakeLaxPolygon(
//    absl::string_view str, std::unique_ptr<S2LaxPolygonShape>* lax_polygon);

//ABSL_DEPRECATED("Use MakeLaxPolygonOrDie.")
//std::unique_ptr<S2LaxPolygonShape> MakeLaxPolygon(absl::string_view str);

// Returns a MutableS2ShapeIndex containing the points, polylines, and loops
// (in the form of a single polygon) described by the following format:
//
//   point1|point2|... # line1|line2|... # polygon1|polygon2|...
//
// Examples:
//   1:2 | 2:3 # #                     // Two points
//   # 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
//   # # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
//   5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each
//   # # empty                         // One empty polygon
//   # # empty | full                  // One empty polygon, one full polygon
//
// Loops should be directed so that the region's interior is on the left.
// Loops can be degenerate (they do not need to meet S2Loop requirements).
//
// CAVEAT: Because whitespace is ignored, empty polygons must be specified
//         as the string "empty" rather than as the empty string ("").
/+ TODO: Resume when S2LaxPolygonShape is implemented.
MutableS2ShapeIndex makeIndexOrDie(string str) {
  auto index = new MutableS2ShapeIndex();
  enforce(makeIndex(str, index), ": str == \"" ~ str ~ "\"");
  return index;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeIndex(string str, MutableS2ShapeIndex index) {
  string[] strs = str.split('#');
  enforce(strs.size() == 3, "Must contain two # characters: " ~ str);

  S2Point[] points;
  foreach (auto point_str; strs[0].split('|')) {
    S2Point point;
    if (!makePoint(point_str, point)) return false;
    points ~= point;
  }
  if (!points.empty()) {
    index.add(new S2PointVectorShape(points));
  }
  foreach (auto line_str; strs[1].split('|')) {
    auto lax_polyline = new S2LaxPolylineShape();
    if (!makeLaxPolyline(line_str, lax_polyline)) return false;
    index.add(new S2Shape(lax_polyline));
  }
  foreach (string polygon_str; strs[2].split('|')) {
    auto lax_polygon = new S2LaxPolygonShape();
    if (!makeLaxPolygon(polygon_str, lax_polygon)) return false;
    index.add(lax_polygon);
  }
  return true;
}

deprecated("Use MakeIndexOrDie.")
MutableS2ShapeIndex makeIndex(string str) {
  return makeIndexOrDie(str);
}
+/
// string ToString(const S2Point& point);
// string ToString(const S2LatLngRect& rect);
// string ToString(const S2LatLng& latlng);
// string ToString(const S2Loop& loop);
// string ToString(const S2Polyline& polyline);
// string ToString(const S2Polygon& polygon);
// string ToString(const std::vector<S2Point>& points);
// string ToString(const std::vector<S2LatLng>& points);
// string ToString(const S2LaxPolylineShape& polyline);
// string ToString(const S2LaxPolygonShape& polygon);

// Convert the contents of an S2ShapeIndex to the format above.  The index may
// contain S2Shapes of any type.  Shapes are reordered if necessary so that
// all point geometry (shapes of dimension 0) are first, followed by all
// polyline geometry, followed by all polygon geometry.
// string ToString(const S2ShapeIndex& index);
