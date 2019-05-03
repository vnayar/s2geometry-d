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

import s2.mutable_s2shape_index;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2lax_polygon_shape;
import s2.s2lax_polyline_shape;
import s2.s2loop;
import s2.s2point;
import s2.s2point_vector_shape;
import s2.s2polygon;
import s2.s2polyline;
import s2.s2shape;
import s2.s2shape_index;
import s2.strings.serialize;

import std.conv;
import std.exception;
import std.format : format;
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
bool makePoint(string str, ref S2Point point) {
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
bool parseLatLngs(string str, ref S2LatLng[] latlngs) {
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
bool parsePoints(string str, ref S2Point[] vertices) {
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

bool makeLatLng(string str, ref S2LatLng latlng) {
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
bool makeLatLngRect(string str, ref S2LatLngRect rect) {
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
S2Loop makeLoopOrDie(string str) {
  S2Loop loop;
  enforce(makeLoop(str, loop), ": str == \"" ~ str ~ "\"");
  return loop;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeLoop(string str, ref S2Loop loop) {
  if (str == "empty") {
    loop = new S2Loop(S2Loop.empty());
    return true;
  }
  if (str == "full") {
    loop = new S2Loop(S2Loop.full());
    return true;
  }
  S2Point[] vertices;
  if (!parsePoints(str, vertices)) return false;
  loop = new S2Loop(vertices);
  return true;
}

deprecated("Use MakeLoopOrDie.")
S2Loop makeLoop(string str) {
  return makeLoopOrDie(str);
}

// Similar to MakeLoop(), but returns an S2Polyline rather than an S2Loop.
S2Polyline makePolylineOrDie(string str) {
  S2Polyline polyline;
  enforce(makePolyline(str, polyline), ": str == \"" ~ str ~ "\"");
  return polyline;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makePolyline(string str, ref S2Polyline polyline) {
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
S2LaxPolylineShape makeLaxPolylineOrDie(string str) {
  auto lax_polyline = new S2LaxPolylineShape();
  enforce(makeLaxPolyline(str, lax_polyline), ": str == \"" ~ str ~ "\"");
  return lax_polyline;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeLaxPolyline(string str, ref S2LaxPolylineShape lax_polyline) {
  S2Point[] vertices;
  if (!parsePoints(str, vertices)) return false;
  lax_polyline = new S2LaxPolylineShape(vertices);
  return true;
}

deprecated("Use MakeLaxPolylineOrDie.")
S2LaxPolylineShape makeLaxPolyline(string str) {
  return makeLaxPolylineOrDie(str);
}

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
S2Polygon makePolygonOrDie(string str) {
  S2Polygon polygon;
  enforce(internalMakePolygon(str, true, polygon), ": str == \"" ~ str ~ "\"");
  return polygon;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
//ABSL_MUST_USE_RESULT bool MakePolygon(absl::string_view str, std::unique_ptr<S2Polygon>* polygon);

bool makePolygon(string str, ref S2Polygon polygon) {
  return internalMakePolygon(str, true, polygon);
}

// TODO: Resume here.
private bool internalMakePolygon(string str, bool normalize_loops, ref S2Polygon polygon) {
  polygon = new S2Polygon();
  if (str == "empty") str = "";
  string[] loop_strs = str.split(';');
  S2Loop[] loops;
  foreach (loop_str; loop_strs) {
    S2Loop loop;
    if (!makeLoop(loop_str, loop)) return false;
    // Don't normalize loops that were explicitly specified as "full".
    if (normalize_loops && !loop.isFull()) loop.normalize();
    loops ~= loop;
  }
  polygon = new S2Polygon(loops);
  return true;
}


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
S2LaxPolygonShape makeLaxPolygonOrDie(string str) {
  S2LaxPolygonShape lax_polygon;
  enforce(makeLaxPolygon(str, lax_polygon), ": str == \"" ~ str ~ "\"");
  return lax_polygon;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeLaxPolygon(string str, ref S2LaxPolygonShape lax_polygon) {
  string[] loop_strs = str.split(";");
  S2Point[][] loops;
  foreach (loop_str; loop_strs) {
    if (loop_str == "full") {
      loops ~= new S2Point[0];
    } else if (loop_str != "empty") {
      S2Point[] points;
      if (!parsePoints(loop_str, points)) return false;
      loops ~= points;
    }
  }
  lax_polygon = new S2LaxPolygonShape(loops);
  return true;
}

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
MutableS2ShapeIndex makeIndexOrDie(string str) {
  auto index = new MutableS2ShapeIndex();
  enforce(makeIndex(str, index), ": str == \"" ~ str ~ "\"");
  return index;
}

// As above, but does not CHECK-fail on invalid input. Returns true if
// conversion is successful.
bool makeIndex(string str, ref MutableS2ShapeIndex index) {
  string[] strs = str.split('#');
  enforce(strs.length == 3, "Must contain two # characters: " ~ str);

  S2Point[] points;
  foreach (point_str; strs[0].strip().split('|')) {
    S2Point point;
    if (!makePoint(point_str, point)) return false;
    points ~= point;
  }
  if (!points.empty()) {
    index.add(new S2PointVectorShape(points));
  }
  foreach (line_str; strs[1].strip().split('|')) {
    auto lax_polyline = new S2LaxPolylineShape();
    if (!makeLaxPolyline(line_str, lax_polyline)) return false;
    index.add(lax_polyline);
  }
  foreach (polygon_str; strs[2].strip().split('|')) {
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

private void appendVertex(in S2LatLng ll, ref string val) {
  val ~= format("%.15g:%.15g", ll.lat().degrees(), ll.lng().degrees());
}

private void appendVertex(in S2Point p, ref string val) {
  auto ll = S2LatLng(p);
  return appendVertex(ll, val);
}

private void appendVertices(in S2Point[] v, ref string val) {
  for (int i = 0; i < v.length; ++i) {
    if (i > 0) val ~= ", ";
    appendVertex(v[i], val);
  }
}

string toString(in S2Point point) {
  string val;
  appendVertex(point, val);
  return val;
}

string toString(in S2LatLngRect rect) {
  string val;
  appendVertex(rect.lo(), val);
  val ~= ", ";
  appendVertex(rect.hi(), val);
  return val;
}

string toString(in S2LatLng latlng) {
  string val;
  appendVertex(latlng, val);
  return val;
}

string toString(in S2Loop loop) {
  if (loop.isEmpty()) {
    return "empty";
  } else if (loop.isFull()) {
    return "full";
  }
  string val;
  if (loop.numVertices() > 0) {
    appendVertices(loop.vertices(), val);
  }
  return val;
}

string toString(in S2Polyline polyline) {
  string val;
  if (polyline.numVertices() > 0) {
    appendVertices(polyline.vertices(), val);
  }
  return val;
}

string toString(in S2Polygon polygon) {
  if (polygon.isEmpty()) {
    return "empty";
  } else if (polygon.isFull()) {
    return "full";
  }
  string val;
  for (int i = 0; i < polygon.numLoops(); ++i) {
    if (i > 0) val ~= ";\n";
    const(S2Loop) loop = polygon.loop(i);
    appendVertices(loop.vertices(), val);
  }
  return val;
}

string toString(in S2Point[] points) {
  string val;
  appendVertices(points, val);
  return val;
}

string toString(in S2LatLng[] latlngs) {
  string val;
  for (int i = 0; i < latlngs.length; ++i) {
    if (i > 0) val ~= ", ";
    appendVertex(latlngs[i], val);
  }
  return val;
}

string toString(in S2LaxPolylineShape polyline) {
  string val;
  if (polyline.numVertices() > 0) {
    appendVertices(polyline.vertices(), val);
  }
  return val;
}

string toString(in S2LaxPolygonShape polygon) {
  string val;
  for (int i = 0; i < polygon.numLoops(); ++i) {
    if (i > 0) val ~= ";\n";
    int n = polygon.numLoopVertices(i);
    if (n > 0) appendVertices(polygon.loopVertices(i), val);
  }
  return val;
}

// Convert the contents of an S2ShapeIndex to the format above.  The index may
// contain S2Shapes of any type.  Shapes are reordered if necessary so that
// all point geometry (shapes of dimension 0) are first, followed by all
// polyline geometry, followed by all polygon geometry.
// string ToString(const S2ShapeIndex& index);
string toString(in S2ShapeIndex index) {
  string val;
  for (int dim = 0; dim < 3; ++dim) {
    if (dim > 0) val ~= "#";
    int count = 0;
    for (int s = 0; s < index.numShapeIds(); ++s) {
      const(S2Shape) shape = index.shape(s);
      if (shape is null || shape.dimension() != dim) continue;
      val ~= (count > 0) ? " | " : (dim > 0) ? " " : "";
      for (int i = 0; i < shape.numChains(); ++i, ++count) {
        if (i > 0) val ~= (dim == 2) ? "; " : " | ";
        S2Shape.Chain chain = shape.chain(i);
        appendVertex(shape.edge(chain.start).v0, val);
        int limit = chain.start + chain.length;
        if (dim != 1) --limit;
        for (int e = chain.start; e < limit; ++e) {
          val ~= ", ";
          appendVertex(shape.edge(e).v1, val);
        }
      }
    }
    // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
    if (dim == 1 || (dim == 0 && count > 0)) val ~= " ";
  }
  return val;
}
