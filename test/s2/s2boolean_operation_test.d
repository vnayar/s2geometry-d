// Copyright 2017 Google Inc. All Rights Reserved.
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

module s2.s2boolean_operation_test;

import s2.s2boolean_operation;

import s2.s2polygon;
import s2.builder.graph;
import s2.builder.layer;
import s2.builder.util.snap_functions;
import s2.s2builder;
import s2.s2error;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2text_format;
import s2.s2text_format;

import std.algorithm;
import std.conv;
import std.range;
import fluent.asserts;

alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;
alias SiblingPairs = GraphOptions.SiblingPairs;

alias OpType = S2BooleanOperation.OpType;
alias PolygonModel = S2BooleanOperation.PolygonModel;
alias PolylineModel = S2BooleanOperation.PolylineModel;

S2Error.Code INDEXES_DO_NOT_MATCH = S2Error.Code.USER_DEFINED_START;

class IndexMatchingLayer : Layer {
public:
  this(in S2ShapeIndex index, int dimension) {
    _index = index;
    _dimension = dimension;
  }

  override
  GraphOptions graphOptions() const {
    return new GraphOptions(
        EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.KEEP, SiblingPairs.KEEP);
  }

  override
  void build(in Graph g, ref S2Error error) {
    S2Shape.Edge[] actual, expected;
    for (int e = 0; e < g.numEdges(); ++e) {
      const(Graph.Edge) edge = g.edge(e);
      actual ~= S2Shape.Edge(g.vertex(edge[0]), g.vertex(edge[1]));
    }
    for (int s = 0; s < _index.numShapeIds(); ++s) {
      const(S2Shape) shape = _index.shape(s);
      if (shape is null || shape.dimension() != _dimension) {
        continue;
      }
      for (int e = shape.numEdges(); --e >= 0; ) {
        expected ~= shape.edge(e);
      }
    }
    sort(actual);
    sort(expected);

    // The edges are a multiset, so we can't use std::set_difference.
    S2Shape.Edge[] missing, extra;
    for (auto ai = 0, ei = 0; ai != actual.length || ei != expected.length; ) {
      if (ei == expected.length || (ai != actual.length && actual[ai] < expected[ei])) {
        extra ~= actual[ai++];
      } else if (ai == actual.length || expected[ei] < actual[ai]) {
        missing ~= expected[ei++];
      } else {
        ++ai;
        ++ei;
      }
    }
    if (!missing.empty() || !extra.empty()) {
      // There may be errors in more than one dimension, so we append to the
      // existing error text.
      error.initialize(INDEXES_DO_NOT_MATCH,
          "%sDimension %d: Missing edges: %s Extra edges: %s\n",
          error.text(), _dimension, toString(missing), toString(extra));
    }
  }

private:
  alias EdgeVector = S2Shape.Edge[];
  static string toString(in EdgeVector edges) {
    string msg;
    foreach (const edge; edges) {
      S2Point[] vertices = [ edge.v0, edge.v1 ];
      msg ~= .toString(vertices);
      msg ~= "; ";
    }
    return msg;
  }

  const(S2ShapeIndex) _index;
  int _dimension;
}

void expectResult(
    S2BooleanOperation.OpType op_type,
    S2BooleanOperation.Options options,
    string a_str, string b_str,
    string expected_str) {
  auto a = makeIndexOrDie(a_str);
  auto b = makeIndexOrDie(b_str);
  auto expected = makeIndexOrDie(expected_str);
  Layer[] layers;
  for (int dim = 0; dim < 3; ++dim) {
    layers ~= new IndexMatchingLayer(expected, dim);
  }
  auto op = new S2BooleanOperation(op_type, layers, options);
  S2Error error;
  Assert.equal(op.build(a, b, error), true,
      op_type.to!string ~ " failed:\n" ~ "Expected result: "
      ~ expected_str ~ "\n" ~ error.toString());

  // Now try the same thing with boolean output.
  Assert.equal(expected.numShapeIds() == 0, S2BooleanOperation.isEmpty(op_type, a, b, options));
}

// The intersections in the "expected" data below were computed in lat-lng
// space (i.e., the rectangular projection), while the actual intersections
// are computed using geodesics.  We can compensate for this by rounding the
// intersection points to a fixed precision in degrees (e.g., 2 decimals).
private S2BooleanOperation.Options roundToE(int exp) {
  S2BooleanOperation.Options options;
  options.setSnapFunction(new IntLatLngSnapFunction(exp));
  return options;
}

// TODO(ericv): Clean up or remove these notes.
//
// Options to test:
//   polygon_model:                   OPEN, SEMI_OPEN, CLOSED
//   polyline_model:                  OPEN, SEMI_OPEN, CLOSED
//   polyline_loops_have_boundaries:  true, false
//   conservative:                    true, false
//
// Geometry combinations to test:
//
// Point/point:
//  - disjoint, coincident
// Point/polyline:
//  - Start vertex, end vertex, interior vertex, degenerate polyline
//  - With polyline_loops_have_boundaries: start/end vertex, degenerate polyline
// Point/polygon:
//  - Polygon interior, exterior, vertex
//  - Vertex of degenerate sibling pair shell, hole
//  - Vertex of degenerate single point shell, hole
// Polyline/polyline:
//  - Vertex intersection:
//    - Start, end, interior, degenerate, loop start/end, degenerate loop
//    - Test cases where vertex is not emitted because an incident edge is.
//  - Edge/edge: interior crossing, duplicate, reversed, degenerate
//  - Test that degenerate edges are ignored unless polyline has a single edge.
//    (For example, AA has one edge but AAA has no edges.)
// Polyline/polygon:
//  - Vertex intersection: polyline vertex cases already covered, but test
//    polygon normal vertex, sibling pair shell/hole, single vertex shell/hole
//    - Also test cases where vertex is not emitted because an edge is.
//  - Edge/edge: interior crossing, duplicate, reversed
//  - Edge/interior: polyline edge in polygon interior, exterior
// Polygon/polygon:
//  - Vertex intersection:
//    - normal vertex, sibling pair shell/hole, single vertex shell/hole
//    - Also test cases where vertex is not emitted because an edge is.
//    - Test that polygons take priority when there is a polygon vertex and
//      also isolated polyline vertices.  (There should not be any points.)
//  - Edge/edge: interior crossing, duplicate, reversed
//  - Interior/interior: polygons in interior/exterior of other polygons

@("S2BooleanOperation.DegeneratePolylines") unittest {
  // Verify that degenerate polylines are preserved under all boundary models.
  S2BooleanOperation.Options options;
  auto a = "# 0:0, 0:0 #";
  auto b = "# #";
  options.setPolylineModel(PolylineModel.OPEN);
  expectResult(OpType.UNION, options, a, b, a);
  options.setPolylineModel(PolylineModel.SEMI_OPEN);
  expectResult(OpType.UNION, options, a, b, a);
  options.setPolylineModel(PolylineModel.CLOSED);
  expectResult(OpType.UNION, options, a, b, a);
}

@("S2BooleanOperation.DegeneratePolygons") unittest {
  // Verify that degenerate polygon features (single-vertex and sibling pair
  // shells and holes) are preserved under all boundary models.
  S2BooleanOperation.Options options;
  auto a = "# # 0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 3:3; 6:6; 7:7, 8:8";
  auto b = "# #";
  options.setPolygonModel(PolygonModel.OPEN);
  expectResult(OpType.UNION, options, a, b, a);
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  expectResult(OpType.UNION, options, a, b, a);
  options.setPolygonModel(PolygonModel.CLOSED);
  expectResult(OpType.UNION, options, a, b, a);
}

@("S2BooleanOperation.PointPoint") unittest {
  S2BooleanOperation.Options options;
  auto a = "0:0 | 1:0 # #";
  auto b = "0:0 | 2:0 # #";
  // Note that these results have duplicates, which is correct.  Clients can
  // eliminated the duplicates with the appropriate GraphOptions.
  expectResult(OpType.UNION, options, a, b, "0:0 | 0:0 | 1:0 | 2:0 # #");
  expectResult(OpType.INTERSECTION, options, a, b, "0:0 | 0:0 # #");
  expectResult(OpType.DIFFERENCE, options, a, b, "1:0 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b, "1:0 | 2:0 # #");
}

@("S2BooleanOperation.PointOpenPolyline") unittest {
  // Tests operations between an open polyline and its vertices.
  //
  // The polyline "3:0, 3:0" consists of a single degenerate edge and contains
  // no points (since polyline_model() is OPEN).  Since S2BooleanOperation
  // preserves degeneracies, this means that the union includes *both* the
  // point 3:0 and the degenerate polyline 3:0, since they do not intersect.
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.OPEN);
  auto a = "0:0 | 1:0 | 2:0 | 3:0 # #";
  auto b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 #";
  expectResult(OpType.UNION, options, a, b,
      "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "1:0 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "0:0 | 2:0 | 3:0 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 #");
}

@("S2BooleanOperation.PointSemiOpenPolyline") unittest {
  // Degenerate polylines are defined not contain any points under the
  // SEMI_OPEN model either, so again the point 3:0 and the degenerate
  // polyline "3:0, 3:0" do not intersect.
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.SEMI_OPEN);
  auto a = "0:0 | 1:0 | 2:0 | 3:0 # #";
  auto b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 #";
  expectResult(OpType.UNION, options, a, b,
      "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "0:0 | 1:0 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "2:0 | 3:0 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 #");
}

@("S2BooleanOperation.PointClosedPolyline") unittest {
  // Under the CLOSED model, the degenerate polyline 3:0 does contain its
  // vertex.  Since polylines take precedence over points, the union of the
  // point 3:0 and the polyline 3:0 is the polyline only.  Similarly, since
  // subtracting a point from a polyline has no effect, the symmetric
  // difference includes only the polyline objects.
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.CLOSED);
  auto a = "0:0 | 1:0 | 2:0 | 3:0 # #";
  auto b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 #";
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 1:0, 2:0 | 3:0, 3:0 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "0:0 | 1:0 | 2:0 | 3:0 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 1:0, 2:0 | 3:0, 3:0 #");
}

@("S2BooleanOperation.PointPolygonInterior") unittest {
  S2BooleanOperation.Options options;  // PolygonModel is irrelevant.
  // One interior point and one exterior point.
  auto a = "1:1 | 4:4 # #";
  auto b = "# # 0:0, 0:3, 3:0";
  expectResult(OpType.UNION, options, a, b,
      "4:4 # # 0:0, 0:3, 3:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "1:1 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "4:4 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "4:4 # # 0:0, 0:3, 3:0");
}

@("S2BooleanOperation.PointOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);
  // See notes about the two vertices below.
  auto a = "0:1 | 1:0 # #";
  auto b = "# # 0:0, 0:1, 1:0";
  expectResult(OpType.UNION, options, a, b,
      "0:1 | 1:0 # # 0:0, 0:1, 1:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "0:1 | 1:0 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "0:1 | 1:0 # # 0:0, 0:1, 1:0");
}

@("S2BooleanOperation.PointSemiOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  // The two vertices are chosen such that the polygon contains one vertex but
  // not the other under PolygonModel::SEMI_OPEN.  (The same vertices are used
  // for all three PolygonModel options.)
  auto polygon = makePolygonOrDie("0:0, 0:1, 1:0");
  Assert.equal(polygon.contains(makePointOrDie("0:1")), true);
  Assert.equal(polygon.contains(makePointOrDie("1:0")), false);
  auto a = "0:1 | 1:0 # #";
  auto b = "# # 0:0, 0:1, 1:0";
  expectResult(OpType.UNION, options, a, b,
      "1:0 # # 0:0, 0:1, 1:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "0:1 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "1:0 # #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "1:0 # # 0:0, 0:1, 1:0");
}

@("S2BooleanOperation.PointClosedPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.CLOSED);
  // See notes about the two vertices above.
  auto a = "0:1 | 1:0 # #";
  auto b = "# # 0:0, 0:1, 1:0";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:1, 1:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "0:1 | 1:0 # #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:1, 1:0");
}

@("S2BooleanOperation.PolylineVertexOpenPolylineVertex") unittest {
  // Test starting, ending, and middle vertices of both polylines.  Degenerate
  // polylines are tested in PolylineEdgePolylineEdgeOverlap below.
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.OPEN);
  auto a = "# 0:0, 0:1, 0:2 #";
  auto b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #";
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");

  // The output consists of the portion of each input polyline that intersects
  // the opposite region, so the intersection vertex is present twice.  This
  // allows reassembling the individual polylins that intersect, if desired.
  // (Otherwise duplicates can be removed using DuplicateEdges::MERGE.)
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 0:1, 0:1 | 0:1, 0:1 #");

  // Note that all operations are defined such that subtracting a
  // lower-dimensional subset of an object has no effect.  In this case,
  // subtracting the middle vertex of a polyline has no effect.
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");
}

@("S2BooleanOperation.PolylineVertexSemiOpenPolylineVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.SEMI_OPEN);
  auto a = "# 0:0, 0:1, 0:2 #";
  auto b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #";
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");
}

@("S2BooleanOperation.PolylineVertexClosedPolylineVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.CLOSED);
  auto a = "# 0:0, 0:1, 0:2 #";
  auto b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #";
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 | 0:2, 0:2 | 0:2, 0:2 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 #");
}

// The polygon used in the polyline/polygon vertex tests below.
private string kVertexTestPolygonStr() {
  return "0:0, 0:1, 0:2, 0:3, 0:4, 0:5, 5:5, 5:4, 5:3, 5:2, 5:1, 5:0";
}

@("S2BooleanOperation.TestSemiOpenPolygonVerticesContained") unittest {
  // Verify whether certain vertices of the test polygon are contained under
  // the semi-open boundary model (for use in the tests below).
  auto polygon = makePolygonOrDie(kVertexTestPolygonStr());
  Assert.equal(polygon.contains(makePointOrDie("0:1")), true);
  Assert.equal(polygon.contains(makePointOrDie("0:2")), true);
  Assert.equal(polygon.contains(makePointOrDie("0:3")), true);
  Assert.equal(polygon.contains(makePointOrDie("0:4")), true);
  Assert.equal(polygon.contains(makePointOrDie("5:1")), false);
  Assert.equal(polygon.contains(makePointOrDie("5:2")), false);
  Assert.equal(polygon.contains(makePointOrDie("5:3")), false);
  Assert.equal(polygon.contains(makePointOrDie("5:4")), false);
}

// Don't bother testing every PolylineModel with every PolygonModel for vertex
// intersection, since we have already tested the PolylineModels individually
// above.  It is sufficient to use PolylineModel::CLOSED with the various
// PolygonModel options.
@("S2BooleanOperation.PolylineVertexOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);

  // Define some constants to reduce code duplication.
  // Test all combinations of polylines that start or end on a polygon vertex,
  // where the polygon vertex is open or closed using semi-open boundaries,
  // and where the incident edge is inside or outside the polygon.
  auto a = "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
      ~ "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
  auto b = "# # " ~ kVertexTestPolygonStr();

  const string kDifferenceResult =
      "# 0:1, 0:1 | 0:2, 0:2 | -1:3, 0:3 | 0:4, -1:4"
      ~ "| 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
  expectResult(OpType.UNION, options, a, b,
               kDifferenceResult ~ kVertexTestPolygonStr());
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 0:1 | 0:2, 1:2 | 4:3, 5:3 | 5:4, 4:4 #");
  expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      kDifferenceResult ~ kVertexTestPolygonStr());
}

// Like the test above, except that every polygon vertex is also incident to a
// closed polyline vertex.  This tests that when an open vertex and a closed
// vertex coincide with each other, the result is considered closed.
// TODO: Fix me.
@("S2BooleanOperation.PolylineVertexOpenPolygonClosedPolylineVertex") unittest {
  const string kTestGeometrySuffix =
      "-2:0, 0:1 | -2:1, 0:2 | -2:2, 0:3 | -2:3, 0:4 | "
      ~ "7:0, 5:1 | 7:1, 5:2 | 7:2, 5:3 | 7:3, 5:4 # " ~ kVertexTestPolygonStr();

  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);
  auto a = "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
      ~ "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
  auto b = ("# " ~ kTestGeometrySuffix);

  const string kDifferencePrefix =
      "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2";
  expectResult(OpType.UNION, options, a, b,
      kDifferencePrefix
      ~ " | 0:1, 0:1 | 0:2, 0:2 | 5:3, 5:3 | 5:4, 5:4 | "
      ~ kTestGeometrySuffix);
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4"
      ~ "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4"
      ~ "| 0:1, 0:1 | 0:2, 0:2 | 0:3, 0:3 | 0:4, 0:4"
      ~ "| 5:1, 5:1 | 5:2, 5:2 | 5:3, 5:3 | 5:4, 5:4 #");
  expectResult(OpType.DIFFERENCE, options, a, b, kDifferencePrefix ~ " #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      kDifferencePrefix ~ " | " ~ kTestGeometrySuffix);
}

@("S2BooleanOperation.PolylineVertexSemiOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  // Test all combinations of polylines that start or end on a polygon vertex,
  // where the polygon vertex is open or closed using semi-open boundaries,
  // and where the incident edge is inside or outside the polygon.
  //
  // The vertices at latitude 0 used below are all closed while the vertices
  // at latitude 5 are all open (see TestSemiOpenPolygonVerticesContained).
  auto a = "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
      ~ "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
  auto b = "# # " ~ kVertexTestPolygonStr();
  const string kDifferenceResult =
      "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
  expectResult(OpType.UNION, options, a, b,
               kDifferenceResult ~ kVertexTestPolygonStr());
  expectResult(OpType.INTERSECTION, options, a, b,
               "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4 "
               ~ "| 4:3, 5:3 | 5:4, 4:4 #");
  expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
               kDifferenceResult ~ kVertexTestPolygonStr());
}

@("S2BooleanOperation.PolylineVertexClosedPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.CLOSED);
  // Test all combinations of polylines that start or end on a polygon vertex,
  // where the polygon vertex is open or closed using semi-open boundaries,
  // and where the incident edge is inside or outside the polygon.
  auto a = "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
      ~ "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
  auto b = "# # " ~ kVertexTestPolygonStr();
  const string kDifferenceResult =
      "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 #";
  expectResult(OpType.UNION, options, a, b,
      kDifferenceResult ~ kVertexTestPolygonStr());
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4"
      ~ "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4 #");
  expectResult(OpType.DIFFERENCE, options, a, b, kDifferenceResult);
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      kDifferenceResult ~ kVertexTestPolygonStr());
}

@("S2BooleanOperation.PolylineEdgePolylineEdgeCrossing") unittest {
  // Two polyline edges that cross at a point interior to both edges.
  S2BooleanOperation.Options options = roundToE(1);
  auto a = "# 0:0, 2:2 #";
  auto b = "# 2:0, 0:2 #";
  expectResult(OpType.UNION, options, a, b,
     "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 1:1 | 1:1, 1:1 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:0, 2:2 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
}

@("S2BooleanOperation.PolylineEdgePolylineEdgeOverlap") unittest {
  // The PolylineModel does not affect this calculation.  In particular the
  // intersection of a degenerate polyline edge with itself is non-empty, even
  // though the edge contains no points in the OPEN and SEMI_OPEN models.
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);
  // Test edges in the same and reverse directions, and degenerate edges.
  auto a = "# 0:0, 1:0, 2:0, 2:5 | 3:0, 3:0 | 6:0, 5:0, 4:0 #";
  auto b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0 #";
  // As usual, the expected output includes the relevant portions of *both*
  // input polylines.  Duplicates can be removed using GraphOptions.
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 1:0, 2:0, 2:5 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 "
      ~ "| 6:0, 5:0, 4:0 | 4:0, 5:0 #");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 0:0, 1:0, 2:0 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 "
      ~ "| 5:0, 4:0 | 4:0, 5:0 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 2:0, 2:5 | 6:0, 5:0 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 2:0, 2:5 | 6:0, 5:0 #");
}

@("S2BooleanOperation.PolylineEdgeOpenPolygonEdgeOverlap") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);
  // A polygon and two polyline edges that coincide with the polygon boundary,
  // one in the same direction and one in the reverse direction.
  auto a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
  auto b = "# # 1:1, 1:3, 3:3, 3:1";
  expectResult(OpType.UNION, options, a, b,
      "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
}

@("S2BooleanOperation.PolylineEdgeSemiOpenPolygonEdgeOverlap") unittest {
  auto polygon = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");
  Assert.equal(polygon.contains(makePointOrDie("1:1")), false);
  Assert.equal(polygon.contains(makePointOrDie("1:3")), true);
  Assert.equal(polygon.contains(makePointOrDie("3:3")), false);
  Assert.equal(polygon.contains(makePointOrDie("3:1")), false);
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  auto a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
  auto b = "# # 1:1, 1:3, 3:3, 3:1";
  expectResult(OpType.UNION, options, a, b,
      "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:3, 1:3 | 1:1, 1:3, 3:3 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
}

@("S2BooleanOperation.PolylineEdgeClosedPolygonEdgeOverlap") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.CLOSED);
  auto a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
  auto b = "# # 1:1, 1:3, 3:3, 3:1";
  expectResult(OpType.UNION, options, a, b,
      "# # 1:1, 1:3, 3:3, 3:1");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 1:1, 1:3, 3:3, 3:1");
}

@("S2BooleanOperation.PolygonVertexMatching") unittest {
  // This test shows that CrossingProcessor::ProcessEdgeCrossings() must set
  // a0_matches_polygon and a1_matches_polygon correctly even when (a0, a1)
  // itself is a polygon edge (or its sibling).  (It requires degenerate
  // polygon geometry to demonstrate this.)
  S2BooleanOperation.Options options;
  options.setPolylineModel(PolylineModel.CLOSED);
  options.setPolygonModel(PolygonModel.CLOSED);
  auto a = "# 0:0, 1:1 # ";
  auto b = "# # 0:0, 1:1";
  expectResult(OpType.UNION, options, a, b, "# # 0:0, 1:1");
}

@("S2BooleanOperation.PolylineEdgePolygonInterior") unittest {
  S2BooleanOperation.Options options;  // PolygonModel is irrelevant.
  // One normal and one degenerate polyline edge in the polygon interior, and
  // similarly for the polygon exterior.
  auto a = "# 1:1, 2:2 | 3:3, 3:3 | 6:6, 7:7 | 8:8, 8:8 # ";
  auto b = "# # 0:0, 0:5, 5:5, 5:0";
  expectResult(OpType.UNION, options, a, b,
      "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 2:2 | 3:3, 3:3 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 6:6, 7:7 | 8:8, 8:8 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
}

@("S2BooleanOperation.PolygonVertexOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.OPEN);
  auto a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
  auto b = "# # 0:0, 5:3, 5:2";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
}

@("S2BooleanOperation.PolygonVertexSemiOpenPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  auto a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
  auto b = "# # 0:0, 5:3, 5:2";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
}

@("S2BooleanOperation.PolygonVertexClosedPolygonVertex") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.CLOSED);
  auto a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
  auto b = "# # 0:0, 5:3, 5:2";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 0:0");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
  expectResult(OpType.DIFFERENCE, options, b, a,
      "# # 0:0, 5:3, 5:2");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
}

/+ TODO: Fix test.
@("S2BooleanOperation.PolygonEdgePolygonEdgeCrossing") unittest {
  // Two polygons whose edges cross at points interior to both edges.
  S2BooleanOperation.Options options = roundToE(2);
  auto a = "# # 0:0, 0:2, 2:2, 2:0";
  auto b = "# # 1:1, 1:3, 3:3, 3:1";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:2, 1:2, 1:3, 3:3, 3:1, 2:1, 2:0");
  // expectResult(OpType.INTERSECTION, options, a, b,
  //     "# # 1:1, 1:2, 2:2, 2:1");
  // expectResult(OpType.DIFFERENCE, options, a, b,
  //     "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0");
  // expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
  //     "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0; "
  //     ~ "1:2, 1:3, 3:3, 3:1, 2:1, 2:2");
}
+/

@("S2BooleanOperation.PolygonEdgeOpenPolygonEdgeOverlap") unittest {
  S2BooleanOperation.Options options;
  // One shape is a rectangle, the other consists of one triangle inside the
  // rectangle and one triangle outside the rectangle, where each triangle
  // shares one edge with the rectangle.  This implies that the edges are in
  // the same direction in one case and opposite directions in the other case.
  options.setPolygonModel(PolygonModel.OPEN);
  auto a = "# # 0:0, 0:4, 2:4, 2:0";
  auto b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0; 0:4, 1:5, 2:4");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 0:0, 1:1, 2:0");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
}

@("S2BooleanOperation.PolygonEdgeSemiOpenPolygonEdgeOverlap") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.SEMI_OPEN);
  auto a = "# # 0:0, 0:4, 2:4, 2:0";
  auto b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:4, 1:5, 2:4, 2:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 0:0, 1:1, 2:0");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1");
  // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
  // normalized, i.e. the output could contain siblings pairs (which can be
  // discarded using S2Builder::GraphOptions).
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
}

@("S2BooleanOperation.PolygonEdgeClosedPolygonEdgeOverlap") unittest {
  S2BooleanOperation.Options options;
  options.setPolygonModel(PolygonModel.CLOSED);
  auto a = "# # 0:0, 0:4, 2:4, 2:0";
  auto b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:4, 1:5, 2:4, 2:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 0:0, 1:1, 2:0; 0:4, 2:4");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1");
  // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are
  // normalized, i.e. the output could contain siblings pairs.
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
}

@("S2BooleanOperation.PolygonPolygonInterior") unittest {
  S2BooleanOperation.Options options;  // PolygonModel is irrelevant.
  // One loop in the interior of another polygon and one loop in the exterior.
  auto a = "# # 0:0, 0:4, 4:4, 4:0";
  auto b = "# # 1:1, 1:2, 2:2, 2:1; 5:5, 5:6, 6:6, 6:5";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:4, 4:4, 4:0; 5:5, 5:6, 6:6, 6:5");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 1:1, 1:2, 2:2, 2:1");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1; "
      ~ "5:5, 5:6, 6:6, 6:5");
}

@("S2BooleanOperation.PolygonEdgesDegenerateAfterSnapping") unittest {
  S2BooleanOperation.Options options = roundToE(0);
  auto a = "# # 0:-1, 0:1, 0.1:1, 0.1:-1";
  auto b = "# # -1:0.1, 1:0.1, 1:0, -1:0";
  // When snapping causes an output edge to become degenerate, it is still
  // emitted (since otherwise loops that contract to a single point would be
  // lost).  If the output layer doesn't want such edges, they can be removed
  // via DegenerateEdges::DISCARD or DISCARD_EXCESS.
  expectResult(OpType.UNION, options, a, b,
      "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | "
      ~ "-1:0, -1:0, 0:0, 1:0, 1:0, 0:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 0:0, 0:0, 0:0, 0:0");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | 0:0, 0:0");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:-1, 0:-1, 0:0, 0:1, 0:1, 0:0 | "
      ~ "-1:0, -1:0, 0:0, 1:0, 1:0, 0:0 | 0:0, 0:0, 0:0, 0:0");
}

///////////////////////////////////////////////////////////////////////////
// The remaining tests are intended to cover combinations of features or
// interesting special cases.

/+ TODO: Fix test.
@("S2BooleanOperation.ThreeOverlappingBars") unittest {
  // Two vertical bars and a horizontal bar that overlaps both of the other
  // bars and connects them.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BooleanOperation.Options options = roundToE(2);
  auto a = "# # 0:0, 0:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3";
  auto b = "# # 1:1, 1:4, 2:4, 2:1";
  expectResult(OpType.UNION, options, a, b,
      "# # 0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 1:1, 1:2, 2:2, 2:1; 1:3, 1:4, 2:4, 2:3");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; "
      ~ "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; "
      ~ "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3; "
      ~ "1:2, 1:3, 2:3, 2:2");
}
+/

/+ TODO: Fix test.
@("S2BooleanOperation.FourOverlappingBars") unittest {
  // Two vertical bars and two horizontal bars.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BooleanOperation.Options options = roundToE(2);
  auto a = "# # 1:88, 1:93, 2:93, 2:88; -1:88, -1:93, 0:93, 0:88";
  auto b = "# # -2:89, -2:90, 3:90, 3:89; -2:91, -2:92, 3:92, 3:91";
  expectResult(OpType.UNION, options, a, b,
      "# # -1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, "
      ~ "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, "
      ~ "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88; "
      ~ "0:90, 1:90, 1:91, 0:91" /*CW*/ );
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # 1:89, 1:90, 2:90, 2:89; 1:91, 1:92, 2:92, 2:91; "
      ~ "-1:89, -1:90, 0:90, 0:89; -1:91, -1:92, 0:92, 0:91");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # 1:88, 1:89, 2:89, 2:88; 1:90, 1:91, 2:91, 2:90; "
      ~ "1:92, 1:93, 2:93, 2:92; -1:88, -1:89, 0:89, 0:88; "
      ~ "-1:90, -1:91, 0:91, 0:90; -1:92, -1:93, 0:93, 0:92");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # 1:88, 1:89, 2:89, 2:88; -1:88, -1:89, 0:89, 0:88; "
      ~ "1:90, 1:91, 2:91, 2:90; -1:90, -1:91, 0:91, 0:90; "
      ~ "1:92, 1:93, 2:93, 2:92; -1:92, -1:93, 0:93, 0:92; "
      ~ "-2:89, -2:90, -1:90, -1:89; -2:91, -2:92, -1:92, -1:91; "
      ~ "0:89, 0:90, 1:90, 1:89; 0:91, 0:92, 1:92, 1:91; "
      ~ "2:89, 2:90, 3:90, 3:89; 2:91, 2:92, 3:92, 3:91");
}
+/

/+ TODO: Resume here.
@("S2BooleanOperation.OverlappingDoughnuts") unittest {
  // Two overlapping square doughnuts whose holes do not overlap.
  // This means that the union polygon has only two holes rather than three.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BooleanOperation.Options options = roundToE(1);
  auto a = "# # -1:-93, -1:-89, 3:-89, 3:-93; "
      ~ "0:-92, 2:-92, 2:-90, 0:-90" /*CW*/ ;
  auto b = "# # -3:-91, -3:-87, 1:-87, 1:-91; "
      ~ "-2:-90, 0:-90, 0:-88, -2:-88" /*CW*/ ;
  expectResult(OpType.UNION, options, a, b,
      "# # -1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93; "
      ~ "0:-92, 2:-92, 2:-90, 1:-90, 1:-91, 0:-91; " /*CW */
      ~ "-2:-90, -1:-90, -1:-89, 0:-89, 0:-88, -2:-88" /* CW */ );
  expectResult(OpType.INTERSECTION, options, a, b,
      "# # -1:-91, -1:-90, 0:-90, 0:-91; "
      ~ "0:-90, 0:-89, 1:-89, 1:-90");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, "
      ~ "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; "
      ~ "-1:-90, -1:-89, 0:-89, 0:-90");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, "
      ~ "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; "
      ~ "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88,-2:-88,-2:-90,-1:-90,-1:-91; "
      ~ "-1:-90, -1:-89, 0:-89, 0:-90; "
      ~ "1:-91, 0:-91, 0:-90, 1:-90");
}
+/

@("S2BooleanOperation.PolylineEnteringRectangle") unittest {
  // A polyline that enters a rectangle very close to one of its vertices.
  S2BooleanOperation.Options options = roundToE(1);
  auto a = "# 0:0, 2:2 #";
  auto b = "# # 1:1, 1:3, 3:3, 3:1";
  expectResult(OpType.UNION, options, a, b,
      "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 1:1, 2:2 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:0, 1:1 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
}

/+ TODO: Fix test.
@("S2BooleanOperation.PolylineCrossingRectangleTwice") unittest {
  // A polyline that crosses a rectangle in one direction, then moves to a
  // different side and crosses the rectangle in the other direction.  Note
  // that an extra vertex is added where the two polyline edges cross.
  S2BooleanOperation.Options options = roundToE(1);
  auto a = "# 0:-5, 0:5, 5:0, -5:0 #";
  auto b = "# # 1:1, 1:-1, -1:-1, -1:1";
  expectResult(OpType.UNION, options, a, b,
      "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 "
      ~ "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
  expectResult(OpType.INTERSECTION, options, a, b,
      "# 0:-1, 0:0, 0:1 | 1:0, 0:0, -1:0 #");
  expectResult(OpType.DIFFERENCE, options, a, b,
      "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 #");
  expectResult(OpType.SYMMETRIC_DIFFERENCE, options, a, b,
      "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 "
      ~ "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
}
+/
