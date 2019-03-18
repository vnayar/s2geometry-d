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

module s2.s2builder_test;

import s2.s2builder;

import s2.builder.graph;
import s2.builder.layer;
import s2.builder.util.s2polygon_layer;
import s2.builder.util.s2polyline_layer;
import s2.builder.util.s2polyline_vector_layer;
import s2.builder.util.snap_functions;
import s2.builder.util.testing;
import s2.logger;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell_id;
import s2.s2debug;
import s2.s2edge_crossings;
import s2.s2edge_distances;
import s2.s2error;
import s2.s2latlng;
import s2.s2loop;
import s2.s2point;
import s2.s2polygon;
import s2.s2polyline;
import s2.s2predicates;
import s2.s2testing;
import s2.s2text_format;

import fluent.asserts;

import std.algorithm;
import std.array;
import std.conv;
import std.format;
import std.math;
import std.typecons;

import std.stdio;

alias EdgeType = S2Builder.EdgeType;
alias InputEdgeId = Graph.InputEdgeId;


// Iteration multiplier for randomized tests.
int iterationMultiplier = 1;

void expectPolygonsEqual(in S2Polygon expected, in S2Polygon actual) {
  Assert.equal(actual, expected,
      "\nExpected:\n'" ~ .toString(expected) ~ "'"
      ~ "\nActual:\n'" ~ .toString(actual) ~ "'");
}

void expectPolygonsApproxEqual(in S2Polygon expected, in S2Polygon actual, S1Angle tolerance) {
  Assert.equal(expected.boundaryApproxEquals(actual, tolerance), true,
      "\nExpected:  " ~ .toString(expected)
      ~ "\nActual:    " ~ .toString(actual)
      ~ "\nTolerance: " ~ tolerance.degrees().to!string);
}

void expectPolylinesEqual(in S2Polyline expected, in S2Polyline actual) {
  Assert.equal(actual, expected,
      "\nExpected:\n" ~ .toString(expected)
      ~ "\nActual:\n" ~ .toString(actual));
}

@("S2Builder.AddShape") unittest {
  auto builder = new S2Builder(new S2Builder.Options());
  S2Polygon output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  auto input = makePolygonOrDie("0:0, 0:5, 5:5, 5:0; 1:1, 1:4, 4:4, 4:1");
  builder.addShape(input.index().shape(0));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  expectPolygonsEqual(input, output);
}

@("S2Builder.SimpleVertexMerging") unittest {
  // When IdentitySnapFunction is used (i.e., no special requirements on
  // vertex locations), check that vertices closer together than the snap
  // radius are merged together.

  S1Angle snap_radius = S1Angle.fromDegrees(0.5);
  auto builder = new S2Builder(new S2Builder.Options(new IdentitySnapFunction(snap_radius)));
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  auto input = makePolygonOrDie("0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9");
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  auto expected = makePolygonOrDie("0:0, 0:1, 1:0.9");
  expectPolygonsApproxEqual(expected, output, snap_radius);
}

@("S2Builder.SimpleS2CellIdSnapping") unittest {
  // When S2CellIdSnapFunction is used, check that all output vertices are the
  // centers of S2CellIds at the specified level level.

  int level = S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.fromDegrees(1.0));
  auto snap_function = new S2CellIdSnapFunction(level);
  auto builder = new S2Builder(new S2Builder.Options(snap_function));
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  auto input = makePolygonOrDie("2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(output.numLoops(), 1);
  const S2Loop loop = output.loop(0);
  for (int i = 0; i < loop.numVertices(); ++i) {
    Assert.equal(loop.vertex(i), new S2CellId(loop.vertex(i)).parent(level).toS2Point());
  }
  expectPolygonsApproxEqual(input, output, snap_function.snapRadius());
}

@("S2Builder.SimpleIntLatLngSnapping") unittest {
  auto builder = new S2Builder(new S2Builder.Options(new IntLatLngSnapFunction(0)));  // E0 coords
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  S2Polygon input = makePolygonOrDie(
      "2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, "
      ~ "5.22:3.88, 5.55:2.49, 4.49:2.51");
  S2Polygon expected = makePolygonOrDie(
      "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(output.numLoops(), 1);
  expectPolygonsEqual(expected, output);
}

@("S2Builder.VerticesMoveLessThanSnapRadius") unittest {
  // Check that chains of closely spaced vertices do not collapse into a
  // single vertex.

  S1Angle snap_radius = S1Angle.fromDegrees(1.0);
  auto builder = new S2Builder(new S2Builder.Options(new IdentitySnapFunction(snap_radius)));
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  // The spacing between input vertices is about 2*pi*20/1000 = 0.125 degrees.
  // The output vertices are spaced between 1 and 2 degrees apart; the average
  // spacing is about 1.33 degrees.
  auto input = new S2Polygon(
      S2Loop.makeRegularLoop(S2Point(1, 0, 0), S1Angle.fromDegrees(20.0), 1000));
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(output.numLoops(), 1);
  Assert.notLessThan(output.loop(0).numVertices(), 90);
  Assert.notGreaterThan(output.loop(0).numVertices(), 100);
  Assert.equal(output.boundaryNear(input, snap_radius), true);
}

@("S2Builder.MinEdgeVertexSeparation") unittest {
  // Check that edges are separted from non-incident vertices by at least
  // min_edge_vertex_separation().  This requires adding new vertices (not
  // present in the input) in some cases.

  // The input is a skinny right triangle with two legs of length 10 and 1,
  // and whose diagonal is subdivided into 10 short edges.  Using a snap
  // radius of 0.5, about half of the long leg is snapped onto the diagonal
  // (which causes that part of the polygon to be removed).  But the real
  // problem is that the remaining part of the long leg gets too close to the
  // remaining vertices on the diagonal, i.e. it would violate the minimum
  // edge-vertex separation guarantee.  S2Builder handles this by creating at
  // least one vertex along the original long leg, to keep the snapped edge
  // far enough away from the diagonal.
  S2Polygon input = makePolygonOrDie(
      "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0");
  S2Polygon expected = makePolygonOrDie(
      "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 4.00021862252687:0");
  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(0.5)));
  auto builder = new S2Builder(options);
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  expectPolygonsApproxEqual(expected, output, S1Angle.fromRadians(1e-15));
}

@("S2Builder.IdempotencySnapsInadequatelySeparatedVertices") unittest {
  // This test checks that when vertices are closer together than
  // min_vertex_separation() then they are snapped together even when
  // options.idempotent() is true.
  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(1.0)));
  auto builder = new S2Builder(options);
  auto output = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output));
  builder.addPolyline(makePolylineOrDie("0:0, 0:0.9, 0:2"));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  string expected = "0:0, 0:2";
  Assert.equal(expected, .toString(output));
}

@("S2Builder.IdempotencySnapsIdenticalVerticesWithZeroSnapRadius") unittest {
  // This test checks that even when the snap radius is zero, identical
  // vertices are snapped together.
  auto builder = new S2Builder(new S2Builder.Options());
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  builder.addPolyline(makePolylineOrDie("0:1, 1:0"));
  builder.addPolyline(makePolylineOrDie("0:0, 0:1"));
  builder.addEdge(makePointOrDie("0:1"), makePointOrDie("0:1"));
  builder.addPolyline(makePolylineOrDie("1:0, 0:0"));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  string expected = "0:0, 0:1, 1:0";
  Assert.equal(expected, .toString(output));
}

@("S2Builder.IdempotencySnapsIdenticalVerticesWithZeroSnapRadiusEdgeSplitting") unittest {
  // This test checks that identical vertices are snapped together even when
  // the snap radius is zero and options.split_crossing_edges() is true.
  auto options = new S2Builder.Options();
  options.setSplitCrossingEdges(true);
  auto builder = new S2Builder(options);
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  builder.addPolyline(makePolylineOrDie("0:1, 1:0"));
  builder.addPolyline(makePolylineOrDie("0:0, 0:1"));
  builder.addEdge(makePointOrDie("0:1"), makePointOrDie("0:1"));
  builder.addPolyline(makePolylineOrDie("1:0, 0:0"));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  string expected = "0:0, 0:1, 1:0";
  Assert.equal(expected, .toString(output));
}

@("S2Builder.IdempotencySnapsUnsnappedVertices") unittest {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that S2Builder detects vertices
  // that are not at a valid location returned by the given snap function.

  // In this example we snap two vertices to integer lat/lng coordinates.  The
  // two vertices are far enough apart (more than min_vertex_separation) so
  // that they might be the result of a previous snapping operation, but one
  // of the two vertices does not have integer lat/lng coordinates.  We use
  // internal knowledge of how snap sites are chosen (namely, that candidates
  // are considered in S2CellId order) to construct two different cases, one
  // where the snapped vertex is processed first and one where the unsnapped
  // vertex is processed first.  This exercises two different code paths.
  auto snap_function = new IntLatLngSnapFunction(0);
  Assert.notLessThan(snap_function.snapRadius(), S1Angle.fromDegrees(0.7));
  Assert.notGreaterThan(snap_function.minVertexSeparation(), S1Angle.fromDegrees(0.35));
  auto builder = new S2Builder(new S2Builder.Options(snap_function));

  // In this example, the snapped vertex (0, 0) is processed first and is
  // selected as a Voronoi site (i.e., output vertex).  The second vertex is
  // closer than min_(), therefore it is snapped to the first vertex
  // and the polyline becomes degenerate.
  S2Point a = S2LatLng.fromDegrees(0, 0).toS2Point();
  S2Point b = S2LatLng.fromDegrees(0.01, 0.6).toS2Point();
  Assert.lessThan(S2CellId(a), S2CellId(b));
  auto input1 = new S2Polyline([a, b]);
  auto output1 = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output1));
  builder.addPolyline(input1);
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal("0:0, 0:1", .toString(output1));

  // In this example the unsnapped vertex is processed first and is snapped to
  // (0, 0).  The second vertex is further than snap_radius() away, so it is
  // also snapped (which does nothing) and is left at (0, 1).
  S2Point c = S2LatLng.fromDegrees(0.01, 0.4).toS2Point();
  S2Point d = S2LatLng.fromDegrees(0, 1).toS2Point();
  Assert.lessThan(S2CellId(c), S2CellId(d));
  auto input2 = new S2Polyline([c, d]);
  auto output2 = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output2));
  builder.addPolyline(input2);
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal("0:0, 0:1", .toString(output2));
}

@("S2Builder.IdempotencySnapsEdgesWithTinySnapRadius") unittest {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that S2Builder detects edges that
  // are too close to vertices even when the snap radius is very small
  // (e.g., S2::kIntersectionError).
  //
  // Previously S2Builder used a conservative approximation to decide whether
  // edges were too close to vertices; unfortunately this meant that when the
  // snap radius was very small then no snapping would be done at all, because
  // even an edge/vertex distance of zero was considered far enough apart.
  //
  // This tests that the current code (which uses exact predicates) handles
  // this situation correctly (i.e., that an edge separated from a
  // non-incident vertex by a distance of zero cannot be the output of a
  // previous snapping operation).
  auto options = new S2Builder.Options();
  options.setSnapFunction(new IdentitySnapFunction(INTERSECTION_ERROR));
  auto layer_options = S2PolylineVectorLayer.Options();
  layer_options.setDuplicateEdges(S2PolylineVectorLayer.Options.DuplicateEdges.MERGE);
  auto builder = new S2Builder(options);
  S2Polyline[] output;
  builder.startLayer(new S2PolylineVectorLayer(&output, layer_options));
  builder.addPolyline(makePolylineOrDie("0:0, 0:10"));
  builder.addPolyline(makePolylineOrDie("0:5, 0:7"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal(output.length, 1);
  Assert.equal("0:0, 0:5, 0:7, 0:10", .toString(output[0]));
}

@("S2Builder.IdempotencyDoesNotSnapAdequatelySeparatedEdges") unittest {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that when an edge is further away
  // than min_edge_vertex_separation() then no snapping is done.
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));
  options.setIdempotent(true);  // Test fails if this is "false".
  auto builder = new S2Builder(options);
  S2Polygon output1 = new S2Polygon();
  S2Polygon output2 = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output1));
  builder.addPolygon(makePolygonOrDie("1.49:0, 0:2, 0.49:3"));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  string expected = "1:0, 0:2, 0:3";
  Assert.equal(.toString(output1), expected);
  builder.startLayer(new S2PolygonLayer(output2));
  builder.addPolygon(output1);
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(expected, .toString(output2));
}

@("S2Builder.kMaxSnapRadiusCanSnapAtLevel0") unittest {
  // Verify that kMaxSnapRadius will allow snapping at S2CellId level 0.
  Assert.notGreaterThan(S2CellIdSnapFunction.minSnapRadiusForLevel(0),
            S2Builder.SnapFunction.kMaxSnapRadius());
}

@("S2Builder.S2CellIdSnappingAtAllLevels") unittest {
  auto input = makePolygonOrDie("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");
  for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
    auto snap_function = new S2CellIdSnapFunction(level);
    auto builder = new S2Builder(new S2Builder.Options(snap_function));
    auto output = new S2Polygon();
    builder.startLayer(new S2PolygonLayer(output));
    builder.addPolygon(input);
    S2Error error;
    Assert.equal(builder.build(error), true, error.toString());
    Assert.equal(output.isValid(), true);
    // The ApproxContains calls below are not guaranteed to succeed in general
    // because ApproxContains works by snapping both polygons together using
    // the given tolerance and then checking for containment.  Since
    // ApproxContains snaps to an arbitrary subset of the input vertices
    // rather than to S2CellId centers at the current level, this means that
    // corresponding vertices in "input" and "output" can snap to different
    // sites, which causes the containment test to fail.  Nevertheless, by
    // using a larger tolerance of 2 * snap_radius, all calls in this test
    // succeed (and would be likely to succeed in other similar tests).
    // (To guarantee correctness we would need to use S2CellIdSnapFunction
    // within the ApproxContains implementation.)
    S1Angle tolerance = min(2 * snap_function.snapRadius(),
                            snap_function.kMaxSnapRadius());
    Assert.equal(output.approxContains(input, tolerance), true);
    Assert.equal(input.approxContains(output, tolerance), true);
  }
}

@("S2Builder.SnappingDoesNotRotateVertices") unittest {
  // This is already tested extensively elsewhere.
  auto input = makePolygonOrDie(
      "49.9305505:-124.8345463, 49.9307448:-124.8299657, "
      ~ "49.9332101:-124.8301996, 49.9331224:-124.8341368; "
      ~ "49.9311087:-124.8327042, 49.9318176:-124.8312621, "
      ~ "49.9318866:-124.8334451");
  auto options = new S2Builder.Options(new S2CellIdSnapFunction());
  auto builder = new S2Builder(options);
  auto output1 = new S2Polygon();
  auto output2 = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output1));
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  // This checks that the vertices are in the same cyclic order, and that
  // vertices have not moved by more than "snap_radius".
  expectPolygonsApproxEqual(input, output1, options.snapFunction().snapRadius());

  // Check that snapping twice doesn't rotate the vertices.  This also
  // verifies that S2Builder can be used again after Build() is called.
  builder.startLayer(new S2PolygonLayer(output2));
  builder.addPolygon(output1);
  Assert.equal(builder.build(error), true, error.toString());
  expectPolygonsEqual(output1, output2);
}

@("S2Builder.SelfIntersectingPolyline") unittest {
  // Check that when two edges of a polyline cross, the intersection point is
  // added to both edges.

  auto options = new S2Builder.Options();
  auto snap_function = new IntLatLngSnapFunction(1);  // Snap to E1 coordinates
  options.setSnapFunction(snap_function);
  options.setSplitCrossingEdges(true);
  auto builder = new S2Builder(options);
  auto output = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output));
  auto input = makePolylineOrDie("3:1, 1:3, 1:1, 3:3");
  auto expected = makePolylineOrDie("3:1, 2:2, 1:3, 1:1, 2:2, 3:3");
  builder.addPolyline(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  expectPolylinesEqual(output, expected);
}

@("S2Builder.SelfIntersectingPolygon") unittest {
  // Check that when two edge of a polygon cross, the intersection point is
  // added to both edges, and that the resulting (undirected) edges can be
  // assembled into a valid polygon.

  auto snap_function = new IntLatLngSnapFunction(1);  // Snap to E1 coordinates
  auto options = new S2Builder.Options();
  options.setSnapFunction(snap_function);
  options.setSplitCrossingEdges(true);
  auto builder = new S2Builder(options);
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output, S2PolygonLayer.Options(EdgeType.UNDIRECTED)));
  auto input = makePolylineOrDie("3:1, 1:3, 1:1, 3:3, 3:1");
  auto expected = makePolygonOrDie("1:1, 1:3, 2:2; 3:3, 3:1, 2:2");
  builder.addPolyline(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  expectPolygonsEqual(expected, output);
}

@("S2Builder.TieBreakingIsConsistent") unittest {
  // Check that when an edge passes between two equally distant vertices, that
  // the choice of which one to snap to does not depend on the edge direction.

  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(2.0)));
  options.setIdempotent(false);
  auto builder = new S2Builder(options);
  builder.forceVertex(S2LatLng.fromDegrees(1, 0).toS2Point());
  builder.forceVertex(S2LatLng.fromDegrees(-1, 0).toS2Point());
  auto output1 = new S2Polyline();
  auto output2 = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output1));
  builder.addPolyline(makePolylineOrDie("0:-5, 0:5"));
  builder.startLayer(new S2PolylineLayer(output2));
  builder.addPolyline(makePolylineOrDie("0:5, 0:-5"));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(output1.numVertices(), 3);
  Assert.equal(output2.numVertices(), 3);
  for (int i = 0; i < 3; ++i) {
    Assert.equal(output1.vertex(i), output2.vertex(2 - i));
  }
}

// Verifies that two graphs have the same vertices and edges.
void expectGraphsEqual(in Graph expected, in Graph actual) {
  Assert.equal(actual.vertices(), expected.vertices());
  Assert.equal(actual.edges(), expected.edges());
  Assert.equal(actual.inputEdgeIdSetIds(), expected.inputEdgeIdSetIds());
}

// This layer makes both a shallow and a deep copy of the Graph object passed
// to its Build() method and appends them to two vectors.  Furthermore, it
// verifies that the shallow and deep copies of any graphs previously appended
// to those vectors are still identical.
class GraphPersistenceLayer : Layer {
public:
  this(
      GraphOptions graph_options,
      Graph[]* graphs,
      GraphClone[]* clones) {
    _graphOptions = graph_options;
    _graphs = graphs;
    _clones = clones;
  }

  override
  GraphOptions graphOptions() {
    return _graphOptions;
  }

  override
  void build(Graph g, ref S2Error error) {
    // Verify that all graphs built so far are unchanged.
    for (int i = 0; i < _graphs.length; ++i) {
      expectGraphsEqual((*_clones)[i].graph(), (*_graphs)[i]);
    }
    *_graphs ~= g;
    *_clones ~= new GraphClone(g);
  }

private:
  GraphOptions _graphOptions;
  Graph[]* _graphs;  // Shallow copies.
  GraphClone[]* _clones;    // Deep copies.
}

@("S2Builder.GraphPersistence") unittest {
  // Ensure that the Graph objects passed to S2Builder::Layer::Build() methods
  // remain valid until all layers have been built.
  Graph[] graphs;
  GraphClone[] clones;
  auto builder = new S2Builder(new S2Builder.Options());
  for (int i = 0; i < 20; ++i) {
    builder.startLayer(new GraphPersistenceLayer(new GraphOptions(), &graphs, &clones));
    for (int n = S2Testing.rnd.uniform(10); n > 0; --n) {
      builder.addEdge(S2Testing.randomPoint(), S2Testing.randomPoint());
    }
  }
  S2Error error;
  Assert.equal(builder.build(error), true);
}

void checkPolylineLayers(
    in string[] input_strs,
    in string[] expected_strs,
    in S2PolylineLayer.Options layer_options,
    S2Builder.Options builder_options = null) {
  //SCOPED_TRACE(layer_options.edge_type() == EdgeType::DIRECTED ?
  //             "DIRECTED" : "UNDIRECTED");
  if (builder_options is null) builder_options = new S2Builder.Options();
  auto builder = new S2Builder(builder_options);
  S2Polyline[] output;
  foreach (input_str; input_strs) {
    output ~= new S2Polyline();
    builder.startLayer(new S2PolylineLayer(output.back(), layer_options));
    builder.addPolyline(makePolylineOrDie(input_str));
  }
  S2Error error;
  Assert.equal(builder.build(error), true);
  string[] output_strs;
  foreach (polyline; output) {
    output_strs ~= .toString(polyline);
  }
  Assert.equal(output_strs.join("; "), expected_strs.join("; "));
}

void checkPolylineVector(
    in string[] input_strs,
    in string[] expected_strs,
    in S2PolylineVectorLayer.Options layer_options,
    S2Builder.Options builder_options = null) {
  if (builder_options is null) builder_options = new S2Builder.Options();
  auto builder = new S2Builder(builder_options);
  S2Polyline[] output;
  builder.startLayer(new S2PolylineVectorLayer(&output, layer_options));
  foreach (input_str; input_strs) {
    builder.addPolyline(makePolylineOrDie(input_str));
  }
  S2Error error;
  Assert.equal(builder.build(error), true);
  string[] output_strs;
  foreach (polyline; output) {
    output_strs ~= .toString(polyline);
  }
  Assert.equal(output_strs.join("; "), expected_strs.join("; "));
}

void checkPolylineLayersBothEdgeTypes(
    in string[] input_strs,
    in string[] expected_strs,
    S2PolylineLayer.Options layer_options,  // by value
    S2Builder.Options builder_options = null) {
  if (builder_options is null) builder_options = new S2Builder.Options();
  layer_options.setEdgeType(EdgeType.DIRECTED);
  checkPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
  layer_options.setEdgeType(EdgeType.UNDIRECTED);
  checkPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
}

@("S2Builder.SimplifyOneEdge") unittest {
  // Simplify a perturbed edge chain into a single edge.

  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(1.0)));
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["0:0, 1:0.5, 2:-0.5, 3:0.5, 4:-0.5, 5:0"],
      ["0:0, 5:0"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyTwoLayers") unittest {
  // Construct two layers, each containing a polyline that could be simplified
  // to a single edge on its own.  However the two polylines actually cross,
  // so make sure that the output still contains the intersection vertex.

  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(0.5)));
  options.setSplitCrossingEdges(true);
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["-2:-1, -1:0, 1:0, 2:1", "1:-2, 0:-1, 0:1, -1:2"],
      ["-2:-1, 0:0, 2:1", "1:-2, 0:0, -1:2"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyOneLoop") unittest {
  // Simplify a regular loop with 1000 vertices and a radius of 20 degrees.
  // Turning on edge chain simplification yields a dramatically smaller number
  // of vertices than snapping alone (10 vertices vs 95 vertices using a snap
  // radius of 1 degree).  This is because snapping alone yields vertices that
  // stay within 1 degree of the input *vertices*, while simplifying edge
  // chains yields edges that stay within 1 degree of the input *edges*.

  for (int i = 0; i < 2; ++i) {
    EdgeType edge_type = cast(EdgeType) i;
    S1Angle snap_radius = S1Angle.fromDegrees(1.0);
    auto options = new S2Builder.Options(new IdentitySnapFunction(snap_radius));
    options.setSimplifyEdgeChains(true);
    auto builder = new S2Builder(options);
    auto output = new S2Polygon();
    builder.startLayer(new S2PolygonLayer(output, S2PolygonLayer.Options(edge_type)));
    // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
    auto input = new S2Polygon(
        S2Loop.makeRegularLoop(S2Point(1, 0, 0), S1Angle.fromDegrees(20.0), 1000));
    builder.addPolygon(input);
    S2Error error;
    Assert.equal(builder.build(error), true, error.toString());
    Assert.equal(output.numLoops(), 1);
    Assert.notLessThan(output.loop(0).numVertices(), 10);
    Assert.notGreaterThan(output.loop(0).numVertices(), 12);
    Assert.equal(output.boundaryNear(input, snap_radius), true);
  }
}

@("S2Builder.SimplifyOppositeDirections") unittest {
  // We build two layers with two polylines that follow the same circular arc
  // in opposite directions, and verify that they are snapped identically.
  // (The snap radius is adjusted so that the arc is simplified into a long
  // edge and a short edge, and therefore we would get a different result if
  // the two layers followed the edge chain in different directions.)

  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(0.5)));
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["-4:0.83, -3:0.46, -2:0.2, -1:0.05, 0:0, 1:0.5, 2:0.2, 3:0.46, 4:0.83",
            "4:.83, 3:.46, 2:.2, 1:.05, 0:0, -1:.5, -2:.2, -3:.46, -4:.83"],
      ["-4:0.83, -2:0.2, 4:0.83", "4:0.83, -2:0.2, -4:0.83"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyKeepsEdgeVertexSeparation") unittest {
  // We build two layers each containing a polyline, such that the polyline in
  // the first layer could be simplified to a straight line except that then
  // it would create an intersection with the second polyline.

  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(1.0)));
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"],
      ["0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyBacktrackingEdgeChain") unittest {
  // Test simplifying an edge chain that backtracks on itself.
  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromDegrees(0.5)));
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 4:0, 3:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0"],
      ["0:0, 2:0, 5:0, 2:0, 5:0, 7:0"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyLimitsEdgeDeviation") unittest {
  // Make sure that simplification does not create long edges such that the
  // midpoint of the edge might be further than max_edge_deviation() from an
  // input edge.  In the example below, vertices are snapped to integer
  // lat/lng coordinates, and the snap radius is approximately 0.707 degrees.
  // Snapping moves the input vertices perpendicular to the input edge by just
  // slightly less than the snap radius (0.693 degrees).  Now the midpoint of
  // the snapped edge is about 0.98 degrees from the input edge, which causes
  // an extra site to be added at the midpoint of the original edge.
  //
  // When simplify_edge_chains() is enabled, then usually an extra site like
  // this would be simplified away (because the simplified edge would still be
  // within snap_radius() of all the input vertices) except that there is an
  // explicit check in S2Builder that prevents this.  (If the check is removed
  // then this test fails.)

  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));  // E0 coordinates
  options.setSimplifyEdgeChains(true);
  checkPolylineLayersBothEdgeTypes(
      ["-30.49:-29.51, 29.51:30.49"], ["-30:-30, -1:1, 30:30"],
      S2PolylineLayer.Options(), options);
}

@("S2Builder.SimplifyPreservesTopology") unittest {
  // Crate several nested concentric loops, and verify that the loops are
  // still nested after simplification.

  // FIXME: Test breaks once kNumLoops is greater than 4.
  //const int kNumLoops = 20;
  const int kNumLoops = 4;
  //const int kNumVerticesPerLoop = 1000;
  const int kNumVerticesPerLoop = 5;
  const S1Angle kBaseRadius = S1Angle.fromDegrees(5.0);
  const S1Angle kSnapRadius = S1Angle.fromDegrees(0.1);
  auto options = new S2Builder.Options(new IdentitySnapFunction(kSnapRadius));
  options.setSimplifyEdgeChains(true);
  auto builder = new S2Builder(options);
  S2Polygon[] input, output;
  for (int j = 0; j < kNumLoops; ++j) {
    // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
    S1Angle radius = kBaseRadius + 0.7 * j * j / kNumLoops * kSnapRadius;
    input ~= new S2Polygon(S2Loop.makeRegularLoop(
        S2Point(1, 0, 0), radius, kNumVerticesPerLoop));
    output ~= new S2Polygon();
    builder.startLayer(new S2PolygonLayer(output.back()));
    builder.addPolygon(input.back());
  }
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  for (int j = 0; j < kNumLoops; ++j) {
    Assert.equal(output[j].boundaryNear(input[j], kSnapRadius), true);
    if (j > 0) Assert.equal(output[j].contains(output[j - 1]), true);
  }
}

@("S2Builder.SimplifyRemovesSiblingPairs") unittest {
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));  // E0 coords
  auto layer_options = S2PolylineVectorLayer.Options();
  layer_options.setSiblingPairs(GraphOptions.SiblingPairs.DISCARD);

  // Check that there is no sibling pair without simplification.
  checkPolylineVector(
      ["0:0, 0:10", "0:10, 0.6:5, 0:0"],
      ["0:0, 0:10, 1:5, 0:0"], layer_options, options);

  // Now check that (1) simplification produces a sibling pair,
  // and (2) the sibling pair is removed (since we requested it).
  options.setSimplifyEdgeChains(true);
  checkPolylineVector(
      ["0:0, 0:10", "0:10, 0.6:5, 0:0"],
      [], layer_options, options);
}

@("S2Builder.SimplifyMergesDuplicateEdges") unittest {
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));  // E0 coords
  S2PolylineVectorLayer.Options layer_options;
  layer_options.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);

  // Check that there are no duplicate edges without simplification.
  checkPolylineVector(
      ["0:0, 0:10", "0:0, 0.6:5, 0:10"],
      ["0:0, 0:10", "0:0, 1:5, 0:10"], layer_options, options);

  // Now check that (1) simplification produces a duplicate edge pair,
  // and (2) the duplicate pair is merged (since we requested it).
  options.setSimplifyEdgeChains(true);
  checkPolylineVector(
      ["0:0, 0:10", "0:0, 0.6:5, 0:10"],
      ["0:0, 0:10"], layer_options, options);
}

@("S2Builder.SimplifyKeepsForcedVertices") unittest {
  auto options = new S2Builder.Options(new IdentitySnapFunction(S1Angle.fromRadians(1e-15)));
  options.setSimplifyEdgeChains(true);
  auto builder = new S2Builder(options);
  auto output = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output));
  builder.addPolyline(makePolylineOrDie("0:0, 0:1, 0:2, 0:3"));
  builder.forceVertex(makePointOrDie("0:1"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal(.toString(output), "0:0, 0:1, 0:3");
}

// A set of (edge string, vector<InputEdgeId>) pairs representing the
// InputEdgeIds attached to the edges of a graph.  Edges are in
// s2textformat::ToString() format, such as "1:3, 4:5".
alias EdgeInputEdgeId = Tuple!(string, const(int[]));
alias EdgeInputEdgeIds = EdgeInputEdgeId[];

S2Error.Code INPUT_EDGE_ID_MISMATCH = S2Error.Code.USER_DEFINED_START;

class InputEdgeIdCheckingLayer : Layer {
public:
  this(EdgeInputEdgeIds expected, GraphOptions graph_options) {
    _expected = expected;
    _graphOptions = graph_options;
  }

  override
  GraphOptions graphOptions() {
    return _graphOptions;
  }

  override
  void build(Graph g, ref S2Error error) {
    EdgeInputEdgeIds actual;
    S2Point[] vertices;
    for (Graph.EdgeId e = 0; e < g.numEdges(); ++e) {
      vertices.length = 0;
      vertices ~= g.vertex(g.edge(e)[0]);
      vertices ~= g.vertex(g.edge(e)[1]);
      string edge = .toString([g.vertex(g.edge(e)[0]), g.vertex(g.edge(e)[1])]);
      const(int[]) ids = g.inputEdgeIds(e);
      actual ~= EdgeInputEdgeId(edge, ids);
    }
    // This comparison doesn't consider multiplicity, but that's fine.
    string missing, extra;
    foreach (p; _expected) {
      if (count(actual, p) > 0) continue;
      missing ~= toString(p);
    }
    foreach (p; actual) {
      if (count(_expected, p) > 0) continue;
      extra ~= toString(p);
    }
    if (!missing.empty() || !extra.empty()) {
      error.initialize(INPUT_EDGE_ID_MISMATCH, "Missing:\n%sExtra:\n%s\n", missing, extra);
    }
  }

private:
  string toString(EdgeInputEdgeId p) {
    string r = "  (" ~ p[0] ~ ")={";
    if (!p[1].empty()) {
      foreach (int id; p[1]) {
        r ~= id.to!string ~ ", ";
      }
      r = r[0 .. $-2];
    }
    r ~= "}\n";
    return r;
  }

  EdgeInputEdgeIds _expected;
  GraphOptions _graphOptions;
}

void checkInputEdgeIds(
    string[] input_strs, EdgeInputEdgeIds expected,
    GraphOptions graph_options, S2Builder.Options options) {
  auto builder = new S2Builder(options);
  builder.startLayer(new InputEdgeIdCheckingLayer(expected, graph_options));
  foreach (input_str; input_strs) {
    builder.addPolyline(makePolylineOrDie(input_str));
  }
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
}

@("S2Builder.InputEdgeIdAssignment") unittest {
  // Check that input edge ids are assigned in order.
  checkInputEdgeIds(["0:0, 0:1, 0:2"], [EdgeInputEdgeId("0:0, 0:1", [0]), EdgeInputEdgeId("0:1, 0:2", [1])],
                   new GraphOptions(), new S2Builder.Options());
}

@("S2Builder.UndirectedSiblingsDontHaveInputEdgeIds") unittest {
  // Check that the siblings of undirected edges do not have InputEdgeIds.
  auto graph_options = new GraphOptions();
  graph_options.setEdgeType(EdgeType.UNDIRECTED);
  checkInputEdgeIds(["0:0, 0:1, 0:2"], [
          EdgeInputEdgeId("0:0, 0:1", [0]),
          EdgeInputEdgeId("0:1, 0:2", [1]),
          EdgeInputEdgeId("0:1, 0:0", []),
          EdgeInputEdgeId("0:2, 0:1", [])
      ], graph_options, new S2Builder.Options());
}

@("S2Builder.CreatedSiblingsDontHaveInputEdgeIds") unittest {
  // Check that edges created by SiblingPairs::CREATE do not have
  // InputEdgeIds.
  auto graph_options = new GraphOptions();
  graph_options.setSiblingPairs(GraphOptions.SiblingPairs.CREATE);
  checkInputEdgeIds(["0:0, 0:1, 0:2"],
      [EdgeInputEdgeId("0:0, 0:1", [0]), EdgeInputEdgeId("0:1, 0:2", [1])],
      new GraphOptions(), new S2Builder.Options());
}

@("S2Builder.EdgeMergingDirected") unittest {
  // Tests that input edge ids are merged when directed edges are merged.
  auto graph_options = new GraphOptions();
  graph_options.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);
  checkInputEdgeIds(["0:0, 0:1", "0:0, 0:1"], [EdgeInputEdgeId("0:0, 0:1", [0, 1])],
      graph_options, new S2Builder.Options());
}

@("S2Builder.EdgeMergingUndirected") unittest {
  // Tests that input edge ids are merged when undirected edges are merged.
  auto graph_options = new GraphOptions();
  graph_options.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);
  graph_options.setSiblingPairs(GraphOptions.SiblingPairs.KEEP);
  checkInputEdgeIds(["0:0, 0:1, 0:2", "0:0, 0:1", "0:2, 0:1"], [
          EdgeInputEdgeId("0:0, 0:1", [0, 2]),
          EdgeInputEdgeId("0:1, 0:2", [1]),
          EdgeInputEdgeId("0:2, 0:1", [3])
      ], graph_options, new S2Builder.Options());
}

@("S2Builder.SimplifyDegenerateEdgeMergingEasy") unittest {
  // Check that when an input edge is snapped to a chain that includes
  // degenerate edges, and the edge chain is simplified, that the InputEdgeIds
  // attached to those degenerate edges are transferred to the simplified
  // edge.  For example (using integers for vertices), an edge chain 1->2,
  // 2->2, 2->3 that is simplified to 1->3 should get the InputEdgeIds
  // associated with all three original edges.  (This ensures that the labels
  // attached to those edges are also transferred.)
  //
  // This also tests that degenerate edges at the start and end of the
  // simplified chain are *not* merged.  (It's up to the output layer to
  // decide what to do with these edges.  The only reason we merge degenerate
  // edges in the interior of the interior of the simplified edge is because
  // those edges are being removed from the graph.)
  auto graph_options = new GraphOptions();
  graph_options.setDegenerateEdges(GraphOptions.DegenerateEdges.KEEP);
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));
  options.setSimplifyEdgeChains(true);
  checkInputEdgeIds(["0:0, 0:0.1, 0:1.1, 0:1, 0:0.9, 0:2, 0:2.1"], [
          EdgeInputEdgeId("0:0, 0:0", [0]),
          EdgeInputEdgeId("0:0, 0:2", [1, 2, 3, 4]),
          EdgeInputEdgeId("0:2, 0:2", [5])
    ], graph_options, options);
}

@("S2Builder.SimplifyDegenerateEdgeMergingHard") unittest {
  // This is a harder version of the test above.  Now there are several edge
  // chains that overlap each other in both directions, and several degenerate
  // edges at that middle vertex.  This tests that if exactly one edge chain
  // contains a degenerate edge in input edge order (e.g., the input order was
  // AB, BB, BC), then the degenerate edge is assigned to that edge chain.
  // Otherwise the edge is assigned to an arbitrary chain.
  auto graph_options = new GraphOptions();  // Default options keep everything.
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));
  options.setSimplifyEdgeChains(true);
  string[] input = [
    "0:1, 0:1.1", "0:0, 0:1, 0:2",  // Degenerate edge defined before chain
    "0:0, 0:0.9, 0:1, 0:1.1, 0:2",  // Degenerate edge defined in chain
    "0:2, 0:1, 0:0.9, 0:0",         // Defined in chain, chain reversed
    "0:2, 0:1, 0:0", "0:1.1, 0:1", "0:1, 0:1.1",  // Defined after chain
  ];
  EdgeInputEdgeIds expected = [
      EdgeInputEdgeId("0:0, 0:2", [0, 1, 2]),
      EdgeInputEdgeId("0:0, 0:2", [3, 4, 5, 6]),
      EdgeInputEdgeId("0:2, 0:0", [7, 8, 9]),
      EdgeInputEdgeId("0:2, 0:0", [10, 11, 12, 13])
  ];
  checkInputEdgeIds(input, expected, graph_options, options);

  // Now try the same test with undirected edges.  This results in four more
  // simplified edges that are not labelled with any input edge ids.
  expected ~= [
      EdgeInputEdgeId("0:0, 0:2", []),
      EdgeInputEdgeId("0:0, 0:2", []),
      EdgeInputEdgeId("0:2, 0:0", []),
      EdgeInputEdgeId("0:2, 0:0", [])
  ];
  graph_options.setEdgeType(EdgeType.UNDIRECTED);
  checkInputEdgeIds(input, expected, graph_options, options);
}

@("S2Builder.SimplifyDegenerateEdgeMergingMultipleLayers") unittest {
  // Check that degenerate edges are assigned to an edge in the correct layer
  // when multiple edge chains in different layers are simplified in the same
  // way (i.e., yielding a set of identical or reversed edges in different
  // layers).
  auto graph_options = new GraphOptions();  // Default options keep everything.
  auto options = new S2Builder.Options(new IntLatLngSnapFunction(0));
  options.setSimplifyEdgeChains(true);

  // Note below that the edge chains in different layers have different vertex
  // locations, different number of interior vertices, different degenerate
  // edges, etc, and yet they can all be simplified together.
  string[][] input = [ [
      "0.1:5, 0:5.2", "0.1:0, 0:9.9",   // Defined before chain
      "0:10.1, 0:0.1", "0:3.1, 0:2.9",  // Defined after chain
    ], [
      "0.1:3, 0:3.2", "-0.1:0, 0:4.1, 0:9.9",  // Defined before chain
      "0.1:9.9, 0:7, 0.1:6.9, 0.1:0.2",        // Defined inside chain
    ], [
      "0.2:0.3, 0.1:6, 0:5.9, 0.1:10.2",       // Defined inside chain
      "0.1:0.1, 0:9.8", "0.1:2, 0:2.1",        // Defined after chain
    ]
  ];
  EdgeInputEdgeIds[] expected = [ [
      EdgeInputEdgeId("0:0, 0:10", [0, 1]), EdgeInputEdgeId("0:10, 0:0", [2, 3])
    ], [
      EdgeInputEdgeId("0:0, 0:10", [4, 5, 6]), EdgeInputEdgeId("0:10, 0:0", [7, 8, 9])
    ], [
      EdgeInputEdgeId("0:0, 0:10", [10, 11, 12]), EdgeInputEdgeId("0:0, 0:10", [13, 14])
    ]
  ];
  auto builder = new S2Builder(options);
  for (int i = 0; i < input.length; ++i) {
    builder.startLayer(new InputEdgeIdCheckingLayer(expected[i], graph_options));
    foreach (input_str; input[i]) {
      builder.addPolyline(makePolylineOrDie(input_str));
    }
  }
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
}

@("S2Builder.HighPrecisionPredicates") unittest {
  // To produce correct output in this example, the algorithm needs fall back
  // to high precision predicates when the output of the normal predicates is
  // uncertain.
  S2Point[] vertices = [
      S2Point(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
      S2Point(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),
      S2Point(-0.10531192039116592, -0.80522217309701472, 0.58354661457028933),
  ];
  auto input = new S2Polyline(vertices);
  S1Angle snap_radius = INTERSECTION_MERGE_RADIUS;
  auto options = new S2Builder.Options(new IdentitySnapFunction(snap_radius));
  options.setIdempotent(false);
  auto builder = new S2Builder(options);
  auto output = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output));
  builder.forceVertex(S2Point(-0.10531192039134191, -0.80522217309705857, 0.58354661457019719));
  builder.addPolyline(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
}

// Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
S2Point choosePoint() {
  S2Point x = S2Testing.randomPoint();
  for (int i = 0; i < 3; ++i) {
    if (S2Testing.rnd.oneIn(3)) {
      x[i] *= pow(1e-50, S2Testing.rnd.randDouble());
    }
  }
  return x.normalize();
}

@("S2Builder.HighPrecisionStressTest") unittest {
  // This test constructs many small, random inputs such that the output is
  // likely to be inconsistent unless high-precision predicates are used.

  S1Angle snap_radius = INTERSECTION_MERGE_RADIUS;
  // Some S2Builder calculations use an upper bound that takes into account
  // S1ChordAngle errors.  We sometimes try perturbing points by very close to
  // that distance in an attempt to expose errors.
  auto ca = S1ChordAngle(snap_radius);
  S1Angle snap_radius_with_error = ca.plusError(
      ca.getS1AngleConstructorMaxError() +
      getUpdateMinDistanceMaxError(ca)).toS1Angle();

  auto rnd = &S2Testing.rnd;
  int non_degenerate = 0;
  const int kIters = 8000 * iterationMultiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    // TODO(ericv): This test fails with a random seed of 96.  Change this
    // back to "iter + 1" once all the exact predicates are implemented.
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.

    // We construct a nearly degenerate triangle where one of the edges is
    // sometimes very short.  Then we add a forced vertex somewhere near the
    // shortest edge.  Then after snapping, we check that (1) the edges still
    // form a loop, and (2) if the loop is non-degenerate, then it has the
    // same orientation as the original triangle.
    //
    // v1 is located randomly.  (v0,v1) is the longest of the three edges.
    // v2 is located along (v0,v1) but is perturbed by up to 2 * snap_radius.
    S2Point v1 = choosePoint(), v0_dir = choosePoint();
    double d0 = pow(1e-16, rnd.randDouble());
    S2Point v0 = interpolateAtDistance(S1Angle.fromRadians(d0), v1, v0_dir);
    double d2 = 0.5 * d0 * pow(1e-16, pow(rnd.randDouble(), 2));
    S2Point v2 = interpolateAtDistance(S1Angle.fromRadians(d2), v1, v0_dir);
    v2 = S2Testing.samplePoint(new S2Cap(v2, 2 * snap_radius));
    // Vary the edge directions by randomly swapping v0 and v2.
    if (rnd.oneIn(2)) swap(v0, v2);

    // The forced vertex (v3) is either located near the (v1, v2) edge.
    // We perturb it either in a random direction from v1 or v2, or
    // perpendicular to (v1, v2) starting from an interior edge point.
    S1Angle d3 = rnd.oneIn(2) ? snap_radius : snap_radius_with_error;
    if (rnd.oneIn(3)) d3 = 1.5 * rnd.randDouble() * d3;
    S2Point v3;
    if (rnd.oneIn(5)) {
      v3 = rnd.oneIn(2) ? v1 : v2;
      v3 = interpolateAtDistance(d3, v3, choosePoint());
    } else {
      v3 = interpolate(pow(1e-16, rnd.randDouble()), v1, v2);
      v3 = interpolateAtDistance(d3, v3, v1.crossProd(v2).normalize());
    }
    auto options = new S2Builder.Options(new IdentitySnapFunction(snap_radius));
    options.setIdempotent(false);
    auto builder = new S2Builder(options);
    auto output = new S2Polygon();
    output.setS2debugOverride(S2Debug.DISABLE);
    builder.startLayer(new S2PolygonLayer(output));
    builder.forceVertex(v3);
    builder.addEdge(v0, v1);
    builder.addEdge(v1, v2);
    builder.addEdge(v2, v0);
    S2Error error;
    if (!builder.build(error)) {
      logger.logError("d0=", d0, ", d2=", d2, ", d3=", d3);
    }
    if (error.ok() && !output.isEmpty()) {
      Assert.equal(output.numLoops(), 1);
      if (output.numLoops() == 1) {
        Assert.equal(output.isValid(), true);
        Assert.equal(sign(v0, v1, v2) > 0, output.loop(0).isNormalized(),
            "d0=" ~ d0.to!string ~ ", d2=" ~ d2.to!string ~ ", d3=" ~ d3.to!string);
        ++non_degenerate;
      }
    }
  }
  logger.logInfo(non_degenerate, " non-degenerate out of ", kIters);
  Assert.notLessThan(non_degenerate, kIters / 10);
}

@("S2Builder.SelfIntersectionStressTest") unittest {
  const int kIters = 50 * iterationMultiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Testing.rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    // TODO: Implement timer.
    //CycleTimer timer;
    //timer.Start();

    // The minimum radius is about 36cm on the Earth's surface.  The
    // performance is reduced for radii much smaller than this because
    // S2ShapeIndex only indexes regions down to about 1cm across.
    auto cap = S2Testing.getRandomCap(1e-14, 1e-2);

    auto options = new S2Builder.Options();
    options.setSplitCrossingEdges(true);
    if (S2Testing.rnd.oneIn(2)) {
      S1Angle radius = cap.getRadius();
      int min_exp = IntLatLngSnapFunction.exponentForMaxSnapRadius(radius);
      int exponent = min(IntLatLngSnapFunction.MAX_EXPONENT,
          min_exp + S2Testing.rnd.uniform(5));
      options.setSnapFunction(new IntLatLngSnapFunction(exponent));
    }
    auto builder = new S2Builder(options);

    // Note that the number of intersections (and the running time) is
    // quadratic in the number of vertices.  With 200 input vertices, the
    // output consists of about 2300 loops and 9000 vertices.
    auto output = new S2Polygon();
    builder.startLayer(new S2PolygonLayer(output, S2PolygonLayer.Options(EdgeType.UNDIRECTED)));
    S2Point[50] vertices;
    foreach (S2Point vertex; vertices) {
      vertex = S2Testing.samplePoint(cap);
    }
    vertices.back() = vertices.front();
    auto input = new S2Polyline(vertices);
    builder.addPolyline(input);
    S2Error error;
    Assert.equal(builder.build(error), true, error.toString());
    Assert.equal(output.findValidationError(error), false, error.toString());
    if (iter == -1) {
      writeln("S2Polyline: ", toString(input));
      writeln("S2Polygon: ", toString(output));
    }
    if (iter < 50) {
      writefln("iter=%4d: ms=%4d, radius=%8.3g, loops=%d, vertices=%d\n",
          iter, 0 /*static_cast<int64_t>(timer.GetInMs())*/,
          cap.getRadius().radians(), output.numLoops(),
          output.numVertices());
    }
  }
}

@("S2Builder.FractalStressTest") unittest {
  // TODO: This test fails on iteration 6 in S2ClosestEdgeQueryBase.addInitialRange().
  //const int kIters = 100 * iterationMultiplier;
  const int kIters = 5 * iterationMultiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Testing.rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    auto fractal = new S2Testing.Fractal();
    fractal.setLevelForApproxMaxEdges(800);
    fractal.setLevelForApproxMinEdges(12);
    fractal.setFractalDimension(1.5 + 0.5 * S2Testing.rnd.randDouble());
    auto input = new S2Polygon(
        fractal.makeLoop(S2Testing.getRandomFrame(), S1Angle.fromDegrees(20.0)));
    auto options = new S2Builder.Options();
    if (S2Testing.rnd.oneIn(3)) {
      int exponent = S2Testing.rnd.uniform(11);
      options.setSnapFunction(new IntLatLngSnapFunction(exponent));
    } else if (S2Testing.rnd.oneIn(2)) {
      int level = S2Testing.rnd.uniform(20);
      options.setSnapFunction(new S2CellIdSnapFunction(level));
    } else {
      options.setSnapFunction(new IdentitySnapFunction(
          S1Angle.fromDegrees(10 * pow(1e-4, S2Testing.rnd.randDouble()))));
    }
    auto builder = new S2Builder(options);
    auto output = new S2Polygon();
    builder.startLayer(new S2PolygonLayer(output));
    builder.addPolygon(input);
    S2Error error;
    Assert.equal(builder.build(error), true, error.toString());
    Assert.equal(output.findValidationError(error), false, error.toString());
    if (iter == -1) {
      writeln("S2Polygon: ", toString(input));
      writeln("S2Polygon: ", toString(output));
    }
    if (iter < 50) {
      writefln("iter=%4d: in_vertices=%d, out_vertices=%d\n",
             iter, input.numVertices(), output.numVertices());
    }
  }
}

void checkSnappingWithForcedVertices(
    string input_str, S1Angle snap_radius, string vertices_str, string expected_str) {
  auto builder = new S2Builder(new S2Builder.Options(new IdentitySnapFunction(snap_radius)));
  S2Point[] vertices = parsePoints(vertices_str);
  foreach (vertex; vertices) {
    builder.forceVertex(vertex);
  }
  auto output = new S2Polyline();
  builder.startLayer(new S2PolylineLayer(output));
  builder.addPolyline(makePolylineOrDie(input_str));
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(expected_str, toString(output));
}

@("S2Builder.AdjacentCoverageIntervalsSpanMoreThan90Degrees") unittest {
  // The test for whether one Voronoi site excludes another along a given
  // input edge boils down to a test of whether two angle intervals "a" and
  // "b" overlap.  Let "ra" and "rb" be the semi-widths of the two intervals,
  // and let "d" be the angle between their centers.  Then "a" contains "b" if
  // (rb + d <= ra), and "b" contains "a" if (rb - d >= ra).  However the
  // actual code uses the sines of the angles, e.g.  sin(rb + d) <= sin(ra).
  // This works fine most of the time, but the first condition (rb + d <= ra)
  // also needs to check that rb + d < 90 degrees.  This test verifies that
  // case.

  // The following 3 tests have d < 90, d = 90, and d > 90 degrees, but in all
  // 3 cases rb + d > 90 degrees.
  checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.fromDegrees(60.0),
      "0:0, 0:70", "0:0, 0:70");
  checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.fromDegrees(60.0),
      "0:0, 0:90", "0:0, 0:90");
  checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.fromDegrees(60.0),
      "0:0, 0:110", "0:0, 0:110");

  // This test has d = 180 degrees, i.e. the two sites project to points that
  // are 180 degrees apart along the input edge.  The snapped edge doesn't
  // stay within max_edge_deviation() of the input edge, so an extra site is
  // added and it is snapped again (yielding two edges).  The case we are
  // testing here is the first call to SnapEdge() before adding the site.
  checkSnappingWithForcedVertices("0:10, 0:170", S1Angle.fromDegrees(50.0),
      "47:0, 49:180", "47:0, 0:90, 49:180");

  // This test has d = 220 degrees, i.e. when the input edge is snapped it
  // goes the "wrong way" around the sphere.  Again, the snapped edge is too
  // far from the input edge so an extra site is added and it is resnapped.
  checkSnappingWithForcedVertices("0:10, 0:170", S1Angle.fromDegrees(70.0),
      "0:-20, 0:-160", "0:-20, 0:90, 0:-160");

  // Without using forced vertices, the maximum angle between the coverage
  // interval centers is d = 300 degrees.  This would use an edge 180 degrees
  // long, and then place two sites 60 degrees past either endpoint.  With
  // forced vertices we can increase the snap radius to 70 degrees and get an
  // angle of up to d = 320 degrees, but the sites are only 40 degrees apart
  // (which is why it requires forced vertices).  The test below is an
  // approximation of this situation with d = 319.6 degrees.
  checkSnappingWithForcedVertices("0:0.1, 0:179.9", S1Angle.fromDegrees(70.0),
      "0:-69.8, 0:-110.2",
      "0:-69.8, 0:90, 0:-110.2");
}

@("S2Builder.OldS2PolygonBuilderBug") unittest {
  // This is a polygon that caused the obsolete S2PolygonBuilder class to
  // generate an invalid output polygon (duplicate edges).
  auto input = makePolygonOrDie(
      "32.2983095:72.3416582, 32.2986281:72.3423059, "
      ~ "32.2985238:72.3423743, 32.2987176:72.3427807, "
      ~ "32.2988174:72.3427056, 32.2991269:72.3433480, "
      ~ "32.2991881:72.3433077, 32.2990668:72.3430462, "
      ~ "32.2991745:72.3429778, 32.2995078:72.3436725, "
      ~ "32.2996075:72.3436269, 32.2985465:72.3413832, "
      ~ "32.2984558:72.3414530, 32.2988015:72.3421839, "
      ~ "32.2991552:72.3429416, 32.2990498:72.3430073, "
      ~ "32.2983764:72.3416059");
  Assert.equal(input.isValid(), true);

  S1Angle snap_radius = S2Testing.metersToAngle(20.0/0.866);
  auto builder = new S2Builder(new S2Builder.Options(new IdentitySnapFunction(snap_radius)));
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(output));
  builder.addPolygon(input);
  S2Error error;
  Assert.equal(builder.build(error), true, error.toString());
  Assert.equal(output.isValid(), true);
  auto expected = makePolygonOrDie(
      "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; "
      ~ "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, "
      ~ "32.2985238:72.3423743, 32.2987176:72.3427807");
  expectPolygonsEqual(expected, output);
}
