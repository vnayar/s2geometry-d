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

module s2.builder.util.s2polygon_layer_test;

import s2.builder.util.s2polygon_layer;

import s2.logger;
import s2.mutable_s2shape_index;
import s2.s2builder;
import s2.s2debug;
import s2.s2error;
import s2.s2loop;
import s2.s2point;
import s2.id_set_lexicon;
import s2.s2polygon;
import s2.s2polyline;
import s2.s2text_format;

import std.algorithm;

import fluent.asserts;

alias EdgeType = S2Builder.EdgeType;

void checkS2Polygon(string[] input_strs, string expected_str, EdgeType edge_type) {
  logger.logTrace(edge_type == EdgeType.DIRECTED ? "DIRECTED" : "UNDIRECTED");
  auto builder = new S2Builder(new S2Builder.Options());
  auto output = new S2Polygon();
  builder.startLayer(new S2PolygonLayer(
      output, S2PolygonLayer.Options(edge_type)));
  foreach (input_str; input_strs) {
    builder.addPolygon(makeVerbatimPolygon(input_str));
  }
  S2Error error;
  Assert.equal(builder.build(error), true);
  // The input strings in tests may not be in normalized form, so we build an
  // S2Polygon and convert it back to a string.
  auto expected = makePolygonOrDie(expected_str);
  Assert.equal(toString(output), toString(expected));
}

// Unlike the methods above, the input consists of a set of *polylines*.
void checkS2PolygonError(string[]  input_strs, S2Error.Code expected_error, EdgeType edge_type) {
  logger.logTrace(edge_type == EdgeType.DIRECTED ? "DIRECTED" : "UNDIRECTED");
  auto builder = new S2Builder(new S2Builder.Options());
  auto output = new S2Polygon();
  auto options = S2PolygonLayer.Options(edge_type);
  options.setValidate(true);
  builder.startLayer(new S2PolygonLayer(output, options));
  foreach (input_str; input_strs) {
    builder.addPolyline(makePolyline(input_str));
  }
  S2Error error;
  Assert.equal(builder.build(error), false);
  Assert.equal(error.code(), expected_error);
}

void checkS2PolygonError(string[] input_strs, S2Error.Code expected_error) {
  checkS2PolygonError(input_strs, expected_error, EdgeType.DIRECTED);
  checkS2PolygonError(input_strs, expected_error, EdgeType.UNDIRECTED);
}

void checkS2Polygon(string[] input_strs, string expected_str) {
  checkS2Polygon(input_strs, expected_str, EdgeType.DIRECTED);
  checkS2Polygon(input_strs, expected_str, EdgeType.UNDIRECTED);
}

void checkS2PolygonUnchanged(string input_str) {
  checkS2Polygon([input_str], input_str);
}

@("S2PolygonLayer.NoLoops") unittest {
  checkS2PolygonUnchanged("");
}

@("S2PolygonLayer.SmallLoop") unittest {
  checkS2PolygonUnchanged("0:0, 0:1, 1:1");
}

@("S2PolygonLayer.ThreeLoops") unittest {
  // The second two loops are nested.
  checkS2PolygonUnchanged("0:1, 1:1, 0:0; "
      ~ "3:3, 3:6, 6:6, 6:3; "
      ~ "4:4, 4:5, 5:5, 5:4");
}

@("S2PolygonLayer.PartialLoop") unittest {
  checkS2PolygonError(["0:1, 2:3, 4:5"],
      S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS);
}

@("S2PolygonLayer.InvalidPolygon") unittest {
  checkS2PolygonError(["0:0, 0:10, 10:0, 10:10, 0:0"],
      S2Error.Code.LOOP_SELF_INTERSECTION);
}

@("S2PolygonLayer.DuplicateInputEdges") unittest {
  // Check that S2PolygonLayer can assemble polygons even when there are
  // duplicate edges (after sibling pairs are removed).
  auto builder = new S2Builder(new S2Builder.Options());
  auto output = new S2Polygon();
  S2PolygonLayer.Options options;
  options.setValidate(true);
  builder.startLayer(new S2PolygonLayer(output, options));
  builder.addPolyline(makePolylineOrDie(
      "0:0, 0:2, 2:2, 1:1, 0:2, 2:2, 2:0, 0:0"));
  S2Error error;
  Assert.equal(builder.build(error), false);
  Assert.equal(error.code(), S2Error.Code.POLYGON_LOOPS_SHARE_EDGE);
  Assert.equal(output.numLoops(), 2);
  auto loop0 = new S2Loop(makeLoopOrDie("0:0, 0:2, 2:2, 2:0"));
  auto loop1 = new S2Loop(makeLoopOrDie("0:2, 2:2, 1:1"));
  Assert.equal(loop0.equals(output.loop(0)), true);
  Assert.equal(loop1.equals(output.loop(1)), true);
}

// Since we don't expect to have any crossing edges, the key for each edge is
// simply the sum of its endpoints.  This key has the advantage of being
// unchanged when the endpoints of an edge are swapped.
alias EdgeLabelMap =  int[][S2Point];

void addPolylineWithLabels(S2Polyline polyline, EdgeType edge_type,
                           int label_begin, S2Builder builder,
                           ref EdgeLabelMap edge_label_map) {
  for (int i = 0; i + 1 < polyline.numVertices(); ++i) {
    int label = label_begin + i;
    builder.setLabel(label);
    // With undirected edges, reverse the direction of every other input edge.
    int dir = edge_type == EdgeType.DIRECTED ? 1 : (i & 1);
    builder.addEdge(polyline.vertex(i + (1 - dir)), polyline.vertex(i + dir));
    S2Point key = polyline.vertex(i) + polyline.vertex(i + 1);
    edge_label_map[key] ~= label;
  }
}

private void checkEdgeLabels(EdgeType edge_type) {
  auto builder = new S2Builder(new S2Builder.Options());
  auto output = new S2Polygon();
  S2PolygonLayer.LabelSetIds label_set_ids = new S2PolygonLayer.LabelSetIds(0);
  IdSetLexicon label_set_lexicon = new IdSetLexicon();
  builder.startLayer(new S2PolygonLayer(
      output, label_set_ids, label_set_lexicon,
      S2PolygonLayer.Options(edge_type)));

  // We use a polygon consisting of 3 loops.  The loops are reordered and
  // some of the loops are inverted during S2Polygon construction.
  EdgeLabelMap edge_label_map;
  addPolylineWithLabels(makePolylineOrDie("0:0, 9:1, 1:9, 0:0, 2:8, 8:2, "
          ~ "0:0, 0:10, 10:10, 10:0, 0:0"),
      edge_type, 0, builder, edge_label_map);
  S2Error error;
  Assert.equal(builder.build(error), true);
  int[] expected_loop_sizes = [ 4, 3, 3 ];
  Assert.equal(label_set_ids.length, expected_loop_sizes.length);
  for (int i = 0; i < expected_loop_sizes.length; ++i) {
    Assert.equal(label_set_ids[i].length, expected_loop_sizes[i]);
    for (int j = 0; j < label_set_ids[i].length; ++j) {
      S2Point key = output.loop(i).vertex(j) + output.loop(i).vertex(j + 1);
      int[] expected_labels = edge_label_map[key];
      Assert.equal(label_set_lexicon.idSet(label_set_ids[i][j]).length, expected_labels.length);
      Assert.equal(isPermutation(
              expected_labels, label_set_lexicon.idSet(label_set_ids[i][j])), true);
    }
  }
}

@("S2PolygonLayer.DirectedEdgeLabels") unittest {
  checkEdgeLabels(EdgeType.DIRECTED);
}

@("S2PolygonLayer.UndirectedEdgeLabels") unittest {
  checkEdgeLabels(EdgeType.UNDIRECTED);
}

@("S2PolygonLayer.ThreeLoopsIntoOne") unittest {
  // Three loops (two shells and one hole) that combine into one.
  checkS2Polygon([
          "10:0, 0:0, 0:10, 5:10, 10:10, 10:5",
          "0:10, 0:15, 5:15, 5:10",
          "10:10, 5:10, 5:5, 10:5"],
      "10:5, 10:0, 0:0, 0:10, 0:15, 5:15, 5:10, 5:5");
}

@("S2PolygonLayer.TrianglePyramid") unittest {
  // A big CCW triangle containing 3 CW triangular holes.  The whole thing
  // looks like a pyramid of nine triangles.  The output consists of 6
  // positive triangles with no holes.
  checkS2Polygon([
          "0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1",
          "0:2, 1:1, 1:3",
          "0:4, 1:3, 1:5",
          "1:3, 2:2, 2:4"],
      "0:4, 0:6, 1:5; 2:4, 3:3, 2:2; 2:2, 1:1, 1:3; "
      ~ "1:1, 0:0, 0:2; 1:3, 0:2, 0:4; 1:3, 1:5, 2:4");
}

@("S2PolygonLayer.ComplexNesting") unittest {
  // A complex set of nested polygons, with the loops in random order and the
  // vertices in random cyclic order within each loop.  This test checks that
  // the order (after S2Polygon::InitNested is called) is preserved exactly,
  // whether directed or undirected edges are used.
  checkS2PolygonUnchanged(
      "47:15, 47:5, 5:5, 5:15; "
      ~ "35:12, 35:7, 27:7, 27:12; "
      ~ "1:50, 50:50, 50:1, 1:1; "
      ~ "42:22, 10:22, 10:25, 42:25; "
      ~ "47:30, 47:17, 5:17, 5:30; "
      ~ "7:27, 45:27, 45:20, 7:20; "
      ~ "37:7, 37:12, 45:12, 45:7; "
      ~ "47:47, 47:32, 5:32, 5:47; "
      ~ "50:60, 50:55, 1:55, 1:60; "
      ~ "25:7, 17:7, 17:12, 25:12; "
      ~ "7:7, 7:12, 15:12, 15:7");
}

@("S2PolygonLayer.FiveLoopsTouchingAtOneCommonPoint") unittest {
  // Five nested loops that touch at one common point.
  checkS2PolygonUnchanged("0:0, 0:10, 10:10, 10:0; "
      ~ "0:0, 1:9, 9:9, 9:1; "
      ~ "0:0, 2:8, 8:8, 8:2; "
      ~ "0:0, 3:7, 7:7, 7:3; "
      ~ "0:0, 4:6, 6:6, 6:4");
}

@("S2PolygonLayer.FourNestedDiamondsTouchingAtTwoPointsPerPair") unittest {
  // Four diamonds nested inside each other, where each diamond shares two
  // vertices with the diamond inside it and shares its other two vertices
  // with the diamond that contains it.  The resulting shape looks vaguely
  // like an eye made out of chevrons.
  checkS2Polygon([
          "0:10, -10:0, 0:-10, 10:0",
          "0:-20, -10:0, 0:20, 10:0",
          "0:-10, -5:0, 0:10, 5:0",
          "0:5, -5:0, 0:-5, 5:0"],
      "10:0, 0:10, -10:0, 0:20; "
      ~ "0:-20, -10:0, 0:-10, 10:0; "
      ~ "5:0, 0:-10, -5:0, 0:-5; "
      ~ "0:5, -5:0, 0:10, 5:0");
}

@("S2PolygonLayer.SevenDiamondsTouchingAtOnePointPerPair") unittest {
  // Seven diamonds nested within each other touching at one
  // point between each nested pair.
  checkS2PolygonUnchanged("0:-70, -70:0, 0:70, 70:0; "
      ~ "0:-70, -60:0, 0:60, 60:0; "
      ~ "0:-50, -60:0, 0:50, 50:0; "
      ~ "0:-40, -40:0, 0:50, 40:0; "
      ~ "0:-30, -30:0, 0:30, 40:0; "
      ~ "0:-20, -20:0, 0:30, 20:0; "
      ~ "0:-10, -20:0, 0:10, 10:0");
}

@("IndexedS2PolygonLayer.AddsShape") unittest {
  auto builder = new S2Builder(new S2Builder.Options());
  auto index = new MutableS2ShapeIndex();
  builder.startLayer(new IndexedS2PolygonLayer(index));
  string polygon_str = "0:0, 0:10, 10:0";
  builder.addPolygon(makePolygonOrDie(polygon_str));
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal(index.numShapeIds(), 1);
  const S2Polygon polygon = (cast(S2Polygon.Shape) index.shape(0)).polygon();
  Assert.equal(polygon_str, toString(polygon));
}

@("IndexedS2PolygonLayer.AddsEmptyShape") unittest {
  auto builder = new S2Builder(new S2Builder.Options());
  auto index = new MutableS2ShapeIndex();
  builder.startLayer(new IndexedS2PolygonLayer(index));
  auto polygon = new S2Polygon();
  builder.addPolygon(polygon);
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal(index.numShapeIds(), 0);
}
