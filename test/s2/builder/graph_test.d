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
//
// Most of S2Builder::Graph is tested by the S2Builder::Layer implementations
// rather than here.

module s2.builder.graph_test;

import s2.builder.graph;

import s2.builder.util.testing;
import s2.id_set_lexicon;
import s2.s2builder;
import s2.s2error;
import s2.s2lax_polyline_shape;
import s2.s2text_format;

import fluent.asserts;

import std.conv;

alias EdgeType = S2Builder.EdgeType;

alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;
alias SiblingPairs = GraphOptions.SiblingPairs;

alias DegenerateBoundaries = Graph.DegenerateBoundaries;
alias Edge = Graph.Edge;
alias EdgeId = Graph.EdgeId;
alias InputEdgeId = Graph.InputEdgeId;
alias InputEdgeIdSetId = Graph.InputEdgeIdSetId;
alias LoopType = Graph.LoopType;
alias PolylineType = Graph.PolylineType;
alias VertexId = Graph.VertexId;

@("GetDirectedLoops.DegenerateEdges") unittest {
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
  builder.addShape(makeLaxPolylineOrDie("0:3, 3:3, 0:3"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  const Graph g = gc.graph();
  EdgeId[][] loops;
  Assert.equal(g.getDirectedLoops(LoopType.SIMPLE, loops, error), true);
  Assert.equal(loops.length, 3);
  Assert.equal(loops[0].length, 1);
  Assert.equal(loops[1].length, 4);
  Assert.equal(loops[2].length, 2);
}

@("GetDirectedComponents.DegenerateEdges") unittest {
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.MERGE, SiblingPairs.CREATE);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  const(Graph) g = gc.graph();
  Graph.DirectedComponent[] components;
  Assert.equal(g.getDirectedComponents(DegenerateBoundaries.KEEP, components, error), true);
  Assert.equal(components.length, 2);
  Assert.equal(components[0].length, 1);
  Assert.equal(components[0][0].length, 1);
  Assert.equal(components[1].length, 2);
  Assert.equal(components[1][0].length, 4);
  Assert.equal(components[1][1].length, 4);
}

@("GetUndirectedComponents.DegenerateEdges") unittest {
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  const(Graph) g = gc.graph();
  Graph.UndirectedComponent[] components;
  Assert.equal(g.getUndirectedComponents(LoopType.CIRCUIT, components, error), true);
  // The result consists of two components, each with two complements.  Each
  // complement in this example has exactly one loop.  The loops in both
  // complements of the first component have 1 vertex, while the loops in both
  // complements of the second component have 4 vertices.
  Assert.equal(components.length, 2);
  Assert.equal(components[0][0].length, 1);
  Assert.equal(components[0][0][0].length, 1);
  Assert.equal(components[0][1].length, 1);
  Assert.equal(components[0][1][0].length, 1);
  Assert.equal(components[1][0].length, 1);
  Assert.equal(components[1][0][0].length, 4);
  Assert.equal(components[1][1].length, 1);
  Assert.equal(components[1][1][0].length, 4);
}

@("GetPolylines.UndirectedDegeneratePaths") unittest {
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  builder.addShape(makeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  const(Graph) g = gc.graph();
  auto polylines = g.getPolylines(PolylineType.PATH);
  Assert.equal(polylines.length, 7);
}

@("GetPolylines.UndirectedDegenerateWalks") unittest {
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  builder.addShape(makeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
  builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
  S2Error error;
  Assert.equal(builder.build(error), true);
  const(Graph) g = gc.graph();
  auto polylines = g.getPolylines(PolylineType.WALK);
  Assert.equal(polylines.length, 2);
  Assert.equal(polylines[0].length, 2);
  Assert.equal(polylines[1].length, 5);
}

struct TestEdge {
  VertexId first;
  VertexId second;
  InputEdgeId[] inputIds;
}

void checkProcessEdges(in TestEdge[] input,
                      in TestEdge[] expected,
                      GraphOptions options,
                      S2Error.Code expected_code = S2Error.Code.OK) {
  Edge[] edges;
  InputEdgeIdSetId[] input_id_set_ids;
  auto id_set_lexicon = new IdSetLexicon();
  foreach (const e; input) {
    edges ~= Edge(e.first, e.second);
    input_id_set_ids ~= id_set_lexicon.add(e.inputIds);
  }
  S2Error error;
  Graph.processEdges(options, edges, input_id_set_ids, id_set_lexicon, error);
  Assert.equal(error.code(), expected_code);
  Assert.equal(edges.length, input_id_set_ids.length);
  for (int i = 0; i < expected.length; ++i) {
    const auto e = expected[i];
    Assert.lessThan(i, cast(int) edges.length, "Not enough output edges");
    Assert.equal(edges[i], Edge(e.first, e.second), "(edge " ~ i.to!string ~ ")");
    auto actual_ids = id_set_lexicon.idSet(input_id_set_ids[i]);
    Assert.equal(e.inputIds, actual_ids, "(edge " ~ i.to!string ~ ")");
  }
  Assert.equal(expected.length, edges.length, "Too many output edges");
}

@("ProcessEdges.DiscardDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 0), TestEdge(0, 0)], [], options);
}

@("ProcessEdges.KeepDuplicateDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 0), TestEdge(0, 0)], [TestEdge(0, 0), TestEdge(0, 0)], options);
}

@("ProcessEdges.MergeDuplicateDegenerateEdges") unittest {
  auto options = new GraphOptions(
      EdgeType.DIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.MERGE, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0, [2])], [TestEdge(0, 0, [1, 2])], options);
}

@("ProcessEdges.MergeUndirectedDuplicateDegenerateEdges") unittest {
  // Edge count should be reduced to 2 (i.e., one undirected edge), and all
  // labels should be merged.
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.MERGE, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0), TestEdge(0, 0), TestEdge(0, 0, [2])],
      [TestEdge(0, 0, [1, 2]), TestEdge(0, 0, [1, 2])], options);
}

@("ProcessEdges.ConvertedUndirectedDegenerateEdges") unittest {
  // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
  // merges any edge labels.
  auto options = new GraphOptions(
      EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0), TestEdge(0, 0), TestEdge(0, 0, [2])],
      [TestEdge(0, 0, [1, 2]), TestEdge(0, 0, [1, 2])], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);
}

@("ProcessEdges.MergeConvertedUndirectedDuplicateDegenerateEdges") unittest {
  // Like the test above, except that we also merge duplicates.
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.MERGE, SiblingPairs.REQUIRE);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0), TestEdge(0, 0), TestEdge(0, 0, [2])],
      [TestEdge(0, 0, [1, 2])], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);
}

@("ProcessEdges.DiscardExcessConnectedDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  // Test that degenerate edges are discarded if they are connnected to any
  // non-degenerate edges (whether they are incoming or outgoing, and whether
  // they are lexicographically before or after the degenerate edge).
  checkProcessEdges([TestEdge(0, 0), TestEdge(0, 1)], [TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(0, 0), TestEdge(1, 0)], [TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 1)], [TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(1, 0), TestEdge(1, 1)], [TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardExcessIsolatedDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  // Test that DISCARD_EXCESS does not merge any remaining duplicate
  // degenerate edges together.
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0, [2])],
      [TestEdge(0, 0, [1]), TestEdge(0, 0, [2])], options);
}

@("ProcessEdges.DiscardExcessUndirectedIsolatedDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  // Test that DISCARD_EXCESS does not merge any remaining duplicate
  // undirected degenerate edges together.
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0), TestEdge(0, 0, [2]), TestEdge(0, 0)],
      [TestEdge(0, 0, [1]), TestEdge(0, 0), TestEdge(0, 0, [2]), TestEdge(0, 0)], options);
}

@("ProcessEdges.DiscardExcessConvertedUndirectedIsolatedDegenerateEdges") unittest {
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD_EXCESS,
      DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
  // Converting from UNDIRECTED to DIRECTED cuts the edge count in half and
  // merges edge labels.
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0, [2]), TestEdge(0, 0, [3]), TestEdge(0, 0)],
      [TestEdge(0, 0, [1, 2, 3]), TestEdge(0, 0, [1, 2, 3])], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);
}

@("ProcessEdges.SiblingPairsDiscardMergesDegenerateEdgeLabels") unittest {
  // Test that when SiblingPairs::DISCARD or SiblingPairs::DISCARD_EXCESS
  // is specified, the edge labels of degenerate edges are merged together
  // (for consistency, since these options merge the labels of all
  // non-degenerate edges as well).
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.KEEP,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0, [2]), TestEdge(0, 0, [3])],
      [TestEdge(0, 0, [1, 2, 3]), TestEdge(0, 0, [1, 2, 3]), TestEdge(0, 0, [1, 2, 3])],
      options);
  options.setSiblingPairs(SiblingPairs.DISCARD_EXCESS);
  checkProcessEdges([TestEdge(0, 0, [1]), TestEdge(0, 0, [2]), TestEdge(0, 0, [3])],
      [TestEdge(0, 0, [1, 2, 3]), TestEdge(0, 0, [1, 2, 3]), TestEdge(0, 0, [1, 2, 3])],
      options);
}

@("ProcessEdges.KeepSiblingPairs") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [TestEdge(0, 1), TestEdge(1, 0)], options);
}

@("ProcessEdges.MergeDuplicateSiblingPairs") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.MERGE, SiblingPairs.KEEP);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardSiblingPairs") unittest {
  // Check that matched pairs are discarded, leaving behind any excess edges.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], [], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(1, 0), TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardSiblingPairsMergeDuplicates") unittest {
  // Check that matched pairs are discarded, and then any remaining edges
  // are merged.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.MERGE, SiblingPairs.DISCARD);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], [], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardUndirectedSiblingPairs") unittest {
  // An undirected sibling pair consists of four edges, two in each direction
  // (see s2builder.h).  Since undirected edges always come in pairs, this
  // means that the result always consists of either 0 or 2 edges.
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], [], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0),
          TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardExcessSiblingPairs") unittest {
  // Like SiblingPairs::DISCARD, except that one sibling pair is kept if the
  // result would otherwise be empty.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(1, 0), TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardExcessSiblingPairsMergeDuplicates") unittest {
  // Like SiblingPairs::DISCARD, except that one sibling pair is kept if the
  // result would otherwise be empty.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.MERGE, SiblingPairs.DISCARD_EXCESS);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(1, 0)], options);
}

@("ProcessEdges.DiscardExcessUndirectedSiblingPairs") unittest {
  // Like SiblingPairs::DISCARD, except that one undirected sibling pair
  // (4 edges) is kept if the result would otherwise be empty.
  auto options = new GraphOptions(EdgeType.UNDIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.DISCARD_EXCESS);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0),
          TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
}

@("ProcessEdges.CreateSiblingPairs") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.CREATE);
  checkProcessEdges([TestEdge(0, 1)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1)],
      [TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], options);
}

@("ProcessEdges.RequireSiblingPairs") unittest {
  // Like SiblingPairs::CREATE, but generates an error.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.REQUIRE);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1)], [TestEdge(0, 1), TestEdge(1, 0)], options,
      S2Error.Code.BUILDER_MISSING_EXPECTED_SIBLING_EDGES);
}

@("ProcessEdges.CreateUndirectedSiblingPairs") unittest {
  // An undirected sibling pair consists of 4 edges, but SiblingPairs::CREATE
  // also converts the graph to EdgeType::DIRECTED and cuts the number of
  // edges in half.
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.KEEP, SiblingPairs.CREATE);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);

  options.setEdgeType(EdgeType.UNDIRECTED);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);

  options.setEdgeType(EdgeType.UNDIRECTED);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0),
          TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0), TestEdge(1, 0)], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);
}

@("ProcessEdges.CreateSiblingPairsMergeDuplicates") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.MERGE, SiblingPairs.CREATE);
  checkProcessEdges([TestEdge(0, 1)], [TestEdge(0, 1), TestEdge(1, 0)], options);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1)], [TestEdge(0, 1), TestEdge(1, 0)], options);
}

@("ProcessEdges.CreateUndirectedSiblingPairsMergeDuplicates") unittest {
  auto options = new GraphOptions(EdgeType.DIRECTED, DegenerateEdges.DISCARD,
      DuplicateEdges.MERGE, SiblingPairs.CREATE);
  checkProcessEdges([TestEdge(0, 1), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);

  options.setEdgeType(EdgeType.UNDIRECTED);
  checkProcessEdges([TestEdge(0, 1), TestEdge(0, 1), TestEdge(0, 1), TestEdge(1, 0),
          TestEdge(1, 0), TestEdge(1, 0)],
      [TestEdge(0, 1), TestEdge(1, 0)], options);
  Assert.equal(options.edgeType(), EdgeType.DIRECTED);
}
