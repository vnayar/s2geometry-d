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

module s2.builder.util.testing_test;

import s2.builder.graph;
import s2.builder.util.testing;
import s2.s2builder;
import s2.s2error;
import s2.s2point;

import fluent.asserts;
import std.range;

alias EdgeType = S2Builder.EdgeType;

alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;
alias SiblingPairs = GraphOptions.SiblingPairs;

@("GraphCloningLayer.MakeIndependentCopy") unittest {
  // Also tests GraphClone.
  auto gc = new GraphClone();
  auto builder = new S2Builder(new S2Builder.Options());
  auto graph_options = new GraphOptions(
      EdgeType.DIRECTED, DegenerateEdges.DISCARD, DuplicateEdges.MERGE, SiblingPairs.KEEP);
  builder.startLayer(new GraphCloningLayer(graph_options, gc));
  auto v0 = S2Point(1, 0, 0);
  auto v1 = S2Point(0, 1, 0);
  builder.setLabel(14);
  builder.addEdge(v0, v1);
  S2Error error;
  Assert.equal(builder.build(error), true);
  const Graph g = gc.graph();
  Assert.equal(graph_options == g.options(), true);
  Assert.equal([v0, v1], g.vertices());
  Assert.equal(g.edges(), [Graph.Edge(0, 1)]);
  Assert.equal(g.inputEdgeIds(0).length, 1);
  Assert.equal(g.inputEdgeIds(0).front(), 0);
  auto fetcher = new Graph.LabelFetcher(g, g.options().edgeType());
  S2Builder.Label[] labels;
  fetcher.fetch(0, labels);
  Assert.equal(labels, [14]);
  // S2Builder sets a default IsFullPolygonPredicate that returns an error.
  Assert.equal(g.isFullPolygonPredicate() !is null, true);
}

@("GraphAppendingLayer.AppendsTwoGraphs") unittest {
  Graph[] graphs;
  GraphClone[] clones;
  auto builder = new S2Builder(new S2Builder.Options());
  builder.startLayer(new GraphAppendingLayer(new GraphOptions(), &graphs, &clones));
  auto v0 = S2Point(1, 0, 0);
  auto v1 = S2Point(0, 1, 0);
  builder.addEdge(v0, v1);
  builder.startLayer(new GraphAppendingLayer(new GraphOptions(), &graphs, &clones));
  builder.addEdge(v1, v0);
  S2Error error;
  Assert.equal(builder.build(error), true);
  Assert.equal(graphs.length, 2);
  Assert.equal(clones.length, 2);
  Assert.equal(graphs[0].vertex(graphs[0].edge(0)[0]), v0);
  Assert.equal(graphs[1].vertex(graphs[1].edge(0)[0]), v1);
}
