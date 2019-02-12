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

module s2.builder.util.testing;

import s2.builder.graph;
import s2.builder.layer;
import s2.id_set_lexicon;
import s2.s2builder;
import s2.s2error;
import s2.s2point;

import std.range;
import std.typecons;

// A class that copies an S2Builder::Graph and owns the underlying data
// (unlike S2Builder::Graph, which is just a view).
class GraphClone {
public:
  this() {}  // Must call Init().

  this(Graph g) { initialize(g); }

  void initialize(in Graph g) {
    _options = g.options();
    _vertices = g.vertices();
    _edges = g.edges();
    _inputEdgeIdSetIds = g.inputEdgeIdSetIds();
    _inputEdgeIdSetLexicon = g.inputEdgeIdSetLexicon();
    _labelSetIds = g.labelSetIds();
    _labelSetLexicon = g.labelSetLexicon();
    _isFullPolygonPredicate = g.isFullPolygonPredicate();
    _g = new Graph(
        _options, _vertices, _edges, _inputEdgeIdSetIds,
        _inputEdgeIdSetLexicon, _labelSetIds, _labelSetLexicon,
        _isFullPolygonPredicate);
  }

  Graph graph() {
    return _g;
  }

private:
  Rebindable!(const(GraphOptions)) _options;
  Rebindable!(const(S2Point[])) _vertices;
  Rebindable!(const(Graph.Edge[])) _edges;
  Rebindable!(const(Graph.InputEdgeIdSetId[])) _inputEdgeIdSetIds;
  Rebindable!(const(IdSetLexicon)) _inputEdgeIdSetLexicon;
  Rebindable!(const(Graph.LabelSetId[])) _labelSetIds;
  Rebindable!(const(IdSetLexicon)) _labelSetLexicon;
  S2Builder.IsFullPolygonPredicate _isFullPolygonPredicate;
  Graph _g;
}

/// A layer type that copies an S2Builder::Graph into a GraphClone object
/// (which owns the underlying data, unlike S2Builder::Graph itself).
class GraphCloningLayer : Layer {
public:
  this(GraphOptions graph_options, GraphClone gc) {
    _graphOptions = graph_options;
    _gc = gc;
  }

  override
  GraphOptions graphOptions() {
    return _graphOptions;
  }

  override
  void build(in Graph g, ref S2Error error) {
    _gc.initialize(g);
  }

private:
  GraphOptions _graphOptions;
  GraphClone _gc;
}

/// A layer type that copies an S2Builder::Graph and appends it to a vector,
/// and appends the corresponding GraphClone object (which owns the Graph data)
/// to a separate vector.
class GraphAppendingLayer : Layer {
public:
  this(
      GraphOptions graph_options,
      Graph[] graphs,
      GraphClone[] clones) {
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
    _clones ~= new GraphClone(g);
    _graphs ~= _clones.back().graph();
  }

private:
  GraphOptions _graphOptions;
  Graph[] _graphs;
  GraphClone[] _clones;
}
