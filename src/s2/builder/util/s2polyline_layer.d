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

module s2.builder.util.s2polyline_layer;

import s2.id_set_lexicon;
import s2.mutable_s2shape_index;
import s2.s2point;
import s2.s2builder;
import s2.builder.graph;
import s2.builder.layer;
import s2.s2debug;
import s2.s2error;
import s2.s2polyline;

alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;
alias SiblingPairs = GraphOptions.SiblingPairs;

/**
 * A layer type that assembles edges (directed or undirected) into an
 * S2Polyline.  Returns an error if the edges cannot be assembled into a
 * single unbroken polyline.
 *
 * Duplicate edges are handled correctly (e.g., if a polyline backtracks on
 * itself, or loops around and retraces some of its previous edges.)  The
 * implementation attempts to preserve the order of directed input edges
 * whenever possible, so that if the input is a polyline and it is not
 * modified by S2Builder, then the output will be the same polyline (even if
 * the polyline backtracks on itself or forms a loop).  With undirected edges,
 * there are no such guarantees; for example, even if the input consists of a
 * single undirected edge, then either directed edge may be returned.
 *
 * S2PolylineLayer does not support options such as discarding sibling pairs
 * or merging duplicate edges because these options can split the polyline
 * into several pieces.  Use S2PolylineVectorLayer if you need these features.
 */
class S2PolylineLayer : Layer {
public:
  static struct Options {
  public:
    /**
     * Indicates whether the input edges provided to S2Builder are directed or
     * undirected.  Directed edges should be used whenever possible to avoid
     * ambiguity.
     *
     * DEFAULT: S2Builder::EdgeType::DIRECTED
     */
    S2Builder.EdgeType edgeType() const {
      return _edgeType;
    }

    void setEdgeType(S2Builder.EdgeType edge_type) {
      _edgeType = edge_type;
    }

    /**
     * If true, calls FindValidationError() on the output polyline.  If any
     * error is found, it will be returned by S2Builder::Build().
     *
     * Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
     * order to turn off the default error checking in debug builds.
     *
     * DEFAULT: false
     */
    bool validate() const {
      return _validate;
    }

    void setValidate(bool validate) {
      _validate = validate;
    }

  private:
    S2Builder.EdgeType _edgeType = S2Builder.EdgeType.DIRECTED;
    bool _validate = false;
  }

  /// Specifies that a polyline should be constructed using the given options.
  this(S2Polyline polyline, Options options = Options()) {
    initialize(polyline, null, null, options);
  }

  alias LabelSetIds = LabelSetId[];

  /**
   * Specifies that a polyline should be constructed using the given options,
   * and that any labels attached to the input edges should be returned in
   * "label_set_ids" and "label_set_lexicion".
   *
   * The labels associated with the edge "polyline.vertex({j, j+1})" can be
   * retrieved as follows:
   *
   *   for (int32 label : label_set_lexicon.id_set(label_set_ids[j])) {...}
   */
  this(S2Polyline polyline, LabelSetIds* label_set_ids,
      IdSetLexicon* label_set_lexicon,
      Options options = Options()) {
    initialize(polyline, label_set_ids, label_set_lexicon, options);
  }

  /// Layer interface:
  override
  GraphOptions graphOptions() const {
    // Remove edges that collapse to a single vertex, but keep duplicate and
    // sibling edges, since merging duplicates or discarding siblings can make
    // it impossible to assemble the edges into a single polyline.
    return new GraphOptions(_options.edgeType(), DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP, SiblingPairs.KEEP);
  }

  override
  void build(in Graph g, ref S2Error error) {
    if (g.numEdges() == 0) {
      _polyline.initialize(new S2Point[0]);
      return;
    }
    Graph.EdgePolyline[] edge_polylines = g.getPolylines(Graph.PolylineType.WALK);
    if (edge_polylines.length != 1) {
      error.initialize(S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_POLYLINE,
          "Input edges cannot be assembled into polyline");
      return;
    }
    const(Graph.EdgePolyline) edge_polyline = edge_polylines[0];
    S2Point[] vertices;  // Temporary storage for vertices.
    vertices.reserve(edge_polyline.length);
    vertices ~= g.vertex(g.edge(edge_polyline[0])[0]);
    foreach (Graph.EdgeId e; edge_polyline) {
      vertices ~= g.vertex(g.edge(e)[1]);
    }
    if (_labelSetIds) {
      auto fetcher = new Graph.LabelFetcher(g, _options.edgeType());
      Label[] labels;  // Temporary storage for labels.
      (*_labelSetIds).reserve(edge_polyline.length);
      foreach (Graph.EdgeId e; edge_polyline) {
        fetcher.fetch(e, labels);
        (*_labelSetIds) ~= _labelSetLexicon.add(labels);
      }
    }
    _polyline.initialize(vertices);
    if (_options.validate()) {
      _polyline.findValidationError(error);
    }
  }

private:
  void initialize(S2Polyline polyline, LabelSetIds* label_set_ids,
      IdSetLexicon* label_set_lexicon, Options options)
  in {
    assert((label_set_ids is null) == (label_set_lexicon is null));
  } do {
    _polyline = polyline;
    _labelSetIds = label_set_ids;
    _labelSetLexicon = label_set_lexicon;
    _options = options;

    if (_options.validate()) {
      _polyline.setS2debugOverride(S2Debug.DISABLE);
    }
  }

  S2Polyline _polyline;
  LabelSetIds* _labelSetIds;
  IdSetLexicon* _labelSetLexicon;
  Options _options;
}

/// Like S2PolylineLayer, but adds the polyline to a MutableS2ShapeIndex (if the
/// polyline is non-empty).
class IndexedS2PolylineLayer : Layer {
public:
  alias Options = S2PolylineLayer.Options;

  this(MutableS2ShapeIndex index, Options options = Options()) {
    _index = index;
    _polyline = new S2Polyline();
    _layer = new S2PolylineLayer(_polyline, options);
  }

  override
  GraphOptions graphOptions() const {
    return _layer.graphOptions();
  }

  override
  void build(in Graph g, ref S2Error error) {
    _layer.build(g, error);
    if (error.ok() && _polyline.numVertices() > 0) {
      _index.add(new S2Polyline.Shape(_polyline));
    }
  }

private:
  MutableS2ShapeIndex _index;
  S2Polyline _polyline;
  S2PolylineLayer _layer;
}
