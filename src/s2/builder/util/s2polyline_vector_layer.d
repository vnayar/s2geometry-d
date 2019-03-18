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

module s2.builder.util.s2polyline_vector_layer;

import s2.id_set_lexicon;
import s2.mutable_s2shape_index;
import s2.s2point;
import s2.s2builder;
import s2.builder.graph;
import s2.builder.layer;
import s2.s2debug;
import s2.s2error;
import s2.s2polyline;

import std.range;

alias DegenerateEdges = GraphOptions.DegenerateEdges;

/**
 * A layer type that assembles edges (directed or undirected) into multiple
 * S2Polylines.  Returns an error if S2Builder found any problem with the
 * input edges; this layer type does not generate any errors of its own.
 *
 * Duplicate edges are handled correctly (e.g., if a polyline backtracks on
 * itself, or loops around and retraces some of its previous edges.)  The
 * implementation attempts to preserve the order of the input edges whenever
 * possible, so that if the input is a polyline and it is not modified by
 * S2Builder, then the output will be the same polyline even if the polyline
 * forms a loop.  However, note that this is not guaranteed when undirected
 * edges are used: for example, if the input consists of a single undirected
 * edge, then either directed edge may be returned.
 */
class S2PolylineVectorLayer : Layer {
public:
  static struct Options {
  public:
    /**
     * Indicates whether the input edges provided to S2Builder are directed or
     * undirected.
     *
     * Directed edges should be used whenever possible to avoid ambiguity.
     * The implementation attempts to preserve the structure of directed input
     * edges whenever possible, so that if the input is a vector of disjoint
     * polylines and none of them need to be modified then the output will be
     * the same polylines in the same order.  With undirected edges, there are
     * no such guarantees.
     *
     * DEFAULT: S2Builder::EdgeType::DIRECTED
     */
    S2Builder.EdgeType edgeType() const {
      return _edgeType;
    }

    void setEdgeType(S2Builder.EdgeType edge_type) {
      _edgeType = edge_type;
    }


    alias PolylineType = Graph.PolylineType;

    /**
     * Indicates whether polylines should be "paths" (which don't allow
     * duplicate vertices, except possibly the first and last vertex) or
     * "walks" (which allow duplicate vertices and edges).
     *
     * If your input consists of polylines, and you want to split them into
     * separate pieces whenever they self-intersect or cross each other, then
     * use PolylineType::PATH (and probably use split_crossing_edges()).  If
     * you don't mind if your polylines backtrack or contain loops, then use
     * PolylineType::WALK.
     *
     * DEFAULT: PolylineType::PATH
     */
    PolylineType polylineType() const {
      return _polylineType;
    }

    void setPolylineType(PolylineType polyline_type) {
      _polylineType = polyline_type;
    }

    alias DuplicateEdges = GraphOptions.DuplicateEdges;

    /**
     * Indicates whether duplicate edges in the input should be kept (KEEP) or
     * merged together (MERGE).  Note you can use edge labels to determine
     * which input edges were merged into a given output edge.
     *
     * DEFAULT: DuplicateEdges::KEEP
     */
    DuplicateEdges duplicateEdges() const {
      return _duplicateEdges;
    }

    void setDuplicateEdges(DuplicateEdges duplicate_edges) {
      _duplicateEdges = duplicate_edges;
    }

    alias SiblingPairs = GraphOptions.SiblingPairs;

    /**
     * Indicates whether sibling edge pairs (i.e., pairs consisting of an edge
     * and its reverse edge) should be kept (KEEP) or discarded (DISCARD).
     * For example, if a polyline backtracks on itself, the DISCARD option
     * would cause this section of the polyline to be removed.  Note that this
     * option may cause a single polyline to split into several pieces (e.g.,
     * if a polyline has a "lollipop" shape).
     *
     * REQUIRES: sibling_pairs == { DISCARD, KEEP }
     *           (the CREATE and REQUIRE options are not allowed)
     *
     * DEFAULT: SiblingPairs::KEEP
     */
    SiblingPairs siblingPairs() const {
      return _siblingPairs;
    }

    void setSiblingPairs(SiblingPairs sibling_pairs)
    in {
      assert(sibling_pairs == SiblingPairs.KEEP || sibling_pairs == SiblingPairs.DISCARD);
    } do {
      _siblingPairs = sibling_pairs;
    }

    /**
     * If true, calls FindValidationError() on each output polyline.  If any
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
      setS2debugOverride(S2Debug.DISABLE);
    }

    // This method can turn off the automatic validity checks triggered by the
    // --s2debug flag (which is on by default in debug builds).  The main
    // reason to do this is if your code already does its own error checking,
    // or if you need to work with invalid geometry for some reason.
    //
    // In any case, polylines have very few restrictions so they are unlikely
    // to have errors.  Errors include vertices that aren't unit length (which
    // can only happen if they are present in the input data), or adjacent
    // vertices that are at antipodal points on the sphere (unlikely with real
    // data).  The other possible error is adjacent identical vertices, but
    // this can't happen because S2Builder does not generate such polylines.
    //
    // DEFAULT: S2Debug::ALLOW
    S2Debug s2debugOverride() const {
      return _s2debugOverride;
    }

    void setS2debugOverride(S2Debug s2debug_override) {
      _s2debugOverride = s2debug_override;
    }

  private:
    S2Builder.EdgeType _edgeType = S2Builder.EdgeType.DIRECTED;
    PolylineType _polylineType = PolylineType.PATH;
    DuplicateEdges _duplicateEdges = DuplicateEdges.KEEP;
    SiblingPairs _siblingPairs = SiblingPairs.KEEP;
    bool _validate = false;
    S2Debug _s2debugOverride = S2Debug.ALLOW;
  }

  /// Specifies that a vector of polylines should be constructed using the
  /// given options.
  this(S2Polyline[]* polylines, Options options = Options()) {
    initialize(polylines, null, null, options);
  }

  alias LabelSetIds = LabelSetId[][];

  /**
   * Specifies that a vector of polylines should be constructed using the
   * given options, and that any labels attached to the input edges should be
   * returned in "label_set_ids" and "label_set_lexicion".
   *
   * The labels associated with the edge "polyline[i].vertex({j, j+1})" can be
   * retrieved as follows:
   *
   *   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
   */
  this(S2Polyline[]* polylines,
      LabelSetIds* label_set_ids,
      IdSetLexicon* label_set_lexicon,
      Options options = Options()) {
    initialize(polylines, label_set_ids, label_set_lexicon, options);
  }

  // Layer interface:
  override
  GraphOptions graphOptions() const {
    return new GraphOptions(_options.edgeType(), DegenerateEdges.DISCARD,
        _options.duplicateEdges(), _options.siblingPairs());
  }

  override
  void build(in Graph g, ref S2Error error) {
    Graph.EdgePolyline[] edge_polylines = g.getPolylines(_options.polylineType());
    (*_polylines).reserve(edge_polylines.length);
    if (_labelSetIds) (*_labelSetIds).reserve(edge_polylines.length);
    S2Point[] vertices;  // Temporary storage for vertices.
    Label[] labels;  // Temporary storage for labels.
    foreach (edge_polyline; edge_polylines) {
      vertices ~= g.vertex(g.edge(edge_polyline[0])[0]);
      foreach (Graph.EdgeId e; edge_polyline) {
        vertices ~= g.vertex(g.edge(e)[1]);
      }
      S2Polyline polyline = new S2Polyline(vertices, _options.s2debugOverride());
      vertices.length = 0;
      if (_options.validate()) {
        polyline.findValidationError(error);
      }
      (*_polylines) ~= polyline;
      if (_labelSetIds) {
        auto fetcher = new Graph.LabelFetcher(g, _options.edgeType());
        LabelSetId[] polyline_labels;
        polyline_labels.reserve(edge_polyline.length);
        foreach (Graph.EdgeId e; edge_polyline) {
          fetcher.fetch(e, labels);
          polyline_labels ~= _labelSetLexicon.add(labels);
        }
        (*_labelSetIds) ~= polyline_labels;
      }
    }
  }

private:
  void initialize(S2Polyline[]* polylines, LabelSetIds* label_set_ids,
      IdSetLexicon* label_set_lexicon, Options options)
  in {
    assert((label_set_ids is null) == (label_set_lexicon is null));
  } do {
    _polylines = polylines;
    _labelSetIds = label_set_ids;
    _labelSetLexicon = label_set_lexicon;
    _options = options;
  }

  S2Polyline[]* _polylines;
  LabelSetIds* _labelSetIds;
  IdSetLexicon* _labelSetLexicon;
  Options _options;
}

/// Like S2PolylineVectorLayer, but adds the polylines to a MutableS2ShapeIndex.
class IndexedS2PolylineVectorLayer : Layer {
public:
  alias Options = S2PolylineVectorLayer.Options;

  this(MutableS2ShapeIndex index, Options options = Options()) {
    _index = index;
    _layer = new S2PolylineVectorLayer(&_polylines, options);
  }

  override
  GraphOptions graphOptions() const {
    return _layer.graphOptions();
  }

  override
  void build(in Graph g, ref S2Error error) {
    _layer.build(g, error);
    if (error.ok()) {
      foreach (polyline; _polylines) {
        _index.add(new S2Polyline.Shape(polyline));
      }
    }
  }

private:
  MutableS2ShapeIndex _index;
  S2Polyline[] _polylines;
  S2PolylineVectorLayer _layer;
}
