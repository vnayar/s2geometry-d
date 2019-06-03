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
//
// Note that there are two supported output types for polygons: S2Polygon and
// S2LaxPolygonShape.  Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use S2LaxPolygonShape if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from a S2LaxPolygonShape to an S2Polygon you will need to use
// S2Builder again.
//
// Similarly, there are two supported output formats for polygon meshes:
// S2LaxPolygonShapeVector and S2PolygonMesh.  Use S2PolygonMesh if you need
// to be able to determine which polygons are adjacent to each edge or vertex;
// otherwise use S2LaxPolygonShapeVector, which uses less memory and is faster
// to construct.

module s2.builder.util.s2polygon_layer;

import s2.builder.graph;
import s2.builder.layer;
import s2.s2point;
import s2.id_set_lexicon;
import s2.mutable_s2shape_index;
import s2.s2builder;
import s2.s2debug;
import s2.s2error;
import s2.s2loop;
import s2.s2polygon;

import std.algorithm;
import std.range;
import std.typecons;

alias LoopType = Graph.LoopType;
alias SiblingPairs = GraphOptions.SiblingPairs;
alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;

/**
 * A layer type that assembles edges (directed or undirected) into an
 * S2Polygon.  Returns an error if the edges cannot be assembled into loops.
 *
 * If the input edges are directed, they must be oriented such that the
 * polygon interior is to the left of all edges.  Directed edges are always
 * preferred (see S2Builder::EdgeType).
 *
 * Before the edges are assembled into loops, "sibling pairs" consisting of an
 * edge and its reverse edge are automatically removed.  Such edge pairs
 * represent zero-area degenerate regions, which S2Polygon does not allow.
 * (If you need to build polygons with degeneracies, use LaxPolygonLayer
 * instead.)
 *
 * S2PolygonLayer is implemented such that if the input to S2Builder is a
 * polygon and is not modified, then the output has the same cyclic ordering
 * of loop vertices and the same loop ordering as the input polygon.
 *
 * CAVEAT: Because polygons are constructed from their boundaries, this method
 * cannot distinguish between the empty and full polygons.  An empty boundary
 * always yields an empty polygon.  If the result should sometimes be the full
 * polygon, such logic must be implemented outside of this class (and will
 * need to consider factors other than the polygon's boundary).
 */
class S2PolygonLayer : Layer {
public:

  struct Options {
  public:
    /**
     * Indicates whether the input edges provided to S2Builder are directed or
     * undirected.  Directed edges should be used whenever possible (see
     * S2Builder::EdgeType for details).
     *
     * If the input edges are directed, they should be oriented so that the
     * polygon interior is to the left of all edges.  This means that for a
     * polygon with holes, the outer loops ("shells") should be directed
     * counter-clockwise while the inner loops ("holes") should be directed
     * clockwise.  Note that S2Builder::AddPolygon() does this automatically.
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
     * If true, calls FindValidationError() on the output polygon.  If any
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

  /// Specifies that a polygon should be constructed using the given options.
  this(S2Polygon polygon, Options options = Options()) {
    initialize(polygon, null, null, options);
  }

  alias LabelSetIds = LabelSetId[][];

  /**
   * Specifies that a polygon should be constructed using the given options,
   * and that any labels attached to the input edges should be returned in
   * "label_set_ids" and "label_set_lexicion".
   *
   * The labels associated with the edge "polygon.loop(i).vertex({j, j+1})"
   * can be retrieved as follows:
   *
   *   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
   */
  this(S2Polygon polygon, ref LabelSetIds label_set_ids,
      IdSetLexicon label_set_lexicon,
      in Options options = Options()) {
    initialize(polygon, &label_set_ids, label_set_lexicon, options);
  }


  // Layer interface:
  override
  GraphOptions graphOptions() const {
    // Prevent degenerate edges and sibling edge pairs.  There should not be any
    // duplicate edges if the input is valid, but if there are then we keep them
    // since this tends to produce more comprehensible errors.
    return new GraphOptions(_options.edgeType(), DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP, SiblingPairs.DISCARD);
  }

  override
  void build(in Graph g, ref S2Error error) {
    if (_labelSetIds) _labelSetIds.length = 0;

    // It's tricky to compute the edge labels for S2Polygons because the
    // S2Polygon::Init methods can reorder and/or invert the loops.  We handle
    // this by remembering the original vector index of each loop and whether or
    // not the loop contained S2::Origin().  By comparing this with the final
    // S2Polygon loops we can fix up the edge labels appropriately.
    LoopMap loop_map;
    if (g.options().edgeType() == EdgeType.DIRECTED) {
      Graph.EdgeLoop[] edge_loops;
      if (!g.getDirectedLoops(LoopType.SIMPLE, edge_loops, error)) {
        return;
      }
      S2Loop[] loops;
      appendS2Loops(g, edge_loops, loops);
      appendEdgeLabels(g, edge_loops);
      initializeLoopMap(loops, loop_map);
      _polygon.initializeOriented(loops);
    } else {
      Graph.UndirectedComponent[] components;
      if (!g.getUndirectedComponents(LoopType.SIMPLE, components, error)) {
        return;
      }
      // It doesn't really matter which complement of each component we use,
      // since below we normalize all the loops so that they enclose at most
      // half of the sphere (to ensure that the loops can always be nested).
      //
      // The only reason to prefer one over the other is that when there are
      // multiple loops that touch, only one of the two complements matches the
      // structure of the input loops.  GetUndirectedComponents() tries to
      // ensure that this is always complement 0 of each component.
      S2Loop[] loops;
      foreach (component; components) {
        appendS2Loops(g, component[0], loops);
        appendEdgeLabels(g, component[0]);
      }
      initializeLoopMap(loops, loop_map);
      foreach (loop; loops) loop.normalize();
      _polygon.initializeNested(loops);
    }
    reorderEdgeLabels(loop_map);
    if (_options.validate()) {
      _polygon.findValidationError(error);
    }
  }

private:
  void initialize(S2Polygon polygon, LabelSetIds* label_set_ids,
      IdSetLexicon label_set_lexicon, in Options options)
  in {
    assert((label_set_ids is null) == (label_set_lexicon is null));
  } do {
    _polygon = polygon;
    _labelSetIds = label_set_ids;
    _labelSetLexicon = label_set_lexicon;
    _options = options;

    if (_options.validate()) {
      _polygon.setS2debugOverride(S2Debug.DISABLE);
    }
  }

  void appendS2Loops(in Graph g,
      in Graph.EdgeLoop[] edge_loops,
      ref S2Loop[] loops) const {
    S2Point[] vertices;
    foreach (edge_loop; edge_loops) {
      vertices.reserve(edge_loop.length);
      foreach (edge_id; edge_loop) {
        vertices ~= g.vertex(g.edge(edge_id)[0]);
      }
      loops ~= new S2Loop(vertices, _polygon.s2debugOverride());
      vertices.length = 0;
    }
  }

  void appendEdgeLabels(in Graph g, in Graph.EdgeLoop[] edge_loops) {
    if (!_labelSetIds) return;

    Label[] labels;  // Temporary storage for labels.
    auto fetcher = new Graph.LabelFetcher(g, _options.edgeType());
    foreach (edge_loop; edge_loops) {
      LabelSetId[] loop_label_set_ids;
      loop_label_set_ids.reserve(edge_loop.length);
      foreach (edge_id; edge_loop) {
        fetcher.fetch(edge_id, labels);
        loop_label_set_ids ~= _labelSetLexicon.add(labels);
      }
      *_labelSetIds ~= loop_label_set_ids;
    }
  }

  alias LoopMap = Tuple!(int, bool)[const(S2Loop)];

  void initializeLoopMap(in S2Loop[] loops, ref LoopMap loop_map) const {
    if (!_labelSetIds) return;
    foreach (i, loop; loops) {
      loop_map[loop] = tuple(cast(int) i, loop.containsOrigin());
    }
  }

  void reorderEdgeLabels(in LoopMap loop_map) {
    if (!_labelSetIds) return;
    auto new_ids = new LabelSetIds(_labelSetIds.length);
    for (int i = 0; i < _polygon.numLoops(); ++i) {
      S2Loop loop = _polygon.loop(i);
      Tuple!(int, bool) old = loop_map[loop];
      new_ids[i] = (*_labelSetIds)[old[0]];
      if (loop.containsOrigin() != old[1]) {
        // S2Loop::Invert() reverses the order of the vertices, which leaves
        // the last edge unchanged.  For example, the loop ABCD (with edges
        // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
        reverse(new_ids[i][0 .. $-1]);
      }
    }
    *_labelSetIds = new_ids;
  }

  S2Polygon _polygon;
  LabelSetIds* _labelSetIds;
  IdSetLexicon _labelSetLexicon;
  Options _options;
}

/**
 * Like S2PolygonLayer, but adds the polygon to a MutableS2ShapeIndex (if the
 * polygon is non-empty).
 */
class IndexedS2PolygonLayer : Layer {
public:
  alias Options = S2PolygonLayer.Options;

  this(MutableS2ShapeIndex index, Options options = Options()) {
    _index = index;
    _polygon = new S2Polygon();
    _layer = new S2PolygonLayer(_polygon, options);
  }

  override
  GraphOptions graphOptions() const {
    return _layer.graphOptions();
  }

  override
  void build(in Graph g, ref S2Error error) {
    _layer.build(g, error);
    if (error.ok() && !_polygon.isEmpty()) {
      _index.add(new S2Polygon.Shape(_polygon));
    }
  }

private:
  MutableS2ShapeIndex _index;
  S2Polygon _polygon;
  S2PolygonLayer _layer;
}
