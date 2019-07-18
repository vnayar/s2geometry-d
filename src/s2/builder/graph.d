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

module s2.builder.graph;

import s2.s2builder;
import s2.s2error;
import s2.s2point;
import s2.id_set_lexicon;

import std.algorithm;
import std.exception;
import std.range;
import std.typecons;

/**
 * An S2Builder::Graph represents a collection of snapped edges that is passed
 * to a Layer for assembly.  (Example layers include polygons, polylines, and
 * polygon meshes.)  The Graph object does not own any of its underlying data;
 * it is simply a view of data that is stored elsewhere.  You will only
 * need this interface if you want to implement a new Layer subtype.
 *
 * The graph consists of vertices and directed edges.  Vertices are numbered
 * sequentially starting from zero.  An edge is represented as a pair of
 * vertex ids.  The edges are sorted in lexicographic order, therefore all of
 * the outgoing edges from a particular vertex form a contiguous range.
 *
 * S2Builder::Graph is movable and copyable.  Note that although this class
 * does not own the underlying vertex and edge data, S2Builder guarantees that
 * all Graph objects passed to S2Builder::Layer::Build() methods will remain
 * valid until all layers have been built.
 *
 * TODO(ericv): Consider pulling out the methods that are helper functions for
 * Layer implementations (such as GetDirectedLoops) into s2builderutil_graph.h.
 */
class Graph {
public:
  // Identifies a vertex in the graph.  Vertices are numbered sequentially
  // starting from zero.
  alias VertexId = int;

  // Defines an edge as an (origin, destination) vertex pair.
  alias Edge = Tuple!(VertexId, VertexId);

  // Identifies an edge in the graph.  Edges are numbered sequentially
  // starting from zero.
  alias EdgeId = int;

  alias EdgeType = S2Builder.EdgeType;

  // Identifies an S2Builder *input* edge (before snapping).
  alias InputEdgeId = S2Builder.InputEdgeId;

  // Identifies a set of S2Builder input edges.
  alias InputEdgeIdSetId = S2Builder.InputEdgeIdSetId;

  // Identifies a set of edge labels.
  alias LabelSetId = S2Builder.LabelSetId;

  // Determines whether a degenerate polygon is empty or full.
  alias IsFullPolygonPredicate = S2Builder.IsFullPolygonPredicate;

  alias SiblingPairs = GraphOptions.SiblingPairs;
  alias DegenerateEdges = GraphOptions.DegenerateEdges;
  alias DuplicateEdges = GraphOptions.DuplicateEdges;

  /// The default constructor exists only for the benefit of STL containers.
  /// The graph must be initialized (using the assignment operator) before it
  /// is used.
  this() {
    _options = new GraphOptions();
    _numVertices = -1;
    _vertices = null;
    _edges = null;
    _inputEdgeIdSetIds = null;
    _inputEdgeIdSetLexicon = null;
    _labelSetIds = null;
    _labelSetLexicon = null;
  }

  /**
   * Note that most of the parameters are passed by const reference and must
   * exist for the duration of the Graph object.  Notes on parameters:
   * "options":
   *    - the GraphOptions used to build the Graph.  In some cases these
   *      can be different than the options provided by the Layer.
   * "vertices":
   *   - a vector of S2Points indexed by VertexId.
   * "edges":
   *   - a vector of VertexId pairs (sorted in lexicographic order)
   *     indexed by EdgeId.
   * "input_edge_id_set_ids":
   *   - a vector indexed by EdgeId that allows access to the set of
   *     InputEdgeIds that were mapped to the given edge, by looking up the
   *     returned value (an InputEdgeIdSetId) in "input_edge_id_set_lexicon".
   * "input_edge_id_set_lexicon":
   *   - a class that maps an InputEdgeIdSetId to a set of InputEdgeIds.
   * "label_set_ids":
   *   - a vector indexed by InputEdgeId that allows access to the set of
   *     labels that were attached to the given input edge, by looking up the
   *     returned value (a LabelSetId) in the "label_set_lexicon".
   * "label_set_lexicon":
   *   - a class that maps a LabelSetId to a set of S2Builder::Labels.
   * "is_full_polygon_predicate":
   *   - a predicate called to determine whether a graph consisting only of
   *     polygon degeneracies represents the empty polygon or the full polygon
   *     (see s2builder.h for details).
   */
  this(
      in GraphOptions options,
      in S2Point[] vertices,
      in Edge[] edges,
      in InputEdgeIdSetId[] input_edge_id_set_ids,
      in IdSetLexicon input_edge_id_set_lexicon,
      in LabelSetId[] label_set_ids,
      in IdSetLexicon label_set_lexicon,
      // TODO(ericv/hagzonal): Fix st_lib and remove default parameter.
      IsFullPolygonPredicate is_full_polygon_predicate = IsFullPolygonPredicate()) {
    _options = options;
    _numVertices = cast(int) vertices.length;
    _vertices = vertices;
    _edges = edges;
    _inputEdgeIdSetIds = input_edge_id_set_ids;
    _inputEdgeIdSetLexicon = input_edge_id_set_lexicon;
    _labelSetIds = label_set_ids;
    _labelSetLexicon = label_set_lexicon;
    _isFullPolygonPredicate = is_full_polygon_predicate;
    debug enforce(isSorted(edges));
  }

  const(GraphOptions) options() const {
    return _options;
  }

  /// Returns the number of vertices in the graph.
  VertexId numVertices() const {
    return _numVertices;  // vertices_.size() requires division by 24.
  }

  /// Returns the vertex at the given index.
  const(S2Point) vertex(VertexId v) const {
    return vertices()[v];
  }

  /// Returns the entire set of vertices.
  const(S2Point[]) vertices() const {
    return _vertices;
  }

  /// Returns the total number of edges in the graph.
  EdgeId numEdges() const {
    return cast(EdgeId)(edges().length);
  }

  /// Returns the endpoints of the given edge (as vertex indices).
  Edge edge(EdgeId e) const {
    return edges()[e];
  }

  /// Returns the entire set of edges.
  const(Edge[]) edges() const {
    return _edges;
  }

  /// Given an edge (src, dst), returns the reverse edge (dst, src).
  static Edge reverse(in Edge e) {
    return Edge(e[1], e[0]);
  }

  /**
   * Returns a vector of edge ids sorted in lexicographic order by
   * (destination, origin).  All of the incoming edges to each vertex form a
   * contiguous subrange of this ordering.
   */
  EdgeId[] getInEdgeIds() const {
    EdgeId[] in_edge_ids = new EdgeId[](numEdges());
    for (EdgeId i = 0; i < numEdges(); i++) in_edge_ids[i] = i;
    sort!((ai, bi) {
          return stableLessThan(reverse(edge(ai)), reverse(edge(bi)), ai, bi);
        })(in_edge_ids);
    return in_edge_ids;
  }

  /**
   * Given a graph such that every directed edge has a sibling, returns a map
   * from EdgeId to the sibling EdgeId.  This method is identical to
   * GetInEdgeIds() except that (1) it requires edges to have siblings, and
   * (2) undirected degenerate edges are grouped together in pairs such that
   * one edge is the sibling of the other.  Handles duplicate edges correctly
   * and is also consistent with GetLeftTurnMap().
   *
   * REQUIRES: An option is chosen that guarantees sibling pairs:
   *     (options.sibling_pairs() == { REQUIRE, CREATE } ||
   *      options.edge_type() == UNDIRECTED)
   */
  EdgeId[] getSiblingMap() const {
    EdgeId[] in_edge_ids = getInEdgeIds();
    makeSiblingMap(in_edge_ids);
    return in_edge_ids;
  }

  /**
   * Like GetSiblingMap(), but constructs the map starting from the vector of
   * incoming edge ids returned by GetInEdgeIds().  (This operation is a no-op
   * except unless undirected degenerate edges are present, in which case such
   * edges are grouped together in pairs to satisfy the requirement that every
   * edge must have a sibling edge.)
   */
  void makeSiblingMap(ref EdgeId[] in_edge_ids) const {
    debug enforce(
        _options.siblingPairs() == SiblingPairs.REQUIRE
        || _options.siblingPairs() == SiblingPairs.CREATE
        || _options.edgeType() == EdgeType.UNDIRECTED);
    for (EdgeId e = 0; e < numEdges(); ++e) {
      debug enforce(edge(e) == reverse(edge(in_edge_ids[e])));
    }
    if (_options.edgeType() == EdgeType.DIRECTED) return;
    if (_options.degenerateEdges() == DegenerateEdges.DISCARD) return;

    for (EdgeId e = 0; e < numEdges(); ++e) {
      VertexId v = edge(e)[0];
      if (edge(e)[1] == v) {
        debug enforce(e + 1 < numEdges());
        debug enforce(edge(e + 1)[0] == v);
        debug enforce(edge(e + 1)[1] == v);
        debug enforce(in_edge_ids[e] == e);
        debug enforce(in_edge_ids[e + 1] == e + 1);
        in_edge_ids[e] = e + 1;
        in_edge_ids[e + 1] = e;
        ++e;
      }
    }
  }

  /// A helper class for VertexOutMap that represents the outgoing edges
  /// from a given vertex.
  // class VertexOutEdges {
  // public:
  //   const Edge* begin() const { return begin_; }
  //   const Edge* end() const { return end_; }
  //   size_t size() const { return end_ - begin_; }
  // private:
  //   friend class VertexOutMap;
  //   VertexOutEdges(const Edge* begin, const Edge* end);
  //   const Edge* begin_;
  //   const Edge* end_;
  // };

  /// A helper class for VertexOutMap that represents the outgoing edge *ids*
  /// from a given vertex.
  // class VertexOutEdgeIds {
  // public:
  //   // An iterator over a range of edge ids (like boost::counting_iterator).
  //   class Iterator {
  //   public:
  //     explicit Iterator(EdgeId id) : id_(id) {}
  //     const EdgeId& operator*() const { return id_; }
  //     Iterator& operator++() { ++id_; return *this; }
  //     Iterator operator++(int) { return Iterator(id_++); }
  //     size_t operator-(const Iterator& x) const { return id_ - x.id_; }
  //     bool operator==(const Iterator& x) const { return id_ == x.id_; }
  //     bool operator!=(const Iterator& x) const { return id_ != x.id_; }
  //   private:
  //     EdgeId id_;
  //   };
  //   Iterator begin() const { return Iterator(begin_); }
  //   Iterator end() const { return Iterator(end_); }
  //   size_t size() const { return end_ - begin_; }
  // private:
  //   friend class VertexOutMap;
  //   VertexOutEdgeIds(EdgeId begin, EdgeId end);
  //   EdgeId begin_, end_;
  // };

  /**
   * A class that maps vertices to their outgoing edge ids.  Example usage:
   *   VertexOutMap out(g);
   *   for (Graph::EdgeId e : out.edge_ids(v)) { ... }
   *   for (const Graph::Edge& edge : out.edges(v)) { ... }
   */
  static class VertexOutMap {
  public:
    this(const(Graph) g) {
      _edges = g.edges();
      _edgeBegins = new EdgeId[](g.numVertices() + 1);
      EdgeId e = 0;
      for (VertexId v = 0; v <= g.numVertices(); ++v) {
        while (e < g.numEdges() && g.edge(e)[0] < v) ++e;
        _edgeBegins[v] = e;
      }
    }

    int degree(VertexId v) const {
      return cast(int) edgeIds(v).length;
    }

    const(Edge[]) edges(VertexId v) const {
      return _edges[_edgeBegins[v] .. _edgeBegins[v + 1]];
    }

    // Returns a range of EdgeId from the given VertexId to the next.
    auto edgeIds(VertexId v) const {
      return iota(_edgeBegins[v], _edgeBegins[v + 1]);
    }

    /// Returns a range of Edge (or edge ids) between a specific pair of vertices.
    auto edges(VertexId v0, VertexId v1) const {
      return assumeSorted(_edges[_edgeBegins[v0] .. _edgeBegins[v0 + 1]])
          .equalRange(Edge(v0, v1));
    }

    EdgeId[] edgeIds(VertexId v0, VertexId v1) const {
      auto trisectRanges = assumeSorted(_edges[_edgeBegins[v0] .. _edgeBegins[v0 + 1]])
          .trisect(Edge(v0, v1));

      // Start past the ids that are less than those that match Edge(v0, v1).
      EdgeId startId = _edgeBegins[v0] + cast(int) trisectRanges[0].length;
      return iota(startId, startId + cast(int) trisectRanges[1].length).array;
    }

  private:
    const(Edge[]) _edges;
    EdgeId[] _edgeBegins;
  }

  // A helper class for VertexInMap that represents the incoming edge *ids*
  // to a given vertex.
  // class VertexInEdgeIds {
  //  public:
  //   const EdgeId* begin() const { return begin_; }
  //   const EdgeId* end() const { return end_; }
  //   size_t size() const { return end_ - begin_; }

  //  private:
  //   friend class VertexInMap;
  //   VertexInEdgeIds(const EdgeId* begin, const EdgeId* end);
  //   const EdgeId* begin_;
  //   const EdgeId* end_;
  // }

  // A class that maps vertices to their incoming edge ids.  Example usage:
  //   VertexInMap in(g);
  //   for (Graph::EdgeId e : in.edge_ids(v)) { ... }
  static class VertexInMap {
  public:
    this(in Graph g) {
      _inEdgeIds = g.getInEdgeIds();
      _inEdgeBegins = new EdgeId[](g.numVertices() + 1);
      EdgeId e = 0;
      for (VertexId v = 0; v <= g.numVertices(); ++v) {
        while (e < g.numEdges() && g.edge(_inEdgeIds[e])[1] < v) ++e;
        _inEdgeBegins[v] = e;
      }
    }

    int degree(VertexId v) const {
      return cast(int)(edgeIds(v).length);
    }

    const(EdgeId[]) edgeIds(VertexId v) const {
      return _inEdgeIds[_inEdgeBegins[v] .. _inEdgeBegins[v + 1]];
    }

    // Returns a sorted vector of all incoming edges (see GetInEdgeIds).
    const(EdgeId[]) inEdgeIds() const {
      return _inEdgeIds;
    }

  private:
    EdgeId[] _inEdgeIds;
    EdgeId[] _inEdgeBegins;
  }

  /// Defines a value larger than any valid InputEdgeId.
  static immutable InputEdgeId MAX_INPUT_EDGE_ID = InputEdgeId.max;

  /// The following value of InputEdgeId means that an edge does not
  /// corresponds to any input edge.
  static immutable InputEdgeId NO_INPUT_EDGE_ID = MAX_INPUT_EDGE_ID - 1;

  /**
   * Returns the set of input edge ids that were snapped to the given
   * edge.  ("Input edge ids" are assigned to input edges sequentially in
   * the order they are added to the builder.)  For example, if input
   * edges 2 and 17 were snapped to edge 12, then input_edge_ids(12)
   * returns a set containing the numbers 2 and 17.  Example usage:
   *
   *   for (InputEdgeId input_edge_id : g.input_edge_ids(e)) { ... }
   *
   * Please note the following:
   *
   *  - When edge chains are simplified, the simplified edge is assigned all
   *    the input edge ids associated with edges of the chain.
   *
   *  - Edges can also have multiple input edge ids due to edge merging
   *    (if DuplicateEdges::MERGE is specified).
   *
   *  - Siblings edges automatically created by EdgeType::UNDIRECTED or
   *    SiblingPairs::CREATE have an empty set of input edge ids.  (However
   *    you can use a LabelFetcher to retrieve the set of labels associated
   *    with both edges of a given sibling pair.)
   */
  const(IdSetLexicon.IdSet) inputEdgeIds(EdgeId e) const {
    return inputEdgeIdSetLexicon().idSet(inputEdgeIdSetIds()[e]);
  }

  /**
   * Low-level method that returns an integer representing the entire set of
   * input edge ids that were snapped to the given edge.  The elements of the
   * IdSet can be accessed using input_edge_id_set_lexicon().
   */
  InputEdgeIdSetId inputEdgeIdSetId(EdgeId e) const {
    return inputEdgeIdSetIds()[e];
  }

  /// Low-level method that returns a vector where each element represents the
  /// set of input edge ids that were snapped to a particular output edge.
  const(InputEdgeIdSetId[]) inputEdgeIdSetIds() const {
    return _inputEdgeIdSetIds;
  }

  /// Returns a mapping from an InputEdgeIdSetId to a set of input edge ids.
  const(IdSetLexicon) inputEdgeIdSetLexicon() const {
    return _inputEdgeIdSetLexicon;
  }

  /**
   * Returns the minimum input edge id that was snapped to this edge, or -1 if
   * no input edges were snapped (see SiblingPairs::CREATE).  This is
   * useful for layers that wish to preserve the input edge ordering as much
   * as possible (e.g., to ensure idempotency).
   */
  InputEdgeId minInputEdgeId(EdgeId e) const {
    const(IdSetLexicon.IdSet) id_set = inputEdgeIds(e);
    return (id_set.length == 0) ? NO_INPUT_EDGE_ID : id_set[0];
  }

  /// Returns a vector containing the minimum input edge id for every edge.
  /// If an edge has no input ids, kNoInputEdgeId is used.
  InputEdgeId[] getMinInputEdgeIds() const {
    auto min_input_ids = new InputEdgeId[numEdges()];
    for (EdgeId e = 0; e < numEdges(); ++e) {
      min_input_ids[e] = minInputEdgeId(e);
    }
    return min_input_ids;
  }

  /// Returns a vector of EdgeIds sorted by minimum input edge id.  This is an
  /// approximation of the input edge ordering.
  EdgeId[] getInputEdgeOrder(in InputEdgeId[] input_ids) const {
    auto order = new EdgeId[input_ids.length];
    for (EdgeId i = 0; i < order.length; i++) order[i] = i;
    // Comparison function ensures sort is stable.
    sort!((EdgeId a, EdgeId b) => Edge(input_ids[a], a) < Edge(input_ids[b], b))(order);
    return order;
  }

  /**
   * Convenience class to return the set of labels associated with a given
   * graph edge.  Note that due to snapping, one graph edge may correspond to
   * several different input edges and will have all of their labels.
   * This class is the preferred way to retrieve edge labels.
   *
   * The reason this is a class rather than a graph method is because for
   * undirected edges, we need to fetch the labels associated with both
   * siblings.  This is because only the original edge of the sibling pair has
   * labels; the automatically generated sibling edge does not.
   */
  static class LabelFetcher {
  public:
    /**
     * Prepares to fetch labels associated with the given edge type.  For
     * EdgeType::UNDIRECTED, labels associated with both edges of the sibling
     * pair will be returned.  "edge_type" is a parameter (rather than using
     * g.options().edge_type()) so that clients can explicitly control whether
     * labels from one or both siblings are returned.
     */
    this(in Graph g, EdgeType edge_type) {
      _g = g;
      _edgeType = edge_type;
      if (edge_type == EdgeType.UNDIRECTED) _siblingMap = g.getSiblingMap();
    }

    /**
     * Returns the set of labels associated with edge "e" (and also the labels
     * associated with the sibling of "e" if edge_type() is UNDIRECTED).
     * Labels are sorted and duplicate labels are automatically removed.
     *
     * This method uses an output parameter rather than returning by value in
     * order to avoid allocating a new vector on every call to this method.
     */
    void fetch(EdgeId e, ref S2Builder.Label[] labels) {
      labels.length = 0;
      foreach (InputEdgeId input_edge_id; _g.inputEdgeIds(e)) {
        foreach (S2Builder.Label label; _g.labels(input_edge_id)) {
          labels ~= label;
        }
      }
      if (_edgeType == EdgeType.UNDIRECTED) {
        foreach (InputEdgeId input_edge_id; _g.inputEdgeIds(_siblingMap[e])) {
          foreach (S2Builder.Label label; _g.labels(input_edge_id)) {
            labels ~= label;
          }
        }
      }
      if (labels.length > 1) {
        labels = sort(labels).uniq.array;
      }
    }

  private:
    const(Graph) _g;
    EdgeType _edgeType;
    EdgeId[] _siblingMap;
  }

  /// Returns the set of labels associated with a given input edge.  Example:
  ///   for (Label label : g.labels(input_edge_id)) { ... }
  const(IdSetLexicon.IdSet) labels(InputEdgeId id) const {
    return labelSetLexicon().idSet(labelSetIds()[id]);
  }

  /**
   * Low-level method that returns an integer representing the set of
   * labels associated with a given input edge.  The elements of
   * the IdSet can be accessed using label_set_lexicon().
   */
  LabelSetId labelSetId(InputEdgeId e) const {
    return _labelSetIds[e];
  }

  /// Low-level method that returns a vector where each element represents the
  /// set of labels associated with a particular output edge.
  const(LabelSetId[]) labelSetIds() const {
    return _labelSetIds;
  }

  /// Returns a mapping from a LabelSetId to a set of labels.
  const(IdSetLexicon) labelSetLexicon() const {
    return _labelSetLexicon;
  }

  /**
   * Convenience method that calls is_full_polygon_predicate() to determine
   * whether a graph that consists only of polygon degeneracies represents the
   * empty polygon or the full polygon (see s2builder.h for details).
   */
  bool isFullPolygon(ref S2Error error) const {
    return _isFullPolygonPredicate(this, error);
  }

  /**
   * Returns a method that determines whether a graph that consists only of
   * polygon degeneracies represents the empty polygon or the full polygon
   * (see s2builder.h for details).
   */
  const(IsFullPolygonPredicate) isFullPolygonPredicate() const {
    return _isFullPolygonPredicate;
  }

  /**
   * Returns a map "m" that maps each edge e=(v0,v1) to the following outgoing
   * edge around "v1" in clockwise order.  \(This corresponds to making a "left
   * turn" at the vertex.\)  By starting at a given edge and making only left
   * turns, you can construct a loop whose interior does not contain any edges
   * in the same connected component.
   *
   * If the incoming and outgoing edges around a vertex do not alternate
   * perfectly \(e.g., there are two incoming edges in a row\), then adjacent
   * \(incoming, outgoing\) pairs are repeatedly matched and removed.  This is
   * similar to finding matching parentheses in a string such as "\(\(\)\(\)\)\(\)".
   *
   * For sibling edge pairs, the incoming edge is assumed to immediately
   * follow the outgoing edge in clockwise order.  Thus a left turn is made
   * from an edge to its sibling only if there are no other outgoing edges.
   * With respect to the parentheses analogy, a sibling pair is "\(\)".
   * Similarly, if there are multiple copies of a sibling edge pair then the
   * duplicate incoming and outgoing edges are sorted in alternating order
   * \(e.g., "\)\(\)\("\).
   *
   * Degenerate edges \(edges from a vertex to itself\) are treated as loops
   * consisting of a single edge.  This avoids the problem of deciding the
   * connectivity and ordering of such edges when they share a vertex with
   * other edges \(possibly including other degenerate edges\).
   *
   * If it is not possible to make a left turn from every input edge, this
   * method returns false and sets "error" appropriately.  In this situation
   * the left turn map is still valid except that any incoming edge where it
   * is not possible to make a left turn will have its entry set to -1.
   *
   * "in_edge_ids" should be equal to GetInEdgeIds() or GetSiblingMap().
   */
  bool getLeftTurnMap(
      in EdgeId[] in_edge_ids, ref EdgeId[] left_turn_map, ref S2Error error) const {
    left_turn_map.length = numEdges();
    left_turn_map[] = -1;
    if (numEdges() == 0) return true;

    // Declare vectors outside the loop to avoid reallocating them each time.
    VertexEdge[] v0_edges;
    EdgeId[] e_in;
    EdgeId[] e_out;

    // Walk through the two sorted arrays of edges (outgoing and incoming) and
    // gather all the edges incident to each vertex.  Then we sort those edges
    // and add an entry to the left turn map from each incoming edge to the
    // immediately following outgoing edge in clockwise order.
    int outId = 0;
    int inId = 0;
    Edge out_edge = edge(outId);
    Edge in_edge = edge(in_edge_ids[inId]);
    auto sentinel = Edge(numVertices(), numVertices());
    Edge min_edge = min(out_edge, reverse(in_edge));
    while (min_edge != sentinel) {
      // Gather all incoming and outgoing edges around vertex "v0".
      VertexId v0 = min_edge[0];
      for (; min_edge[0] == v0; min_edge = min(out_edge, reverse(in_edge))) {
        VertexId v1 = min_edge[1];
        // Count the number of copies of "min_edge" in each direction.
        int out_begin = outId;
        int in_begin = inId;
        while (out_edge == min_edge) {
          out_edge = (++outId == numEdges()) ? sentinel : Edge(edge(outId));
        }
        while (reverse(in_edge) == min_edge) {
          in_edge = (++inId == numEdges()) ? sentinel : Edge(edge(in_edge_ids[inId]));
        }
        if (v0 != v1) {
          addVertexEdges(out_begin, outId, in_begin, inId, v1, v0_edges);
        } else {
          // Each degenerate edge becomes its own loop.
          for (; in_begin < inId; ++in_begin) {
            left_turn_map[in_begin] = in_begin;
          }
        }
      }
      if (v0_edges.empty()) continue;

      // Sort the edges in clockwise order around "v0".
      VertexId min_endpoint = v0_edges.front().endpoint;
      sort!((in VertexEdge a, in VertexEdge b) {
            import s2.s2predicates : orderedCCW;
            if (a.endpoint == b.endpoint) return a.rank < b.rank;
            if (a.endpoint == min_endpoint) return true;
            if (b.endpoint == min_endpoint) return false;
            return !orderedCCW(
                vertex(a.endpoint), vertex(b.endpoint), vertex(min_endpoint), vertex(v0));
          })(v0_edges);
      // Match incoming with outgoing edges.  We do this by keeping a stack of
      // unmatched incoming edges.  We also keep a stack of outgoing edges with
      // no previous incoming edge, and match these at the end by wrapping
      // around circularly to the start of the edge ordering.
      foreach (const VertexEdge e; v0_edges) {
        if (e.incoming) {
          e_in ~= in_edge_ids[e.index];
        } else if (!e_in.empty()) {
          left_turn_map[e_in.back()] = e.index;
          e_in.popBack();
        } else {
          e_out ~= e.index;  // Matched below.
        }
      }
      // Pair up additional edges using the fact that the ordering is circular.
      .reverse(e_out);
      for (; !e_out.empty() && !e_in.empty(); e_out.popBack(), e_in.popBack()) {
        left_turn_map[e_in.back()] = e_out.back();
      }
      // We only need to process unmatched incoming edges, since we are only
      // responsible for creating left turn map entries for those edges.
      if (!e_in.empty() && error.ok()) {
        error.initialize(S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS,
            "Given edges do not form loops (indegree != outdegree)");
      }
      e_in.length = 0;
      e_out.length = 0;
      v0_edges.length = 0;
    }
    return error.ok();
  }

  /**
   * Rotates the edges of "loop" if necessary so that the edge(s) with the
   * largest input edge ids are last.  This ensures that when an output loop
   * is equivalent to an input loop, their cyclic edge orders are the same.
   * "min_input_ids" is the output of GetMinInputEdgeIds().
   */
  static void canonicalizeLoopOrder(in InputEdgeId[] min_input_ids, ref EdgeId[] loop) {
    if (loop.empty()) return;
    // Find the position of the element with the highest input edge id.  If
    // there are multiple such elements together (i.e., the edge was split
    // into several pieces by snapping it to several vertices), then we choose
    // the last such position in cyclic order (this attempts to preserve the
    // original loop order even when new vertices are added).  For example, if
    // the input edge id sequence is (7, 7, 4, 5, 6, 7) then we would rotate
    // it to obtain (4, 5, 6, 7, 7, 7).

    // The reason that we put the highest-numbered edge last, rather than the
    // lowest-numbered edge first, is that S2Loop::Invert() reverses the loop
    // edge order *except* for the last edge.  For example, the loop ABCD (with
    // edges AB, BC, CD, DA) becomes DCBA (with edges DC, CB, BA, AD).  Note
    // that the last edge is the same except for its direction (DA vs. AD).
    // This has the advantage that if an undirected loop is assembled with the
    // wrong orientation and later inverted (e.g. by S2Polygon::InitOriented),
    // we still end up preserving the original cyclic vertex order.
    int pos = 0;
    bool saw_gap = false;
    for (int i = 1; i < loop.length; ++i) {
      int cmp = min_input_ids[loop[i]] - min_input_ids[loop[pos]];
      if (cmp < 0) {
        saw_gap = true;
      } else if (cmp > 0 || !saw_gap) {
        pos = i;
        saw_gap = false;
      }
    }
    if (++pos == loop.length) pos = 0;  // Convert loop end to loop start.
    bringToFront(loop[0 .. pos], loop[pos .. $]);
  }

  /**
   * Sorts the given edge chains (i.e., loops or polylines) by the minimum
   * input edge id of each chains's first edge.  This ensures that when the
   * output consists of multiple loops or polylines, they are sorted in the
   * same order as they were provided in the input.
   */
  static void canonicalizeVectorOrder(in InputEdgeId[] min_input_ids, ref EdgeId[][] chains) {
    sort!((in EdgeId[] a, in EdgeId[] b) {
          return min_input_ids[a[0]] < min_input_ids[b[0]];
        })(chains);
  }

  /// A loop consisting of a sequence of edges.
  alias EdgeLoop = EdgeId[];

  /**
   * Indicates whether loops should be simple cycles (no repeated vertices) or
   * circuits (which allow repeated vertices but not repeated edges).  In
   * terms of how the loops are built, this corresponds to closing off a loop
   * at the first repeated vertex vs. the first repeated edge.
   */
  enum LoopType { SIMPLE, CIRCUIT }

  /**
   * Builds loops from a set of directed edges, turning left at each vertex
   * until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
   * (for LoopType::CIRCUIT) is found.  (Use LoopType::SIMPLE if you intend to
   * construct an S2Loop.)
   *
   * Each loop is represented as a sequence of edges.  The edge ordering and
   * loop ordering are automatically canonicalized in order to preserve the
   * input ordering as much as possible.  Loops are non-crossing provided that
   * the graph contains no crossing edges.  If some edges cannot be turned
   * into loops, returns false and sets "error" appropriately.
   *
   * If any degenerate edges are present, then each such edge is treated as a
   * separate loop.  This is mainly useful in conjunction with
   * options.degenerate_edges() == DISCARD_EXCESS, in order to build polygons
   * that preserve degenerate geometry.
   *
   * REQUIRES: options.degenerate_edges() == {DISCARD, DISCARD_EXCESS}
   * REQUIRES: options.edge_type() == DIRECTED
   */
  bool getDirectedLoops(LoopType loop_type, ref EdgeLoop[] loops, ref S2Error error) const
  in {
    assert(_options.degenerateEdges() == DegenerateEdges.DISCARD
        || _options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);
    assert(_options.edgeType() == EdgeType.DIRECTED);
  } do {
    EdgeId[] left_turn_map;
    if (!getLeftTurnMap(getInEdgeIds(), left_turn_map, error)) return false;
    InputEdgeId[] min_input_ids = getMinInputEdgeIds();

    // If we are breaking loops at repeated vertices, we maintain a map from
    // VertexId to its position in "path".
    int[] path_index;
    if (loop_type == LoopType.SIMPLE) {
      path_index.length = numVertices();
      fill(path_index, -1);
    }

    // Visit edges in arbitrary order, and try to build a loop from each edge.
    EdgeId[] path;
    for (EdgeId start = 0; start < numEdges(); ++start) {
      if (left_turn_map[start] < 0) continue;

      // Build a loop by making left turns at each vertex until we return to
      // "start".  We use "left_turn_map" to keep track of which edges have
      // already been visited by setting its entries to -1 as we go along.  If
      // we are building vertex cycles, then whenever we encounter a vertex that
      // is already part of the path, we "peel off" a loop by removing those
      // edges from the path so far.
      for (EdgeId e = start, next; left_turn_map[e] >= 0; e = next) {
        path ~= e;
        next = left_turn_map[e];
        left_turn_map[e] = -1;
        if (loop_type == LoopType.SIMPLE) {
          path_index[edge(e)[0]] = cast(int) path.length - 1;
          int loop_start = path_index[edge(e)[1]];
          if (loop_start < 0) continue;
          // Peel off a loop from the path.
          EdgeId[] loop = path[loop_start .. $];
          path.length = loop_start;
          foreach (EdgeId e2; loop) path_index[edge(e2)[0]] = -1;
          canonicalizeLoopOrder(min_input_ids, loop);
          loops ~= loop;
        }
      }
      if (loop_type == LoopType.SIMPLE) {
        debug enforce(path.empty());  // Invariant.
      } else {
        canonicalizeLoopOrder(min_input_ids, path);
        loops ~= path;
        path.length = 0;
      }
    }
    canonicalizeVectorOrder(min_input_ids, loops);
    return true;
  }

  /**
   * Builds loops from a set of directed edges, turning left at each vertex
   * until a repeated edge is found (i.e., LoopType::CIRCUIT).  The loops are
   * further grouped into connected components, where each component consists
   * of one or more loops connected by shared vertices.
   *
   * This method is used to build polygon meshes from directed or undirected
   * input edges.  To convert the output of this method into a mesh, the
   * client must determine how the loops in different components are related
   * to each other: for example, several loops from different components may
   * bound the same region on the sphere, in which case all of those loops are
   * combined into a single polygon.  (See s2shapeutil::BuildPolygonBoundaries
   * and s2builderutil::LaxPolygonVectorLayer for details.)
   *
   * Note that loops may include both edges of a sibling pair.  When several
   * such edges are connected in a chain or a spanning tree, they form a
   * zero-area "filament".  The entire loop may be a filament (i.e., a
   * degenerate loop with an empty interior), or the loop may have have
   * non-empty interior with several filaments that extend inside it, or the
   * loop may consist of several "holes" connected by filaments.  These
   * filaments do not change the interior of any loop, so if you are only
   * interested in point containment then they can safely be removed by
   * setting the "degenerate_boundaries" parameter to DISCARD.  (They can't be
   * removed by setting (options.sibling_pairs() == DISCARD) because the two
   * siblings might belong to different polygons of the mesh.)  Note that you
   * can prevent multiple copies of sibling pairs by specifying
   * options.duplicate_edges() == MERGE.
   *
   * Each loop is represented as a sequence of edges.  The edge ordering and
   * loop ordering are automatically canonicalized in order to preserve the
   * input ordering as much as possible.  Loops are non-crossing provided that
   * the graph contains no crossing edges.  If some edges cannot be turned
   * into loops, returns false and sets "error" appropriately.
   *
   * REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
   *           (but requires DISCARD if degenerate_boundaries == DISCARD)
   * REQUIRES: options.sibling_pairs() == { REQUIRE, CREATE }
   *           [i.e., every edge must have a sibling edge]
   */
  enum DegenerateBoundaries { DISCARD, KEEP }

  alias DirectedComponent = EdgeLoop[];

  bool getDirectedComponents(
      DegenerateBoundaries degenerate_boundaries,
      ref DirectedComponent[] components,
      ref S2Error error) const
  in {
    assert(_options.degenerateEdges() == DegenerateEdges.DISCARD
        || (_options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS
            && degenerate_boundaries == DegenerateBoundaries.KEEP));
    assert(_options.siblingPairs() == SiblingPairs.REQUIRE
        || _options.siblingPairs() == SiblingPairs.CREATE);
    assert(_options.edgeType() == EdgeType.DIRECTED);  // Implied by above.
  } do {
    EdgeId[] sibling_map = getInEdgeIds();
    EdgeId[] left_turn_map;
    if (!getLeftTurnMap(sibling_map, left_turn_map, error)) return false;
    makeSiblingMap(sibling_map);
    InputEdgeId[] min_input_ids = getMinInputEdgeIds();
    EdgeId[] frontier;  // Unexplored sibling edges.

    // A map from EdgeId to the position of that edge in "path".  Only needed if
    // degenerate boundaries are being discarded.
    int[] path_index;
    if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
      path_index.length = numEdges();
      path_index[] = -1;
    }
    for (EdgeId min_start = 0; min_start < numEdges(); ++min_start) {
      if (left_turn_map[min_start] < 0) continue;  // Already used.

      // Build a connected component by keeping a stack of unexplored siblings
      // of the edges used so far.
      DirectedComponent component;
      frontier ~= min_start;
      while (!frontier.empty()) {
        EdgeId start = frontier.back();
        frontier.popBack();
        if (left_turn_map[start] < 0) continue;  // Already used.

        // Build a path by making left turns at each vertex until we return to
        // "start".  Whenever we encounter an edge that is a sibling of an edge
        // that is already on the path, we "peel off" a loop consisting of any
        // edges that were between these two edges.
        EdgeId[] path;
        for (EdgeId e = start, next; left_turn_map[e] >= 0; e = next) {
          path ~= e;
          next = left_turn_map[e];
          left_turn_map[e] = -1;
          // If the sibling hasn't been visited yet, add it to the frontier.
          EdgeId sibling = sibling_map[e];
          if (left_turn_map[sibling] >= 0) {
            frontier ~= sibling;
          }
          if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
            path_index[e] = cast(int) path.length - 1;
            int sibling_index = path_index[sibling];
            if (sibling_index < 0) continue;

            // Common special case: the edge and its sibling are adjacent, in
            // which case we can simply remove them from the path and continue.
            if (sibling_index == path.length - 2) {
              path.length = sibling_index;
              // We don't need to update "path_index" for these two edges
              // because both edges of the sibling pair have now been used.
              continue;
            }
            // Peel off a loop from the path.
            // TODO: Resume here.
            EdgeId[] loop = path[sibling_index + 1 .. $ - 1];
            path.length = sibling_index;
            // Mark the edges that are no longer part of the path.
            foreach (EdgeId e2; loop) path_index[e2] = -1;
            canonicalizeLoopOrder(min_input_ids, loop);
            component ~= loop;
          }
        }
        // Mark the edges that are no longer part of the path.
        if (degenerate_boundaries == DegenerateBoundaries.DISCARD) {
          foreach (EdgeId e2; path) path_index[e2] = -1;
        }
        canonicalizeLoopOrder(min_input_ids, path);
        component ~= path;
      }
      canonicalizeVectorOrder(min_input_ids, component);
      components ~= component;
    }
    // Sort the components to correspond to the input edge ordering.
    sort!((in DirectedComponent a, in DirectedComponent b) {
          return min_input_ids[a[0][0]] < min_input_ids[b[0][0]];
        })(components);
    return true;
  }

  alias UndirectedComponent = EdgeLoop[][2];

  /**
   * Builds loops from a set of undirected edges, turning left at each vertex
   * until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
   * (for LoopType::CIRCUIT) is found.  The loops are further grouped into
   * "components" such that all the loops in a component are connected by
   * shared vertices.  Finally, the loops in each component are divided into
   * two "complements" such that every edge in one complement is the sibling
   * of an edge in the other complement.  This corresponds to the fact that
   * given any set of non-crossing undirected loops, there are exactly two
   * possible interpretations of the region that those loops represent (where
   * one possibility is the complement of the other).  This method does not
   * attempt to resolve this ambiguity, but instead returns both possibilities
   * for each connected component and lets the client choose among them.
   *
   * This method is used to build single polygons.  (Use GetDirectedComponents
   * to build polygon meshes, even when the input edges are undirected.)  To
   * convert the output of this method into a polygon, the client must choose
   * one complement from each component such that the entire set of loops is
   * oriented consistently (i.e., they define a region such that the interior
   * of the region is always on the left).  The non-chosen complements form
   * another set of loops that are also oriented consistently but represent
   * the complementary region on the sphere.  Finally, the client needs to
   * choose one of these two sets of loops based on heuristics (e.g., the area
   * of each region), since both sets of loops are equally valid
   * interpretations of the input.
   *
   * Each loop is represented as a sequence of edges.  The edge ordering and
   * loop ordering are automatically canonicalized in order to preserve the
   * input ordering as much as possible.  Loops are non-crossing provided that
   * the graph contains no crossing edges.  If some edges cannot be turned
   * into loops, returns false and sets "error" appropriately.
   *
   * REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
   * REQUIRES: options.edge_type() == UNDIRECTED
   * REQUIRES: options.siblings_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
   *           [since REQUIRE, CREATE convert the edge_type() to DIRECTED]
   */
  bool getUndirectedComponents(
      LoopType loop_type, ref UndirectedComponent[] components, ref S2Error error) const
  in {
    assert(_options.degenerateEdges() == DegenerateEdges.DISCARD
        || _options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);
    assert(_options.edgeType() == EdgeType.UNDIRECTED);
  } do {
    import std.typecons : Tuple, tuple;
    EdgeId[] sibling_map = getInEdgeIds();
    EdgeId[] left_turn_map;
    if (!getLeftTurnMap(sibling_map, left_turn_map, error)) return false;
    makeSiblingMap(sibling_map);
    InputEdgeId[] min_input_ids = getMinInputEdgeIds();

    // A stack of unexplored sibling edges.  Each sibling edge has a "slot"
    // (0 or 1) that indicates which of the two complements it belongs to.
    Tuple!(EdgeId, int)[] frontier;

    // If we are breaking loops at repeated vertices, we maintain a map from
    // VertexId to its position in "path".
    int[] path_index;
    if (loop_type == LoopType.SIMPLE) {
      path_index.length = numVertices();
      path_index[] = -1;
    }

    for (EdgeId min_start = 0; min_start < numEdges(); ++min_start) {
      if (left_turn_map[min_start] < 0) continue;  // Already used.

      // Build a connected component by keeping a stack of unexplored siblings
      // of the edges used so far.
      UndirectedComponent component;
      frontier ~= tuple(min_start, 0);
      while (!frontier.empty()) {
        EdgeId start = frontier.back()[0];
        int slot = frontier.back()[1];
        frontier.popBack();
        if (left_turn_map[start] < 0) continue;  // Already used.

        // Build a path by making left turns at each vertex until we return to
        // "start".  We use "left_turn_map" to keep track of which edges have
        // already been visited, and which complement they were assigned to, by
        // setting its entries to negative values as we go along.
        EdgeId[] path;
        for (EdgeId e = start, next; left_turn_map[e] >= 0; e = next) {
          path ~= e;
          next = left_turn_map[e];
          left_turn_map[e] = markEdgeUsed(slot);
          // If the sibling hasn't been visited yet, add it to the frontier.
          EdgeId sibling = sibling_map[e];
          if (left_turn_map[sibling] >= 0) {
            frontier ~= tuple(sibling, 1 - slot);
          } else if (left_turn_map[sibling] != markEdgeUsed(1 - slot)) {
            // Two siblings edges can only belong the same complement if the
            // given undirected edges do not form loops.
            error.initialize(S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS,
                "Given undirected edges do not form loops");
            return false;
          }
          if (loop_type == LoopType.SIMPLE) {
            // Whenever we encounter a vertex that is already part of the path,
            // we "peel off" a loop by removing those edges from the path.
            path_index[edge(e)[0]] = cast(int) path.length - 1;
            int loop_start = path_index[edge(e)[1]];
            if (loop_start < 0) continue;
            EdgeId[] loop = path[loop_start .. $];
            path.length = loop_start;
            // Mark the vertices that are no longer part of the path.
            foreach (EdgeId e2; loop) path_index[edge(e2)[0]] = -1;
            canonicalizeLoopOrder(min_input_ids, loop);
            component[slot] ~= loop;
          }
        }
        if (loop_type == LoopType.SIMPLE) {
          enforce(path.empty());  // Invariant.
        } else {
          canonicalizeLoopOrder(min_input_ids, path);
          component[slot] ~= path;
        }
      }
      canonicalizeVectorOrder(min_input_ids, component[0]);
      canonicalizeVectorOrder(min_input_ids, component[1]);
      // To save some work in S2PolygonLayer, we swap the two loop sets of the
      // component so that the loop set whose first loop most closely follows
      // the input edge ordering is first.  (If the input was a valid S2Polygon,
      // then this component will contain normalized loops.)
      if (min_input_ids[component[0][0][0]] > min_input_ids[component[1][0][0]]) {
        swap(component[0], component[1]);
      }
      components ~= component;
    }
    // Sort the components to correspond to the input edge ordering.
    sort!((in UndirectedComponent a, in UndirectedComponent b) {
          return min_input_ids[a[0][0][0]] < min_input_ids[b[0][0][0]];
        })(components);
    return true;
  }

  /// Encodes the index of one of the two complements of each component
  /// (a.k.a. the "slot", either 0 or 1) as a negative EdgeId.
  private static EdgeId markEdgeUsed(int slot) {
    return -1 - slot;
  }

  /**
   * Indicates whether polylines should be "paths" (which don't allow
   * duplicate vertices, except possibly the first and last vertex) or
   * "walks" (which allow duplicate vertices and edges).
   */
  enum PolylineType { PATH, WALK }

  /**
   * Builds polylines from a set of edges.  If "polyline_type" is PATH, then
   * only vertices of indegree and outdegree 1 (or degree 2 in the case of
   * undirected edges) will appear in the interior of polylines.  This
   * essentially generates one polyline for each edge chain in the graph.  If
   * "polyline_type" is WALK, then polylines may pass through the same vertex
   * or even the same edge multiple times (if duplicate edges are present),
   * and each polyline will be as long as possible.  This option is useful for
   * reconstructing a polyline that has been snapped to a lower resolution,
   * since snapping can cause edges to become identical.
   *
   * This method attempts to preserve the input edge ordering in order to
   * implement idempotency, even when there are repeated edges or loops.  This
   * is true whether directed or undirected edges are used.  Degenerate edges
   * are also handled appropriately.
   *
   * REQUIRES: options.sibling_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
   */
  alias EdgePolyline = EdgeId[];
  EdgePolyline[] getPolylines(PolylineType polyline_type) const
  in {
    assert(_options.siblingPairs() == SiblingPairs.DISCARD
        || _options.siblingPairs() == SiblingPairs.DISCARD_EXCESS
        || _options.siblingPairs() == SiblingPairs.KEEP);
  } do {
    auto builder = new PolylineBuilder(this);
    if (polyline_type == PolylineType.PATH) {
      return builder.buildPaths();
    } else {
      return builder.buildWalks();
    }
  }

  ////////////////////////////////////////////////////////////////////////
  //////////////// Helper Functions for Creating Graphs //////////////////

  /**
   * Given an unsorted collection of edges, transform them according to the
   * given set of GraphOptions.  This includes actions such as discarding
   * degenerate edges; merging duplicate edges; and canonicalizing sibling
   * edge pairs in several possible ways (e.g. discarding or creating them).
   * The output is suitable for passing to the Graph constructor.
   *
   * If options.edge_type() == EdgeType::UNDIRECTED, then all input edges
   * should already have been transformed into a pair of directed edges.
   *
   * "input_ids" is a vector of the same length as "edges" that indicates
   * which input edges were snapped to each edge.  This vector is also updated
   * appropriately as edges are discarded, merged, etc.
   *
   * Note that "options" may be modified by this method: in particular, the
   * edge_type() can be changed if sibling_pairs() is CREATE or REQUIRE (see
   * the description of S2Builder::GraphOptions).
   */
  static void processEdges(
      GraphOptions options, ref Edge[] edges, ref InputEdgeIdSetId[] input_ids,
      IdSetLexicon id_set_lexicon, ref S2Error error) {
    auto processor = new EdgeProcessor(options, &edges, &input_ids, id_set_lexicon);
    processor.run(error);
    // Certain values of sibling_pairs() discard half of the edges and change
    // the edge_type() to DIRECTED (see the description of GraphOptions).
    if (options.siblingPairs() == SiblingPairs.REQUIRE
        || options.siblingPairs() == SiblingPairs.CREATE) {
      options.setEdgeType(EdgeType.DIRECTED);
    }
  }

  /**
   * Given a set of vertices and edges, removes all vertices that do not have
   * any edges and returned the new, minimal set of vertices.  Also updates
   * each edge in "edges" to correspond to the new vertex numbering.  (Note
   * that this method does *not* merge duplicate vertices, it simply removes
   * vertices of degree zero.)
   *
   * The new vertex ordering is a subsequence of the original ordering,
   * therefore if the edges were lexicographically sorted before calling this
   * method then they will still be sorted after calling this method.
   *
   * The extra argument "tmp" points to temporary storage used by this method.
   * All calls to this method from a single thread can reuse the same
   * temporary storage.  It should initially point to an empty vector.  This
   * can make a big difference to efficiency when this method is called many
   * times (e.g. to extract the vertices for different layers), since the
   * incremental running time for each layer becomes O(edges.size()) rather
   * than O(vertices.size() + edges.size()).
   */
  static S2Point[] filterVertices(in S2Point[] vertices, ref Edge[] edges, ref VertexId[] tmp) {
    // Gather the vertices that are actually used.
    VertexId[] used;
    used.reserve(2 * edges.length);
    foreach (e; edges) {
      used ~= e[0];
      used ~= e[1];
    }
    // Sort the vertices and find the distinct ones.
    used = sort(used).uniq().array;

    // Build the list of new vertices, and generate a map from old vertex id to
    // new vertex id.
    VertexId[] vmap = tmp;
    vmap.length = vertices.length;
    auto new_vertices = new S2Point[used.length];
    for (int i = 0; i < used.length; ++i) {
      new_vertices[i] = vertices[used[i]];
      vmap[used[i]] = i;
    }
    // Update the edges.
    foreach (e; edges) {
      e[0] = vmap[e[0]];
      e[1] = vmap[e[1]];
    }
    return new_vertices;
  }

  // A comparison function that allows stable sorting with std::sort (which is
  // fast but not stable).  It breaks ties between equal edges by comparing
  // their edge ids.
  static bool stableLessThan(in Edge a, in Edge b, EdgeId ai, EdgeId bi) {
    // The following is simpler but the compiler (2016) doesn't optimize it as
    // well as it should:
    //   return make_pair(a, ai) < make_pair(b, bi);
    if (a[0] < b[0]) return true;
    if (b[0] < a[0]) return false;
    if (a[1] < b[1]) return true;
    if (b[1] < a[1]) return false;
    return ai < bi;  // Stable sort.
  }

private:

  static class EdgeProcessor {
  public:
    this(
        in GraphOptions options,
        Edge[]* edges,
        InputEdgeIdSetId[]* input_ids,
        IdSetLexicon id_set_lexicon) {
      _options = options;
      _edges = edges;
      _inputIds = input_ids;
      _idSetLexicon = id_set_lexicon;

      // TODO: Resume here.

      // Sort the outgoing and incoming edges in lexigraphic order.  We use a
      // stable sort to ensure that each undirected edge becomes a sibling pair,
      // even if there are multiple identical input edges.
      _outEdges = iota(0, cast(int) (*_edges).length).array;
      sort!((EdgeId a, EdgeId b) {
            return stableLessThan((*_edges)[a], (*_edges)[b], a, b);
          })(_outEdges);
      _inEdges = iota(0, cast(int) (*_edges).length).array;
      sort!((EdgeId a, EdgeId b) {
            return stableLessThan(reverse((*_edges)[a]), reverse((*_edges)[b]), a, b);
          })(_inEdges);
      _newEdges.reserve((*_edges).length);
      _newInputIds.reserve((*_edges).length);
    }

    void run(ref S2Error error) {
      int num_edges = cast(int) (*_edges).length;
      if (num_edges == 0) return;

      // Walk through the two sorted arrays performing a merge join.  For each
      // edge, gather all the duplicate copies of the edge in both directions
      // (outgoing and incoming).  Then decide what to do based on "options_" and
      // how many copies of the edge there are in each direction.
      int outId = 0;
      int inId = 0;
      Edge out_edge = (*_edges)[_outEdges[outId]];
      Edge in_edge = (*_edges)[_inEdges[inId]];
      Edge sentinel = tuple(VertexId.max, VertexId.max);
      for (;;) {
        Edge edge = min(out_edge, reverse(in_edge));
        if (edge == sentinel) break;

        int out_begin = outId;
        int in_begin = inId;
        while (out_edge == edge) {
          out_edge = (++outId == num_edges) ? sentinel : (*_edges)[_outEdges[outId]];
        }
        while (reverse(in_edge) == edge) {
          in_edge = (++inId == num_edges) ? sentinel : (*_edges)[_inEdges[inId]];
        }
        int n_out = outId - out_begin;
        int n_in = inId - in_begin;
        if (edge[0] == edge[1]) {
          debug enforce(n_in == n_out);
          if (_options.degenerateEdges() == DegenerateEdges.DISCARD) {
            continue;
          }
          if (_options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS
              && ((out_begin > 0 && (*_edges)[_outEdges[out_begin - 1]][0] == edge[0])
                  || (outId < num_edges && (*_edges)[_outEdges[outId]][0] == edge[0])
                  || (in_begin > 0 && (*_edges)[_inEdges[in_begin - 1]][1] == edge[0])
                  || (inId < num_edges && (*_edges)[_inEdges[inId]][1] == edge[0]))) {
            continue;  // There were non-degenerate incident edges, so discard.
          }
          if (_options.edgeType() == EdgeType.UNDIRECTED
              && (_options.siblingPairs() == SiblingPairs.REQUIRE
                  || _options.siblingPairs() == SiblingPairs.CREATE)) {
            // When we have undirected edges and are guaranteed to have siblings,
            // we cut the number of edges in half (see s2builder.h).
            debug enforce((n_out & 1) == 0);  // Number of edges is always even.
            addEdges(_options.duplicateEdges() == DuplicateEdges.MERGE ? 1 : (n_out / 2),
                edge, mergeInputIds(out_begin, outId));
          } else if (_options.duplicateEdges() == DuplicateEdges.MERGE) {
            addEdges(_options.edgeType() == EdgeType.UNDIRECTED ? 2 : 1,
                edge, mergeInputIds(out_begin, outId));
          } else if (_options.siblingPairs() == SiblingPairs.DISCARD ||
              _options.siblingPairs() == SiblingPairs.DISCARD_EXCESS) {
            // Any SiblingPair option that discards edges causes the labels of all
            // duplicate edges to be merged together (see s2builder.h).
            addEdges(n_out, edge, mergeInputIds(out_begin, outId));
          } else {
            copyEdges(out_begin, outId);
          }
        } else if (_options.siblingPairs() == SiblingPairs.KEEP) {
          if (n_out > 1 && _options.duplicateEdges() == DuplicateEdges.MERGE) {
            addEdge(edge, mergeInputIds(out_begin, outId));
          } else {
            copyEdges(out_begin, outId);
          }
        } else if (_options.siblingPairs() == SiblingPairs.DISCARD) {
          if (_options.edgeType() == EdgeType.DIRECTED) {
            // If n_out == n_in: balanced sibling pairs
            // If n_out < n_in:  unbalanced siblings, in the form AB, BA, BA
            // If n_out > n_in:  unbalanced siblings, in the form AB, AB, BA
            if (n_out <= n_in) continue;
            // Any option that discards edges causes the labels of all duplicate
            // edges to be merged together (see s2builder.h).
            addEdges(_options.duplicateEdges() == DuplicateEdges.MERGE
                ? 1 : (n_out - n_in), edge, mergeInputIds(out_begin, outId));
          } else {
            if ((n_out & 1) == 0) continue;
            addEdge(edge, mergeInputIds(out_begin, outId));
          }
        } else if (_options.siblingPairs() == SiblingPairs.DISCARD_EXCESS) {
          if (_options.edgeType() == EdgeType.DIRECTED) {
            // See comments above.  The only difference is that if there are
            // balanced sibling pairs, we want to keep one such pair.
            if (n_out < n_in) continue;
            addEdges(_options.duplicateEdges() == DuplicateEdges.MERGE
                ? 1 : max(1, n_out - n_in), edge, mergeInputIds(out_begin, outId));
          } else {
            addEdges((n_out & 1) ? 1 : 2, edge, mergeInputIds(out_begin, outId));
          }
        } else {
          debug enforce(_options.siblingPairs() == SiblingPairs.REQUIRE
              || _options.siblingPairs() == SiblingPairs.CREATE);
          if (error.ok() && _options.siblingPairs() == SiblingPairs.REQUIRE
              && (_options.edgeType() == EdgeType.DIRECTED
                  ? (n_out != n_in) : ((n_out & 1) != 0))) {
            error.initialize(S2Error.Code.BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
                "Expected all input edges to have siblings, but some were missing");
          }
          if (_options.duplicateEdges() == DuplicateEdges.MERGE) {
            addEdge(edge, mergeInputIds(out_begin, outId));
          } else if (_options.edgeType() == EdgeType.UNDIRECTED) {
            addEdges((n_out + 1) / 2, edge, mergeInputIds(out_begin, outId));
          } else {
            copyEdges(out_begin, outId);
            if (n_in > n_out) {
              addEdges(n_in - n_out, edge, mergeInputIds(outId, outId));
            }
          }
        }
      }
      *_edges = _newEdges;
      *_inputIds = _newInputIds;
    }

  private:
    void addEdge(Edge edge, InputEdgeIdSetId input_edge_id_set_id) {
      _newEdges ~= edge;
      _newInputIds ~= input_edge_id_set_id;
    }

    void addEdges(int num_edges, Edge edge, InputEdgeIdSetId input_edge_id_set_id) {
      for (int i = 0; i < num_edges; ++i) {
        addEdge(edge, input_edge_id_set_id);
      }
    }

    void copyEdges(int out_begin, int out_end) {
      for (int i = out_begin; i < out_end; ++i) {
        addEdge((*_edges)[_outEdges[i]], (*_inputIds)[_outEdges[i]]);
      }
    }

    InputEdgeIdSetId mergeInputIds(int out_begin, int out_end) {
      if (out_end - out_begin == 1) {
        return (*_inputIds)[_outEdges[out_begin]];
      }
      _tmpIds.length = 0;
      for (int i = out_begin; i < out_end; ++i) {
        foreach (id; _idSetLexicon.idSet((*_inputIds)[_outEdges[i]])) {
          _tmpIds ~= id;
        }
      }
      return _idSetLexicon.add(_tmpIds);
    }

    const(GraphOptions) _options;
    Edge[]* _edges;
    InputEdgeIdSetId[]* _inputIds;
    IdSetLexicon _idSetLexicon;
    EdgeId[] _outEdges;
    EdgeId[] _inEdges;

    Edge[] _newEdges;
    InputEdgeIdSetId[] _newInputIds;

    InputEdgeId[] _tmpIds;
  }

  static class PolylineBuilder {
  public:
    this(in Graph g) {
      _g = g;
      _in = new VertexInMap(g);
      _out = new VertexOutMap(g);
      _minInputIds = g.getMinInputEdgeIds();
      _directed = _g.options().edgeType() == Graph.EdgeType.DIRECTED;
      _edgesLeft = g.numEdges() / (_directed ? 1 : 2);
      _used.length = g.numEdges();
      _used[] = false;
      if (!_directed) {
        _siblingMap = _in.inEdgeIds().dup;
        g.makeSiblingMap(_siblingMap);
      }
    }

    EdgePolyline[] buildPaths() {
      // First build polylines starting at all the vertices that cannot be in the
      // polyline interior (i.e., indegree != 1 or outdegree != 1 for directed
      // edges, or degree != 2 for undirected edges).  We consider the possible
      // starting edges in input edge id order so that we preserve the input path
      // direction even when undirected edges are used.  (Undirected edges are
      // represented by sibling pairs where only the edge in the input direction
      // is labeled with an input edge id.)
      EdgePolyline[] polylines;
      EdgeId[] edges = _g.getInputEdgeOrder(_minInputIds);
      foreach (i; 0 .. edges.length) {
        EdgeId e = edges[i];
        if (!_used[e] && !isInterior(_g.edge(e)[0])) {
          polylines ~= buildPath(e);
        }
      }
      // If there are any edges left, they form non-intersecting loops.  We build
      // each loop and then canonicalize its edge order.  We consider candidate
      // starting edges in input edge id order in order to preserve the input
      // direction of undirected loops.  Even so, we still need to canonicalize
      // the edge order to ensure that when an input edge is split into an edge
      // chain, the loop does not start in the middle of such a chain.
      for (int i = 0; i < edges.length && _edgesLeft > 0; ++i) {
        EdgeId e = edges[i];
        if (_used[e]) continue;
        EdgePolyline polyline = buildPath(e);
        canonicalizeLoopOrder(_minInputIds, polyline);
        polylines ~= polyline;
      }
      debug enforce(_edgesLeft == 0);

      // Sort the polylines to correspond to the input order (if possible).
      canonicalizeVectorOrder(_minInputIds, polylines);
      return polylines;
    }

    EdgePolyline[] buildWalks() {
      // Note that some of this code is worst-case quadratic in the maximum vertex
      // degree.  This could be fixed with a few extra arrays, but it should not
      // be a problem in practice.

      // First, build polylines from all vertices where outdegree > indegree (or
      // for undirected edges, vertices whose degree is odd).  We consider the
      // possible starting edges in input edge id order, for idempotency in the
      // case where multiple input polylines share vertices or edges.
      EdgePolyline[] polylines;
      EdgeId[] edges = _g.getInputEdgeOrder(_minInputIds);
      for (int i = 0; i < edges.length; ++i) {
        EdgeId e = edges[i];
        if (_used[e]) continue;
        VertexId v = _g.edge(e)[0];
        int excess = excessDegree(v);
        if (excess <= 0) continue;
        if (v !in _excessUsed) _excessUsed[v] = 0;
        excess -= _excessUsed[v];
        if (_directed ? (excess <= 0) : (excess % 2 == 0)) continue;
        ++_excessUsed[v];
        polylines ~= buildWalk(v);
        --_excessUsed[_g.edge(polylines.back().back())[1]];
      }
      // Now all vertices have outdegree == indegree (or even degree if undirected
      // edges are being used).  Therefore all remaining edges can be assembled
      // into loops.  We first try to expand the existing polylines if possible by
      // adding loops to them.
      if (_edgesLeft > 0) {
        foreach (ref Graph.EdgePolyline polyline; polylines) {
          maximizeWalk(polyline);
        }
      }
      // Finally, if there are still unused edges then we build loops.  If the
      // input is a polyline that forms a loop, then for idempotency we need to
      // start from the edge with minimum input edge id.  If the minimal input
      // edge was split into several edges, then we start from the first edge of
      // the chain.
      for (int i = 0; i < edges.length && _edgesLeft > 0; ++i) {
        EdgeId e = edges[i];
        if (_used[e]) continue;

        // Determine whether the origin of this edge is the start of an edge
        // chain.  To do this, we test whether (outdegree - indegree == 1) for the
        // origin, considering only unused edges with the same minimum input edge
        // id.  (Undirected edges have input edge ids in one direction only.)
        VertexId v = _g.edge(e)[0];
        InputEdgeId id = _minInputIds[e];
        int excess = 0;
        for (int j = i; j < edges.length && _minInputIds[edges[j]] == id; ++j) {
          EdgeId e2 = edges[j];
          if (_used[e2]) continue;
          if (_g.edge(e2)[0] == v) ++excess;
          if (_g.edge(e2)[1] == v) --excess;
        }
        // It is also acceptable to start a polyline from any degenerate edge.
        if (excess == 1 || _g.edge(e)[1] == v) {
          Graph.EdgePolyline polyline = buildWalk(v);
          maximizeWalk(polyline);
          polylines ~= polyline;
        }
      }
      debug enforce(_edgesLeft == 0);

      // Sort the polylines to correspond to the input order (if possible).
      canonicalizeVectorOrder(_minInputIds, polylines);
      return polylines;
    }

  private:
    bool isInterior(VertexId v) {
      if (_directed) {
        return _in.degree(v) == 1 && _out.degree(v) == 1;
      } else {
        return _out.degree(v) == 2;
      }
    }

    int excessDegree(VertexId v) {
      return _directed ? _out.degree(v) - _in.degree(v) : _out.degree(v) % 2;
    }

    EdgePolyline buildPath(EdgeId e) {
      // We simply follow edges until either we reach a vertex where there is a
      // choice about which way to go (where is_interior(v) is false), or we
      // return to the starting vertex (if the polyline is actually a loop).
      EdgePolyline polyline;
      VertexId start = _g.edge(e)[0];
      for (;;) {
        polyline ~= e;
        debug enforce(!_used[e]);
        _used[e] = true;
        if (!_directed) _used[_siblingMap[e]] = true;
        --_edgesLeft;
        VertexId v = _g.edge(e)[1];
        if (!isInterior(v) || v == start) break;
        if (_directed) {
          debug enforce(_out.degree(v) == 1);
          e = _out.edgeIds(v).front();
        } else {
          debug enforce(_out.degree(v) == 2);
          foreach (EdgeId e2; _out.edgeIds(v)) if (!_used[e2]) e = e2;
        }
      }
      return polyline;
    }

    EdgePolyline buildWalk(VertexId v) {
      EdgePolyline polyline;
      for (;;) {
        // Follow the edge with the smallest input edge id.
        EdgeId best_edge = -1;
        InputEdgeId best_out_id = InputEdgeId.max;
        foreach (EdgeId e; _out.edgeIds(v)) {
          if (_used[e] || _minInputIds[e] >= best_out_id) continue;
          best_out_id = _minInputIds[e];
          best_edge = e;
        }
        if (best_edge < 0) return polyline;
        // For idempotency when there are multiple input polylines, we stop the
        // walk early if "best_edge" might be a continuation of a different
        // incoming edge.
        if (v !in _excessUsed) _excessUsed[v] = 0;
        int excess = excessDegree(v) - _excessUsed[v];
        if (_directed ? (excess < 0) : (excess % 2) == 1) {
          foreach (EdgeId e; _in.edgeIds(v)) {
            if (!_used[e] && _minInputIds[e] <= best_out_id) {
              return polyline;
            }
          }
        }
        polyline ~= best_edge;
        _used[best_edge] = true;
        if (!_directed) _used[_siblingMap[best_edge]] = true;
        --_edgesLeft;
        v = _g.edge(best_edge)[1];
      }
    }

    void maximizeWalk(ref EdgePolyline polyline) {
      // Examine all vertices of the polyline and check whether there are any
      // unused outgoing edges.  If so, then build a loop starting at that vertex
      // and insert it into the polyline.  (The walk is guaranteed to be a loop
      // because this method is only called when all vertices have equal numbers
      // of unused incoming and outgoing edges.)
      for (int i = 0; i <= polyline.length; ++i) {
        VertexId v = (i == 0 ? _g.edge(polyline[i])[0] : _g.edge(polyline[i - 1])[1]);
        foreach (EdgeId e; _out.edgeIds(v)) {
          if (!_used[e]) {
            EdgePolyline loop = buildWalk(v);
            debug enforce(_g.edge(loop.back())[1] == v);
            polyline = polyline[0 .. i] ~ loop ~ polyline[i .. $];
            debug enforce(_used[e]);  // All outgoing edges from "v" are now used.
            break;
          }
        }
      }
    }

    const(Graph) _g;
    VertexInMap _in;
    VertexOutMap _out;
    EdgeId[] _siblingMap;
    InputEdgeId[] _minInputIds;
    bool _directed;
    int _edgesLeft;
    bool[] _used;

    // A map of (outdegree(v) - indegree(v)) considering used edges only.
    int[VertexId] _excessUsed;
  }

  const(GraphOptions) _options;
  VertexId _numVertices;  // Cached to avoid division by 24.

  const(S2Point[]) _vertices;
  const(Edge[]) _edges;
  const(InputEdgeIdSetId[]) _inputEdgeIdSetIds;
  const(IdSetLexicon) _inputEdgeIdSetLexicon;
  const(LabelSetId[]) _labelSetIds;
  const(IdSetLexicon) _labelSetLexicon;
  IsFullPolygonPredicate _isFullPolygonPredicate;
}

/// A struct for sorting the incoming and outgoing edges around a vertex "v0".
private struct VertexEdge {
  bool incoming;            // Is this an incoming edge to "v0"?
  Graph.EdgeId index;       // Index of this edge in "edges_" or "in_edge_ids"
  Graph.VertexId endpoint;  // The other (not "v0") endpoint of this edge
  int rank;                 // Secondary key for edges with the same endpoint
}

/**
 * Given a set of duplicate outgoing edges (v0, v1) and a set of duplicate
 * incoming edges (v1, v0), this method assigns each edge an integer "rank" so
 * that the edges are sorted in a consistent order with respect to their
 * orderings around "v0" and "v1".  Usually there is just one edge, in which
 * case this is easy.  Sometimes there is one edge in each direction, in which
 * case the outgoing edge is always ordered before the incoming edge.
 *
 * In general, we allow any number of duplicate edges in each direction, in
 * which case outgoing edges are interleaved with incoming edges so as to
 * create as many degenerate (two-edge) loops as possible.  In order to get a
 * consistent ordering around "v0" and "v1", we move forwards through the list
 * of outgoing edges and backwards through the list of incoming edges.  If
 * there are more incoming edges, they go at the beginning of the ordering,
 * while if there are more outgoing edges then they go at the end.
 *
 * For example, suppose there are 2 edges "a,b" from "v0" to "v1", and 4 edges
 * "w,x,y,z" from "v1" to "v0".  Using lower/upper case letters to represent
 * incoming/outgoing edges, the clockwise ordering around v0 would be zyAxBw,
 * and the clockwise ordering around v1 would be WbXaYZ.  (Try making a
 * diagram with each edge as a separate arc.)
 */
private void addVertexEdges(
    Graph.EdgeId out_begin, Graph.EdgeId out_end,
    Graph.EdgeId in_begin, Graph.EdgeId in_end,
    Graph.VertexId v1, ref VertexEdge[] v0_edges) {
  int rank = 0;
  // Any extra incoming edges go at the beginning of the ordering.
  while (in_end - in_begin > out_end - out_begin) {
    v0_edges ~= VertexEdge(true, --in_end, v1, rank++);
  }
  // Next we interleave as many outgoing and incoming edges as possible.
  while (in_end > in_begin) {
    v0_edges ~= VertexEdge(false, out_begin++, v1, rank++);
    v0_edges ~= VertexEdge(true, --in_end, v1, rank++);
  }
  // Any extra outgoing edges to at the end of the ordering.
  while (out_end > out_begin) {
    v0_edges ~= VertexEdge(false, out_begin++, v1, rank++);
  }
}

