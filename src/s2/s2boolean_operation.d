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

module s2.s2boolean_operation;

import s2.builder.graph;
import s2.builder.layer;
import s2.builder.util.snap_functions : IdentitySnapFunction;
import s2.id_set_lexicon;
import s2.logger;
import s2.s1angle;
import s2.s2builder;
import s2.s2contains_point_query;
import s2.s2error;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.shapeutil.shape_edge;
import s2.shapeutil.shape_edge_id;
import s2.util.container.btree_map;
import s2.value_lexicon;

import std.algorithm;
import std.bitmanip;
import std.exception;
import std.range;
import std.stdio;
import std.typecons;

alias EdgeType = S2Builder.EdgeType;
alias SnapFunction = S2Builder.SnapFunction;
alias DegenerateEdges = GraphOptions.DegenerateEdges;
alias DuplicateEdges = GraphOptions.DuplicateEdges;
alias SiblingPairs = GraphOptions.SiblingPairs;

alias EdgeId = Graph.EdgeId;
alias VertexId = Graph.VertexId;
alias InputEdgeId = Graph.InputEdgeId;
alias InputEdgeIdSetId = Graph.InputEdgeIdSetId;

alias PolygonModel = S2BooleanOperation.PolygonModel;
alias PolylineModel = S2BooleanOperation.PolylineModel;

// A collection of special InputEdgeIds that allow the GraphEdgeClipper state
// modifications to be inserted into the list of edge crossings.
private enum InputEdgeId SET_INSIDE = -1;
private enum InputEdgeId SET_INVERT_B = -2;
private enum InputEdgeId SET_REVERSE_A = -3;

/**
 * This class implements boolean operations (intersection, union, difference,
 * and symmetric difference) for regions whose boundaries are defined by
 * geodesic edges.
 *
 * S2BooleanOperation operates on exactly two input regions at a time.  Each
 * region is represented as an S2ShapeIndex and may contain any number of
 * points, polylines, and polygons.  The region is essentially the union of
 * these objects, except that polygon interiors must be disjoint from all
 * other geometry (including other polygon interiors).  If the input geometry
 * for a region does not meet this condition, it can be normalized by
 * computing its union first.  Note that points or polylines are allowed to
 * coincide with the boundaries of polygons.
 *
 * Degeneracies are supported.  A polygon loop or polyline may consist of a
 * single edge from a vertex to itself, and polygons may contain "sibling
 * pairs" consisting of an edge and its corresponding reverse edge.  Polygons
 * must not have any duplicate edges (due to the requirement that polygon
 * interiors are disjoint), but polylines may have duplicate edges or can even
 * be self-intersecting.
 *
 * Points and polyline edges are treated as multisets: if the same point or
 * polyline edge appears multiple times in the input, it will appear multiple
 * times in the output.  For example, the union of a point with an identical
 * point consists of two points.  This feature is useful for modeling large
 * sets of points or polylines as a single region while maintaining their
 * distinct identities, even when the points or polylines intersect each
 * other.  It is also useful for reconstructing polylines that loop back on
 * themselves.  If duplicate geometry is not desired, it can be merged by
 * GraphOptions::DuplicateEdges::MERGE in the S2Builder output layer.
 *
 * Polylines are always considered to be directed.  Polyline edges between the
 * same pair of vertices are defined to intersect even if the two edges are in
 * opposite directions.  (Undirected polylines can be modeled by specifying
 * GraphOptions::EdgeType::UNDIRECTED in the S2Builder output layer.)
 *
 * The output of each operation is sent to an S2Builder::Layer provided by the
 * client.  This allows clients to build any representation of the geometry
 * they choose.  It also allows the client to do additional postprocessing of
 * the output before building data structures; for example, the client can
 * easily discard degeneracies or convert them to another data type.
 *
 * The boundaries of polygons and polylines can be modeled as open, semi-open,
 * or closed.  Polyline boundaries are controlled by the PolylineModel class,
 * whose options are as follows:
 *
 *  - In the OPEN model, polylines do not contain their first or last vertex.
 *
 *  - In the SEMI_OPEN model, polylines contain vertices except the last.
 *    Therefore if one polyline starts where another polyline stops, the two
 *    polylines do not intersect.
 *
 *  - In the CLOSED model, polylines contain all of their vertices.
 *
 * When multiple polylines are present, they are processed independently and
 * have no effect on each other.  For example, in the OPEN boundary model the
 * polyline ABC contains the vertex B, while set of polylines {AB, BC} does
 * not.  (If you want to treat the polylines as a union instead, with
 * boundaries merged according to the "mod 2" rule, this can be achieved by
 * reassembling the edges into maximal polylines using S2PolylineVectorLayer
 * with EdgeType::UNDIRECTED, DuplicateEdges::MERGE, and PolylineType::WALK.)
 *
 * Polygon boundaries are controlled by the PolygonModel class, which has the
 * following options:
 *
 *  - In the OPEN model, polygons do not contain their vertices or edges.
 *    This implies that a polyline that follows the boundary of a polygon will
 *    not intersect it.
 *
 *  - In the SEMI_OPEN model, polygon point containment is defined such that
 *    if several polygons tile the region around a vertex, then exactly one of
 *    those polygons contains that vertex.  Similarly polygons contain all of
 *    their edges, but none of their reversed edges.  This implies that a
 *    polyline and polygon edge with the same endpoints intersect if and only
 *    if they are in the same direction.  (This rule ensures that if a
 *    polyline is intersected with a polygon and its complement, the two
 *    resulting polylines do not have any edges in common.)
 *
 *  - In the CLOSED model, polygons contain all their vertices, edges, and
 *    reversed edges.  This implies that a polyline that shares an edge (in
 *    either direction) with a polygon is defined to intersect it.  Similarly,
 *    this is the only model where polygons that touch at a vertex or along an
 *    edge intersect.
 *
 * Note that PolylineModel and PolygonModel are defined as separate classes in
 * order to allow for possible future extensions.
 *
 * Operations between geometry of different dimensions are defined as follows:
 *
 *  - For UNION, the higher-dimensional shape always wins.  For example the
 *    union of a closed polygon A with a polyline B that coincides with the
 *    boundary of A consists only of the polygon A.
 *
 *  - For INTERSECTION, the lower-dimensional shape always wins.  For example,
 *    the intersection of a closed polygon A with a point B that coincides
 *    with a vertex of A consists only of the point B.
 *
 *  - For DIFFERENCE, higher-dimensional shapes are not affected by
 *    subtracting lower-dimensional shapes.  For example, subtracting a point
 *    or polyline from a polygon A yields the original polygon A.  This rule
 *    exists because in general, it is impossible to represent the output
 *    using the specified boundary model(s).  (Consider subtracting one vertex
 *    from a PolylineModel::CLOSED polyline, or subtracting one edge from a
 *    PolygonModel::CLOSED polygon.)  If you want to perform operations like
 *    this, consider representing all boundaries explicitly (topological
 *    boundaries) using OPEN boundary models.  Another option for polygons is
 *    to subtract a degenerate loop, which yields a polygon with a degenerate
 *    hole (see S2LaxPolygonShape).
 *
 * Note that in the case of Precision::EXACT operations, the above remarks
 * only apply to the output before snapping.  Snapping may cause nearby
 * distinct edges to become coincident, e.g. a polyline may become coincident
 * with a polygon boundary.  However also note that S2BooleanOperation is
 * perfectly happy to accept such geometry as input.
 *
 * Note the following differences between S2BooleanOperation and the similar
 * S2MultiBooleanOperation class:
 *
 *  - S2BooleanOperation operates on exactly two regions at a time, whereas
 *    S2MultiBooleanOperation operates on any number of regions.
 *
 *  - S2BooleanOperation is potentially much faster when the input is already
 *    represented as S2ShapeIndexes.  The algorithm is output sensitive and is
 *    often sublinear in the input size.  This can be a big advantage if, say,
 *
 *  - S2BooleanOperation supports exact predicates and the corresponding
 *    exact operations (i.e., operations that are equivalent to computing the
 *    exact result and then snap rounding it).
 *
 *  - S2MultiBooleanOperation has better error guarantees when there are many
 *    regions, since it requires only one snapping operation for any number of
 *    input regions.
 *
 * Example usage:
 *   S2ShapeIndex a, b;  // Input geometry, e.g. containing polygons.
 *   S2Polygon polygon;  // Output geometry.
 *   S2BooleanOperation::Options options;
 *   options.set_snap_function(snap_function);
 *   S2BooleanOperation op(S2BooleanOperation::OpType::INTERSECTION,
 *                         absl::make_unique<S2PolygonLayer>(&polygon),
 *                         options);
 *   S2Error error;
 *   if (!op.Build(a, b, &error)) {
 *     LOG(ERROR) << error;
 *     ...
 *   }
 *
 * If the output includes objects of different dimensions, they can be
 * assembled into different layers with code like this:
 *
 *   vector<S2Point> points;
 *   vector<unique_ptr<S2Polyline>> polylines;
 *   S2Polygon polygon;
 *   S2BooleanOperation op(
 *       S2BooleanOperation::OpType::UNION,
 *       absl::make_unique<s2builderutil::PointVectorLayer>(&points),
 *       absl::make_unique<s2builderutil::S2PolylineVectorLayer>(&polylines),
 *       absl::make_unique<S2PolygonLayer>(&polygon));
 */
class S2BooleanOperation {
public:
  /// The supported operation types.
  enum OpType {
    UNION,                // Contained by either region.
    INTERSECTION,         // Contained by both regions.
    DIFFERENCE,           // Contained by the first region but not the second.
    SYMMETRIC_DIFFERENCE  // Contained by one region but not the other.
  }

  /// Translates OpType to one of the strings above.
  // Use std.conv.to!string intead.
  //static const char* OpTypeToString(OpType op_type);

  /// Defines whether polygons are considered to contain their vertices and/or
  /// edges (see definitions above).
  enum PolygonModel { OPEN, SEMI_OPEN, CLOSED }

  /// Defines whether polylines are considered to contain their endpoints
  /// (see definitions above).
  enum PolylineModel { OPEN, SEMI_OPEN, CLOSED }

  /**
   * With Precision::EXACT, the operation is evaluated using the exact input
   * geometry.  Predicates that use this option will produce exact results;
   * for example, they can distinguish between a polyline that barely
   * intersects a polygon from one that barely misses it.  Constructive
   * operations (ones that yield new geometry, as opposed to predicates) are
   * implemented by computing the exact result and then snap rounding it
   * according to the given snap_function() (see below).  This is as close as
   * it is possible to get to the exact result while requiring that vertex
   * coordinates have type "double".
   *
   * With Precision::SNAPPED, the input regions are snapped together *before*
   * the operation is evaluated.  So for example, two polygons that overlap
   * slightly will be treated as though they share a common boundary, and
   * similarly two polygons that are slightly separated from each other will
   * be treated as though they share a common boundary.  Snapped results are
   * useful for dealing with points, since in S2 the only points that lie
   * exactly on a polyline or polygon edge are the endpoints of that edge.
   *
   * Conceptually, the difference between these two options is that with
   * Precision::SNAPPED, the inputs are snap rounded (together), whereas with
   * Precision::EXACT only the result is snap rounded.
   */
  enum Precision { EXACT, SNAPPED }

  /**
   * SourceId identifies an edge from one of the two input S2ShapeIndexes.
   * It consists of a region id (0 or 1), a shape id within that region's
   * S2ShapeIndex, and an edge id within that shape.
   */
  struct SourceId {
  public:
    this(this) {
      _regionId = 0;
      _shapeId = 0;
      _edgeId = -1;
    }

    this(int region_id, int shape_id, int edge_id) {
      _regionId = region_id;
      _shapeId = shape_id;
      _edgeId = edge_id;
    }

    this(int special_edge_id) {
      _regionId = 0;
      _shapeId = 0;
      _edgeId = special_edge_id;
    }

    int regionId() const {
      return _regionId;
    }

    int shapeId() const {
      return _shapeId;
    }

    int edgeId() const {
      return _edgeId;
    }

    bool opEquals(in SourceId other) const {
      return _regionId == other._regionId
          && _shapeId == other._shapeId
          && _edgeId == other._edgeId;
    }

    int opCmp(in SourceId other) const {
      import std.typecons;
      return tuple(_regionId, _shapeId, _edgeId)
          .opCmp(tuple(other._regionId, other._shapeId, other._edgeId));
    }

  private:
    mixin(bitfields!(
            uint, "_regionId", 1,
            uint, "_shapeId", 31));
    int _edgeId;
  }

  struct Options {
  public:
    this(this) {
      _sourceIdLexicon = new ValueLexicon!SourceId();
      _snapFunction = new IdentitySnapFunction(S1Angle.zero());
      _polygonModel = PolygonModel.SEMI_OPEN;
      _polylineModel = PolylineModel.CLOSED;
    }

    /// Convenience constructor that calls set_snap_function().
    this(in S2Builder.SnapFunction snap_function) {
      _sourceIdLexicon = new ValueLexicon!SourceId();
      _snapFunction = snap_function.clone();
      _polygonModel = PolygonModel.SEMI_OPEN;
      _polylineModel = PolylineModel.CLOSED;
    }

    // Options may be assigned and copied.
    this(in Options options) {
      _sourceIdLexicon = new ValueLexicon!SourceId();
      _snapFunction = options._snapFunction.clone();
      _polygonModel = options._polygonModel;
      _polylineModel = options._polylineModel;
    }

    /**
     * Specifies the function to be used for snap rounding.
     *
     * DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
     * [This does no snapping and preserves all input vertices exactly unless
     *  there are crossing edges, in which case the snap radius is increased
     *  to the maximum intersection point error (S2::kIntersectionError.]
     */
    const(S2Builder.SnapFunction) snapFunction() const {
      return _snapFunction;
    }

    void setSnapFunction(in S2Builder.SnapFunction snap_function) {
      _snapFunction = snap_function.clone();
    }

    /**
     * Defines whether polygons are considered to contain their vertices
     * and/or edges.
     *
     * DEFAULT: PolygonModel::SEMI_OPEN
     */
    PolygonModel polygonModel() const {
      return _polygonModel;
    }

    void setPolygonModel(PolygonModel model) {
      _polygonModel = model;
    }

    /**
     * Defines whether polylines are considered to contain their vertices.
     *
     * DEFAULT: PolylineModel::CLOSED
     */
    PolylineModel polylineModel() const {
        return _polylineModel;
    }

    void setPolylineModel(PolylineModel model) {
      _polylineModel = model;
    }

    /**
     * Specifies whether the operation should use the exact input geometry
     * (Precision::EXACT), or whether the two input regions should be snapped
     * together first (Precision::SNAPPED).
     *
     * DEFAULT: Precision::EXACT
     */
    Precision precision() const {
      return _precision;
    }

    void setPrecision(Precision precision) {
      _precision = precision;
    }

    /**
     * If true, the input geometry is interpreted as representing nearby
     * geometry that has been snapped or simplified.  It then outputs a
     * conservative result based on the value of polygon_model() and
     * polyline_model().  For the most part, this only affects the handling of
     * degeneracies.
     *
     * - If the model is OPEN, the result is as open as possible.  For
     *   example, the intersection of two identical degenerate shells is empty
     *   under PolygonModel::OPEN because they could have been disjoint before
     *   snapping.  Similarly, two identical degenerate polylines have an
     *   empty intersection under PolylineModel::OPEN.
     *
     * - If the model is CLOSED, the result is as closed as possible.  In the
     *   case of the DIFFERENCE operation, this is equivalent to evaluating
     *   A - B as Closure(A) - Interior(B).  For other operations, it affects
     *   only the handling of degeneracies.  For example, the union of two
     *   identical degenerate holes is empty under PolygonModel::CLOSED
     *   (i.e., the hole disappears) because the holes could have been
     *   disjoint before snapping.
     *
     * - If the model is SEMI_OPEN, the result is as degenerate as possible.
     *   New degeneracies will not be created, but all degeneracies that
     *   coincide with the opposite region's boundary are retained unless this
     *   would cause a duplicate polygon edge to be created.  This model is
     *   is very useful for working with input data that has both positive and
     *   negative degeneracies (i.e., degenerate shells and holes).
     *
     * DEFAULT: false
     */
    bool conservativeOutput() const {
      return _conservative;
    }

    void setConservativeOutput(bool conservative) {
      _conservative = conservative;
    }

    /**
     * If specified, then each output edge will be labelled with one or more
     * SourceIds indicating which input edge(s) it corresponds to.  This
     * can be useful if your input geometry has additional data that needs to
     * be propagated from the input to the output (e.g., elevations).
     *
     * You can access the labels by using an S2Builder::Layer type that
     * supports labels, such as S2PolygonLayer.  The layer outputs a
     * "label_set_lexicon" and an "label_set_id" for each edge.  You can then
     * look up the source information for each edge like this:
     *
     * for (int32 label : label_set_lexicon.id_set(label_set_id)) {
     *   const SourceId& src = source_id_lexicon.value(label);
     *   // region_id() specifies which S2ShapeIndex the edge is from (0 or 1).
     *   DoSomething(src.region_id(), src.shape_id(), src.edge_id());
     * }
     *
     * DEFAULT: nullptr
     */
    const(ValueLexicon!SourceId) sourceIdLexicon() const {
      return _sourceIdLexicon;
    }

    void setSourceIdLexicon(ValueLexicon!SourceId source_id_lexicon) {
      _sourceIdLexicon = source_id_lexicon;
    }

  private:
    S2Builder.SnapFunction _snapFunction;
    PolygonModel _polygonModel;
    PolylineModel _polylineModel;
    Precision _precision;
    bool _conservative;
    ValueLexicon!SourceId _sourceIdLexicon;
  }

  this(OpType op_type, Layer layer, Options options = Options()) {
    _opType = op_type;
    _options = options;
    _resultEmpty = null;
    _layers ~= layer;
  }

  /**
   * Specifies that the output boundary edges should be sent to three
   * different layers according to their dimension.  Points (represented by
   * degenerate edges) are sent to layer 0, polyline edges are sent to
   * layer 1, and polygon edges are sent to layer 2.
   *
   * The dimension of an edge is defined as the minimum dimension of the two
   * input edges that produced it.  For example, the intersection of two
   * crossing polyline edges is a considered to be a degenerate polyline
   * rather than a point, so it is sent to layer 1.  Clients can easily
   * reclassify such polylines as points if desired, but this rule makes it
   * easier for clients that want to process point, polyline, and polygon
   * inputs differently.
   *
   * The layers are always built in the order 0, 1, 2, and all arguments to
   * the Build() calls are guaranteed to be valid until the last call returns.
   * All Graph objects have the same set of vertices and the same lexicon
   * objects, in order to make it easier to write classes that process all the
   * edges in parallel.
   */
  this(OpType op_type, Layer[] layers, Options options = Options()) {
    _opType = op_type;
    _options = options;
    _layers = layers;
    _resultEmpty = null;
  }

  OpType opType() const {
    return _opType;
  }

  /**
   * Executes the given operation.  Returns true on success, and otherwise
   * sets "error" appropriately.  (This class does not generate any errors
   * itself, but the S2Builder::Layer might.)
   */
  bool build(S2ShapeIndex a, S2ShapeIndex b, ref S2Error error) {
    _regions[0] = a;
    _regions[1] = b;
    return new Impl(this).build(error);
  }

  /// Convenience method that returns true if the result of the given operation is empty.
  static bool isEmpty(
      OpType op_type, S2ShapeIndex a, S2ShapeIndex b, Options options = Options()) {
    bool result_empty;
    auto op = new S2BooleanOperation(op_type, &result_empty, options);
    S2Error error;
    op.build(a, b, error);
    debug enforce(error.ok());
    return result_empty;
  }

  /// Convenience method that returns true if A intersects B.
  static bool intersects(S2ShapeIndex a, S2ShapeIndex b, Options options = Options()) {
    return !isEmpty(OpType.INTERSECTION, b, a, options);
  }

  /// Convenience method that returns true if A contains B, i.e., if the
  /// difference (B - A) is empty.
  static bool contains(S2ShapeIndex a, S2ShapeIndex b, Options options = Options()) {
    return isEmpty(OpType.DIFFERENCE, b, a, options);
  }

  /**
   * Convenience method that returns true if the symmetric difference of A and
   * B is empty.  (Note that A and B may still not be identical, e.g. A may
   * contain two copies of a polyline while B contains one.)
   */
  static bool equals(S2ShapeIndex a, S2ShapeIndex b, Options options = Options()) {
    return isEmpty(OpType.SYMMETRIC_DIFFERENCE, b, a, options);
  }

private:

  /// The actual implementation.
  class Impl {
  public:
    this(S2BooleanOperation op) {
      _op = op;
      _indexCrossingsFirstRegionId = -1;
    }

    bool build(ref S2Error error) {
      error.clear();
      if (isBooleanOutput()) {
        // BuildOpType() returns true if and only if the result is empty.
        *(_op._resultEmpty) = buildOpType(_op.opType());
        return true;
      }

      // TODO(ericv): Rather than having S2Builder split the edges, it would be
      // faster to call AddVertex() in this class and have a new S2Builder
      // option that increases the edge_snap_radius_ to account for errors in
      // the intersection point (the way that split_crossing_edges does).
      auto options = new S2Builder.Options(_op._options.snapFunction());
      options.setSplitCrossingEdges(true);

      // TODO(ericv): Ideally idempotent() should be true, but existing clients
      // expect vertices closer than the full "snap_radius" to be snapped.
      options.setIdempotent(false);
      _builder = new S2Builder(options);
      _builder.startLayer(new EdgeClippingLayer(_op._layers, &_inputDimensions, &_inputCrossings));

      // Polygons with no edges are assumed to be empty.  It is the responsibility
      // of clients to fix this if desired (e.g. S2Polygon has code for this).
      //
      // TODO(ericv): Implement a predicate that can determine whether a
      // degenerate polygon is empty or full based on the input S2ShapeIndexes.
      // (It is possible to do this 100% robustly, but tricky.)
      _builder.addIsFullPolygonPredicate(&isFullPolygonNever);
      buildOpType(_op.opType());
      return _builder.build(error);
    }

  private:
    struct CrossingIterator {
    public:
      /**
       * Creates an iterator over crossing edge pairs (a, b) where "b" is an edge
       * from "b_index".  "crossings_complete" indicates that "crossings" contains
       * all edge crossings between the two regions (rather than a subset).
       */
      this(S2ShapeIndex b_index, IndexCrossings crossings, bool crossings_complete) {
        _bIndex = b_index;
        _crossings = crossings;
        _bShapeId = -1;
        _crossingsComplete = crossings_complete;
        update();
      }

      void next() {
        _crossings.popFront();
        update();
      }

      bool done(ShapeEdgeId id) const {
        return aId() != id;
      }

      /// True if all edge crossings are available (see above).
      bool crossingsComplete() const {
        return _crossingsComplete;
      }

      /// True if this crossing occurs at a point interior to both edges.
      bool isInteriorCrossing() const {
        return cast(bool) _crossings[0].isInteriorCrossing;
      }

      /// Equal to S2::VertexCrossing(a_edge, b_edge), provided that a_edge and
      /// b_edge have exactly one vertex in common and neither edge is degenerate.
      bool isVertexCrossing() const {
        return cast(bool) _crossings[0].isVertexCrossing;
      }

      /// True if a_edge crosses b_edge from left to right (for interior crossings).
      bool leftToRight() const {
        return cast(bool) _crossings[0].leftToRight;
      }

      ShapeEdgeId aId() const {
        return _crossings[0].a;
      }

      ShapeEdgeId bId() const {
        return _crossings[0].b;
      }

      inout(S2ShapeIndex) bIndex() inout {
        return _bIndex;
      }

      inout(S2Shape) bShape() inout {
        return _bShape;
      }

      int bDimension() const {
        return _bDimension;
      }

      int bShapeId() const {
        return _bShapeId;
      }

      int bEdgeId() const {
        return bId().edgeId;
      }

      S2Shape.Edge bEdge() const {
        return _bShape.edge(bEdgeId());  // Opportunity to cache this.
      }

      /// Information about the chain to which an edge belongs.
      struct ChainInfo {
        int chainId;   // chain id
        int start;     // starting edge id
        int limit;     // limit edge id
      }

      /// Returns a description of the chain to which the current B edge belongs.
      ChainInfo bChainInfo() {
        if (_bInfo.chainId < 0) {
          _bInfo.chainId = bShape().chainPosition(bEdgeId()).chainId;
          auto chain = bShape().chain(_bInfo.chainId);
          _bInfo.start = chain.start;
          _bInfo.limit = chain.start + chain.length;
        }
        return _bInfo;
      }

    private:
      /// Updates information about the B shape whenever it changes.
      void update() {
        if (aId() != SENTINEL && bId().shapeId != _bShapeId) {
          _bShapeId = bId().shapeId;
          _bShape = _bIndex.shape(_bShapeId);
          _bDimension = _bShape.dimension();
          _bInfo.chainId = -1;  // Computed on demand.
        }
      }

      S2ShapeIndex _bIndex;
      IndexCrossings _crossings;
      S2Shape _bShape;
      int _bShapeId;
      int _bDimension;
      ChainInfo _bInfo;  // Computed on demand.
      bool _crossingsComplete;
    }

    /**
     * CrossingProcessor is a helper class that processes all the edges from one
     * region that cross a specific edge of the other region.  It outputs the
     * appropriate edges to an S2Builder, and outputs other information required
     * by GraphEdgeClipper to the given vectors.
     */
    class CrossingProcessor {
    public:

      /**
       * Prepares to build output for the given polygon and polyline boundary
       * models.  Edges are emitted to "builder", while other auxiliary data is
       * appended to the given vectors.
       *
       * If a predicate is being evaluated (i.e., we do not need to construct the
       * actual result), then "builder" and the various output vectors should all
       * be nullptr.
       */
      this(
          in PolygonModel polygon_model, in PolylineModel polyline_model, S2Builder builder,
          byte[]* input_dimensions, InputEdgeCrossings* input_crossings) {
        _polygonModel = polygon_model;
        _polylineModel = polyline_model;
        _builder = builder;
        _inputDimensions = input_dimensions;
        _inputCrossings = input_crossings;
        _prevInside = false;
        _sourceIdMap = new SourceIdMap();
      }

      /**
       * Starts processing edges from the given region.  "invert_a", "invert_b",
       * and "invert_result" indicate whether region A, region B, and/or the
       * result should be inverted, which allows operations such as union and
       * difference to be implemented.  For example, union is ~(~A & ~B).
       *
       * This method should be called in pairs, once to process the edges from
       * region A and once to process the edges from region B.
       */
      void startBoundary(int a_region_id, bool invert_a, bool invert_b, bool invert_result) {
        _aRegionId = a_region_id;
        _bRegionId = 1 - a_region_id;
        _invertA = invert_a;
        _invertB = invert_b;
        _invertResult = invert_result;
        _isUnion = invert_b && invert_result;

        // Specify to GraphEdgeClipper how these edges should be clipped.
        setClippingState(SET_REVERSE_A, invert_a != invert_result);
        setClippingState(SET_INVERT_B, invert_b);
      }

      /// Starts processing edges from the given shape.
      void startShape(S2Shape a_shape) {
        _aShape = a_shape;
        _aDimension = a_shape.dimension();
      }

      /// Starts processing edges from the given chain.
      void startChain(int chain_id, S2Shape.Chain chain, bool inside) {
        _chainId = chain_id;
        _chainStart = chain.start;
        _chainLimit = chain.start + chain.length;
        _inside = inside;
        _v0EmittedMaxEdgeId = chain.start - 1;  // No edges emitted yet.
        _chainV0Emitted = false;
      }

      /**
       * Processes the given edge "a_id".  "it" should be positioned to the set of
       * edges from the other region that cross "a_id" (if any).
       *
       * Supports "early exit" in the case of boolean results by returning false
       * as soon as the result is known to be non-empty.
       */
      bool processEdge(ShapeEdgeId a_id, ref CrossingIterator it) {
        // chain_edge() is faster than edge() when there are multiple chains.
        auto a = _aShape.chainEdge(_chainId, a_id.edgeId - _chainStart);
        if (_aDimension == 0) {
          return processEdge0(a_id, a, it);
        } else if (_aDimension == 1) {
          return processEdge1(a_id, a, it);
        } else {
          debug enforce(_aDimension == 2);
          return processEdge2(a_id, a, it);
        }
      }

      /**
       * This method should be called after each pair of calls to StartBoundary.
       * (The only operation that processes more than one pair of boundaries is
       * SYMMETRIC_DIFFERENCE, which computes the union of A-B and B-A.)
       *
       * Resets the state of the CrossingProcessor.
       *
       * Translates the temporary representation of crossing edges (SourceId) into
       * the format expected by EdgeClippingLayer (InputEdgeId).
       */
      void doneBoundaryPair() {
        // Add entries that translate the "special" crossings.
        _sourceIdMap[SourceId(SET_INSIDE)] = SET_INSIDE;
        _sourceIdMap[SourceId(SET_INVERT_B)] = SET_INVERT_B;
        _sourceIdMap[SourceId(SET_REVERSE_A)] = SET_REVERSE_A;
        (*_inputCrossings).reserve(_inputCrossings.length + _sourceEdgeCrossings.length);
        foreach (tmp; _sourceEdgeCrossings) {
          auto eqRange = _sourceIdMap.equalRange(tmp[1][0]);
          debug enforce(!eqRange.empty());
          *_inputCrossings ~= tuple(tmp[0], CrossingInputEdge(eqRange.front.value, tmp[1][1]));
        }
        _sourceEdgeCrossings.length = 0;
        _sourceIdMap.clear();
      }

      /**
       * Indicates whether the point being processed along the current edge chain
       * is in the polygonal interior of the opposite region, using semi-open
       * boundaries.  If "invert_b_" is true then this field is inverted.
       *
       * This value along with the set of incident edges can be used to compute
       * whether the opposite region contains this point under any of the
       * supported boundary models (PolylineModel::CLOSED, etc).
       */
      bool inside() const {
        return _inside;
      }

    private:
      // SourceEdgeCrossing represents an input edge that crosses some other
      // edge; it crosses the edge from left to right iff the second parameter
      // is "true".
      alias SourceEdgeCrossing = Tuple!(SourceId, bool);

      /**
       * PointCrossingResult describes the relationship between a point from region A
       * and a set of crossing edges from region B.  For example, "matches_polygon"
       * indicates whether a polygon vertex from region B matches the given point.
       */
      struct PointCrossingResult {
        // Note that "matches_polyline" is true only if the point matches a polyline
        // vertex of B *and* the polyline contains that vertex, whereas
        // "matches_polygon" is true if the point matches any polygon vertex.
        bool matchesPoint;     // Matches point.
        bool matchesPolyline;  // Matches contained polyline vertex.
        bool matchesPolygon;   // Matches polygon vertex.
      }

      /**
       * EdgeCrossingResult describes the relationship between an edge from region A
       * ("a_edge") and a set of crossing edges from region B.  For example,
       * "matches_polygon" indicates whether "a_edge" matches a polygon edge from
       * region B.
       */
      struct EdgeCrossingResult {
        // These fields indicate that "a_edge" exactly matches an edge of B.
        bool matchesPolyline;     // Matches polyline edge (either direction).
        bool matchesPolygon;      // Matches polygon edge (same direction).
        bool matchesSibling;      // Matches polygon edge (reverse direction).

        // These fields indicate that a vertex of "a_edge" matches a polyline vertex
        // of B *and* the polyline contains that vertex.
        bool a0MatchesPolyline;  // Start vertex matches contained polyline vertex.
        bool a1MatchesPolyline;  // End vertex matches contained polyline vertex.

        // These fields indicate that a vertex of "a_edge" matches a polygon vertex
        // of B.  (Unlike with polylines, the polygon may not contain that vertex.)
        bool a0MatchesPolygon;   // Start vertex matches polygon vertex.
        bool a1MatchesPolygon;   // End vertex matches polygon vertex.

        // These fields count the number of edge crossings at the start vertex, end
        // vertex, and interior of "a_edge".
        int a0Crossings;          // Count of polygon crossings at start vertex.
        int a1Crossings;          // Count of polygon crossings at end vertex.
        int interiorCrossings;    // Count of polygon crossings in edge interior.
      }

      InputEdgeId inputEdgeId() const {
        return cast(int) _inputDimensions.length;
      }

      /**
       * Returns true if the edges on either side of the first vertex of the
       * current edge have not been emitted.
       *
       * REQUIRES: This method is called just after updating "inside_" for "v0".
       */
      bool isV0Isolated(ShapeEdgeId a_id) const {
        return !_inside && _v0EmittedMaxEdgeId < a_id.edgeId;
      }

      /**
       * Returns true if "a_id" is the last edge of the current chain, and the
       * edges on either side of the last vertex have not been emitted (including
       * the possibility that the chain forms a loop).
       */
      bool isChainLastVertexIsolated(ShapeEdgeId a_id) const {
        return (a_id.edgeId == _chainLimit - 1 && !_chainV0Emitted
            && _v0EmittedMaxEdgeId <= a_id.edgeId);
      }

      /// Returns true if the given polyline edge contains "v0", taking into
      /// account the specified PolylineModel.
      bool polylineContainsV0(int edge_id, int chain_start) const {
        return (_polylineModel != PolylineModel.OPEN || edge_id > chain_start);
      }

      void addCrossing(SourceEdgeCrossing crossing) {
        _sourceEdgeCrossings ~= tuple(inputEdgeId(), crossing);
      }

      void setClippingState(InputEdgeId parameter, bool state) {
        addCrossing(SourceEdgeCrossing(SourceId(parameter), state));
      }

      /// Supports "early exit" in the case of boolean results by returning false
      /// as soon as the result is known to be non-empty.
      bool addEdge(ShapeEdgeId a_id, in S2Shape.Edge a, int dimension, int interior_crossings) {
        if (_builder is null) return false;  // Boolean output.
        if (interior_crossings > 0) {
          // Build a map that translates temporary edge ids (SourceId) to
          // the representation used by EdgeClippingLayer (InputEdgeId).
          auto src_id = SourceId(_aRegionId, a_id.shapeId, a_id.edgeId);
          _sourceIdMap[src_id] = inputEdgeId();
        }
        // Set the GraphEdgeClipper's "inside" state to match ours.
        if (_inside != _prevInside) setClippingState(SET_INSIDE, _inside);
        *_inputDimensions ~= cast(byte) dimension;
        _builder.addEdge(a.v0, a.v1);
        _inside ^= (interior_crossings & 1);
        _prevInside = _inside;
        return true;
      }

      /// Supports "early exit" in the case of boolean results by returning false
      /// as soon as the result is known to be non-empty.
      bool addPointEdge(in S2Point p, int dimension) {
        if (_builder is null) return false;  // Boolean output.
        if (!_prevInside) setClippingState(SET_INSIDE, true);
        *_inputDimensions ~= cast(byte) dimension;
        _builder.addEdge(p, p);
        _prevInside = true;
        return true;
      }

      /**
       * Processes an edge of dimension 0 (i.e., a point) from region A.
       *
       * Supports "early exit" in the case of boolean results by returning false
       * as soon as the result is known to be non-empty.
       */
      bool processEdge0(ShapeEdgeId a_id, in S2Shape.Edge a, ref CrossingIterator it)
      in {
        assert(a.v0 == a.v1);
      } do {
        // When a region is inverted, all points and polylines are discarded.
        if (_invertA != _invertResult) {
          skipCrossings(a_id, it);
          return true;
        }
        PointCrossingResult r = processPointCrossings(a_id, a.v0, it);

        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region, using semi-open boundaries.
        bool contained = _inside ^ _invertB;
        if (r.matchesPolygon && _polygonModel != PolygonModel.SEMI_OPEN) {
          contained = (_polygonModel == PolygonModel.CLOSED);
        }
        if (r.matchesPolyline) contained = true;

        // The output of UNION includes duplicate values, so ensure that points are
        // not suppressed by other points.
        if (r.matchesPoint && !_isUnion) contained = true;

        // Test whether the point is contained after region B is inverted.
        if (contained == _invertB) return true;  // Don't exit early.
        return addPointEdge(a.v0, 0);
      }

      /**
       * Processes an edge of dimension 1 (i.e., a polyline edge) from region A.
       *
       * Supports "early exit" in the case of boolean results by returning false
       * as soon as the result is known to be non-empty.
       */
      bool processEdge1(ShapeEdgeId a_id, in S2Shape.Edge a, ref CrossingIterator it) {
        // When a region is inverted, all points and polylines are discarded.
        if (_invertA != _invertResult) {
          skipCrossings(a_id, it);
          return true;
        }
        // Evaluate whether the start vertex should belong to the output, in case it
        // needs to be emitted as an isolated vertex.
        EdgeCrossingResult r = processEdgeCrossings(a_id, a, it);
        bool a0_inside = isPolylineVertexInside(r.a0MatchesPolyline, r.a0MatchesPolygon);

        // Test whether the entire polyline edge should be emitted (or not emitted)
        // because it matches a polyline or polygon edge.
        _inside ^= (r.a0Crossings & 1);
        if (_inside != isPolylineEdgeInside(r)) {
          _inside ^= true;   // Invert the inside_ state.
          ++r.a1Crossings;  // Restore the correct (semi-open) state later.
        }

        // If neither edge adjacent to v0 was emitted, and this polyline contains
        // v0, and the other region contains v0, then emit an isolated vertex.
        if (isV0Isolated(a_id) && polylineContainsV0(a_id.edgeId, _chainStart) && a0_inside) {
          if (!addPointEdge(a.v0, 1)) return false;
        }

        // Test whether the entire edge or any part of it belongs to the output.
        if (_inside || r.interiorCrossings > 0) {
          // Note: updates "inside_" to correspond to the state just before a1.
          if (!addEdge(a_id, a, 1 /*dimension*/, r.interiorCrossings)) {
            return false;
          }
        }
        // Remember whether the edge portion just before "a1" was emitted, so that
        // we can decide whether "a1" need to be emitted as an isolated vertex.
        if (_inside) _v0EmittedMaxEdgeId = a_id.edgeId + 1;

        // Verify that edge crossings are being counted correctly.
        _inside ^= (r.a1Crossings & 1);
        if (it.crossingsComplete()) {
          debug enforce(
              makeS2ContainsPointQuery(it.bIndex()).contains(a.v1) == (_inside ^ _invertB));
        }

        // Special case to test whether the last vertex of a polyline should be
        // emitted as an isolated vertex.
        if (_polylineModel == PolylineModel.CLOSED && it.crossingsComplete()
            && isChainLastVertexIsolated(a_id)
            && isPolylineVertexInside(r.a1MatchesPolyline, r.a1MatchesPolygon)) {
          if (!addPointEdge(a.v1, 1)) return false;
        }
        return true;
      }

      /**
       * Processes an edge of dimension 2 (i.e., a polygon edge) from region A.
       *
       * Supports "early exit" in the case of boolean results by returning false
       * as soon as the result is known to be non-empty.
       */
      bool processEdge2(ShapeEdgeId a_id, in S2Shape.Edge a, ref CrossingIterator it) {

        // In order to keep only one copy of any shared polygon edges, we only
        // output shared edges when processing the second region.
        bool emit_shared = (_aRegionId == 1);

        // Degeneracies such as isolated vertices and sibling pairs can only be
        // created by intersecting CLOSED polygons or unioning OPEN polygons.
        bool emit_degenerate =
            (_polygonModel == PolygonModel.CLOSED && !_invertA && !_invertB)
            || (_polygonModel == PolygonModel.OPEN && _invertA && _invertB);

        EdgeCrossingResult r = processEdgeCrossings(a_id, a, it);
        debug enforce(!r.matchesPolyline);
        _inside ^= (r.a0Crossings & 1);

        // If only one region is inverted, matching/sibling relations are reversed.
        // TODO(ericv): Update the following code to handle degenerate loops.
        debug enforce(!r.matchesPolygon || !r.matchesSibling);
        if (_invertA != _invertB) swap(r.matchesPolygon, r.matchesSibling);

        // Test whether the entire polygon edge should be emitted (or not emitted)
        // because it matches a polygon edge or its sibling.
        bool new_inside = _inside;

        // Shared edge are emitted only while processing the second region.
        if (r.matchesPolygon) new_inside = emit_shared;

        // Sibling pairs are emitted only when degeneracies are desired.
        if (r.matchesSibling) new_inside = emit_degenerate;
        if (_inside != new_inside) {
          _inside ^= true;   // Invert the inside_ state.
          ++r.a1Crossings;  // Restore the correct (semi-open) state later.
        }

        // Test whether the first vertex of this edge should be emitted as an
        // isolated degenerate vertex.
        if (a_id.edgeId == _chainStart) {
          _chainV0Emitted = _inside;
        } else if (emit_shared && emit_degenerate && r.a0MatchesPolygon && isV0Isolated(a_id)) {
          if (!addPointEdge(a.v0, 2)) return false;
        }

        // Test whether the entire edge or any part of it belongs to the output.
        if (_inside || r.interiorCrossings > 0) {
          // Note: updates "inside_" to correspond to the state just before a1.
          if (!addEdge(a_id, a, 2 /*dimension*/, r.interiorCrossings)) {
            return false;
          }
        }

        // Remember whether the edge portion just before "a1" was emitted, so that
        // we can decide whether "a1" need to be emitted as an isolated vertex.
        if (_inside) _v0EmittedMaxEdgeId = a_id.edgeId + 1;
        _inside ^= (r.a1Crossings & 1);

        // Verify that edge crossings are being counted correctly.
        if (it.crossingsComplete()) {
          debug enforce(
              makeS2ContainsPointQuery(it.bIndex()).contains(a.v1) == (_inside ^ _invertB));
        }

        // Special case to test whether the last vertex of a loop should be emitted
        // as an isolated degenerate vertex.
        if (emit_shared && emit_degenerate && r.a1MatchesPolygon
            && isChainLastVertexIsolated(a_id)) {
          if (!addPointEdge(a.v1, 2)) return false;
        }
        return true;
      }

      /// Skip any crossings that were not needed to determine the result.
      void skipCrossings(ShapeEdgeId a_id, ref CrossingIterator it) {
        while (!it.done(a_id)) it.next();
      }

      /// Returns a summary of the relationship between a point from region A and
      /// a set of crossing edges from region B (see PointCrossingResult).
      PointCrossingResult processPointCrossings(
          ShapeEdgeId a_id, in S2Point a0, ref CrossingIterator it) const {
        PointCrossingResult r;
        for (; !it.done(a_id); it.next()) {
          if (it.bDimension() == 0) {
            r.matchesPoint = true;
          } else if (it.bDimension() == 1) {
            if (polylineEdgeContainsVertex(a0, it)) {
              r.matchesPolyline = true;
            }
          } else {
            r.matchesPolygon = true;
          }
        }
        return r;
      }

      /**
       * Returns a summary of the relationship between a test edge from region A and
       * a set of crossing edges from region B (see EdgeCrossingResult).
       *
       * NOTE(ericv): We could save a bit of work when matching polygon vertices by
       * passing in a flag saying whether this information is needed.  For example
       * if is only needed in ProcessEdge2 when (emit_shared && emit_degenerate).
       */
      EdgeCrossingResult processEdgeCrossings(
          ShapeEdgeId a_id, in S2Shape.Edge a, ref CrossingIterator it) {
        EdgeCrossingResult r;
        if (it.done(a_id)) return r;

        // TODO(ericv): bool a_degenerate = (a.v0 == a.v1);
        for (; !it.done(a_id); it.next()) {
          // Polylines and polygons are not affected by point geometry.
          if (it.bDimension() == 0) continue;
          S2Shape.Edge b = it.bEdge();
          if (it.isInteriorCrossing()) {
            // The crossing occurs in the edge interior.  The condition below says
            // that (1) polyline crossings don't affect polygon output, and (2)
            // subtracting a crossing polyline from a polyline has no effect.
            if (_aDimension <= it.bDimension()
                && !(_invertB != _invertResult && it.bDimension() == 1)) {
              auto src_id = SourceId(_bRegionId, it.bShapeId(), it.bEdgeId());
              addCrossing(tuple(src_id, it.leftToRight()));
            }
            r.interiorCrossings += (it.bDimension() == 1) ? 2 : 1;
          } else if (it.bDimension() == 1) {
            // Polygons are not affected by polyline geometry.
            if (_aDimension == 2) continue;
            if ((a.v0 == b.v0 && a.v1 == b.v1) || (a.v0 == b.v1 && a.v1 == b.v0)) {
              r.matchesPolyline = true;
            }
            if ((a.v0 == b.v0 || a.v0 == b.v1) && polylineEdgeContainsVertex(a.v0, it)) {
              r.a0MatchesPolyline = true;
            }
            if ((a.v1 == b.v0 || a.v1 == b.v1) && polylineEdgeContainsVertex(a.v1, it)) {
              r.a1MatchesPolyline = true;
            }
          } else {
            debug enforce(it.bDimension() == 2);
            if (a.v0 == b.v0 && a.v1 == b.v1) {
              ++r.a0Crossings;
              r.matchesPolygon = true;
            } else if (a.v0 == b.v1 && a.v1 == b.v0) {
              ++r.a0Crossings;
              r.matchesSibling = true;
            } else if (it.isVertexCrossing()) {
              if (a.v0 == b.v0 || a.v0 == b.v1) {
                ++r.a0Crossings;
              } else {
                ++r.a1Crossings;
              }
            }
            if (a.v0 == b.v0 || a.v0 == b.v1) {
              r.a0MatchesPolygon = true;
            }
            if (a.v1 == b.v0 || a.v1 == b.v1) {
              r.a1MatchesPolygon = true;
            }
          }
        }
        return r;
      }

      /**
       * Returns true if the current point being processed (which must be a polyline
       * vertex) is contained by the opposite region (after inversion if "invert_b_"
       * is true).  "matches_polyline" and "matches_polygon" indicate whether the
       * vertex matches a polyline/polygon vertex of the opposite region.
       */
      bool isPolylineVertexInside(bool matches_polyline, bool matches_polygon) const {
        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region using semi-open boundaries.
        bool contained = _inside ^ _invertB;

        // For UNION the output includes duplicate polylines.  The test below
        // ensures that isolated polyline vertices are not suppressed by other
        // polyline vertices in the output.
        if (matches_polyline && !_isUnion) {
          contained = true;
        } else if (matches_polygon && _polygonModel != PolygonModel.SEMI_OPEN) {
          contained = (_polygonModel == PolygonModel.CLOSED);
        }
        // Finally, invert the result if the opposite region should be inverted.
        return contained ^ _invertB;
      }

      /**
       * Returns true if the current polyline edge is contained by the opposite
       * region (after inversion if "invert_b_" is true).
       */
      bool isPolylineEdgeInside(in EdgeCrossingResult r) const {
        // "contained" indicates whether the current point is inside the polygonal
        // interior of the opposite region using semi-open boundaries.
        bool contained = _inside ^ _invertB;
        if (r.matchesPolyline && !_isUnion) {
          contained = true;
        } else if (r.matchesPolygon) {
          // In the SEMI_OPEN model, polygon sibling pairs cancel each other and
          // have no effect on point or edge containment.
          if (!(r.matchesSibling && _polygonModel == PolygonModel.SEMI_OPEN)) {
            contained = (_polygonModel != PolygonModel.OPEN);
          }
        } else if (r.matchesSibling) {
          contained = (_polygonModel == PolygonModel.CLOSED);
        }
        // Finally, invert the result if the opposite region should be inverted.
        return contained ^ _invertB;
      }

      /**
       * Returns true if the vertex "v" is contained by the polyline edge referred
       * to by the CrossingIterator "it", taking into account the PolylineModel.
       *
       * REQUIRES: it.b_dimension() == 1
       * REQUIRES: "v" is an endpoint of it.b_edge()
       */
      bool polylineEdgeContainsVertex(in S2Point v, ref CrossingIterator it) const
      in {
        assert(it.bDimension() == 1);
        assert(it.bEdge().v0 == v || it.bEdge().v1 == v);
      } do {

        // Closed polylines contain all their vertices.
        if (_polylineModel == PolylineModel.CLOSED) return true;

        // All interior polyline vertices are contained.
        // The last polyline vertex is contained iff the model is CLOSED.
        // The first polyline vertex is contained iff the model is not OPEN.
        // The test below is written such that usually b_edge() is not needed.
        const auto b_chain = it.bChainInfo();
        int b_edge_id = it.bEdgeId();
        return (b_edge_id < b_chain.limit - 1 || it.bEdge().v1 != v)
            && (polylineContainsV0(b_edge_id, b_chain.start) || it.bEdge().v0 != v);
      }

      // Constructor parameters:

      PolygonModel _polygonModel;
      PolylineModel _polylineModel;

      // The output of the CrossingProcessor consists of a subset of the input
      // edges that are emitted to "builder_", and some auxiliary information
      // that allows GraphEdgeClipper to determine which segments of those input
      // edges belong to the output.  The auxiliary information consists of the
      // dimension of each input edge, and set of input edges from the other
      // region that cross each input input edge.
      S2Builder _builder;
      byte[]* _inputDimensions;
      InputEdgeCrossings* _inputCrossings;

      // Fields set by StartBoundary:

      int _aRegionId, _bRegionId;
      bool _invertA, _invertB, _invertResult;
      bool _isUnion;  // True if this is a UNION operation.

      // Fields set by StartShape:

      S2Shape _aShape;
      int _aDimension;

      // Fields set by StartChain:

      int _chainId;
      int _chainStart;
      int _chainLimit;

      // Fields updated by ProcessEdge:

      alias SourceEdgeCrossings = Tuple!(InputEdgeId, SourceEdgeCrossing)[];

      /**
       * A temporary representation of input_crossings_ that is used internally
       * until all necessary edges from *both* polygons have been emitted to the
       * S2Builder.  This field is then converted by DoneBoundaryPair() into
       * the InputEdgeCrossings format expected by GraphEdgeClipper.
       *
       * The reason that we can't construct input_crossings_ directly is that it
       * uses InputEdgeIds to identify the edges from both polygons, and when we
       * are processing edges from the first polygon, InputEdgeIds have not yet
       * been assigned to the second polygon.  So instead this field identifies
       * edges from the first polygon using an InputEdgeId, and edges from the
       * second polygon using a (region_id, shape_id, edge_id) tuple (i.e., a
       * SourceId).
       *
       * All crossings are represented twice, once to indicate that an edge from
       * polygon 0 is crossed by an edge from polygon 1, and once to indicate that
       * an edge from polygon 1 is crossed by an edge from polygon 0.
       */
      SourceEdgeCrossings _sourceEdgeCrossings;

      alias SourceIdMap = BTreeMap!(SourceId, InputEdgeId);

      /**
       * A map that translates from SourceId (the (region_id, shape_id,
       * edge_id) triple that identifies an S2ShapeIndex edge) to InputEdgeId (the
       * sequentially increasing numbers assigned to input edges by S2Builder).
       */
      SourceIdMap _sourceIdMap;

      /**
       * Indicates whether the point being processed along the current edge chain
       * is in the polygonal interior of the opposite region, using semi-open
       * boundaries.  If "invert_b_" is true then this field is inverted.
       *
       * Equal to: b_index_.Contains(current point) ^ invert_b_
       */
      bool _inside;

      /**
       * The value of that "inside_" would have just before the end of the
       * previous edge added to S2Builder.  This value is used to determine
       * whether the GraphEdgeClipper state needs to be updated when jumping from
       * one edge chain to another.
       */
      bool _prevInside;

      /**
       * The maximum edge id of any edge in the current chain whose v0 vertex has
       * already been emitted.  This is used to determine when an isolated vertex
       * needs to be emitted, e.g. when two closed polygons share only a vertex.
       */
      int _v0EmittedMaxEdgeId;

      /**
       * True if the first vertex of the current chain has been emitted.  This is
       * used when processing loops in order to determine whether the first/last
       * vertex of the loop should be emitted as an isolated vertex.
       */
      bool _chainV0Emitted;
    }

    /**
     * An IndexCrossing represents a pair of intersecting S2ShapeIndex edges
     * ("a_edge" and "b_edge").  We store all such intersections because the
     * algorithm needs them twice, once when processing the boundary of region A
     * and once when processing the boundary of region B.
     */
    struct IndexCrossing {
      ShapeEdgeId a, b;

      mixin(bitfields!(
              /// True if S2::CrossingSign(a_edge, b_edge) > 0.
              uint,  "isInteriorCrossing", 1,

              // True if "a_edge" crosses "b_edge" from left to right.  Undefined if
              // is_interior_crossing is false.
              uint, "leftToRight", 1,

              // Equal to S2::VertexCrossing(a_edge, b_edge).  Undefined if "a_edge" and
              // "b_edge" do not share exactly one vertex or either edge is degenerate.
              uint, "isVertexCrossing", 1,

              // Unused bits.
              uint, "", 5));

      // All flags are "false" by default.
      this(ShapeEdgeId a, ShapeEdgeId b) {
        this.a = a;
        this.b = b;
        this.isInteriorCrossing = false;
        this.leftToRight = false;
        this.isVertexCrossing = false;
      }

      bool opEquals(in IndexCrossing x) const {
        return a == x.a && b == x.b;
      }

      int opCmp(in IndexCrossing x) const {
        // The compiler (2017) doesn't optimize the following as well:
        // return x.a < y.a || (x.a == y.a && x.b < y.b);
        if (a.shapeId != x.a.shapeId) return a.shapeId - x.a.shapeId;
        if (a.edgeId != x.a.edgeId) return a.edgeId - x.a.edgeId;
        if (b.shapeId != x.b.shapeId) return b.shapeId - x.b.shapeId;
        if (b.edgeId != x.b.edgeId) return b.edgeId - x.b.edgeId;
        return 0;
      }
    }

    alias IndexCrossings = IndexCrossing[];

    bool isBooleanOutput() const {
      return _op._resultEmpty !is null;
    }

    // All of the methods below support "early exit" in the case of boolean
    // results by returning "false" as soon as the result is known to be
    // non-empty.

    /**
     * Clips the boundary of A to the interior of the opposite region B and adds
     * the resulting edges to the output.  Optionally, any combination of region
     * A, region B, and the result may be inverted, which allows operations such
     * as union and difference to be implemented.
     *
     * Note that when an input region is inverted with respect to the output
     * (e.g., invert_a != invert_result), all polygon edges are reversed and all
     * points and polylines are discarded, since the complement of such objects
     * cannot be represented.  (If you want to compute the complement of points
     * or polylines, you can use S2LaxPolygonShape to represent your geometry as
     * degenerate polygons instead.)
     *
     * This method must be called an even number of times (first to clip A to B
     * and then to clip B to A), calling DoneBoundaryPair() after each pair.
     *
     * Supports "early exit" in the case of boolean results by returning false
     * as soon as the result is known to be non-empty.
     */
    bool addBoundary(int a_region_id, bool invert_a, bool invert_b,
        bool invert_result,
        in ShapeEdgeId[] a_chain_starts,
        CrossingProcessor cp) {
      S2ShapeIndex a_index = _op._regions[a_region_id];
      S2ShapeIndex b_index = _op._regions[1 - a_region_id];
      if (!getIndexCrossings(a_region_id)) return false;
      cp.startBoundary(a_region_id, invert_a, invert_b, invert_result);

      // Walk the boundary of region A and build a list of all edge crossings.
      // We also keep track of whether the current vertex is inside region B.
      auto next_start_pos = 0; // a_chain_starts.begin();
      auto next_crossing =
          CrossingIterator(b_index, _indexCrossings, true /*crossings_complete*/);
      ShapeEdgeId next_id = min(a_chain_starts[next_start_pos], next_crossing.aId());
      while (next_id != SENTINEL) {
        int a_shape_id = next_id.shapeId;
        S2Shape a_shape = a_index.shape(a_shape_id);
        cp.startShape(a_shape);
        while (next_id.shapeId == a_shape_id) {
          // TODO(ericv): Special handling of dimension 0?  Can omit most of this
          // code, including the loop, since all chains are of length 1.
          int edge_id = next_id.edgeId;
          S2Shape.ChainPosition chain_position = a_shape.chainPosition(edge_id);
          int chain_id = chain_position.chainId;
          S2Shape.Chain chain = a_shape.chain(chain_id);
          bool start_inside = (next_id == a_chain_starts[next_start_pos]);
          if (start_inside) ++next_start_pos;
          cp.startChain(chain_id, chain, start_inside);
          int chain_limit = chain.start + chain.length;
          while (edge_id < chain_limit) {
            auto a_id = ShapeEdgeId(a_shape_id, edge_id);
            debug enforce(cp.inside() || next_crossing.aId() == a_id);
            if (!cp.processEdge(a_id, next_crossing)) {
              return false;
            }
            if (cp.inside()) {
              ++edge_id;
            } else if (next_crossing.aId().shapeId == a_shape_id
                && next_crossing.aId().edgeId < chain_limit) {
              edge_id = next_crossing.aId().edgeId;
            } else {
              break;
            }
          }
          next_id = min(a_chain_starts[next_start_pos], next_crossing.aId());
        }
      }
      return true;
    }

    /**
     * Returns the first edge of each edge chain from "a_region_id" whose first
     * vertex is contained by opposite region's polygons (using the semi-open
     * boundary model).  Each input region and the result region are inverted as
     * specified (invert_a, invert_b, and invert_result) before testing for
     * containment.  The algorithm uses these "chain starts" in order to clip the
     * boundary of A to the interior of B in an output-senstive way.
     *
     * This method supports "early exit" in the case where a boolean predicate is
     * being evaluated and the algorithm discovers that the result region will be
     * non-empty.
     */
    bool getChainStarts(int a_region_id, bool invert_a, bool invert_b,
        bool invert_result, CrossingProcessor cp,
        ref ShapeEdgeId[] chain_starts) {
      S2ShapeIndex a_index = _op._regions[a_region_id];
      S2ShapeIndex b_index = _op._regions[1 - a_region_id];

      if (isBooleanOutput()) {
        // If boolean output is requested, then we use the CrossingProcessor to
        // determine whether the first edge of each chain will be emitted to the
        // output region.  This lets us terminate the operation early in many
        // cases.
        cp.startBoundary(a_region_id, invert_a, invert_b, invert_result);
      }

      // If region B has no two-dimensional shapes and is not inverted, then by
      // definition no chain starts are contained.  However if boolean output is
      // requested then we check for containment anyway, since as a side effect we
      // may discover that the result region is non-empty and terminate the entire
      // operation early.
      bool b_has_interior = hasInterior(b_index);
      if (b_has_interior || invert_b || isBooleanOutput()) {
        auto query = makeS2ContainsPointQuery(b_index);
        int num_shape_ids = a_index.numShapeIds();
        for (int shape_id = 0; shape_id < num_shape_ids; ++shape_id) {
          S2Shape a_shape = a_index.shape(shape_id);
          if (a_shape is null) continue;

          // If region A is being subtracted from region B, points and polylines
          // in region A can be ignored since these shapes never contribute to the
          // output (they can only remove edges from region B).
          if (invert_a != invert_result && a_shape.dimension() < 2) continue;

          if (isBooleanOutput()) cp.startShape(a_shape);
          int num_chains = a_shape.numChains();
          for (int chain_id = 0; chain_id < num_chains; ++chain_id) {
            S2Shape.Chain chain = a_shape.chain(chain_id);
            if (chain.length == 0) continue;
            auto a = new ShapeEdge(shape_id, chain.start, a_shape.chainEdge(chain_id, 0));
            bool inside = (b_has_interior && query.contains(a.v0())) != invert_b;
            if (inside) {
              chain_starts ~= ShapeEdgeId(shape_id, chain.start);
            }
            if (isBooleanOutput()) {
              cp.startChain(chain_id, chain, inside);
              if (!processIncidentEdges(a, query, cp)) return false;
            }
          }
        }
      }
      chain_starts ~= SENTINEL;
      return true;
    }

    bool processIncidentEdges(
        in ShapeEdge a, S2ContainsPointQuery!S2ShapeIndex query, CrossingProcessor cp) {
      _tmpCrossings.length = 0;
      query.visitIncidentEdges(a.v0(), (in ShapeEdge b) {
            return addIndexCrossing(a, b, false /*is_interior*/, _tmpCrossings);
          });
      // Fast path for the common case where there are no incident edges.  We
      // return false (terminating early) if the first chain edge will be emitted.
      if (_tmpCrossings.empty()) {
        return !cp.inside();
      }
      // Otherwise we invoke the full CrossingProcessor logic to determine whether
      // the first chain edge will be emitted.
      if (_tmpCrossings.length > 1) {
        sort(_tmpCrossings);
        // VisitIncidentEdges() should not generate any duplicate values.
        debug enforce(findAdjacent(_tmpCrossings[]).empty());
      }
      _tmpCrossings ~= IndexCrossing(SENTINEL, SENTINEL);
      auto next_crossing =
          CrossingIterator(query.index(), _tmpCrossings, false /*crossings_complete*/);
      return cp.processEdge(a.id(), next_crossing);
    }

    static bool hasInterior(in S2ShapeIndex index) {
      for (int s = index.numShapeIds(); --s >= 0; ) {
        const(S2Shape) shape = index.shape(s);
        if (shape && shape.hasInterior()) return true;
      }
      return false;
    }

    static bool addIndexCrossing(
        in ShapeEdge a, in ShapeEdge b, bool is_interior, IndexCrossings crossings) {
      import s2.s2predicates : sign;
      import s2.s2edge_crossings : vertexCrossing;

      crossings ~= IndexCrossing(a.id(), b.id());
      IndexCrossing crossing = crossings.back();
      if (is_interior) {
        crossing.isInteriorCrossing = true;
        if (sign(a.v0(), a.v1(), b.v0()) > 0) {
          crossing.leftToRight = true;
        }
      } else {
        // TODO(ericv): This field isn't used unless one shape is a polygon and
        // the other is a polyline or polygon, but we don't have the shape
        // dimension information readily available here.
        if (vertexCrossing(a.v0(), a.v1(), b.v0(), b.v1())) {
          crossing.isVertexCrossing = true;
        }
      }
      return true;  // Continue visiting.
    }

    /**
     * Initialize index_crossings_ to the set of crossing edge pairs such that the
     * first element of each pair is an edge from "region_id".
     *
     * Supports "early exit" in the case of boolean results by returning false
     * as soon as the result is known to be non-empty.
     */
    bool getIndexCrossings(int region_id) {
      import s2.shapeutil.visit_crossing_edge_pairs;
      import s2.s2crossing_edge_query : CrossingType;

      if (region_id == _indexCrossingsFirstRegionId) return true;
      if (_indexCrossingsFirstRegionId < 0) {
        debug enforce(region_id == 0);  // For efficiency, not correctness.
        if (!visitCrossingEdgePairs(
                _op._regions[0], _op._regions[1],
                CrossingType.ALL,
                (in ShapeEdge a, in ShapeEdge b, bool is_interior) {
                  // For all supported operations (union, intersection, and
                  // difference), if the input edges have an interior crossing
                  // then the output is guaranteed to have at least one edge.
                  if (is_interior && isBooleanOutput()) return false;
                  return addIndexCrossing(a, b, is_interior, _indexCrossings);
                })) {
          return false;
        }
        if (_indexCrossings.length > 1) {
          sort(_indexCrossings);
          _indexCrossings = uniq(_indexCrossings).array;
        }
        // Add a sentinel value to simplify the loop logic.
        _indexCrossings ~= IndexCrossing(SENTINEL, SENTINEL);
        _indexCrossingsFirstRegionId = 0;
      }
      if (region_id != _indexCrossingsFirstRegionId) {
        foreach (crossing; _indexCrossings) {
          swap(crossing.a, crossing.b);
          // The following predicates get inverted when the edges are swapped.
          crossing.leftToRight(crossing.leftToRight ^ true);
          crossing.isVertexCrossing(crossing.isVertexCrossing ^ true);
        }
        sort(_indexCrossings);
        _indexCrossingsFirstRegionId = region_id;
      }
      return true;
    }

    /// Supports "early exit" in the case of boolean results by returning false
    /// as soon as the result is known to be non-empty.
    bool addBoundaryPair(bool invert_a, bool invert_b, bool invert_result, CrossingProcessor cp) {
      // Optimization: if the operation is DIFFERENCE or SYMMETRIC_DIFFERENCE,
      // it is worthwhile checking whether the two regions are identical (in which
      // case the output is empty).
      //
      // TODO(ericv): When boolean output is requested there are other quick
      // checks that could be done here, such as checking whether a full cell from
      // one S2ShapeIndex intersects a non-empty cell of the other S2ShapeIndex.
      auto type = _op.opType();
      if (type == OpType.DIFFERENCE || type == OpType.SYMMETRIC_DIFFERENCE) {
        if (areRegionsIdentical()) return true;
      } else if (!isBooleanOutput()) {
      }
      ShapeEdgeId[] a_starts, b_starts;
      if (!getChainStarts(0, invert_a, invert_b, invert_result, cp, a_starts)
          || !getChainStarts(1, invert_b, invert_a, invert_result, cp, b_starts)
          || !addBoundary(0, invert_a, invert_b, invert_result, a_starts, cp)
          || !addBoundary(1, invert_b, invert_a, invert_result, b_starts, cp)) {
        return false;
      }
      if (!isBooleanOutput()) cp.doneBoundaryPair();
      return true;
    }

    bool areRegionsIdentical() const {
      const(S2ShapeIndex) a = _op._regions[0];
      const(S2ShapeIndex) b = _op._regions[1];
      if (a == b) return true;
      int num_shape_ids = a.numShapeIds();
      if (num_shape_ids != b.numShapeIds()) return false;
      for (int s = 0; s < num_shape_ids; ++s) {
        const(S2Shape) a_shape = a.shape(s);
        const(S2Shape) b_shape = b.shape(s);
        if (a_shape.dimension() != b_shape.dimension()) return false;
        if (a_shape.dimension() == 2) {
          auto a_ref = a_shape.getReferencePoint();
          auto b_ref = b_shape.getReferencePoint();
          if (a_ref.point != b_ref.point) return false;
          if (a_ref.contained != b_ref.contained) return false;
        }
        int num_chains = a_shape.numChains();
        if (num_chains != b_shape.numChains()) return false;
        for (int c = 0; c < num_chains; ++c) {
          S2Shape.Chain a_chain = a_shape.chain(c);
          S2Shape.Chain b_chain = b_shape.chain(c);
          debug enforce(a_chain.start == b_chain.start);
          if (a_chain.length != b_chain.length) return false;
          for (int i = 0; i < a_chain.length; ++i) {
            S2Shape.Edge a_edge = a_shape.chainEdge(c, i);
            S2Shape.Edge b_edge = b_shape.chainEdge(c, i);
            if (a_edge.v0 != b_edge.v0) return false;
            if (a_edge.v1 != b_edge.v1) return false;
          }
        }
      }
      return true;
    }

    /// Supports "early exit" in the case of boolean results by returning false
    /// as soon as the result is known to be non-empty.
    bool buildOpType(OpType op_type) {
      // CrossingProcessor does the real work of emitting the output edges.
      auto cp = new CrossingProcessor(_op._options.polygonModel(),
          _op._options.polylineModel(),
          _builder, &_inputDimensions, &_inputCrossings);
      final switch (op_type) {
        case OpType.UNION:
          // A | B == ~(~A & ~B)
          return addBoundaryPair(true, true, true, cp);

        case OpType.INTERSECTION:
          // A & B
          return addBoundaryPair(false, false, false, cp);

        case OpType.DIFFERENCE:
          // A - B = A & ~B
          return addBoundaryPair(false, true, false, cp);

        case OpType.SYMMETRIC_DIFFERENCE:
          // Compute the union of (A - B) and (B - A).
          return (addBoundaryPair(false, true, false, cp)
              && addBoundaryPair(true, false, false, cp));
      }
    }

    S2BooleanOperation _op;

    /// The S2Builder used to construct the output.
    S2Builder _builder;

    /// A vector specifying the dimension of each edge added to S2Builder.
    byte[] _inputDimensions;

    /// The set of all input edge crossings, which is used by EdgeClippingLayer
    /// to construct the clipped output polygon.
    InputEdgeCrossings _inputCrossings;

    /// kSentinel is a sentinel value used to mark the end of vectors.
    enum ShapeEdgeId SENTINEL = ShapeEdgeId(int.max, 0);

    /**
     * A vector containing all pairs of crossing edges from the two input
     * regions (including edge pairs that share a common vertex).  The first
     * element of each pair is an edge from "index_crossings_first_region_id_",
     * while the second element of each pair is an edge from the other region.
     */
    IndexCrossings _indexCrossings;

    /**
     * Indicates that the first element of each crossing edge pair in
     * "index_crossings_" corresponds to an edge from the given region.
     * This field is negative if index_crossings_ has not been computed yet.
     */
    int _indexCrossingsFirstRegionId;

    /// Temporary storage used in GetChainStarts(), declared here to avoid
    /// repeatedly allocating memory.
    IndexCrossings _tmpCrossings;
  }

  // TODO: Resume here.

  /// Internal constructor to reduce code duplication.
  this(OpType op_type, Options options) {
    _opType = op_type;
    _options = options;
    _resultEmpty = null;
  }

  /**
   * Specifies that "result_empty" should be set to indicate whether the exact
   * result of the operation is empty (contains no edges).  This constructor
   * is used to efficiently test boolean relationships (see IsEmpty above).
   */
  this(OpType op_type, bool* result_empty, Options options = Options()) {
    _opType = op_type;
    _options = options;
    _resultEmpty = result_empty;
  }

  OpType _opType;
  Options _options;

  /// The input regions.
  S2ShapeIndex[2] _regions;

  /// The output consists either of zero layers, one layer, or three layers.
  Layer[] _layers;

  /// The following field is set if and only if there are no output layers.
  bool* _resultEmpty;
}

/**
 * CrossingInputEdge represents an input edge B that crosses some other input
 * edge A.  It stores the input edge id of edge B and also whether it crosses
 * edge A from left to right (or vice versa).
 */
struct CrossingInputEdge {
public:
  /// Indicates that input edge "input_id" crosses another edge (from left to
  /// right if "left_to_right" is true).
  this(InputEdgeId input_id, bool left_to_right) {
    _leftToRight = left_to_right;
    _inputId = input_id;
  }

  InputEdgeId inputId() const {
    return _inputId;
  }

  bool leftToRight() const {
    return _leftToRight;
  }

  bool opCmp(in CrossingInputEdge other) const {
    return _inputId < other._inputId;
  }

  bool opCmp(in InputEdgeId other) const {
    return _inputId < other;
  }

private:
  mixin(bitfields!(
          bool,  "_leftToRight", 1,
          InputEdgeId, "_inputId", 31));
}

/// InputEdgeCrossings represents all pairs of intersecting input edges.
/// It is sorted in lexicographic order.
alias InputEdgeCrossings = Tuple!(InputEdgeId, CrossingInputEdge)[];

/**
 * Given two input edges A and B that intersect, suppose that A maps to a
 * chain of snapped edges A_0, A_1, ..., A_m and B maps to a chain of snapped
 * edges B_0, B_1, ..., B_n.  CrossingGraphEdge represents an edge from chain
 * B that shares a vertex with chain A.  It is used as a temporary data
 * representation while processing chain A.  The arguments are:
 *
 *   "id" - the Graph::EdgeId of an edge from chain B.
 *   "a_index" - the index of the vertex (A_i) that is shared with chain A.
 *   "outgoing" - true if the shared vertex is the first vertex of the B edge.
 *   "dst" - the Graph::VertexId of the vertex that is not shared with chain A.
 *
 * Note that if an edge from the B chain shares both vertices with the A
 * chain, there will be two entries: an outgoing edge that treats its first
 * vertex as being shared, and an incoming edge that treats its second vertex
 * as being shared.
 */
struct CrossingGraphEdge {
  EdgeId id;
  int aIndex;
  bool outgoing;
  VertexId dst;
}
alias CrossingGraphEdgeVector = CrossingGraphEdge[];

/**
 * Returns a vector of EdgeIds sorted by input edge id.  When more than one
 * output edge has the same input edge id (i.e., the input edge snapped to a
 * chain of edges), the edges are sorted so that they form a directed edge
 * chain.
 *
 * This function could possibily be moved to S2Builder::Graph, but note that
 * it has special requirements.  Namely, duplicate edges and sibling pairs
 * must be kept in order to ensure that every output edge corresponds to
 * exactly one input edge.  (See also S2Builder::Graph::GetInputEdgeOrder.)
 */
private EdgeId[] getInputEdgeChainOrder(in Graph g, in InputEdgeId[] input_ids)
in {
  assert(g.options().edgeType() == EdgeType.DIRECTED);
  assert(g.options().duplicateEdges() == DuplicateEdges.KEEP);
  assert(g.options().siblingPairs() == SiblingPairs.KEEP);
} do {
  // First, sort the edges so that the edges corresponding to each input edge
  // are consecutive.  (Each input edge was snapped to a chain of output
  // edges, or two chains in the case of undirected input edges.)
  EdgeId[] order = g.getInputEdgeOrder(input_ids);

  // Now sort the group of edges corresponding to each input edge in edge
  // chain order (e.g.  AB, BC, CD).
  Tuple!(VertexId, EdgeId)[] vmap;     // Map from source vertex to edge id.
  int[] indegree = new int[](g.numVertices());  // Restricted to current input edge.
  for (int end, begin = 0; begin < order.length; begin = end) {
    // Gather the edges that came from a single input edge.
    InputEdgeId input_id = input_ids[order[begin]];
    for (end = begin; end < order.length; ++end) {
      if (input_ids[order[end]] != input_id) break;
    }
    if (end - begin == 1) continue;

    // Build a map from the source vertex of each edge to its edge id,
    // and also compute the indegree at each vertex considering only the edges
    // that came from the current input edge.
    for (int i = begin; i < end; ++i) {
      EdgeId e = order[i];
      vmap ~= tuple(g.edge(e)[0], e);
      indegree[g.edge(e)[1]] += 1;
    }
    sort(vmap);

    // Find the starting edge for building the edge chain.
    EdgeId next = g.numEdges();
    for (int i = begin; i < end; ++i) {
      EdgeId e = order[i];
      if (indegree[g.edge(e)[0]] == 0) next = e;
    }
    // Build the edge chain.
    for (int i = begin; ;) {
      order[i] = next;
      VertexId v = g.edge(next)[1];
      indegree[v] = 0;  // Clear as we go along.
      if (++i == end) break;
      auto triRanges = assumeSorted(vmap).trisect(tuple(v, 0));
      auto output = chain(triRanges[1], triRanges[2]);
      debug enforce(!output.empty() && output.front[0] == v);
      next = output.front[1];
    }
    vmap.length = 0;
  }
  return order;
}

/**
 * Given a set of clipping instructions encoded as a set of InputEdgeCrossings,
 * GraphEdgeClipper determines which graph edges correspond to clipped
 * portions of input edges and removes them.
 *
 * The clipping model is as follows.  The input consists of edge chains.  The
 * clipper maintains an "inside" boolean state as it clips each chain, and
 * toggles this state whenever an input edge is crossed.  Any edges that are
 * deemed to be "outside" after clipping are removed.
 *
 * The "inside" state can be reset when necessary (e.g., when jumping to the
 * start of a new chain) by adding a special crossing marked kSetInside.
 * There are also two other special "crossings" that modify the clipping
 * parameters: kSetInvertB specifies that edges should be clipped to the
 * exterior of the other region, and kSetReverseA specifies that edges should
 * be reversed before emitting them (which is needed to implement difference
 * operations).
 */
class GraphEdgeClipper {
public:
  // "input_dimensions" is a vector specifying the dimension of each input
  // edge (0, 1, or 2).  "input_crossings" is the set of all crossings to be
  // used when clipping the edges of "g", sorted in lexicographic order.
  //
  // The clipped set of edges and their corresponding set of input edge ids
  // are returned in "new_edges" and "new_input_edge_ids".  (These can be used
  // to construct a new S2Builder::Graph.)
  this(
      in Graph g, in byte[] input_dimensions,
      in InputEdgeCrossings input_crossings,
      Graph.Edge[]* new_edges,
      InputEdgeIdSetId[]* new_input_edge_ids) {
    _g = g;
    _in = new Graph.VertexInMap(g);
    _out = new Graph.VertexOutMap(g);
    _inputDimensions = input_dimensions;
    _inputCrossings = input_crossings;
    _newEdges = new_edges;
    _newInputEdgeIds = new_input_edge_ids;
    _inputIds = g.inputEdgeIdSetIds();
    _order = getInputEdgeChainOrder(_g, _inputIds);
    _rank.length = _order.length;
    for (int i = 0; i < _order.length; ++i) {
      _rank[_order[i]] = i;
    }
  }

  void run() {
    // Declare vectors here and reuse them to avoid reallocation.
    VertexId[] a_vertices;
    int[] a_num_crossings;
    bool[] a_isolated;
    CrossingInputEdge[] b_input_edges;
    CrossingGraphEdgeVector[] b_edges;

    bool inside = false;
    bool invert_b = false;
    bool reverse_a = false;
    for (int i = 0; i < _order.length; ++i) {
      // For each input edge (the "A" input edge), gather all the input edges
      // that cross it (the "B" input edges).
      InputEdgeId a_input_id = _inputIds[_order[i]];
      const(Graph.Edge) edge0 = _g.edge(_order[i]);
      b_input_edges.length = 0;
      foreach (next; _inputCrossings) {
        if (next[0] != a_input_id) break;
        if (next[1].inputId() >= 0) {
          b_input_edges ~= next[1];
        } else if (next[1].inputId() == SET_INSIDE) {
          inside = next[1].leftToRight();
        } else if (next[1].inputId() == SET_INVERT_B) {
          invert_b = next[1].leftToRight();
        } else {
          debug enforce(next[1].inputId() == SET_REVERSE_A);
          reverse_a = next[1].leftToRight();
        }
      }
      // Optimization for degenerate edges.
      // TODO(ericv): If the output layer for this edge dimension specifies
      // DegenerateEdges::DISCARD, then remove the edge here.
      if (edge0[0] == edge0[1]) {
        inside ^= (b_input_edges.length & 1);
        addEdge(edge0, a_input_id);
        continue;
      }
      // Optimization for the case where there are no crossings.
      if (b_input_edges.empty()) {
        // In general the caller only passes edges that are part of the output
        // (i.e., we could DCHECK(inside) here).  The one exception is for
        // polyline/polygon operations, where the polygon edges are needed to
        // compute the polyline output but are not emitted themselves.
        if (inside) {
          addEdge(reverse_a ? Graph.reverse(edge0) : edge0, a_input_id);
        }
        continue;
      }
      // Walk along the chain of snapped edges for input edge A, and at each
      // vertex collect all the incident edges that belong to one of the
      // crossing edge chains (the "B" input edges).
      a_vertices.length = 0;
      a_vertices ~= edge0[0];
      b_edges.length = 0;
      b_edges.length = b_input_edges.length;
      gatherIncidentEdges(a_vertices, 0, b_input_edges, b_edges);
      for (; i < _order.length && _inputIds[_order[i]] == a_input_id; ++i) {
        a_vertices ~= _g.edge(_order[i])[1];
        gatherIncidentEdges(a_vertices, cast(int) a_vertices.length - 1, b_input_edges, b_edges);
      }
      --i;
      if (s2builderVerbose) {
        writeln("input edge ", a_input_id, " (inside=", inside, "):");
        foreach (VertexId id; a_vertices) writeln(" ", id);
      }
      // Now for each B edge chain, decide which vertex of the A chain it
      // crosses, and keep track of the number of signed crossings at each A
      // vertex.  The sign of a crossing depends on whether the other edge
      // crosses from left to right or right to left.
      //
      // This would not be necessary if all calculations were done in exact
      // arithmetic, because crossings would have strictly alternating signs.
      // But because we have already snapped the result, some crossing locations
      // are ambiguous, and GetCrossedVertexIndex() handles this by choosing a
      // candidate vertex arbitrarily.  The end result is that rarely, we may
      // see two crossings in a row with the same sign.  We correct for this by
      // adding extra output edges that essentially link up the crossings in the
      // correct (alternating sign) order.  Compared to the "correct" behavior,
      // the only difference is that we have added some extra sibling pairs
      // (consisting of an edge and its corresponding reverse edge) which do not
      // affect the result.
      a_num_crossings.length = a_vertices.length;
      a_isolated.length = 0;
      a_isolated.length = a_vertices.length;
      for (int bi = 0; bi < b_input_edges.length; ++bi) {
        bool left_to_right = b_input_edges[bi].leftToRight();
        int a_index = getCrossedVertexIndex(a_vertices, b_edges[bi], left_to_right);
        if (s2builderVerbose) {
          writeln("  b input edge ", b_input_edges[bi].inputId(), " (l2r=", left_to_right,
              ", crossing=", a_vertices[a_index], ")");
          foreach (x; b_edges[bi]) {
            const Graph.Edge e = _g.edge(x.id);
            writeln(" (", e[0], ", ", e[1], ")");
          }
        }
        // Keep track of the number of signed crossings (see above).
        bool is_line = _inputDimensions[b_input_edges[bi].inputId()] == 1;
        int sign = is_line ? 0 : (left_to_right == invert_b) ? -1 : 1;
        a_num_crossings[a_index] += sign;

        // Any polyline or polygon vertex that has at least one crossing but no
        // adjacent emitted edge may be emitted as an isolated vertex.
        a_isolated[a_index] = true;
      }
      if (s2builderVerbose) writeln();

      // Finally, we iterate through the A edge chain, keeping track of the
      // number of signed crossings as we go along.  The "multiplicity" is
      // defined as the cumulative number of signed crossings, and indicates how
      // many edges should be output (and in which direction) in order to link
      // up the edge crossings in the correct order.  (The multiplicity is
      // almost always either 0 or 1 except in very rare cases.)
      int multiplicity = inside + a_num_crossings[0];
      for (int ai = 1; ai < a_vertices.length; ++ai) {
        if (multiplicity != 0) {
          a_isolated[ai - 1] = a_isolated[ai] = false;
        }
        int edge_count = reverse_a ? -multiplicity : multiplicity;
        // Output any forward edges required.
        for (int j = 0; j < edge_count; ++j) {
          addEdge(Graph.Edge(a_vertices[ai - 1], a_vertices[ai]), a_input_id);
        }
        // Output any reverse edges required.
        for (int j = edge_count; j < 0; ++j) {
          addEdge(Graph.Edge(a_vertices[ai], a_vertices[ai - 1]), a_input_id);
        }
        multiplicity += a_num_crossings[ai];
      }
      // Multiplicities other than 0 or 1 can only occur in the edge interior.
      debug enforce(multiplicity == 0 || multiplicity == 1);
      inside = (multiplicity != 0);

      // Output any isolated polyline vertices.
      // TODO(ericv): Only do this if an output layer wants degenerate edges.
      if (_inputDimensions[a_input_id] != 0) {
        for (int ai = 0; ai < a_vertices.length; ++ai) {
          if (a_isolated[ai]) {
            addEdge(Graph.Edge(a_vertices[ai], a_vertices[ai]), a_input_id);
          }
        }
      }
    }
  }

private:
  void addEdge(Graph.Edge edge, InputEdgeId input_edge_id) {
    *_newEdges ~= edge;
    *_newInputEdgeIds ~= input_edge_id;
  }

  /**
   * Given the vertices of the snapped edge chain for an input edge A and the
   * set of input edges B that cross input edge A, this method gathers all of
   * the snapped edges of B that are incident to a given snapped vertex of A.
   * The incident edges for each input edge of B are appended to a separate
   * output vector.  (A and B can refer to either the input edge or the
   * corresponding snapped edge chain.)
   */
  void gatherIncidentEdges(
      in VertexId[] a, int ai, in CrossingInputEdge[] b_input_edges,
      ref CrossingGraphEdgeVector[] b_edges) const
  in {
    assert(b_input_edges.length == b_edges.length);
  } do {
    // Examine all of the edges incident to the given vertex of A.  If any edge
    // comes from a B input edge, append it to the appropriate vector.
    foreach (EdgeId e; _in.edgeIds(a[ai])) {
      InputEdgeId id = _inputIds[e];
      auto trisectRanges = b_input_edges.assumeSorted().trisect(id);
      if (!trisectRanges[1].empty()) {
        auto edges = &b_edges[trisectRanges[0].length];
        *edges ~= CrossingGraphEdge(e, ai, false, _g.edge(e)[0]);
      }
    }
    foreach (EdgeId e; _out.edgeIds(a[ai])) {
      InputEdgeId id = _inputIds[e];
      auto trisectRanges = b_input_edges.assumeSorted().trisect(id);
      if (!trisectRanges[1].empty()) {
        auto edges = &b_edges[trisectRanges[0].length];
        *edges ~= CrossingGraphEdge(e, ai, true, _g.edge(e)[1]);
      }
    }
  }

  /**
   * Given an edge chain A that is crossed by another edge chain B (where
   * "left_to_right" indicates whether B crosses A from left to right), this
   * method decides which vertex of A the crossing takes place at.  The
   * parameters are the vertices of the A chain ("a") and the set of edges in
   * the B chain ("b") that are incident to vertices of A.  The B chain edges
   * are sorted in increasing order of (a_index, outgoing) tuple.
   */
  int getCrossedVertexIndex(
      in VertexId[] a, in CrossingGraphEdgeVector b, bool left_to_right) const
  in {
    assert(!a.empty());
    assert(!b.empty());
  } do {
    import s2.s2predicates : orderedCCW;

    // The reason this calculation is tricky is that after snapping, the A and B
    // chains may meet and separate several times.  For example, if B crosses A
    // from left to right, then B may touch A, make an excursion to the left of
    // A, come back to A, then make an excursion to the right of A and come back
    // to A again, like this:
    //
    //  *--B--*-\             /-*-\
    //           B-\       /-B     B-\      6     7     8     9
    //  *--A--*--A--*-A,B-*--A--*--A--*-A,B-*--A--*--A--*-A,B-*
    //  0     1     2     3     4     5      \-B     B-/
    //                                          \-*-/
    //
    // (where "*" is a vertex, and "A" and "B" are edge labels).  Note that B
    // may also follow A for one or more edges whenever they touch (e.g. between
    // vertices 2 and 3 ).  In this case the only vertices of A where the
    // crossing could take place are 5 and 6, i.e. after all excursions of B to
    // the left of A, and before all excursions of B to the right of A.
    //
    // Other factors to consider are that the portion of B before and/or after
    // the crossing may be degenerate, and some or all of the B edges may be
    // reversed relative to the A edges.

    // First, check whether edge A is degenerate.
    int n = cast(int) a.length;
    if (n == 1) return 0;

    // If edge chain B is incident to only one vertex of A, we're done.
    if (b[0].aIndex == b.back().aIndex) return b[0].aIndex;

    // Determine whether the B chain visits the first and last vertices that it
    // shares with the A chain in the same order or the reverse order.  This is
    // only needed to implement one special case (see below).
    bool b_reversed = getVertexRank(b[0]) > getVertexRank(b.back());

    // Examine each incident B edge and use it to narrow the range of positions
    // where the crossing could occur in the B chain.  Vertex positions are
    // represented as a range [lo, hi] of vertex ranks in the B chain (see
    // GetVertexRank).
    //
    // Note that if an edge of B is incident to the first or last vertex of A,
    // we can't test which side of the A chain it is on.  There can be up to 4
    // such edges (one incoming and one outgoing edge at each vertex).  Two of
    // these edges logically extend past the end of the A chain and place no
    // restrictions on the crossing vertex.  The other two edges define the ends
    // of the subchain where B shares vertices with A.  We save these edges in
    // order to handle a special case (see below).
    int lo = -1, hi = cast(int) _order.length;   // Vertex ranks of acceptable crossings
    EdgeId b_first = -1, b_last = -1;  // "b" subchain connecting "a" endpoints
    foreach (e; b) {
      int ai = e.aIndex;
      if (ai == 0) {
        if (e.outgoing != b_reversed && e.dst != a[1]) b_first = e.id;
      } else if (ai == n - 1) {
        if (e.outgoing == b_reversed && e.dst != a[n - 2]) b_last = e.id;
      } else {
        // This B edge is incident to an interior vertex of the A chain.  First
        // check whether this edge is identical (or reversed) to an edge in the
        // A chain, in which case it does not create any restrictions.
        if (e.dst == a[ai - 1] || e.dst == a[ai + 1]) continue;

        // Otherwise we can test which side of the A chain the edge lies on.
        bool on_left = orderedCCW(
            _g.vertex(a[ai + 1]), _g.vertex(e.dst), _g.vertex(a[ai - 1]), _g.vertex(a[ai]));

        // Every B edge that is incident to an interior vertex of the A chain
        // places some restriction on where the crossing vertex could be.
        if (left_to_right == on_left) {
          // This is a pre-crossing edge, so the crossing cannot be before the
          // destination vertex of this edge.  (For example, the input B edge
          // crosses the input A edge from left to right and this edge of the B
          // chain is to the left of the A chain.)
          lo = max(lo, _rank[e.id] + 1);
        } else {
          // This is a post-crossing edge, so the crossing cannot be after the
          // source vertex of this edge.
          hi = min(hi, _rank[e.id]);
        }
      }
    }
    // There is one special case.  If there were no B edges incident to interior
    // vertices of A, then we can't reliably test which side of A the B edges
    // are on.  (An s2pred::Sign test doesn't work, since an edge of B can snap to
    // the "wrong" side of A while maintaining topological guarantees.)  So
    // instead we construct a loop consisting of the A edge chain plus the
    // portion of the B chain that connects the endpoints of A.  We can then
    // test the orientation of this loop.
    //
    // Note that it would be possible to avoid this in some situations by
    // testing whether either endpoint of the A chain has two incident B edges,
    // in which case we could check which side of the B chain the A edge is on
    // and use this to limit the possible crossing locations.
    if (lo < 0 && hi >= _order.length && b_first >= 0 && b_last >= 0) {
      // Swap the edges if necessary so that they are in B chain order.
      if (b_reversed) swap(b_first, b_last);
      bool on_left = edgeChainOnLeft(a, b_first, b_last);
      if (left_to_right == on_left) {
        lo = max(lo, _rank[b_last] + 1);
      } else {
        hi = min(hi, _rank[b_first]);
      }
    }

    // Otherwise we choose the smallest shared VertexId in the acceptable range,
    // in order to ensure that both chains choose the same crossing vertex.
    int best = -1;
    debug enforce(lo <= hi);
    foreach (e; b) {
      int ai = e.aIndex;
      int vrank = getVertexRank(e);
      if (vrank >= lo && vrank <= hi && (best < 0 || a[ai] < a[best])) {
        best = ai;
      }
    }
    return best;
  }

  /**
   * Returns the "vertex rank" of the shared vertex associated with the given
   * CrossingGraphEdge.  Recall that graph edges are sorted in input edge order,
   * and that the rank of an edge is its position in this order (rank_[e]).
   * VertexRank(e) is defined such that VertexRank(e.src) == rank_[e] and
   * VertexRank(e.dst) == rank_[e] + 1.  Note that the concept of "vertex rank"
   * is only defined within a single edge chain (since different edge chains can
   * have overlapping vertex ranks).
   */
  int getVertexRank(in CrossingGraphEdge e) const {
    return _rank[e.id] + !e.outgoing;
  }

  /**
   * Given edge chains A and B that form a loop (after possibly reversing the
   * direction of chain B), returns true if chain B is to the left of chain A.
   * Chain A is given as a sequence of vertices, while chain B is specified as
   * the first and last edges of the chain.
   */
  bool edgeChainOnLeft(in VertexId[] a, EdgeId b_first, EdgeId b_last) const {
    import s2.s2measures : turnAngle;

    // Gather all the interior vertices of the B subchain.
    VertexId[] loop;
    for (int i = _rank[b_first]; i < _rank[b_last]; ++i) {
      loop ~= _g.edge(_order[i])[1];
    }
    // Possibly reverse the chain so that it forms a loop when "a" is appended.
    if (_g.edge(b_last)[1] != a[0]) reverse(loop);
    loop ~= a;
    // Duplicate the first two vertices to simplify vertex indexing.
    loop ~= loop[0 .. 2];
    // Now B is to the left of A if and only if the loop is counterclockwise.
    double sum = 0;
    for (int i = 2; i < loop.length; ++i) {
      sum += turnAngle(_g.vertex(loop[i - 2]), _g.vertex(loop[i - 1]), _g.vertex(loop[i]));
    }
    return sum > 0;
  }

  const(Graph) _g;
  Graph.VertexInMap _in;
  Graph.VertexOutMap _out;
  const(byte[]) _inputDimensions;
  const(InputEdgeCrossings) _inputCrossings;
  Graph.Edge[]* _newEdges;
  InputEdgeIdSetId[]* _newInputEdgeIds;

  // Every graph edge is associated with exactly one input edge in our case,
  // which means that we can declare g_.input_edge_id_set_ids() as a vector of
  // InputEdgeIds rather than a vector of InputEdgeIdSetIds.  (This also takes
  // advantage of the fact that IdSetLexicon represents a singleton set as the
  // value of its single element.)
  const(InputEdgeId[]) _inputIds;

  EdgeId[] _order;  // Graph edges sorted in input edge id order.
  int[] _rank;      // The rank of each graph edge within order_.
}

/**
 * Given a set of clipping instructions encoded as a set of intersections
 * between input edges, EdgeClippingLayer determines which graph edges
 * correspond to clipped portions of input edges and removes them.  It
 * assembles the remaining edges into a new S2Builder::Graph and passes the
 * result to the given output layer for assembly.
 */
class EdgeClippingLayer : Layer {
public:
  this(Layer[] layers, in byte[]* input_dimensions, in InputEdgeCrossings* input_crossings) {
    _layers = layers;
    _inputDimensions = input_dimensions;
    _inputCrossings = input_crossings;
  }

  /// Layer interface:
  override
  GraphOptions graphOptions() const {
  // We keep all edges, including degenerate ones, so that we can figure out
  // the correspondence between input edge crossings and output edge
  // crossings.
    return new GraphOptions(
        EdgeType.DIRECTED, DegenerateEdges.KEEP,
        DuplicateEdges.KEEP, SiblingPairs.KEEP);
  }

  override
  void build(Graph g, ref S2Error error) {
    // The bulk of the work is handled by GraphEdgeClipper.
    Graph.Edge[] new_edges;
    InputEdgeIdSetId[] new_input_edge_ids;
    // Destroy the GraphEdgeClipper immediately to save memory.
    new GraphEdgeClipper(g, *_inputDimensions, *_inputCrossings, &new_edges, &new_input_edge_ids)
        .run();
    if (s2builderVerbose) {
      writeln("Edges after clipping: ");
      for (int i = 0; i < new_edges.length; ++i) {
        writeln("  ", new_input_edge_ids[i], " (", new_edges[i][0], ", ", new_edges[i][1], ")");
      }
    }
    // Construct one or more graphs from the clipped edges and pass them to the
    // given output layer(s).
    auto new_input_edge_id_set_lexicon = new IdSetLexicon();
    if (_layers.length == 1) {
      GraphOptions options = _layers[0].graphOptions();
      Graph new_graph = makeGraph(
          g, options, new_edges, new_input_edge_ids, new_input_edge_id_set_lexicon, error);
      _layers[0].build(new_graph, error);
    } else {
      // The Graph objects must be valid until the last Build() call completes,
      // so we store all of the graph data in arrays with 3 elements.
      debug enforce(_layers.length == 3);
      Graph.Edge[][3] layer_edges;
      InputEdgeIdSetId[][3] layer_input_edge_ids;
      GraphOptions[3] layer_options;
      Graph[] layer_graphs;  // No default constructor.
      layer_graphs.reserve(3);
      // Separate the edges according to their dimension.
      for (int i = 0; i < new_edges.length; ++i) {
        int d = (*_inputDimensions)[new_input_edge_ids[i]];
        layer_edges[d] ~= new_edges[i];
        layer_input_edge_ids[d] ~= new_input_edge_ids[i];
      }
      // Clear variables to save space.
      new_edges.length = 0;
      new_input_edge_ids.length = 0;
      for (int d = 0; d < 3; ++d) {
        layer_options[d] = _layers[d].graphOptions();
        layer_graphs ~= makeGraph(
            g, layer_options[d], layer_edges[d], layer_input_edge_ids[d],
            new_input_edge_id_set_lexicon, error);
        _layers[d].build(layer_graphs[d], error);
      }
    }
  }

private:
  // Helper function (in anonymous namespace) to create an S2Builder::Graph from
  // a vector of edges.
  static Graph makeGraph(
      in Graph g, GraphOptions options, ref Graph.Edge[] new_edges,
      ref InputEdgeIdSetId[] new_input_edge_ids,
      ref IdSetLexicon new_input_edge_id_set_lexicon, ref S2Error error) {
    if (options.edgeType() == EdgeType.UNDIRECTED) {
      // Create a reversed edge for every edge.
      size_t n = new_edges.length;
      new_edges.reserve(2 * n);
      new_input_edge_ids.reserve(2 * n);
      for (int i = 0; i < n; ++i) {
        new_edges ~= Graph.reverse(new_edges[i]);
        new_input_edge_ids ~= IdSetLexicon.emptySetId();
      }
    }
    Graph.processEdges(
        options, new_edges, new_input_edge_ids, new_input_edge_id_set_lexicon, error);
    return new Graph(
        options, g.vertices(), new_edges, new_input_edge_ids,
        new_input_edge_id_set_lexicon, g.labelSetIds(),
        g.labelSetLexicon(), g.isFullPolygonPredicate());
  }

  Layer[] _layers;
  const(byte[])* _inputDimensions;
  const(InputEdgeCrossings)* _inputCrossings;
}

/**
 * Given a polygon edge graph containing only degenerate edges and sibling
 * edge pairs, the purpose of this function is to decide whether the polygon
 * is empty or full except for the degeneracies, i.e. whether the degeneracies
 * represent shells or holes.
 *
 * This function always returns false, meaning that the polygon is empty and
 * the degeneracies represent shells.  The main side effect of this is that
 * operations whose result should be the full polygon will instead be the
 * empty polygon.  (Classes such as S2Polygon already have code to correct for
 * this, but if that functionality were moved here then it would be useful for
 * other polygon representations such as S2LaxPolygonShape.)
 */
private bool isFullPolygonNever(in Graph g, out S2Error error) {
  return false;  // Assumes the polygon is empty.
}
