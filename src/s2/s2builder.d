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
// This class is a replacement for S2PolygonBuilder.  Once all clients have
// been updated to use this class, S2PolygonBuilder will be removed.

module s2.s2builder;

import s2.builder.graph;
import s2.builder.layer;
import s2.builder.util.snap_functions;
import s2.id_set_lexicon;
import s2.logger;
import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cell_id;
import s2.s2closest_edge_query;
import s2.s2closest_point_query;
import s2.s2edge_crossings : INTERSECTION_ERROR;
import s2.s2edge_distances : getUpdateMinDistanceMaxError;
import s2.s2error;
import s2.s2loop;
import s2.s2point;
import s2.s2point_index;
import s2.s2polyline;
import s2.s2polyline_simplifier;
import s2.s2predicates;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2text_format;
import s2.util.math.vector;

import std.algorithm;
import std.exception;
import std.math;
import std.range;
import std.stdio;
import std.typecons;

/// Internal flag intended to be set from within a debugger.
bool s2builderVerbose = false;

/**
 * S2Builder is a tool for assembling polygonal geometry from edges.  Here are
 * some of the things it is designed for:
 *
 * 1. Building polygons, polylines, and polygon meshes from unsorted
 *    collections of edges.
 *
 * 2. Snapping geometry to discrete representations (such as S2CellId centers
 *    or E7 lat/lng coordinates) while preserving the input topology and with
 *    guaranteed error bounds.
 *
 * 3. Simplifying geometry (e.g. for indexing, display, or storage).
 *
 * 4. Importing geometry from other formats, including repairing geometry
 *    that has errors.
 *
 * 5. As a tool for implementing more complex operations such as polygon
 *    intersections and unions.
 *
 * The implementation is based on the framework of "snap rounding".  Unlike
 * most snap rounding implementations, S2Builder defines edges as geodesics on
 * the sphere (straight lines) and uses the topology of the sphere (i.e.,
 * there are no "seams" at the poles or 180th meridian).  The algorithm is
 * designed to be 100% robust for arbitrary input geometry.  It offers the
 * following properties:
 *
 *   - Guaranteed bounds on how far input vertices and edges can move during
 *     the snapping process (i.e., at most the given "snap_radius").
 *
 *   - Guaranteed minimum separation between edges and vertices other than
 *     their endpoints (similar to the goals of Iterated Snap Rounding).  In
 *     other words, edges that do not intersect in the output are guaranteed
 *     to have a minimum separation between them.
 *
 *   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
 *     input already meets the output criteria then it will not be modified.
 *
 *   - Preservation of the input topology (up to the creation of
 *     degeneracies).  This means that there exists a continuous deformation
 *     from the input to the output such that no vertex crosses an edge.  In
 *     other words, self-intersections won't be created, loops won't change
 *     orientation, etc.
 *
 *   - The ability to snap to arbitrary discrete point sets (such as S2CellId
 *     centers, E7 lat/lng points on the sphere, or simply a subset of the
 *     input vertices), rather than being limited to an integer grid.
 *
 * Here are some of its other features:
 *
 *  - It can handle both directed and undirected edges.  Undirected edges can
 *    be useful for importing data from other formats, e.g. where loops have
 *    unspecified orientations.
 *
 *  - It can eliminate self-intersections by finding all edge pairs that cross
 *    and adding a new vertex at each intersection point.
 *
 *  - It can simplify polygons to within a specified tolerance.  For example,
 *    if two vertices are close enough they will be merged, and if an edge
 *    passes nearby a vertex then it will be rerouted through that vertex.
 *    Optionally, it can also detect nearly straight chains of short edges and
 *    replace them with a single long edge, while maintaining the same
 *    accuracy, separation, and topology guarantees ("simplify_edge_chains").
 *
 *  - It supports many different output types through the concept of "layers"
 *    (polylines, polygons, polygon meshes, etc).  You can build multiple
 *    layers at once in order to ensure that snapping does not create
 *    intersections between different objects (for example, you can simplify a
 *    set of contour lines without the risk of having them cross each other).
 *
 *  - It supports edge labels, which allow you to attach arbitrary information
 *    to edges and have it preserved during the snapping process.  (This can
 *    also be achieved using layers, at a coarser level of granularity.)
 *
 * Caveats:
 *
 *  - Because S2Builder only works with edges, it cannot distinguish between
 *    the empty and full polygons.  If your application can generate both the
 *    empty and full polygons, you must implement logic outside of this class.
 *
 * Example showing how to snap a polygon to E7 coordinates:
 *
 *  using s2builderutil::IntLatLngSnapFunction;
 *  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(7)));
 *  S2Polygon output;
 *  builder.StartLayer(absl::make_unique<s2builderutil::S2PolygonLayer>(&output));
 *  builder.AddPolygon(input);
 *  S2Error error;
 *  if (!builder.Build(&error)) {
 *    LOG(ERROR) << error;
 *    ...
 *  }
 */
class S2Builder {
public:
  /**
   * Indicates whether the input edges are undirected.  Typically this is
   * specified for each output layer (e.g., s2builderutil::S2PolygonLayer).
   *
   * Directed edges are preferred, since otherwise the output is ambiguous.
   * For example, output polygons may be the *inverse* of the intended result
   * (e.g., a polygon intended to represent the world's oceans may instead
   * represent the world's land masses).  Directed edges are also somewhat
   * more efficient.
   *
   * However even with undirected edges, most S2Builder layer types try to
   * preserve the input edge direction whenever possible.  Generally, edges
   * are reversed only when it would yield a simpler output.  For example,
   * S2PolygonLayer assumes that polygons created from undirected edges should
   * cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
   * assembles edges into as few polylines as possible, even if this means
   * reversing some of the "undirected" input edges.
   *
   * For shapes with interiors, directed edges should be oriented so that the
   * interior is to the left of all edges.  This means that for a polygon with
   * holes, the outer loops ("shells") should be directed counter-clockwise
   * while the inner loops ("holes") should be directed clockwise.  Note that
   * S2Builder::AddPolygon() follows this convention automatically.
   */
  enum EdgeType { DIRECTED, UNDIRECTED }

  /**
   * A SnapFunction restricts the locations of the output vertices.  For
   * example, there are predefined snap functions that require vertices to be
   * located at S2CellId centers or at E5/E6/E7 coordinates.  The SnapFunction
   * can also specify a minimum spacing between vertices (the "snap radius").
   *
   * A SnapFunction defines the following methods:
   *
   * 1. The SnapPoint() method, which snaps a point P to a nearby point (the
   *    "candidate snap site").  Any point may be returned, including P
   *    itself (this is the "identity snap function").
   *
   * 2. "snap_radius", the maximum distance that vertices can move when
   *    snapped.  The snap_radius must be at least as large as the maximum
   *    distance between P and SnapPoint(P) for any point P.
   *
   * 3. "max_edge_deviation", the maximum distance that edges can move when
   *    snapped.  It is slightly larger than "snap_radius" because when a
   *    geodesic edge is snapped, the center of the edge moves further than
   *    its endpoints.  This value is computed automatically by S2Builder.
   *
   * 4. "min_vertex_separation", the guaranteed minimum distance between
   *    vertices in the output.  This is generally a fraction of
   *    "snap_radius" where the fraction depends on the snap function.
   *
   * 5. A "min_edge_vertex_separation", the guaranteed minimum distance
   *    between edges and non-incident vertices in the output.  This is
   *    generally a fraction of "snap_radius" where the fraction depends on
   *    the snap function.
   *
   * It is important to note that SnapPoint() does not define the actual
   * mapping from input vertices to output vertices, since the points it
   * returns (the candidate snap sites) are further filtered to ensure that
   * they are separated by at least the snap radius.  For example, if you
   * specify E7 coordinates (2cm resolution) and a snap radius of 10m, then a
   * subset of points returned by SnapPoint will be chosen (the "snap sites"),
   * and each input vertex will be mapped to the closest site.  Therefore you
   * cannot assume that P is necessarily snapped to SnapPoint(P).
   *
   * S2Builder makes the following guarantees:
   *
   * 1. Every vertex is at a location returned by SnapPoint().
   *
   * 2. Vertices are within "snap_radius" of the corresponding input vertex.
   *
   * 3. Edges are within "max_edge_deviation" of the corresponding input edge
   *    (a distance slightly larger than "snap_radius").
   *
   * 4. Vertices are separated by at least "min_vertex_separation"
   *    (a fraction of "snap_radius" that depends on the snap function).
   *
   * 5. Edges and non-incident vertices are separated by at least
   *    "min_edge_vertex_separation" (a fraction of "snap_radius").
   *
   * 6. Vertex and edge locations do not change unless one of the conditions
   *    above is not already met (idempotency / stability).
   *
   * 7. The topology of the input geometry is preserved (up to the creation
   *    of degeneracies).  This means that there exists a continuous
   *    deformation from the input to the output such that no vertex
   *    crosses an edge.
   */
  abstract static class SnapFunction {
  public:

    /**
     * The maximum distance that vertices can move when snapped.
     *
     * If the snap radius is zero, then vertices are snapped together only if
     * they are identical.  Edges will not be snapped to any vertices other
     * than their endpoints, even if there are vertices whose distance to the
     * edge is zero, unless split_crossing_edges() is true.
     *
     * REQUIRES: snap_radius() <= kMaxSnapRadius
     */
    abstract S1Angle snapRadius() const;

    /**
     * The maximum snap radius is just large enough to support snapping to
     * S2CellId level 0.  It is equivalent to 7800km on the Earth's surface.
     */
    static S1Angle kMaxSnapRadius() {
      // This value can't be larger than 85.7 degrees without changing the code
      // related to min_edge_length_to_split_ca_, and increasing it to 90 degrees
      // or more would most likely require significant changes to the algorithm.
      return S1Angle.fromDegrees(70.0);
    }

    /**
     * The maximum distance that the center of an edge can move when snapped.
     * This is slightly larger than "snap_radius" because when a geodesic edge
     * is snapped, the center of the edge moves further than its endpoints.
     */
    S1Angle maxEdgeDeviation() const
    in {
      assert(snapRadius() <= kMaxSnapRadius());
    } body {
      // We want max_edge_deviation() to be large enough compared to snap_radius()
      // such that edge splitting is rare.
      //
      // Using spherical trigonometry, if the endpoints of an edge of length L
      // move by at most a distance R, the center of the edge moves by at most
      // asin(sin(R) / cos(L / 2)).  Thus the (max_edge_deviation / snap_radius)
      // ratio increases with both the snap radius R and the edge length L.
      //
      // We arbitrarily limit the edge deviation to be at most 10% more than the
      // snap radius.  With the maximum allowed snap radius of 70 degrees, this
      // means that edges up to 30.6 degrees long are never split.  For smaller
      // snap radii, edges up to 49 degrees long are never split.  (Edges of any
      // length are not split unless their endpoints move far enough so that the
      // actual edge deviation exceeds the limit; in practice, splitting is rare
      // even with long edges.)  Note that it is always possible to split edges
      // when max_edge_deviation() is exceeded; see MaybeAddExtraSites().
      enum double kMaxEdgeDeviationRatio = 1.1;
      return kMaxEdgeDeviationRatio * snapRadius();
    }

    /**
     * The guaranteed minimum distance between vertices in the output.
     * This is generally some fraction of "snap_radius".
     */
    abstract S1Angle minVertexSeparation() const;

    /**
     * The guaranteed minimum spacing between edges and non-incident vertices
     * in the output.  This is generally some fraction of "snap_radius".
     */
    abstract S1Angle minEdgeVertexSeparation() const;

    /**
     * Returns a candidate snap site for the given point.  The final vertex
     * locations are a subset of the snap sites returned by this function
     * (spaced at least "min_vertex_separation" apart).
     *
     * The only requirement is that SnapPoint(x) must return a point whose
     * distance from "x" is no greater than "snap_radius".
     */
    abstract S2Point snapPoint(in S2Point point) const;

    // Returns a deep copy of this SnapFunction.
    abstract SnapFunction clone() const;

    override
    string toString() const {
      import std.conv;
      return "SnapFunction[]";
    }
  }

  static class Options {
  public:
    this() {
      _snapFunction = new IdentitySnapFunction(S1Angle.zero());
    }

    /// Convenience constructor that calls set_snap_function().
    this(in SnapFunction snap_function) {
      _snapFunction = snap_function.clone();
    }

    this(in Options options) {
      _snapFunction = options._snapFunction.clone();
      _splitCrossingEdges = options._splitCrossingEdges;
      _simplifyEdgeChains = options._simplifyEdgeChains;
      _idempotent = options._idempotent;
    }

    /**
     * Sets the desired snap function.  The snap function is copied
     * internally, so you can safely pass a temporary object.
     *
     * Note that if your input data includes vertices that were created using
     * S2::GetIntersection(), then you should use a "snap_radius" of
     * at least S2::kIntersectionSnapRadius, e.g. by calling
     *
     *  options.set_snap_function(s2builderutil::IdentitySnapFunction(
     *      S2::kIntersectionSnapRadius));
     *
     * DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
     * [This does no snapping and preserves all input vertices exactly.]
     */
    const(SnapFunction) snapFunction() const {
      return _snapFunction;
    }

    void setSnapFunction(in SnapFunction snap_function) {
      _snapFunction = snap_function.clone();
    }

    /**
     * If true, then detect all pairs of crossing edges and eliminate them by
     * adding a new vertex at their intersection point.
     *
     * When this option is true, the effective snap_radius() for edges is
     * increased by S2::kIntersectionError to take into account the
     * additional error when computing intersection points.  In other words,
     * edges may move by up to snap_radius() + S2::kIntersectionError.
     *
     * Undirected edges should always be used when the output is a polygon,
     * since splitting a directed loop at a self-intersection converts it into
     * two loops that don't define a consistent interior according to the
     * "interior is on the left" rule.  (On the other hand, it is fine to use
     * directed edges when defining a polygon *mesh* because in that case the
     * input consists of sibling edge pairs.)
     *
     * Self-intersections can also arise when importing data from a 2D
     * projection.  You can minimize this problem by subdividing the input
     * edges so that the S2 edges (which are geodesics) stay close to the
     * original projected edges (which are curves on the sphere).  This can
     * be done using s2builderutil::EdgeSplitter(), for example.
     *
     * DEFAULT: false
     */
    bool splitCrossingEdges() const {
      return _splitCrossingEdges;
    }

    void setSplitCrossingEdges(bool split_crossing_edges) {
      _splitCrossingEdges = split_crossing_edges;
    }

    /**
     * If true, then simplify the output geometry by replacing nearly straight
     * chains of short edges with a single long edge.
     *
     * The combined effect of snapping and simplifying will not change the
     * input by more than the guaranteed tolerances (see the list documented
     * with the SnapFunction class).  For example, simplified edges are
     * guaranteed to pass within snap_radius() of the *original* positions of
     * all vertices that were removed from that edge.  This is a much tighter
     * guarantee than can be achieved by snapping and simplifying separately.
     *
     * However, note that this option does not guarantee idempotency.  In
     * other words, simplifying geometry that has already been simplified once
     * may simplify it further.  (This is unavoidable, since tolerances are
     * measured with respect to the original geometry, which is no longer
     * available when the geometry is simplified a second time.)
     *
     * When the output consists of multiple layers, simplification is
     * guaranteed to be consistent: for example, edge chains are simplified in
     * the same way across layers, and simplification preserves topological
     * relationships between layers (e.g., no crossing edges will be created).
     * Note that edge chains in different layers do not need to be identical
     * (or even have the same number of vertices, etc) in order to be
     * simplified together.  All that is required is that they are close
     * enough together so that the same simplified edge can meet all of their
     * individual snapping guarantees.
     *
     * Note that edge chains are approximated as parametric curves rather than
     * point sets.  This means that if an edge chain backtracks on itself (for
     * example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
     * within snap_radius() (for example, if the preceding point were all in a
     * straight line then the edge chain would be simplified to ACFCFH, noting
     * that C and F have degree > 2 and therefore can't be simplified away).
     *
     * Simplified edges are assigned all labels associated with the edges of
     * the simplified chain.
     *
     * For this option to have any effect, a SnapFunction with a non-zero
     * snap_radius() must be specified.  Also note that vertices specified
     * using ForceVertex are never simplified away.
     *
     * DEFAULT: false
     */
    bool simplifyEdgeChains() const {
      return _simplifyEdgeChains;
    }

    void setSimplifyEdgeChains(bool simplify_edge_chains) {
      _simplifyEdgeChains = simplify_edge_chains;

      // Simplification requires a non-zero snap radius, and while it might be
      // possible to do some simplifying without snapping, it is much simpler to
      // always snap (even if the input geometry already meets the other output
      // requirements).  We need to compute edge_sites_ in order to avoid
      // approaching non-incident vertices too closely, for example.
      setIdempotent(false);
    }

    /**
     * If true, then snapping occurs only when the input geometry does not
     * already meet the S2Builder output guarantees (see the SnapFunction
     * class description for details).  This means that if all input vertices
     * are at snapped locations, all vertex pairs are separated by at least
     * min_vertex_separation(), and all edge-vertex pairs are separated by at
     * least min_edge_vertex_separation(), then no snapping is done.
     *
     * If false, then all vertex pairs and edge-vertex pairs closer than
     * "snap_radius" will be considered for snapping.  This can be useful, for
     * example, if you know that your geometry contains errors and you want to
     * make sure that features closer together than "snap_radius" are merged.
     *
     * This option is automatically turned off by simplify_edge_chains(),
     * since simplifying edge chains is never guaranteed to be idempotent.
     *
     * DEFAULT: true
     */
    bool idempotent() const {
      return _idempotent;
    }

    void setIdempotent(bool idempotent) {
      _idempotent = idempotent;
    }

  private:
    SnapFunction _snapFunction;
    bool _splitCrossingEdges;
    bool _simplifyEdgeChains;
    bool _idempotent = true;
  }

  /**
   * For output layers that represent polygons, there is an ambiguity inherent
   * in spherical geometry that does not exist in planar geometry.  Namely, if
   * a polygon has no edges, does it represent the empty polygon (containing
   * no points) or the full polygon (containing all points)?  This ambiguity
   * also occurs for polygons that consist only of degeneracies, e.g. a
   * degenerate loop with only two edges could be either a degenerate shell in
   * the empty polygon or a degenerate hole in the full polygon.
   *
   * To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
   * each input layer (see AddIsFullPolygonPredicate below).  If the layer
   * consists only of polygon degeneracies, the layer implementation may call
   * this method to determine whether the polygon is empty or full except for
   * the given degeneracies.  (Note that under the semi-open boundary model,
   * degeneracies do not affect point containment.)
   *
   * This predicate is only required for layers that will be assembled into
   * polygons.  It is not used by other layer types.
   */
  alias IsFullPolygonPredicate = bool function(in Graph g, out S2Error error);

  /// Default constructor; requires Init() to be called.
  this() {
    _labelSetLexicon = new IdSetLexicon();
  }

  /**
   * Convenience constructor that calls Init().  Note that to use the default
   * options, C++ syntax requires an extra layer of parentheses:
   *
   *   S2Builder builder((S2Builder::Options()));
   */
  this(Options options) {
    this();
    initialize(options);
  }

  /// Initializes an S2Builder with the given options.
  void initialize(Options options) {
    import std.stdio;
    import core.stdc.stdio;

    _options = options;
    const(SnapFunction) snap_function = options.snapFunction();
    S1Angle snap_radius = snap_function.snapRadius();
    debug enforce(snap_radius <= SnapFunction.kMaxSnapRadius());

    // Convert the snap radius to an S1ChordAngle.  This is the "true snap
    // radius" used when evaluating exact predicates (s2predicates.h).
    _siteSnapRadiusCa = S1ChordAngle(snap_radius);

    // When split_crossing_edges() is true, we need to use a larger snap radius
    // for edges than for vertices to ensure that both edges are snapped to the
    // edge intersection location.  This is because the computed intersection
    // point is not exact; it may be up to kIntersectionError away from its true
    // position.  The computed intersection point might then be snapped to some
    // other vertex up to snap_radius away.  So to ensure that both edges are
    // snapped to a common vertex, we need to increase the snap radius for edges
    // to at least the sum of these two values (calculated conservatively).
    S1Angle edge_snap_radius = snap_radius;
    if (!options.splitCrossingEdges()) {
      _edgeSnapRadiusCa = _siteSnapRadiusCa;
    } else {
      edge_snap_radius += INTERSECTION_ERROR;
      _edgeSnapRadiusCa = roundUp(edge_snap_radius);
    }
    _snappingRequested = edge_snap_radius > S1Angle.zero();

    // Compute the maximum distance that a vertex can be separated from an
    // edge while still affecting how that edge is snapped.
    _maxEdgeDeviation = snap_function.maxEdgeDeviation();
    _edgeSiteQueryRadiusCa =
        S1ChordAngle(_maxEdgeDeviation + snap_function.minEdgeVertexSeparation());

    // Compute the maximum edge length such that even if both endpoints move by
    // the maximum distance allowed (i.e., snap_radius), the center of the edge
    // will still move by less than max_edge_deviation().  This saves us a lot
    // of work since then we don't need to check the actual deviation.
    _minEdgeLengthToSplitCa =
        S1ChordAngle.fromRadians(2 * acos(sin(snap_radius) / sin(_maxEdgeDeviation)));

    // If the condition below is violated, then AddExtraSites() needs to be
    // modified to check that snapped edges pass on the same side of each "site
    // to avoid" as the input edge.  Currently it doesn't need to do this
    // because the condition below guarantees that if the snapped edge passes on
    // the wrong side of the site then it is also too close, which will cause a
    // separation site to be added.
    //
    // Currently max_edge_deviation() is at most 1.1 * snap_radius(), whereas
    // min_edge_vertex_separation() is at least 0.219 * snap_radius() (based on
    // S2CellIdSnapFunction, which is currently the worst case).
    enforce(snap_function.maxEdgeDeviation()
        <= snap_function.snapRadius() + snap_function.minEdgeVertexSeparation());

    // To implement idempotency, we check whether the input geometry could
    // possibly be the output of a previous S2Builder invocation.  This involves
    // testing whether any site/site or edge/site pairs are too close together.
    // This is done using exact predicates, which require converting the minimum
    // separation values to an S1ChordAngle.
    _minSiteSeparation = snap_function.minVertexSeparation();
    _minSiteSeparationCa = S1ChordAngle(_minSiteSeparation);
    _minEdgeSiteSeparationCa = S1ChordAngle(snap_function.minEdgeVertexSeparation());

    // This is an upper bound on the distance computed by S2ClosestPointQuery
    // where the true distance might be less than min_edge_site_separation_ca_.
    _minEdgeSiteSeparationCaLimit = addPointToEdgeError(_minEdgeSiteSeparationCa);

    // Compute the maximum possible distance between two sites whose Voronoi
    // regions touch.  (The maximum radius of each Voronoi region is
    // edge_snap_radius_.)  Then increase this bound to account for errors.
    _maxAdjacentSiteSeparationCa = addPointToPointError(roundUp(2 * edge_snap_radius));

    // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
    // squared distance between a vertex and an edge measured perpendicular to
    // the plane containing the edge, and increase this value by the maximum
    // error in the calculation to compare this distance against the bound.
    double d = sin(edge_snap_radius);
    _edgeSnapRadiusSin2 = d * d;
    _edgeSnapRadiusSin2 +=
        ((9.5 * d + 2.5 + 2 * sqrt(3.0)) * d + 9 * double.epsilon) * double.epsilon;

    // Initialize the current label set.
    _labelSetId = _labelSetLexicon.emptySetId();
    _labelSetModified = false;

    // If snapping was requested, we try to determine whether the input geometry
    // already meets the output requirements.  This is necessary for
    // idempotency, and can also save work.  If we discover any reason that the
    // input geometry needs to be modified, snapping_needed_ is set to true.
    _snappingNeeded = false;
  }

  const(Options) options() const {
    return _options;
  }

  /**
   * Starts a new output layer.  This method must be called before adding any
   * edges to the S2Builder.  You may call this method multiple times to build
   * multiple geometric objects that are snapped to the same set of sites.
   *
   * For example, if you have a set of contour lines, then you could put each
   * contour line in a separate layer.  This keeps the contour lines separate
   * from each other, while also ensuring that no crossing edges are created
   * when they are snapped and/or simplified.  (This is not true if the
   * contour lines are snapped or simplified independently.)
   *
   * Similarly, if you have a set of polygons that share common boundaries
   * (e.g., countries), you can snap and/or simplify them at the same time by
   * putting them in different layers, while ensuring that their boundaries
   * remain consistent (i.e., no crossing edges or T-vertices are introduced).
   *
   * Ownership of the layer is transferred to the S2Builder.  Example usage:
   *
   * S2Polyline line1, line2;
   * builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line1)));
   * ... Add edges using builder.AddEdge(), etc ...
   * builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line2)));
   * ... Add edges using builder.AddEdge(), etc ...
   * S2Error error;
   * CHECK(builder.Build(&error)) << error;  // Builds "line1" & "line2"
   */
  void startLayer(Layer layer) {
    _layerOptions ~= layer.graphOptions();
    _layerBegins ~= cast(int) _inputEdges.length;
    _layerIsFullPolygonPredicates ~= &isFullPolygonUnspecified;
    _layers ~= layer;
  }

  /// Adds a degenerate edge (representing a point) to the current layer.
  void addPoint(in S2Point v) {
    addEdge(v, v);
  }

  /// Adds the given edge to the current layer.
  void addEdge(in S2Point v0, in S2Point v1)
  in {
    assert(!_layers.empty(), "Call StartLayer before adding any edges");
  } body {
    if (v0 == v1
        && _layerOptions.back().degenerateEdges() == GraphOptions.DegenerateEdges.DISCARD) {
      return;
    }
    InputVertexId j0 = addVertex(v0);
    InputVertexId j1 = addVertex(v1);
    _inputEdges ~= [j0, j1];

    // If there are any labels, then attach them to this input edge.
    if (_labelSetModified) {
      if (_labelSetIds.empty()) {
        // Populate the missing entries with empty label sets.
        _labelSetIds.length = _inputEdges.length - 1;
        _labelSetIds.fill(_labelSetId);
      }
      _labelSetId = _labelSetLexicon.add(_labelSet);
      _labelSetIds ~= _labelSetId;
      _labelSetModified = false;
    } else if (!_labelSetIds.empty()) {
      _labelSetIds ~= _labelSetId;
    }
  }

  /**
   * Adds the edges in the given polyline.  (Note that if the polyline
   * consists of 0 or 1 vertices, this method does nothing.)
   */
  void addPolyline(in S2Polyline polyline) {
    const int n = polyline.numVertices();
    for (int i = 1; i < n; ++i) {
      addEdge(polyline.vertex(i - 1), polyline.vertex(i));
    }
  }

  /**
   * Adds the edges in the given loop.  If the sign() of the loop is negative
   * (i.e. this loop represents a hole within a polygon), the edge directions
   * are automatically reversed to ensure that the polygon interior is always
   * to the left of every edge.
   */
  void addLoop(in S2Loop loop) {
    // Ignore loops that do not have a boundary.
    if (loop.isEmptyOrFull()) return;

    // For loops that represent holes, we add the edge from vertex n-1 to vertex
    // n-2 first.  This is because these edges will be assembled into a
    // clockwise loop, which will eventually be normalized in S2Polygon by
    // calling S2Loop::Invert().  S2Loop::Invert() reverses the order of the
    // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we
    // need to build a clockwise loop with vertex order (n-1, n-2, ..., 0).
    // This is done by adding the edge (n-1, n-2) first, and then ensuring that
    // Build() assembles loops starting from edges in the order they were added.
    const int n = loop.numVertices();
    for (int i = 0; i < n; ++i) {
      addEdge(loop.orientedVertex(i), loop.orientedVertex(i + 1));
    }
  }

  /**
   * Adds the loops in the given polygon.  Loops representing holes have their
   * edge directions automatically reversed as described for AddLoop().  Note
   * that this method does not distinguish between the empty and full polygons,
   * i.e. adding a full polygon has the same effect as adding an empty one.
   */
  // TODO: Add when S2Polygon is complete.
  // void addPolygon(in S2Polygon polygon) {
  //   for (int i = 0; i < polygon.numLoops(); ++i) {
  //     addLoop(polygon.loop(i));
  //   }
  // }

  /// Adds the edges of the given shape to the current layer.
  void addShape(in S2Shape shape) {
    for (int e = 0, n = shape.numEdges(); e < n; ++e) {
      S2Shape.Edge edge = shape.edge(e);
      addEdge(edge.v0, edge.v1);
    }
  }

  /**
   * For layers that will be assembled into polygons, this method specifies a
   * predicate that will be called to determine whether the polygon is empty
   * or full except for the given degeneracies.  (See IsFullPolygonPredicate
   * above.)
   *
   * This method should be called at most once per layer; additional calls
   * simply overwrite the previous value for the current layer.
   *
   * The default implementation sets an appropriate error and returns false
   * (i.e., degenerate polygons are assumed to be empty).
   */
  void addIsFullPolygonPredicate(IsFullPolygonPredicate predicate) {
    _layerIsFullPolygonPredicates[$ - 1] = predicate;
  }

  /**
   * Forces a vertex to be located at the given position.  This can be used to
   * prevent certain input vertices from moving.  However if you are trying to
   * preserve part of the input boundary, be aware that this option does not
   * prevent edges from being split by new vertices.
   *
   * Forced vertices are never snapped; if this is desired then you need to
   * call options().snap_function().SnapPoint() explicitly.  Forced vertices
   * are also never simplified away (if simplify_edge_chains() is used).
   *
   * Caveat: Since this method can place vertices arbitrarily close together,
   * S2Builder makes no minimum separation guaranteees with forced vertices.
   */
  void forceVertex(in S2Point vertex) {
    _sites ~= vertex;
  }

  /**
   * Every edge can have a set of non-negative integer labels attached to it.
   * When used with an appropriate layer type, you can then retrieve the
   * labels associated with each output edge.  This can be useful when merging
   * or combining data from several sources.  (Note that in many cases it is
   * easier to use separate output layers rather than labels.)
   *
   * Labels are 32-bit non-negative integers.  To support other label types,
   * you can use ValueLexicon to store the set of unique labels seen so far:
   *
   *   ValueLexicon<MyLabel> my_label_lexicon;
   *   builder.set_label(my_label_lexicon.Add(label));
   *
   * The current set of labels is represented as a stack.  This makes it easy
   * to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
   * set_label() and clear_labels() if you need at most one label per edge.
   */
  alias Label = int;

  /// Clear the stack of labels.
  void clearLabels() {
    _labelSet.length = 0;
    _labelSetModified = true;
  }

  /**
   * Add a label to the stack.
   * REQUIRES: label >= 0.
   */
  void pushLabel(Label label)
  in {
    assert(label >= 0);
  } body {
    _labelSet ~= label;
    _labelSetModified = true;
  }

  /// Remove a label from the stack.
  void popLabel() {
    _labelSet.popBack();
    _labelSetModified = true;
  }

  /**
   * Convenience function that clears the stack and adds a single label.
   * REQUIRES: label >= 0.
   */
  void setLabel(Label label)
  in {
    assert(label >= 0);
  } body {
    _labelSet.length = 1;
    _labelSet[0] = label;
    _labelSetModified = true;
  }

  /**
   * Performs the requested edge splitting, snapping, simplification, etc, and
   * then assembles the resulting edges into the requested output layers.
   *
   * Returns true if all edges were assembled; otherwise sets "error"
   * appropriately.  Depending on the error, some or all output layers may
   * have been created.  Automatically resets the S2Builder state so that it
   * can be reused.
   */
  bool build(out S2Error error) {
    error.clear();
    _error = error;

    // Mark the end of the last layer.
    _layerBegins ~= cast(int) _inputEdges.length;

    // See the algorithm overview at the top of this file.
    if (_snappingRequested && !_options.idempotent()) {
      _snappingNeeded = true;
    }
    chooseSites();
    buildLayers();
    reset();
    return _error.ok();
  }

  /**
   * Clears all input data and resets the builder state.  Any options
   * specified are preserved.
   */
  void reset() {
    _inputVertices.length = 0;
    _inputEdges.length = 0;
    _layers.length = 0;
    _layerOptions.length = 0;
    _layerBegins.length = 0;
    _labelSetIds.length = 0;
    _labelSetLexicon.clear();
    _labelSet.length = 0;
    _labelSetModified = false;
    _sites.length = 0;
    _edgeSites.length = 0;
    _snappingNeeded = false;
  }

  /// Identifies an input vertex.
  alias InputVertexId = int;

  /// Identifies an input edge.
  alias InputEdgeId = int;

  /// Identifies the set of input edge ids that were snapped to a given edge.
  alias InputEdgeIdSetId = int;

  alias LabelSetId = int;

private:
  //////////////////////  Input Types  /////////////////////////
  // All types associated with the S2Builder inputs are prefixed with "Input".

  /// Defines an input edge.
  alias InputEdge = InputVertexId[2];

  /**
   * Sort key for prioritizing input vertices.  (Note that keys are *not*
   * compared using std::less; see SortInputVertices for details.)
   */
  struct InputVertexKey {
    S2CellId first;
    InputVertexId second;
  }

  //////////////////////  Output Types  /////////////////////////
  // These types define the output vertices and edges.

  /**
   * Identifies a snapped vertex ("snap site").  If there is only one layer,
   * than SiteId is the same as Graph::VertexId, but if there are many layers
   * then each Graph may contain only a subset of the sites.  Also see
   * GraphOptions::allow_vertex_filtering().
   */
  alias SiteId = int;

  /// Defines an output edge.
  alias Edge = Tuple!(SiteId, SiteId);

  /// Identifies an output edge.
  alias EdgeId = int;

  /// Identifies an output edge in a particular layer.
  struct LayerEdgeId {
    int first;
    EdgeId second;

    int opCmp(LayerEdgeId other) const {
      return (first != other.first)
          ? first - other.first
          : second - other.second;
    }
  }

  /**
   * Input vertices are stored in a vector, with some removal of duplicates.
   * Edges are represented as (VertexId, VertexId) pairs.  All edges are stored
   * in a single vector; each layer corresponds to a contiguous range.
   */
  InputVertexId addVertex(in S2Point v) {
    // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
    // to do anything more sophisticated, either use a ValueLexicon, or sort the
    // vertices once they have all been added, remove duplicates, and update the
    // edges.
    if (_inputVertices.empty() || v != _inputVertices.back()) {
      _inputVertices ~= v;
    }
    return cast(InputVertexId) _inputVertices.length - 1;
  }

  void chooseSites() {
    if (_inputVertices.empty()) return;

    auto input_edge_index = new MutableS2ShapeIndex();
    input_edge_index.add(new VertexIdEdgeVectorShape(_inputEdges, _inputVertices));
    if (_options.splitCrossingEdges()) {
      addEdgeCrossings(input_edge_index);
    }
    if (_snappingRequested) {
      S2PointIndex!SiteId site_index;
      addForcedSites(site_index);
      chooseInitialSites(site_index);
      collectSiteEdges(site_index);
    }
    if (_snappingNeeded) {
      addExtraSites(input_edge_index);
    } else {
      copyInputEdges();
    }
  }

  /**
   * Sort the input vertices, discard duplicates, and update the input edges
   * to refer to the pruned vertex list.  (We sort in the same order used by
   * ChooseInitialSites() to avoid inconsistencies in tests.)
   */
  void copyInputEdges() {
    InputVertexKey[] sorted = sortInputVertices();
    auto vmap = new InputVertexId[_inputVertices.length];
    _sites.length = 0;
    _sites.reserve(_inputVertices.length);
    for (int ind = 0; ind < sorted.length; ) {
      const S2Point site = _inputVertices[sorted[ind].second];
      vmap[sorted[ind].second] = cast(InputVertexId) _sites.length;
      while (++ind < sorted.length && _inputVertices[sorted[ind].second] == site) {
        vmap[sorted[ind].second] = cast(InputVertexId) _sites.length;
      }
      _sites ~= site;
    }
    _inputVertices = _sites;
    for (int i = 0; i < _inputEdges.length; ++i) {
      InputEdge* e = &(_inputEdges[i]);
      (*e)[0] = vmap[(*e)[0]];
      (*e)[1] = vmap[(*e)[1]];
    }
  }

  InputVertexKey[] sortInputVertices() {
    // Sort all the input vertices in the order that we wish to consider them as
    // candidate Voronoi sites.  Any sort order will produce correct output, so
    // we have complete flexibility in choosing the sort key.  We could even
    // leave them unsorted, although this would have the disadvantage that
    // changing the order of the input edges could cause S2Builder to snap to a
    // different set of Voronoi sites.
    //
    // We have chosen to sort them primarily by S2CellId since this improves the
    // performance of many S2Builder phases (due to better spatial locality).
    // It also allows the possibility of replacing the current S2PointIndex
    // approach with a more efficient recursive divide-and-conquer algorithm.
    //
    // However, sorting by leaf S2CellId alone has two small disadvantages in
    // the case where the candidate sites are densely spaced relative to the
    // snap radius (e.g., when using the IdentitySnapFunction, or when snapping
    // to E6/E7 near the poles, or snapping to S2CellId/E6/E7 using a snap
    // radius larger than the minimum value required):
    //
    //  - First, it tends to bias the Voronoi site locations towards points that
    //    are earlier on the S2CellId Hilbert curve.  For example, suppose that
    //    there are two parallel rows of input vertices on opposite sides of the
    //    edge between two large S2Cells, and the rows are separated by less
    //    than the snap radius.  Then only vertices from the cell with the
    //    smaller S2CellId are selected, because they are considered first and
    //    prevent us from selecting the sites from the other cell (because they
    //    are closer than "snap_radius" to an existing site).
    //
    //  - Second, it tends to choose more Voronoi sites than necessary, because
    //    at each step we choose the first site along the Hilbert curve that is
    //    at least "snap_radius" away from all previously selected sites.  This
    //    tends to yield sites whose "coverage discs" overlap quite a bit,
    //    whereas it would be better to cover all the input vertices with a
    //    smaller set of coverage discs that don't overlap as much.  (This is
    //    the "geometric set cover problem", which is NP-hard.)
    //
    // It is not worth going to much trouble to fix these problems, because they
    // really aren't that important (and don't affect the guarantees made by the
    // algorithm), but here are a couple of heuristics that might help:
    //
    // 1. Sort the input vertices by S2CellId at a coarse level (down to cells
    // that are O(snap_radius) in size), and then sort by a fingerprint of the
    // S2Point coordinates (i.e., quasi-randomly).  This would retain most of
    // the advantages of S2CellId sorting, but makes it more likely that we will
    // select sites that are further apart.
    //
    // 2. Rather than choosing the first uncovered input vertex and snapping it
    // to obtain the next Voronoi site, instead look ahead through following
    // candidates in S2CellId order and choose the furthest candidate whose
    // snapped location covers all previous uncovered input vertices.
    //
    // TODO(ericv): Experiment with these approaches.
    import std.algorithm : sort;

    InputVertexKey[] keys;
    keys.reserve(_inputVertices.length);
    for (InputVertexId i = 0; i < _inputVertices.length; ++i) {
      keys ~= InputVertexKey(S2CellId(_inputVertices[i]), i);
    }
    sort!((a, b) {
          if (a.first < b.first) return true;
          if (b.first < a.first) return false;
          return _inputVertices[a.second] < _inputVertices[b.second];
        })(keys);
    return keys;
  }

  /**
   * Check all edge pairs for crossings, and add the corresponding intersection
   * points to input_vertices_.  (The intersection points will be snapped and
   * merged with the other vertices during site selection.)
   */
  void addEdgeCrossings(MutableS2ShapeIndex input_edge_index) {
    import s2.s2crossing_edge_query : CrossingType;
    import s2.s2edge_crossings : getIntersection;
    import s2.shapeutil.shape_edge : ShapeEdge;
    import s2.shapeutil.visit_crossing_edge_pairs;

    // We need to build a list of intersections and add them afterwards so that
    // we don't reallocate vertices_ during the VisitCrossings() call.
    S2Point[] new_vertices;
    visitCrossingEdgePairs(
        input_edge_index, CrossingType.INTERIOR,
        (in ShapeEdge a, in ShapeEdge b, bool) {
          new_vertices ~= getIntersection(a.v0(), a.v1(), b.v0(), b.v1());
          return true;  // Continue visiting.
        });
    if (!new_vertices.empty()) {
      _snappingNeeded = true;
      foreach (const vertex; new_vertices) addVertex(vertex);
    }
  }

  void addForcedSites(S2PointIndex!SiteId site_index) {
    import std.array, std.algorithm;

    // Sort the forced sites and remove duplicates.
    _sites = _sites.sort.uniq.array;
    // Add the forced sites to the index.
    for (SiteId id = 0; id < _sites.length; ++id) {
      site_index.add(_sites[id], id);
    }
    _numForcedSites = cast(SiteId) _sites.length;
  }

  bool isForced(SiteId v) const {
    return v < _numForcedSites;
  }

  void chooseInitialSites(S2PointIndex!SiteId site_index) {
    import s2.s2predicates : compareDistance;
    // Find all points whose distance is <= min_site_separation_ca_.
    auto options = new S2ClosestPointQueryOptions();
    options.setConservativeMaxDistance(_minSiteSeparationCa);
    auto site_query = new S2ClosestPointQuery!SiteId(site_index, options);
    S2ClosestPointQuery!(SiteId).Result[] results;

    // Apply the snap_function() to each input vertex, then check whether any
    // existing site is closer than min_vertex_separation().  If not, then add a
    // new site.
    //
    // NOTE(ericv): There are actually two reasonable algorithms, which we call
    // "snap first" (the one above) and "snap last".  The latter checks for each
    // input vertex whether any existing site is closer than snap_radius(), and
    // only then applies the snap_function() and adds a new site.  "Snap last"
    // can yield slightly fewer sites in some cases, but it is also more
    // expensive and can produce surprising results.  For example, if you snap
    // the polyline "0:0, 0:0.7" using IntLatLngSnapFunction(0), the result is
    // "0:0, 0:0" rather than the expected "0:0, 0:1", because the snap radius
    // is approximately sqrt(2) degrees and therefore it is legal to snap both
    // input points to "0:0".  "Snap first" produces "0:0, 0:1" as expected.
    InputVertexKey[] sorted = sortInputVertices();
    for (int i = 0; i < sorted.length; ++i) {
      const S2Point vertex = _inputVertices[sorted[i].second];
      S2Point site = snapSite(vertex);
      // If any vertex moves when snapped, the output cannot be idempotent.
      _snappingNeeded = _snappingNeeded || site != vertex;

      // FindClosestPoints() measures distances conservatively, so we need to
      // recheck the distances using exact predicates.
      //
      // NOTE(ericv): When the snap radius is large compared to the average
      // vertex spacing, we could possibly avoid the call the FindClosestPoints
      // by checking whether sites_.back() is close enough.
      auto target = new S2ClosestPointQueryPointTarget(site);
      site_query.findClosestPoints(target, results);
      bool add_site = true;
      foreach (const result; results) {
        if (compareDistance(site, result.point(), _minSiteSeparationCa) <= 0) {
          add_site = false;
          // This pair of sites is too close.  If the sites are distinct, then
          // the output cannot be idempotent.
          _snappingNeeded = _snappingNeeded || site != result.point();
        }
      }
      if (add_site) {
        site_index.add(site, cast(SiteId) _sites.length);
        _sites ~= site;
        site_query.reInitialize();
      }
    }
  }

  S2Point snapSite(in S2Point point) {
    if (!_snappingRequested) return point;
    S2Point site = _options.snapFunction().snapPoint(point);
    auto dist_moved = S1ChordAngle(site, point);
    if (dist_moved > _siteSnapRadiusCa) {
      _error.initialize(
          S2Error.Code.BUILDER_SNAP_RADIUS_TOO_SMALL,
          "Snap function moved vertex (%.15g, %.15g, %.15g) "
              ~ "by %.15g, which is more than the specified snap "
              ~ "radius of %.15g",
          point.x(), point.y(), point.z(),
          dist_moved.toS1Angle().radians(),
          _siteSnapRadiusCa.toS1Angle().radians());
    }
    return site;
  }

  /**
   * For each edge, find all sites within min_edge_site_query_radius_ca_ and
   * store them in edge_sites_.  Also, to implement idempotency this method also
   * checks whether the input vertices and edges may already satisfy the output
   * criteria.  If any problems are found then snapping_needed_ is set to true.
   */
  void collectSiteEdges(S2PointIndex!SiteId site_index) {
    import s2.s2predicates : compareEdgeDistance;
    // Find all points whose distance is <= edge_site_query_radius_ca_.
    auto options = new S2ClosestPointQueryOptions();
    options.setConservativeMaxDistance(_edgeSiteQueryRadiusCa);
    auto site_query = new S2ClosestPointQuery!SiteId(site_index, options);
    S2ClosestPointQuery!(SiteId).Result[] results;
    _edgeSites.length = _inputEdges.length;
    for (InputEdgeId e = 0; e < _inputEdges.length; ++e) {
      const InputEdge edge = _inputEdges[e];
      const S2Point v0 = _inputVertices[edge[0]];
      const S2Point v1 = _inputVertices[edge[1]];
      if (s2builderVerbose) {
        writeln("S2Polyline: ", v0, ", ", v1, "\n");
      }
      auto target = new S2ClosestPointQueryEdgeTarget(v0, v1);
      site_query.findClosestPoints(target, results);
      auto sites = _edgeSites[e];
      sites.reserve(results.length);
      foreach (const result; results) {
        sites ~= result.data();
        if (!_snappingNeeded
            && result.distance() < _minEdgeSiteSeparationCaLimit
            && result.point() != v0
            && result.point() != v1
            && compareEdgeDistance(result.point(), v0, v1, _minEdgeSiteSeparationCa) < 0) {
          _snappingNeeded = true;
        }
      }
      sortSitesByDistance(v0, sites);
    }
  }

  void sortSitesByDistance(in S2Point x, ref SiteId[] sites) const {
    import s2.s2predicates : compareDistances;
    // Sort sites in increasing order of distance to X.
    sort!((SiteId i, SiteId j) => compareDistances(x, _sites[i], _sites[j]) < 0)(sites);
  }

  /**
   * There are two situatons where we need to add extra Voronoi sites in order
   * to ensure that the snapped edges meet the output requirements:
   *
   *  (1) If a snapped edge deviates from its input edge by more than
   *      max_edge_deviation(), we add a new site on the input edge near the
   *      middle of the snapped edge.  This causes the snapped edge to split
   *      into two pieces, so that it follows the input edge more closely.
   *
   *  (2) If a snapped edge is closer than min_edge_vertex_separation() to any
   *      nearby site (the "site to avoid"), then we add a new site (the
   *      "separation site") on the input edge near the site to avoid.  This
   *      causes the snapped edge to follow the input edge more closely and is
   *      guaranteed to increase the separation to the required distance.
   *
   * We check these conditions by snapping all the input edges to a chain of
   * Voronoi sites and then testing each edge in the chain.  If a site needs to
   * be added, we mark all nearby edges for re-snapping.
   */
  void addExtraSites(MutableS2ShapeIndex input_edge_index) {
    // When options_.split_crossing_edges() is true, this function may be called
    // even when site_snap_radius_ca_ == 0 (because edge_snap_radius_ca_ > 0).
    // However neither of the conditions above apply in that case.
    if (_siteSnapRadiusCa == S1ChordAngle.zero()) return;

    SiteId[] chain;  // Temporary
    InputEdgeId[] snap_queue;
    for (InputEdgeId max_e = 0; max_e < _inputEdges.length; ++max_e) {
      snap_queue ~= max_e;
      while (!snap_queue.empty()) {
        InputEdgeId e = snap_queue.back();
        snap_queue.popBack();
        snapEdge(e, chain);
        // We could save the snapped chain here in a snapped_chains_ vector, to
        // avoid resnapping it in AddSnappedEdges() below, however currently
        // SnapEdge only accounts for less than 5% of the runtime.
        maybeAddExtraSites(e, max_e, chain, input_edge_index, snap_queue);
      }
    }
  }

  void maybeAddExtraSites(
      in InputEdgeId edge_id,
      in InputEdgeId max_edge_id,
      in SiteId[] chain,
      MutableS2ShapeIndex input_edge_index,
      ref InputEdgeId[] snap_queue) {
    import s2.s2edge_distances : isEdgeBNearEdgeA, project;
    import s2.s2predicates : compareEdgeDistance;
    // The snapped chain is always a *subsequence* of the nearby sites
    // (edge_sites_), so we walk through the two arrays in parallel looking for
    // sites that weren't snapped.  We also keep track of the current snapped
    // edge, since it is the only edge that can be too close.
    int i = 0;
    foreach (SiteId id; _edgeSites[edge_id]) {
      if (id == chain[i]) {
        if (++i == chain.length) break;
        // Check whether this snapped edge deviates too far from its original
        // position.  If so, we split the edge by adding an extra site.
        const S2Point v0 = _sites[chain[i - 1]];
        const S2Point v1 = _sites[chain[i]];
        if (S1ChordAngle(v0, v1) < _minEdgeLengthToSplitCa) continue;

        const InputEdge edge = _inputEdges[edge_id];
        const S2Point a0 = _inputVertices[edge[0]];
        const S2Point a1 = _inputVertices[edge[1]];
        if (!isEdgeBNearEdgeA(a0, a1, v0, v1, _maxEdgeDeviation)) {
          // Add a new site on the input edge, positioned so that it splits the
          // snapped edge into two approximately equal pieces.  Then we find all
          // the edges near the new site (including this one) and add them to
          // the snap queue.
          //
          // Note that with large snap radii, it is possible that the snapped
          // edge wraps around the sphere the "wrong way".  To handle this we
          // find the preferred split location by projecting both endpoints onto
          // the input edge and taking their midpoint.
          S2Point mid = (project(v0, a0, a1) + project(v1, a0, a1)).normalize();
          S2Point new_site = getSeparationSite(mid, v0, v1, edge_id);
          addExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
          return;
        }
      } else if (i > 0 && id >= _numForcedSites) {
        // Check whether this "site to avoid" is closer to the snapped edge than
        // min_edge_vertex_separation().  Note that this is the only edge of the
        // chain that can be too close because its vertices must span the point
        // where "site_to_avoid" projects onto the input edge XY (this claim
        // relies on the fact that all sites are separated by at least the snap
        // radius).  We don't try to avoid sites added using ForceVertex()
        // because we don't guarantee any minimum separation from such sites.
        const S2Point site_to_avoid = _sites[id];
        const S2Point v0 = _sites[chain[i - 1]];
        const S2Point v1 = _sites[chain[i]];
        if (compareEdgeDistance(site_to_avoid, v0, v1, _minEdgeSiteSeparationCa) < 0) {
          // A snapped edge can only approach a site too closely when there are
          // no sites near the input edge near that point.  We fix that by
          // adding a new site along the input edge (a "separation site"), then
          // we find all the edges near the new site (including this one) and
          // add them to the snap queue.
          S2Point new_site = getSeparationSite(site_to_avoid, v0, v1, edge_id);
          debug enforce(site_to_avoid != new_site);
          addExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
          return;
        }
      }
    }
  }

  /**
   * Adds a new site, then updates "edge_sites"_ for all edges near the new site
   * and adds them to "snap_queue" for resnapping (unless their edge id exceeds
   * "max_edge_id", since those edges have not been snapped the first time yet).
   */
  void addExtraSite(
      in S2Point new_site, InputEdgeId max_edge_id,
      MutableS2ShapeIndex input_edge_index, ref InputEdgeId[] snap_queue) {
    SiteId new_site_id = cast(SiteId) _sites.length;
    _sites ~= new_site;
    // Find all edges whose distance is <= edge_site_query_radius_ca_.
    auto options = new S2ClosestEdgeQuery.Options();
    options.setConservativeMaxDistance(_edgeSiteQueryRadiusCa);
    options.setIncludeInteriors(false);
    auto query = new S2ClosestEdgeQuery(input_edge_index, options);
    auto target = new S2ClosestEdgeQuery.PointTarget(new_site);
    foreach (const result; query.findClosestEdges(target)) {
      InputEdgeId e = result.edgeId;
      SiteId[]* site_ids = &_edgeSites[e];  // Use a pointer so we can modify the original.
      *site_ids ~= new_site_id;
      sortSitesByDistance(_inputVertices[_inputEdges[e][0]], *site_ids);
      if (e <= max_edge_id) snap_queue ~= e;
    }
  }

  S2Point getSeparationSite(
      in S2Point site_to_avoid, in S2Point v0, in S2Point v1, InputEdgeId input_edge_id) {
    import s2.s2pointutil : robustCrossProd;
    import s2.s2edge_distances : project;
    // Define the "coverage disc" of a site S to be the disc centered at S with
    // radius "snap_radius".  Similarly, define the "coverage interval" of S for
    // an edge XY to be the intersection of XY with the coverage disc of S.  The
    // SnapFunction implementations guarantee that the only way that a snapped
    // edge can be closer than min_edge_vertex_separation() to a non-snapped
    // site (i.e., site_to_avoid) if is there is a gap in the coverage of XY
    // near this site.  We can fix this problem simply by adding a new site to
    // fill this gap, located as closely as possible to the site to avoid.
    //
    // To calculate the coverage gap, we look at the two snapped sites on
    // either side of site_to_avoid, and find the endpoints of their coverage
    // intervals.  The we place a new site in the gap, located as closely as
    // possible to the site to avoid.  Note that the new site may move when it
    // is snapped by the snap_function, but it is guaranteed not to move by
    // more than snap_radius and therefore its coverage interval will still
    // intersect the gap.
    const InputEdge edge = _inputEdges[input_edge_id];
    const S2Point x = _inputVertices[edge[0]];
    const S2Point y = _inputVertices[edge[1]];
    Vector3_d xy_dir = y - x;
    S2Point n = robustCrossProd(x, y);
    S2Point new_site = project(site_to_avoid, x, y, n);
    S2Point gap_min = getCoverageEndpoint(v0, x, y, n);
    S2Point gap_max = getCoverageEndpoint(v1, y, x, -n);
    if ((new_site - gap_min).dotProd(xy_dir) < 0) {
      new_site = gap_min;
    } else if ((gap_max - new_site).dotProd(xy_dir) < 0) {
      new_site = gap_max;
    }
    new_site = snapSite(new_site);
    debug enforce(v0 != new_site);
    debug enforce(v1 != new_site);
    return new_site;
  }

  /**
   * Given a site P and an edge XY with normal N, intersect XY with the disc of
   * radius snap_radius() around P, and return the intersection point that is
   * further along the edge XY toward Y.
   */
  S2Point getCoverageEndpoint(in S2Point p, in S2Point x, in S2Point y, in S2Point n) const {
    // Consider the plane perpendicular to P that cuts off a spherical cap of
    // radius snap_radius().  This plane intersects the plane through the edge
    // XY (perpendicular to N) along a line, and that line intersects the unit
    // sphere at two points Q and R, and we want to return the point R that is
    // further along the edge XY toward Y.
    //
    // Let M be the midpoint of QR.  This is the point along QR that is closest
    // to P.  We can now express R as the sum of two perpendicular vectors OM
    // and MR in the plane XY.  Vector MR is in the direction N x P, while
    // vector OM is in the direction (N x P) x N, where N = X x Y.
    //
    // The length of OM can be found using the Pythagorean theorem on triangle
    // OPM, and the length of MR can be found using the Pythagorean theorem on
    // triangle OMR.
    //
    // In the calculations below, we save some work by scaling all the vectors
    // by n.CrossProd(p).Norm2(), and normalizing at the end.
    double n2 = n.norm2();
    double nDp = n.dotProd(p);
    S2Point nXp = n.crossProd(p);
    S2Point nXpXn = n2 * p - nDp * n;
    Vector3_d om = sqrt(1 - _edgeSnapRadiusSin2) * nXpXn;
    double mr2 = _edgeSnapRadiusSin2 * n2 - nDp * nDp;

    // MR is constructed so that it points toward Y (rather than X).
    Vector3_d mr = sqrt(max(0.0, mr2)) * nXp;
    return (om + mr).normalize();
  }

  void snapEdge(InputEdgeId e, out SiteId[] chain) const {
    const InputEdge edge = _inputEdges[e];
    if (!_snappingNeeded) {
      chain ~= edge[0];
      chain ~= edge[1];
      return;
    }

    const S2Point x = _inputVertices[edge[0]];
    const S2Point y = _inputVertices[edge[1]];

    // Optimization: if there is only one nearby site, return.
    // Optimization: if there are exactly two nearby sites, and one is close
    // enough to each vertex, then return.

    // Now iterate through the sites.  We keep track of the sequence of sites
    // that are visited.
    const candidates = _edgeSites[e];
    for (int i = 0; i < candidates.length; ++i) {
      const S2Point c = _sites[candidates[i]];
      // Skip any sites that are too far away.  (There will be some of these,
      // because we also keep track of "sites to avoid".)  Note that some sites
      // may be close enough to the line containing the edge, but not to the
      // edge itself, so we can just use the dot product with the edge normal.
      if (compareEdgeDistance(c, x, y, _edgeSnapRadiusCa) > 0) {
        continue;
      }
      // Check whether the new site C excludes the previous site B.  If so,
      // repeat with the previous site, and so on.
      bool add_site_c = true;
      for (; !chain.empty(); chain.popBack()) {
        S2Point b = _sites[chain.back()];

        // First, check whether B and C are so far apart that their clipped
        // Voronoi regions can't intersect.
        auto bc = S1ChordAngle(b, c);
        if (bc >= _maxAdjacentSiteSeparationCa) break;

        // Otherwise, we want to check whether site C prevents the Voronoi
        // region of B from intersecting XY, or vice versa.  This can be
        // determined by computing the "coverage interval" (the segment of XY
        // intersected by the coverage disc of radius snap_radius) for each
        // site.  If the coverage interval of one site contains the coverage
        // interval of the other, then the contained site can be excluded.
        Excluded result = getVoronoiSiteExclusion(b, c, x, y, _edgeSnapRadiusCa);
        if (result == Excluded.FIRST) continue;  // Site B excluded by C
        if (result == Excluded.SECOND) {
          add_site_c = false;  // Site C is excluded by B.
          break;
        }
        debug enforce(Excluded.NEITHER == result);

        // Otherwise check whether the previous site A is close enough to B and
        // C that it might further clip the Voronoi region of B.
        if (chain.length < 2) break;
        S2Point a = _sites[chain[$ - 2]];
        auto ac = S1ChordAngle(a, c);
        if (ac >= _maxAdjacentSiteSeparationCa) break;

        // If triangles ABC and XYB have the same orientation, the circumcenter
        // Z of ABC is guaranteed to be on the same side of XY as B.
        int xyb = sign(x, y, b);
        if (sign(a, b, c) == xyb) {
          break;  // The circumcenter is on the same side as B but further away.
        }
        // Other possible optimizations:
        //  - if AB > max_adjacent_site_separation_ca_ then keep B.
        //  - if d(B, XY) < 0.5 * min(AB, BC) then keep B.

        // If the circumcenter of ABC is on the same side of XY as B, then B is
        // excluded by A and C combined.  Otherwise B is needed and we can exit.
        if (edgeCircumcenterSign(x, y, a, b, c) != xyb) break;
      }
      if (add_site_c) {
        chain ~= candidates[i];
      }
    }
    debug enforce(!chain.empty());
    debug {
      for (int i = 0; i < candidates.length; ++i) {
        if (compareDistances(y, _sites[chain.back()], _sites[candidates[i]]) > 0) {
          logger.logError("Snapping invariant broken!");
        }
      }
    }
    if (s2builderVerbose) {
      writeln("(", edge[0], ",", edge[1], "): ");
      foreach (SiteId id; chain) write(id, " ");
      writeln();
    }
  }

  void buildLayers() {
    // Each output edge has an "input edge id set id" (an int32) representing
    // the set of input edge ids that were snapped to this edge.  The actual
    // InputEdgeIds can be retrieved using "input_edge_id_set_lexicon".
    Edge[][] layer_edges;
    InputEdgeIdSetId[][] layer_input_edge_ids;
    IdSetLexicon input_edge_id_set_lexicon = new IdSetLexicon();
    buildLayerEdges(layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);

    // At this point we have no further need for the input geometry or nearby
    // site data, so we clear those fields to save space.
    _inputVertices.length = 0;
    _inputEdges.length = 0;
    _edgeSites.length = 0;

    // If there are a large number of layers, then we build a minimal subset of
    // vertices for each layer.  This ensures that layer types that iterate over
    // vertices will run in time proportional to the size of that layer rather
    // than the size of all layers combined.
    S2Point[][] layer_vertices;
    enum int kMinLayersForVertexFiltering = 10;
    if (_layers.length >= kMinLayersForVertexFiltering) {
      // Disable vertex filtering if it is disallowed by any layer.  (This could
      // be optimized, but in current applications either all layers allow
      // filtering or none of them do.)
      bool allow_vertex_filtering = false;
      foreach (const options; _layerOptions) {
        allow_vertex_filtering &= options.allowVertexFiltering();
      }
      if (allow_vertex_filtering) {
        Graph.VertexId[] filter_tmp;  // Temporary used by FilterVertices.
        layer_vertices.length = _layers.length;
        for (int i = 0; i < _layers.length; ++i) {
          layer_vertices[i] = Graph.filterVertices(_sites, layer_edges[i], filter_tmp);
        }
        _sites.length = 0;  // Release memory
      }
    }
    for (int i = 0; i < _layers.length; ++i) {
      const S2Point[] vertices = layer_vertices.empty() ? _sites : layer_vertices[i];
      auto graph = new Graph(_layerOptions[i], vertices, layer_edges[i],
          layer_input_edge_ids[i], input_edge_id_set_lexicon,
          _labelSetIds, _labelSetLexicon,
          _layerIsFullPolygonPredicates[i]);
      _layers[i].build(graph, _error);
      // Don't free the layer data until all layers have been built, in order to
      // support building multiple layers at once (e.g. ClosedSetNormalizer).
    }
  }

  /**
   * Snaps and possibly simplifies the edges for each layer, populating the
   * given output arguments.  The resulting edges can be used to construct an
   * S2Builder::Graph directly (no further processing is necessary).
   *
   * This method is not "const" because Graph::ProcessEdges can modify
   * layer_options_ in some cases (changing undirected edges to directed ones).
   */
  void buildLayerEdges(ref Edge[][] layer_edges, ref InputEdgeIdSetId[][] layer_input_edge_ids,
      ref IdSetLexicon input_edge_id_set_lexicon) {
    // Edge chains are simplified only when a non-zero snap radius is specified.
    // If so, we build a map from each site to the set of input vertices that
    // snapped to that site.
    InputVertexId[][] site_vertices;
    bool simplify = _snappingNeeded && _options.simplifyEdgeChains();
    if (simplify) {
      site_vertices.length = _sites.length;
    }

    layer_edges.length = _layers.length;
    layer_input_edge_ids.length = _layers.length;
    for (int i = 0; i < _layers.length; ++i) {
      addSnappedEdges(_layerBegins[i], _layerBegins[i + 1], _layerOptions[i],
          layer_edges[i], layer_input_edge_ids[i],
          input_edge_id_set_lexicon, site_vertices);
    }
    if (simplify) {
      simplifyEdgeChains(
          site_vertices, layer_edges, layer_input_edge_ids, input_edge_id_set_lexicon);
    }
    // We simplify edge chains before processing the per-layer GraphOptions
    // because simplification can create duplicate edges and/or sibling edge
    // pairs which may need to be removed.
    for (int i = 0; i < _layers.length; ++i) {
      // The errors generated by ProcessEdges are really warnings, so we simply
      // record them and continue.
      Graph.processEdges(_layerOptions[i], layer_edges[i],
          layer_input_edge_ids[i], input_edge_id_set_lexicon, _error);
    }
  }

  /**
   * Snaps all the input edges for a given layer, populating the given output
   * arguments.  If (*site_vertices) is non-empty then it is updated so that
   * (*site_vertices)[site] contains a list of all input vertices that were
   * snapped to that site.
   */
  void addSnappedEdges(
      InputEdgeId begin, InputEdgeId end, in GraphOptions options,
      ref Edge[] edges, ref InputEdgeIdSetId[] input_edge_ids,
      IdSetLexicon input_edge_id_set_lexicon,
      ref InputVertexId[][] site_vertices) const {
    bool discard_degenerate_edges =
        options.degenerateEdges() == GraphOptions.DegenerateEdges.DISCARD;
    SiteId[] chain;
    for (InputEdgeId e = begin; e < end; ++e) {
      InputEdgeIdSetId id = input_edge_id_set_lexicon.addSingleton(e);
      snapEdge(e, chain);
      maybeAddInputVertex(_inputEdges[e][0], chain[0], site_vertices);
      if (chain.length == 1) {
        if (discard_degenerate_edges) continue;
        addSnappedEdge(chain[0], chain[0], id, options.edgeType(), edges, input_edge_ids);
      } else {
        maybeAddInputVertex(_inputEdges[e][1], chain.back(), site_vertices);
        for (int i = 1; i < chain.length; ++i) {
          addSnappedEdge(chain[i - 1], chain[i], id, options.edgeType(), edges, input_edge_ids);
        }
      }
    }
    if (s2builderVerbose) dumpEdges(edges, _sites);
  }

  /**
   * If "site_vertices" is non-empty, ensures that (*site_vertices)[id] contains
   * "v".  Duplicate entries are allowed.
   */
  void maybeAddInputVertex(InputVertexId v, SiteId id, ref InputVertexId[][] site_vertices) const {
    if (site_vertices.empty) return;

    // Optimization: check if we just added this vertex.  This is worthwhile
    // because the input edges usually form a continuous chain, i.e. the
    // destination of one edge is the same as the source of the next edge.
    InputVertexId[]* vertices = &site_vertices[id];
    if ((*vertices).empty() || (*vertices).back() != v) {
      (*vertices) ~= v;
    }
  }

  /**
   * Adds the given edge to "edges" and "input_edge_ids".  If undirected edges
   * are being used, also adds an edge in the opposite direction.
   */
  void addSnappedEdge(SiteId src, SiteId dst, InputEdgeIdSetId id,
      EdgeType edge_type, ref Edge[] edges, ref InputEdgeIdSetId[] input_edge_ids) const {
    edges ~= Edge(src, dst);
    input_edge_ids ~= id;
    if (edge_type == EdgeType.UNDIRECTED) {
      edges ~= Edge(dst, src);
      input_edge_ids ~= IdSetLexicon.emptySetId();
    }
  }

  void simplifyEdgeChains(
      in InputVertexId[][] site_vertices,
      ref Edge[][] layer_edges,
      ref InputEdgeIdSetId[][] layer_input_edge_ids,
      IdSetLexicon input_edge_id_set_lexicon) const {
    if (_layers.empty()) return;

    // Merge the edges from all layers (in order to build a single graph).
    Edge[] merged_edges;
    InputEdgeIdSetId[] merged_input_edge_ids;
    int[] merged_edge_layers;
    mergeLayerEdges(layer_edges, layer_input_edge_ids,
        merged_edges, merged_input_edge_ids, merged_edge_layers);

    // The following fields will be reconstructed by EdgeChainSimplifier.
    foreach (edges; layer_edges) edges.length = 0;
    foreach (input_edge_ids; layer_input_edge_ids) input_edge_ids.length = 0;

    // The graph options are irrelevant for edge chain simplification, but we
    // try to set them appropriately anyway.
    auto graph_options = new GraphOptions(
        EdgeType.DIRECTED,
        GraphOptions.DegenerateEdges.KEEP,
        GraphOptions.DuplicateEdges.KEEP,
        GraphOptions.SiblingPairs.KEEP);
    auto graph = new Graph(graph_options, _sites, merged_edges, merged_input_edge_ids,
        input_edge_id_set_lexicon, null, null, IsFullPolygonPredicate());
    auto simplifier = new EdgeChainSimplifier(
        this, graph, merged_edge_layers, site_vertices,
        &layer_edges, &layer_input_edge_ids, input_edge_id_set_lexicon);
    simplifier.run();
  }

  /**
   * Merges the edges from all layers and sorts them in lexicographic order so
   * that we can construct a single graph.  The sort is stable, which means that
   * any duplicate edges within each layer will still be sorted by InputEdgeId.
   */
  void mergeLayerEdges(
      in Edge[][] layer_edges,
      in InputEdgeIdSetId[][] layer_input_edge_ids,
      ref Edge[] edges,
      ref InputEdgeIdSetId[] input_edge_ids,
      ref int[] edge_layers) const {
    LayerEdgeId[] order;
    for (int i = 0; i < layer_edges.length; ++i) {
      for (int e = 0; e < layer_edges[i].length; ++e) {
        order ~= LayerEdgeId(i, e);
      }
    }
    .sort!((a, b) => stableLessThan(layer_edges[a.first][a.second], layer_edges[b.first][b.second], a, b))(
            order);
    edges.reserve(order.length);
    input_edge_ids.reserve(order.length);
    edge_layers.reserve(order.length);
    for (int i = 0; i < order.length; ++i) {
      const LayerEdgeId id = order[i];
      edges ~= layer_edges[id.first][id.second];
      input_edge_ids ~= layer_input_edge_ids[id.first][id.second];
      edge_layers ~= id.first;
    }
  }

  /**
   * A comparision function that allows stable sorting with std::sort (which is
   * fast but not stable).  It breaks ties between equal edges by comparing
   * their LayerEdgeIds.
   */
  static bool stableLessThan(in Edge a, in Edge b, in LayerEdgeId ai, in LayerEdgeId bi) {
    // The compiler doesn't optimize this as well as it should:
    //   return make_pair(a, ai) < make_pair(b, bi);
    if (a[0] < b[0]) return true;
    if (b[0] < a[0]) return false;
    if (a[1] < b[1]) return true;
    if (b[1] < a[1]) return false;
    return ai < bi;  // Stable sort.
  }

  //////////// Parameters /////////////

  /// S2Builder options.
  Options _options;

  /// The maximum distance (inclusive) that a vertex can move when snapped,
  /// equal to S1ChordAngle(options_.snap_function().snap_radius()).
  S1ChordAngle _siteSnapRadiusCa;

  /**
   * The maximum distance (inclusive) that an edge can move when snapping to a
   * snap site.  It can be slightly larger than the site snap radius when
   * edges are being split at crossings.
   */
  S1ChordAngle _edgeSnapRadiusCa;

  S1Angle _maxEdgeDeviation;
  S1ChordAngle _edgeSiteQueryRadiusCa;
  S1ChordAngle _minEdgeLengthToSplitCa;

  S1Angle _minSiteSeparation;
  S1ChordAngle _minSiteSeparationCa;
  S1ChordAngle _minEdgeSiteSeparationCa;
  S1ChordAngle _minEdgeSiteSeparationCaLimit;

  S1ChordAngle _maxAdjacentSiteSeparationCa;

  /**
   * The squared sine of the edge snap radius.  This is equivalent to the snap
   * radius (squared) for distances measured through the interior of the
   * sphere to the plane containing an edge.  This value is used only when
   * interpolating new points along edges (see GetSeparationSite).
   */
  double _edgeSnapRadiusSin2;

  /// A copy of the argument to Build().
  S2Error _error;

  /**
   * True if snapping was requested.  This is true if either snap_radius() is
   * positive, or split_crossing_edges() is true (which implicitly requests
   * snapping to ensure that both crossing edges are snapped to the
   * intersection point).
   */
  bool _snappingRequested;

  /**
   * Initially false, and set to true when it is discovered that at least one
   * input vertex or edge does not meet the output guarantees (e.g., that
   * vertices are separated by at least snap_function.min_vertex_separation).
   */
  bool _snappingNeeded;

  //////////// Input Data /////////////

  /// A flag indicating whether label_set_ has been modified since the last
  /// time label_set_id_ was computed.
  bool _labelSetModified;

  S2Point[] _inputVertices;
  InputEdge[] _inputEdges;

  Layer[] _layers;
  GraphOptions[] _layerOptions;
  InputEdgeId[] _layerBegins;
  IsFullPolygonPredicate[] _layerIsFullPolygonPredicates;

  /**
   * Each input edge has "label set id" (an int32) representing the set of
   * labels attached to that edge.  This vector is populated only if at least
   * one label is used.
   */
  LabelSetId[] _labelSetIds;
  IdSetLexicon _labelSetLexicon;

  /// The current set of labels (represented as a stack).
  Label[] _labelSet;

  /// The LabelSetId corresponding to the current label set, computed on demand
  /// (by adding it to label_set_lexicon()).
  LabelSetId _labelSetId;

  ////////////// Data for Snapping and Simplifying //////////////

  /// The number of sites specified using ForceVertex().  These sites are
  /// always at the beginning of the sites_ vector.
  SiteId _numForcedSites;

  // The set of snapped vertex locations ("sites").
  S2Point[] _sites;

  /**
   * A map from each input edge to the set of sites "nearby" that edge,
   * defined as the set of sites that are candidates for snapping and/or
   * avoidance.  Note that compact_array will inline up to two sites, which
   * usually takes care of the vast majority of edges.  Sites are kept sorted
   * by increasing distance from the origin of the input edge.
   *
   * Once snapping is finished, this field is discarded unless edge chain
   * simplification was requested, in which case instead the sites are
   * filtered by removing the ones that each edge was snapped to, leaving only
   * the "sites to avoid" (needed for simplification).
   */
  SiteId[][] _edgeSites;
}

/**
 * This class is only needed by S2Builder::Layer implementations.  A layer is
 * responsible for assembling an S2Builder::Graph of snapped edges into the
 * desired output format (e.g., an S2Polygon).  The GraphOptions class allows
 * each Layer type to specify requirements on its input graph: for example, if
 * DegenerateEdges::DISCARD is specified, then S2Builder will ensure that all
 * degenerate edges are removed before passing the graph to the S2Layer::Build
 * method.
 */
class GraphOptions {
public:
  alias EdgeType = S2Builder.EdgeType;

  /**
   * All S2Builder::Layer subtypes should specify GraphOptions explicitly
   * using this constructor, rather than relying on default values.
   */
  this(EdgeType edge_type, DegenerateEdges degenerate_edges,
      DuplicateEdges duplicate_edges, SiblingPairs sibling_pairs) {
    _edgeType = edge_type;
    _degenerateEdges = degenerate_edges;
    _duplicateEdges = duplicate_edges;
    _siblingPairs = sibling_pairs;
    _allowVertexFiltering = true;
  }

  /**
   * The default options specify that all edges should be kept, since this
   * produces the least surprising output and makes it easier to diagnose the
   * problem when an option is left unspecified.
   */
  this() {
    _edgeType = EdgeType.DIRECTED;
    _degenerateEdges = DegenerateEdges.KEEP;
    _duplicateEdges = DuplicateEdges.KEEP;
    _siblingPairs = SiblingPairs.KEEP;
    _allowVertexFiltering = true;
  }

  /**
   * Specifies whether the S2Builder input edges should be treated as
   * undirected.  If true, then all input edges are duplicated into pairs
   * consisting of an edge and a sibling (reverse) edge.  The layer
   * implementation is responsible for ensuring that exactly one edge from
   * each pair is used in the output, i.e. *only half* of the graph edges will
   * be used.  (Note that some values of the sibling_pairs() option
   * automatically take care of this issue by removing half of the edges and
   * changing edge_type() to DIRECTED.)
   */
  EdgeType edgeType() const {
    return _edgeType;
  }

  void setEdgeType(EdgeType edge_type) {
    _edgeType = edge_type;
  }

  /**
   * Controls how degenerate edges (i.e., an edge from a vertex to itself) are
   * handled.  Such edges may be present in the input, or they may be created
   * when both endpoints of an edge are snapped to the same output vertex.
   * The options available are:
   *
   * DISCARD: Discards all degenerate edges.  This is useful for layers that
   *          do not support degeneracies, such as S2PolygonLayer.
   *
   * DISCARD_EXCESS: Discards all degenerate edges that are connected to
   *                 non-degenerate edges.  (Any remaining duplicate edges can
   *                 be merged using DuplicateEdges::MERGE.)  This is useful
   *                 for simplifying polygons while ensuring that loops that
   *                 collapse to a single point do not disappear.
   *
   * KEEP: Keeps all degenerate edges.  Be aware that this may create many
   *       redundant edges when simplifying geometry (e.g., a polyline of the
   *       form AABBBBBCCCCCCDDDD).  DegenerateEdges::KEEP is mainly useful
   *       for algorithms that require an output edge for every input edge.
   */
  enum DegenerateEdges { DISCARD, DISCARD_EXCESS, KEEP }

  DegenerateEdges degenerateEdges() const {
    return _degenerateEdges;
  }

  void setDegenerateEdges(DegenerateEdges degenerate_edges) {
    _degenerateEdges = degenerate_edges;
  }

  /**
   * Controls how duplicate edges (i.e., edges that are present multiple
   * times) are handled.  Such edges may be present in the input, or they can
   * be created when vertices are snapped together.  When several edges are
   * merged, the result is a single edge labelled with all of the original
   * input edge ids.
   */
  enum DuplicateEdges { MERGE, KEEP }

  DuplicateEdges duplicateEdges() const {
    return _duplicateEdges;
  }

  void setDuplicateEdges(DuplicateEdges duplicate_edges) {
    _duplicateEdges = duplicate_edges;
  }

  /**
   * Controls how sibling edge pairs (i.e., pairs consisting of an edge and
   * its reverse edge) are handled.  Layer types that define an interior
   * (e.g., polygons) normally discard such edge pairs since they do not
   * affect the result (i.e., they define a "loop" with no interior).  The
   * various options include:
   *
   * DISCARD: Discards all sibling edge pairs.
   *
   * DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept
   *                 if the result would otherwise be empty.  This is useful
   *                 for polygons with degeneracies (S2LaxPolygonShape), and
   *                 for simplifying polylines while ensuring that they are
   *                 not split into multiple disconnected pieces.
   *
   * KEEP: Keeps sibling pairs.  This can be used to create polylines that
   *       double back on themselves, or degenerate loops (with a layer type
   *       such as S2LaxPolygonShape).
   *
   * REQUIRE: Requires that all edges have a sibling (and returns an error
   *          otherwise).  This is useful with layer types that create a
   *          collection of adjacent polygons (a polygon mesh).
   *
   * CREATE: Ensures that all edges have a sibling edge by creating them if
   *         necessary.  This is useful with polygon meshes where the input
   *         polygons do not cover the entire sphere.  Such edges always
   *         have an empty set of labels.
   *
   * If edge_type() is EdgeType::UNDIRECTED, a sibling edge pair is considered
   * to consist of four edges (two duplicate edges and their siblings), since
   * only two of these four edges will be used in the final output.
   *
   * Furthermore, since the options REQUIRE and CREATE guarantee that all
   * edges will have siblings, S2Builder implements these options for
   * undirected edges by discarding half of the edges in each direction and
   * changing the edge_type() to EdgeType::DIRECTED.  For example, two
   * undirected input edges between vertices A and B would first be converted
   * into two directed edges in each direction, and then one edge of each pair
   * would be discarded leaving only one edge in each direction.
   *
   * Degenerate edges are considered not to have siblings.  If such edges are
   * present, they are passed through unchanged by SiblingPairs::DISCARD.  For
   * SiblingPairs::REQUIRE or SiblingPairs::CREATE with undirected edges, the
   * number of copies of each degenerate edge is reduced by a factor of two.
   *
   * Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and
   * REQUIRE/CREATE in the case of undirected edges) have the side effect that
   * when duplicate edges are present, all of the corresponding edge labels
   * are merged together and assigned to the remaining edges.  (This avoids
   * the problem of having to decide which edges are discarded.)  Note that
   * this merging takes place even when all copies of an edge are kept, and
   * that even labels attached to duplicate degenerate edges are merged.  For
   * example, consider the graph {AB1, AB2, BA3, CD4, CD5} (where XYn denotes
   * an edge from X to Y with label "n").  With SiblingPairs::DISCARD, we need
   * to discard one of the copies of AB.  But which one?  Rather than choosing
   * arbitrarily, instead we merge the labels of all duplicate edges (even
   * ones where no sibling pairs were discarded), yielding {AB12, CD45, CD45}
   * (assuming that duplicate edges are being kept).
   */
  enum SiblingPairs { DISCARD, DISCARD_EXCESS, KEEP, REQUIRE, CREATE }

  SiblingPairs siblingPairs() const {
    return _siblingPairs;
  }

  void setSiblingPairs(SiblingPairs sibling_pairs) {
    _siblingPairs = sibling_pairs;
  }

  /**
   * This is a specialized option that is only needed by clients want to work
   * with the graphs for multiple layers at the same time (e.g., in order to
   * check whether the same edge is present in two different graphs).  [Note
   * that if you need to do this, usually it is easier just to build a single
   * graph with suitable edge labels.]
   *
   * When there are a large number of layers, then by default S2Builder builds
   * a minimal subgraph for each layer containing only the vertices needed by
   * the edges in that layer.  This ensures that layer types that iterate over
   * the vertices run in time proportional to the size of that layer rather
   * than the size of all layers combined.  (For example, if there are a
   * million layers with one edge each, then each layer would be passed a
   * graph with 2 vertices rather than 2 million vertices.)
   *
   * If this option is set to false, this optimization is disabled.  Instead
   * the graph passed to this layer will contain the full set of vertices.
   * (This is not recommended when the number of layers could be large.)
   *
   * DEFAULT: true
   */
  bool allowVertexFiltering() const {
    return _allowVertexFiltering;
  }

  void setAllowVertexFiltering(bool allow_vertex_filtering) {
    _allowVertexFiltering = allow_vertex_filtering;
  }

  override
  bool opEquals(Object other) const {
    GraphOptions y = cast(GraphOptions) other;
    if (y is null) return false;
    return (edgeType() == y.edgeType()
        && degenerateEdges() == y.degenerateEdges()
        && duplicateEdges() == y.duplicateEdges()
        && siblingPairs() == y.siblingPairs()
        && allowVertexFiltering() == y.allowVertexFiltering());
  }

private:
  EdgeType _edgeType;
  DegenerateEdges _degenerateEdges;
  DuplicateEdges _duplicateEdges;
  SiblingPairs _siblingPairs;
  bool _allowVertexFiltering;
}

/**
 * An S2Shape used to represent the entire collection of S2Builder input edges.
 * Vertices are specified as indices into a vertex vector to save space.
 */
class VertexIdEdgeVectorShape : S2Shape {
public:
  /// Requires that "edges" is constant for the lifetime of this object.
  this(in int[2][] edges, in S2Point[] vertices) {
    _edges = edges;
    _vertices = vertices;
  }

  const(S2Point) vertex0(int e) const {
    return vertex(_edges[e][0]);
  }

  const(S2Point) vertex1(int e) const {
    return vertex(_edges[e][1]);
  }

  // S2Shape interface:
  final override
  int numEdges() const {
    return cast(int) _edges.length;
  }

  final override
  Edge edge(int e) const {
    return Edge(_vertices[_edges[e][0]], _vertices[_edges[e][1]]);
  }

  final override
  int dimension() const {
    return 1;
  }

  final override
  ReferencePoint getReferencePoint() const {
    return ReferencePoint(false);
  }

  final override
  int numChains() const {
    return cast(int) _edges.length;
  }

  final override
  Chain chain(int i) const {
    return Chain(i, 1);
  }

  final override
  Edge chainEdge(int i, int j) const {
    return edge(i);
  }

  final override
  ChainPosition chainPosition(int e) const {
    return ChainPosition(e, 0);
  }

private:
  const(S2Point) vertex(int i) const {
    return _vertices[i];
  }

  const(int[2][]) _edges;
  const(S2Point[]) _vertices;
}

// A class that encapsulates the state needed for simplifying edge chains.
class EdgeChainSimplifier {
 public:
  // The graph "g" contains all edges from all layers.  "edge_layers"
  // indicates the original layer for each edge.  "site_vertices" is a map
  // from SiteId to the set of InputVertexIds that were snapped to that site.
  // "layer_edges" and "layer_input_edge_ids" are output arguments where the
  // simplified edge chains will be placed.  The input and output edges are
  // not sorted.
  this(
      in S2Builder builder,
      in Graph g,
      in int[] edge_layers,
      in InputVertexId[][] site_vertices,
      Edge[][]* layer_edges,
      InputEdgeIdSetId[][]* layer_input_edge_ids,
      IdSetLexicon input_edge_id_set_lexicon) {
    _builder = builder;
    _g = g;
    _in = new VertexInMap(g);
    _out = new VertexOutMap(g);
    _edgeLayers = edge_layers;
    _siteVertices = site_vertices;
    _layerEdges = layer_edges;
    _layerInputEdgeIds = layer_input_edge_ids;
    _inputEdgeIdSetLexicon = input_edge_id_set_lexicon;
    _layerBegins = _builder._layerBegins;
    _isInterior = new bool[](g.numVertices());
    _used = new bool[](g.numEdges());
    _newEdges.reserve(g.numEdges());
    _newInputEdgeIds.reserve(g.numEdges());
    _newEdgeLayers.reserve(g.numEdges());
  }

  void run() {
    // Determine which vertices can be interior vertices of an edge chain.
    for (VertexId v = 0; v < _g.numVertices(); ++v) {
      _isInterior[v] = isInterior(v);
    }
    // Attempt to simplify all edge chains that start from a non-interior
    // vertex.  (This takes care of all chains except loops.)
    for (EdgeId e = 0; e < _g.numEdges(); ++e) {
      if (_used[e]) continue;
      Edge edge = _g.edge(e);
      if (_isInterior[edge[0]]) continue;
      if (!_isInterior[edge[1]]) {
        outputEdge(e);  // An edge between two non-interior vertices.
      } else {
        simplifyChain(edge[0], edge[1]);
      }
    }
    // If there are any edges left, they form one or more disjoint loops where
    // all vertices are interior vertices.
    //
    // TODO(ericv): It would be better to start from the edge with the smallest
    // min_input_edge_id(), since that would make the output more predictable
    // for testing purposes.  It also means that we won't create an edge that
    // spans the start and end of a polyline if the polyline is snapped into a
    // loop.  (Unfortunately there are pathological examples that prevent us
    // from guaranteeing this in general, e.g. there could be two polylines in
    // different layers that snap to the same loop but start at different
    // positions.  In general we only consider input edge ids to be a hint
    // towards the preferred output ordering.)
    for (EdgeId e = 0; e < _g.numEdges(); ++e) {
      if (_used[e]) continue;
      Edge edge = _g.edge(e);
      if (edge[0] == edge[1]) {
        // Note that it is safe to output degenerate edges as we go along,
        // because this vertex has at least one non-degenerate outgoing edge and
        // therefore we will (or just did) start an edge chain here.
        outputEdge(e);
      } else {
        simplifyChain(edge[0], edge[1]);
      }
    }

    // Finally, copy the output edges into the appropriate layers.  They don't
    // need to be sorted because the input edges were also unsorted.
    for (int e = 0; e < _newEdges.length; ++e) {
      int layer = _newEdgeLayers[e];
      (*_layerEdges)[layer] ~= _newEdges[e];
      (*_layerInputEdgeIds)[layer] ~= _newInputEdgeIds[e];
    }
  }

private:
  alias VertexId = Graph.VertexId;
  alias VertexInMap = Graph.VertexInMap;
  alias VertexOutMap = Graph.VertexOutMap;

  alias InputVertexId = S2Builder.InputVertexId;
  alias Edge = S2Builder.Edge;
  alias EdgeId = S2Builder.EdgeId;
  alias InputEdgeId = S2Builder.InputEdgeId;
  alias InputEdgeIdSetId = S2Builder.InputEdgeIdSetId;

  /**
   * A helper class for determining whether a vertex can be an interior vertex
   * of a simplified edge chain.  Such a vertex must be adjacent to exactly two
   * vertices (across all layers combined), and in each layer the number of
   * incoming edges from one vertex must equal the number of outgoing edges to
   * the other vertex (in both directions).  Furthermore the vertex cannot have
   * any degenerate edges in a given layer unless it has at least one
   * non-degenerate edge in that layer as well.  (Note that usually there will
   * not be any degenerate edges at all, since most layer types discard them.)
   *
   * The last condition is necessary to prevent the following: suppose that one
   * layer has a chain ABC and another layer has a degenerate edge BB (with no
   * other edges incident to B).  Then we can't simplify ABC to AC because there
   * would be no suitable replacement for the edge BB (since the input edge that
   * mapped to BB can't be replaced by any of the edges AA, AC, or CC without
   * moving further than snap_radius).
   */
  static class InteriorVertexMatcher {
  public:
    /// Checks whether "v0" can be an interior vertex of an edge chain.
    this(VertexId v0) {
      _v0 = v0;
      _v1 = -1;
      _v2 = -1;
      _n0 = 0;
      _n1 = 0;
      _n2 = 0;
      _excessOut = 0;
      _tooManyEndpoints = false;
    }

    /// Starts analyzing the edges of a new layer.
    void startLayer() {
      _excessOut = _n0 = _n1 = _n2 = 0;
    }

    /// This method should be called for each edge incident to "v0" in a given
    /// layer.  (For degenerate edges, it should be called twice.)
    void tally(VertexId v, bool outgoing) {
      _excessOut += outgoing ? 1 : -1;  // outdegree - indegree
      if (v == _v0) {
        ++_n0;  // Counts both endpoints of each degenerate edge.
      } else {
        // We keep track of the total number of edges (incoming or outgoing)
        // connecting v0 to up to two adjacent vertices.
        if (_v1 < 0) _v1 = v;
        if (_v1 == v) {
          ++_n1;
        } else {
          if (_v2 < 0) _v2 = v;
          if (_v2 == v) {
            ++_n2;
          } else {
            _tooManyEndpoints = true;
          }
        }
      }
    }

    /// This method should be called after processing the edges for each layer.
    /// It returns true if "v0" is an interior vertex based on the edges so far.
    bool matches() const {
      // We check that there are the same number of incoming and outgoing edges
      // in each direction by verifying that (1) indegree(v0) == outdegree(v0)
      // and (2) the total number of edges (incoming and outgoing) to "v1" and
      // "v2" are equal.  We also check the condition on degenerate edges that
      // is documented above.
      return (!_tooManyEndpoints && _excessOut == 0 && _n1 == _n2 && (_n0 == 0 || _n1 > 0));
    }

  private:
    VertexId _v0, _v1, _v2;
    int _n0, _n1, _n2;
    int _excessOut;           // outdegree(v0) - indegree(v0)
    bool _tooManyEndpoints;  // Have we seen more than two adjacent vertices?
  }

  /// Copies the given edge to the output and marks it as used.
  void outputEdge(EdgeId e) {
    _newEdges ~= _g.edge(e);
    _newInputEdgeIds ~= _g.inputEdgeIdSetId(e);
    _newEdgeLayers ~= _edgeLayers[e];
    _used[e] = true;
  }

  /// Returns the layer that a given graph edge belongs to.
  int graphEdgeLayer(EdgeId e) const {
    return _edgeLayers[e];
  }

  /// Returns the layer than a given input edge belongs to.
  int inputEdgeLayer(InputEdgeId id) const
  in {
    // NOTE(ericv): If this method shows up in profiling, the result could be
    // stored with each edge (i.e., edge_layers_ and new_edge_layers_).
    assert(id >= 0);
  } do {
    auto triResults = _layerBegins.assumeSorted().trisect(id);
    return cast(int) (triResults[0].length + triResults[1].length) - 1;
  }

  /// Returns true if VertexId "v" can be an interior vertex of a simplified edge
  /// chain.  (See the InteriorVertexMatcher class for what this implies.)
  bool isInterior(VertexId v) {
    // Check a few simple prerequisites.
    if (_out.degree(v) == 0) return false;
    if (_out.degree(v) != _in.degree(v)) return false;
    if (v < _builder._numForcedSites) return false;  // Keep forced vertices.

    // Sort the edges so that they are grouped by layer.
    EdgeId[] edges = _tmpEdges;  // Avoid allocating each time.
    edges.length = 0;
    foreach (EdgeId e; _out.edgeIds(v)) edges ~= e;
    foreach (EdgeId e; _in.edgeIds(v)) edges ~= e;
    sort!((EdgeId x, EdgeId y) {
          return graphEdgeLayer(x) < graphEdgeLayer(y);
        })(edges);
    // Now feed the edges in each layer to the InteriorVertexMatcher.
    auto matcher = new InteriorVertexMatcher(v);
    for (auto i = 0; i < edges.length; ) {
      int layer = graphEdgeLayer(edges[i]);
      matcher.startLayer();
      for (; i < edges.length && graphEdgeLayer(edges[i]) == layer; ++i) {
        Edge edge = _g.edge(edges[i]);
        if (edge[0] == v) matcher.tally(edge[0], true /*outgoing*/);
        if (edge[1] == v) matcher.tally(edge[1], false /*outgoing*/);
      }
      if (!matcher.matches()) return false;
    }
    return true;
  }

  /**
   * Follows the edge chain starting with (v0, v1) until either we find a
   * non-interior vertex or we return to the original vertex v0.  At each vertex
   * we simplify a subchain of edges that is as long as possible.
   */
  void simplifyChain(VertexId v0, VertexId v1) {
    // Avoid allocating "chain" each time by reusing it.
    VertexId[] chain = _tmpVertices;
    auto simplifier = new S2PolylineSimplifier();
    VertexId vstart = v0;
    bool done = false;
    do {
      // Simplify a subchain of edges starting (v0, v1).
      simplifier.initialize(_g.vertex(v0));
      avoidSites(v0, v0, v1, simplifier);
      chain ~= v0;
      do {
        chain ~= v1;
        done = !_isInterior[v1] || v1 == vstart;
        if (done) break;

        // Attempt to extend the chain to the next vertex.
        VertexId vprev = v0;
        v0 = v1;
        v1 = followChain(vprev, v0);
      } while (targetInputVertices(v0, simplifier)
          && avoidSites(chain[0], v0, v1, simplifier)
          && simplifier.extend(_g.vertex(v1)));

      if (chain.length == 2) {
        outputAllEdges(chain[0], chain[1]);  // Could not simplify.
      } else {
        mergeChain(chain);
      }
      // Note that any degenerate edges that were not merged into a chain are
      // output by EdgeChainSimplifier::Run().
      chain.length = 0;
    } while (!done);
  }

  /// Given an edge (v0, v1) where v1 is an interior vertex, returns the (unique)
  /// next vertex in the edge chain.
  VertexId followChain(VertexId v0, VertexId v1) const
  in {
    assert(_isInterior[v1]);
  } do {
    foreach (EdgeId e; _out.edgeIds(v1)) {
      VertexId v = _g.edge(e)[1];
      if (v != v0 && v != v1) return v;
    }
    logger.logFatal("Could not find next edge in edge chain");
    return -1;
  }

  /// Copies all input edges between v0 and v1 (in both directions) to the output.
  void outputAllEdges(VertexId v0, VertexId v1) {
    foreach (EdgeId e; _out.edgeIds(v0, v1)) outputEdge(e);
    foreach (EdgeId e; _out.edgeIds(v1, v0)) outputEdge(e);
  }

  /// Ensures that the simplified edge passes within "edge_snap_radius" of all
  /// the *input* vertices that snapped to the given vertex "v".
  bool targetInputVertices(VertexId v, S2PolylineSimplifier simplifier) const {
    foreach (InputVertexId i; _siteVertices[v]) {
      if (!simplifier.targetDisc(_builder._inputVertices[i], _builder._edgeSnapRadiusCa)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Given the starting vertex v0 and last edge (v1, v2) of an edge chain,
   * restricts the allowable range of angles in order to ensure that all sites
   * near the edge (v1, v2) are avoided by at least min_edge_vertex_separation.
   */
  bool avoidSites(VertexId v0, VertexId v1, VertexId v2, S2PolylineSimplifier simplifier) const {
    import s2.s2predicates : orderedCCW, sign;
    const S2Point p0 = _g.vertex(v0);
    const S2Point p1 = _g.vertex(v1);
    const S2Point p2 = _g.vertex(v2);
    auto r1 = S1ChordAngle(p0, p1);
    auto r2 = S1ChordAngle(p0, p2);

    // The distance from the start of the edge chain must increase monotonically
    // for each vertex, since we don't want to simplify chains that backtrack on
    // themselves (we want a parametric approximation, not a geometric one).
    if (r2 < r1) return false;

    // We also limit the maximum edge length in order to guarantee that the
    // simplified edge stays with max_edge_deviation() of all the input edges
    // that snap to it.
    if (r2 >= _builder._minEdgeLengthToSplitCa) return false;

    // Otherwise it is sufficient to consider the nearby sites (edge_sites_) for
    // a single input edge that snapped to (v1, v2) or (v2, v1).  This is
    // because each edge has a list of all sites within (max_edge_deviation +
    // min_edge_vertex_separation), and since the output edge is within
    // max_edge_deviation of all input edges, this list includes all sites
    // within min_edge_vertex_separation of the output edge.
    //
    // Usually there is only one edge to choose from, but it's not much more
    // effort to choose the edge with the shortest list of edge_sites_.
    InputEdgeId best = -1;
    const edge_sites = _builder._edgeSites;
    foreach (EdgeId e; _out.edgeIds(v1, v2)) {
      foreach (InputEdgeId id; _g.inputEdgeIds(e)) {
        if (best < 0 || edge_sites[id].length < edge_sites[best].length)
          best = id;
      }
    }
    foreach (EdgeId e; _out.edgeIds(v2, v1)) {
      foreach (InputEdgeId id; _g.inputEdgeIds(e)) {
        if (best < 0 || edge_sites[id].length < edge_sites[best].length)
          best = id;
      }
    }
    debug enforce(best >= 0);  // Because there is at least one outgoing edge.

    foreach (VertexId v; edge_sites[best]) {
      // This test is optional since these sites are excluded below anyway.
      if (v == v0 || v == v1 || v == v2) continue;

      // We are only interested in sites whose distance from "p0" is in the
      // range (r1, r2).  Sites closer than "r1" have already been processed,
      // and sites further than "r2" aren't relevant yet.
      S2Point p = _g.vertex(v);
      auto r = S1ChordAngle(p0, p);
      if (r <= r1 || r >= r2) continue;

      // We need to figure out whether this site is to the left or right of the
      // edge chain.  For the first edge this is easy.  Otherwise, since we are
      // only considering sites in the radius range (r1, r2), we can do this by
      // checking whether the site is to the left of the wedge (p0, p1, p2).
      bool disc_on_left = (v1 == v0) ? (sign(p1, p2, p) > 0) : orderedCCW(p0, p2, p, p1);
      if (!simplifier.avoidDisc(p, _builder._minEdgeSiteSeparationCa, disc_on_left)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Given the vertices in a simplified edge chain, adds the corresponding
   * simplified edge(s) to the output.  Note that (1) the edge chain may exist
   * in multiple layers, (2) the edge chain may exist in both directions, (3)
   * there may be more than one copy of an edge chain (in either direction)
   * within a single layer.
   */
  void mergeChain(in VertexId[] vertices) {
    // Suppose that all interior vertices have M outgoing edges and N incoming
    // edges.  Our goal is to group the edges into M outgoing chains and N
    // incoming chains, and then replace each chain by a single edge.
    S2Builder.InputEdgeId[][] merged_input_ids;
    InputEdgeId[] degenerate_ids;
    int num_out;  // Edge count in the outgoing direction.
    for (int i = 1; i < vertices.length; ++i) {
      VertexId v0 = vertices[i - 1];
      VertexId v1 = vertices[i];
      auto out_edges = _out.edgeIds(v0, v1);
      auto in_edges = _out.edgeIds(v1, v0);
      if (i == 1) {
        // Allocate space to store the input edge ids associated with each edge.
        num_out = cast(int) out_edges.length;
        merged_input_ids.length = num_out + in_edges.length;
        foreach (InputEdgeId[] ids; merged_input_ids) {
          ids.reserve(vertices.length - 1);
        }
      } else {
        // For each interior vertex, we build a list of input edge ids
        // associated with degenerate edges.  Each input edge ids will be
        // assigned to one of the output edges later.  (Normally there are no
        // degenerate edges at all since most layer types don't want them.)
        debug enforce(_isInterior[v0]);
        foreach (EdgeId e; _out.edgeIds(v0, v0)) {
          foreach (InputEdgeId id; _g.inputEdgeIds(e)) {
            degenerate_ids ~= id;
          }
          _used[e] = true;
        }
      }
      // Because the edges were created in layer order, and all sorts used are
      // stable, the edges are still in layer order.  Therefore we can simply
      // merge together all the edges in the same relative position.
      int j = 0;
      foreach (EdgeId e; out_edges) {
        foreach (InputEdgeId id; _g.inputEdgeIds(e)) {
          merged_input_ids[j] ~= id;
        }
        _used[e] = true;
        ++j;
      }
      foreach (EdgeId e; in_edges) {
        foreach (InputEdgeId id; _g.inputEdgeIds(e)) {
          merged_input_ids[j] ~= id;
        }
        _used[e] = true;
        ++j;
      }
      debug enforce(merged_input_ids.length == j);
    }
    if (!degenerate_ids.empty()) {
      sort(degenerate_ids);
      assignDegenerateEdges(degenerate_ids, merged_input_ids);
    }
    // Output the merged edges.
    VertexId v0 = vertices[0], v1 = vertices[1], vb = vertices.back();
    foreach (EdgeId e; _out.edgeIds(v0, v1)) {
      _newEdges ~= Edge(v0, vb);
      _newEdgeLayers ~= graphEdgeLayer(e);
    }
    foreach (EdgeId e; _out.edgeIds(v1, v0)) {
      _newEdges ~= Edge(vb, v0);
      _newEdgeLayers ~= graphEdgeLayer(e);
    }
    foreach (ids; merged_input_ids) {
      _newInputEdgeIds ~= _inputEdgeIdSetLexicon.add(ids);
    }
  }

  /**
   * Given a list of the input edge ids associated with degenerate edges in the
   * interior of an edge chain, assigns each input edge id to one of the output
   * edges.
   */
  void assignDegenerateEdges(
      in S2Builder.InputEdgeId[] degenerate_ids, ref S2Builder.InputEdgeId[][] merged_ids) const {
    // Each degenerate edge is assigned to an output edge in the appropriate
    // layer.  If there is more than one candidate, we use heuristics so that if
    // the input consists of a chain of edges provided in consecutive order
    // (some of which became degenerate), then all those input edges are
    // assigned to the same output edge.  For example, suppose that one output
    // edge is labeled with input edges 3,4,7,8, while another output edge is
    // labeled with input edges 50,51,54,57.  Then if we encounter degenerate
    // edges labeled with input edges 5 and 6, we would prefer to assign them to
    // the first edge (yielding the continuous range 3,4,5,6,7,8).
    //
    // The heuristic below is only smart enough to handle the case where the
    // candidate output edges have non-overlapping ranges of input edges.
    // (Otherwise there is probably not a good heuristic for assigning the
    // degenerate edges in any case.)

    // Duplicate edge ids don't affect the heuristic below, so we don't bother
    // removing them.  (They will be removed by IdSetLexicon::Add.)
    foreach (ids; merged_ids) sort(ids);

    // Sort the output edges by their minimum input edge id.  This is sufficient
    // for the purpose of determining which layer they belong to.  With
    // EdgeType::UNDIRECTED, some edges might not have any input edge ids (i.e.,
    // if they consist entirely of siblings of input edges).  We simply remove
    // such edges from the lists of candidates.
    uint[] order;
    order.reserve(merged_ids.length);
    for (int i = 0; i < merged_ids.length; ++i) {
      if (!merged_ids[i].empty()) order ~= i;
    }
    sort!((int i, int j) {
          return merged_ids[i][0] < merged_ids[j][0];
        })(order);

    // Now determine where each degenerate edge should be assigned.
    for (int i = 0; i < degenerate_ids.length; ++i) {
      InputEdgeId degenerate_id = degenerate_ids[i];
      int layer = inputEdgeLayer(degenerate_id);

      // Find the first output edge whose range of input edge ids starts after
      // "degenerate_id".  If the previous edge belongs to the correct layer,
      // then we assign the degenerate edge to it.
      auto it = countUntil!((uint y) {
            return degenerate_id < merged_ids[y][0];
          })(order);
      if (it != 0) {
        if (merged_ids[order[it - 1]][0] >= _layerBegins[layer]) --it;
      }
      debug enforce(inputEdgeLayer(merged_ids[order[it]][0]) == layer);
      merged_ids[order[it]] ~= degenerate_id;
    }
  }

  const(S2Builder) _builder;
  const(Graph) _g;
  const(Graph.VertexInMap) _in;
  const(Graph.VertexOutMap) _out;
  const(int[]) _edgeLayers;
  const(S2Builder.InputVertexId[][]) _siteVertices;
  S2Builder.Edge[][]* _layerEdges;
  S2Builder.InputEdgeIdSetId[][]* _layerInputEdgeIds;
  IdSetLexicon _inputEdgeIdSetLexicon;

  /// Convenience member copied from builder_.
  const(S2Builder.InputEdgeId[]) _layerBegins;

  /**
   * is_interior_[v] indicates that VertexId "v" is eligible to be an interior
   * vertex of a simplified edge chain.  You can think of it as vertex whose
   * indegree and outdegree are both 1 (although the actual definition is a
   * bit more complicated because of duplicate edges and layers).
   */
  bool[] _isInterior;

  /// used_[e] indicates that EdgeId "e" has already been processed.
  bool[] _used;

  // Temporary vectors, declared here to avoid repeated allocation.
  VertexId[] _tmpVertices;
  S2Builder.EdgeId[] _tmpEdges;

  // The output edges after simplification.
  S2Builder.Edge[] _newEdges;
  S2Builder.InputEdgeIdSetId[] _newInputEdgeIds;
  int[] _newEdgeLayers;
}

// Helper functions for computing error bounds:

private S1ChordAngle roundUp(S1Angle a) {
  auto ca = S1ChordAngle(a);
  return ca.plusError(ca.getS1AngleConstructorMaxError());
}

private S1ChordAngle addPointToPointError(S1ChordAngle ca) {
  return ca.plusError(ca.getS2PointConstructorMaxError());
}

private S1ChordAngle addPointToEdgeError(S1ChordAngle ca) {
  return ca.plusError(getUpdateMinDistanceMaxError(ca));
}

private bool isFullPolygonUnspecified(in Graph g, out S2Error error) {
  error.initialize(S2Error.Code.BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED,
      "A degenerate polygon was found, but no predicate was specified "
      ~ "to determine whether the polygon is empty or full.  Call "
      ~ "S2Builder::AddIsFullPolygonPredicate() to fix this problem.");
  return false;  // Assumes the polygon is empty.
}

private void dumpEdges(in Graph.Edge[] edges, in S2Point[] vertices) {
  foreach (e; edges) {
    S2Point[] v;
    v ~= vertices[e[0]];
    v ~= vertices[e[1]];
    writeln("S2Polyline: ", toString(v), "(", e[0], ",", e[1], ")");
  }
}
