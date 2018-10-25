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

module s2.s2contains_point_query;

import s2.s2edge_crosser;
import s2.s2edge_crossings : vertexCrossing;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.shapeutil.shape_edge;

import std.typecons : Rebindable;

/**
 * Defines whether shapes are considered to contain their vertices.  Note that
 * these definitions differ from the ones used by S2BooleanOperation.
 *
 *  - In the OPEN model, no shapes contain their vertices (not even points).
 *    Therefore Contains(S2Point) returns true if and only if the point is
 *    in the interior of some polygon.
 *
 *  - In the SEMI_OPEN model, polygon point containment is defined such that
 *    if several polygons tile the region around a vertex, then exactly one of
 *    those polygons contains that vertex.  Points and polylines still do not
 *    contain any vertices.
 *
 *  - In the CLOSED model, all shapes contain their vertices (including points
 *    and polylines).
 *
 * Note that points other than vertices are never contained by polylines.
 * If you want need this behavior, use S2ClosestEdgeQuery::IsDistanceLess()
 * with a suitable distance threshold instead.
 */
enum S2VertexModel { OPEN, SEMI_OPEN, CLOSED }

/// This class defines the options supported by S2ContainsPointQuery.
class S2ContainsPointQueryOptions {
public:
  this() {}

  /// Convenience constructor that sets the vertex_model() option.
  this(S2VertexModel vertex_model) {
    _vertexModel = vertex_model;
  }

  /**
   * Controls whether shapes are considered to contain their vertices (see
   * definitions above).  By default the SEMI_OPEN model is used.
   *
   * DEFAULT: S2VertexModel::SEMI_OPEN
   */
  S2VertexModel vertexModel() const {
    return _vertexModel;
  }

  void setVertexModel(S2VertexModel model) {
    _vertexModel = model;
  }

private:
  S2VertexModel _vertexModel = S2VertexModel.SEMI_OPEN;
}

// S2ContainsPointQuery determines whether one or more shapes in an
// S2ShapeIndex contain a given S2Point.  The S2ShapeIndex may contain any
// number of points, polylines, and/or polygons (possibly overlapping).
// Shape boundaries may be modeled as OPEN, SEMI_OPEN, or CLOSED (this affects
// whether or not shapes are considered to contain their vertices).
//
// Example usage:
//   auto query = MakeS2ContainsPointQuery(&index, S2VertexModel::CLOSED);
//   return query.Contains(point);
//
// This class is not thread-safe.  To use it in parallel, each thread should
// construct its own instance (this is not expensive).
//
// However, note that if you need to do a large number of point containment
// tests, it is more efficient to re-use the S2ContainsPointQuery object
// rather than constructing a new one each time.
class S2ContainsPointQuery(IndexT) {
private:
  alias Iterator = IndexT.Iterator;

public:
  /// Default constructor; requires Init() to be called.
  this() {
    _index = null;
  }

  /**
   * Rather than calling this constructor, which requires specifying the
   * IndexT template argument explicitly, the preferred idiom is to call
   * MakeS2ContainsPointQuery() instead.  For example:
   *
   *   return MakeS2ContainsPointQuery(&index).Contains(p);
   */
  alias Options = S2ContainsPointQueryOptions;

  this(IndexT index, Options options = null) {
    _index = index;
    _options = options !is null ? options : new Options();
    _it = new Iterator(_index);
  }


  // Convenience constructor that accepts the S2VertexModel directly.
  this(IndexT index, S2VertexModel vertex_model) {
    this(index, new Options(vertex_model));
  }

  const(IndexT) index() const {
    return _index;
  }

  const(Options) options() const {
    return _options;
  }

  /// Equivalent to the two-argument constructor above.
  void initialize(IndexT index, Options options = null) {
    _index = index;
    _options = options !is null ? options : new Options();
    _it.initialize(index);
  }

  /**
   * Returns true if any shape in the given index() contains the point "p"
   * under the vertex model specified (OPEN, SEMI_OPEN, or CLOSED).
   */
  bool contains(in S2Point p) {
    if (!_it.locate(p)) return false;

    const S2ShapeIndexCell cell = _it.cell();
    int num_clipped = cell.numClipped();
    for (int s = 0; s < num_clipped; ++s) {
      if (shapeContains(_it, cell.clipped(s), p)) return true;
    }
    return false;
  }

  /**
   * Returns true if the given shape contains the point "p" under the vertex
   * model specified (OPEN, SEMI_OPEN, or CLOSED).
   *
   * REQUIRES: "shape" belongs to index().
   */
  bool shapeContains(in S2Shape shape, in S2Point p) {
    if (!_it.locate(p)) return false;
    const(S2ClippedShape)* clipped = _it.cell().findClipped(shape.id());
    if (clipped is null) return false;
    return shapeContains(_it, *clipped, p);
  }

  // Visits all shapes in the given index() that contain the given point "p",
  // terminating early if the given ShapeVisitor function returns false (in
  // which case VisitContainingShapes returns false as well).  Each shape is
  // visited at most once.
  //
  // Note that the API allows non-const access to the visited shapes.
  alias ShapeVisitor = bool delegate(S2Shape);

  bool visitContainingShapes(in S2Point p, ShapeVisitor visitor) {
    // This function returns "false" only if the algorithm terminates early
    // because the "visitor" function returned false.
    if (!_it.locate(p)) return true;

    const(S2ShapeIndexCell) cell = _it.cell();
    int num_clipped = cell.numClipped();
    for (int s = 0; s < num_clipped; ++s) {
      const(S2ClippedShape) clipped = cell.clipped(s);
      if (shapeContains(_it, clipped, p) && !visitor(_index.shape(clipped.shapeId()))) {
        return false;
      }
    }
    return true;
  }

  /// Convenience function that returns all the shapes that contain the given point "p".
  S2Shape[] getContainingShapes(in S2Point p) {
    S2Shape[] results;
    visitContainingShapes(p, (S2Shape shape) {
          results ~= shape;
          return true;
        });
    return results;
  }

  // Visits all edges in the given index() that are incident to the point "p"
  // (i.e., "p" is one of the edge endpoints), terminating early if the given
  // EdgeVisitor function returns false (in which case VisitIncidentEdges
  // returns false as well).  Each edge is visited at most once.
  alias EdgeVisitor = bool delegate(in ShapeEdge);

  bool visitIncidentEdges(in S2Point p, in EdgeVisitor visitor) {
    // This function returns "false" only if the algorithm terminates early
    // because the "visitor" function returned false.
    if (!_it.locate(p)) return true;

    const S2ShapeIndexCell cell = _it.cell();
    int num_clipped = cell.numClipped();
    for (int s = 0; s < num_clipped; ++s) {
      const S2ClippedShape clipped = cell.clipped(s);
      int num_edges = clipped.numEdges();
      if (num_edges == 0) continue;
      const S2Shape shape = _index.shape(clipped.shapeId());
      for (int i = 0; i < num_edges; ++i) {
        int edge_id = clipped.edge(i);
        auto edge = shape.edge(edge_id);
        if ((edge.v0 == p || edge.v1 == p) && !visitor(new ShapeEdge(shape.id(), edge_id, edge))) {
          return false;
        }
      }
    }
    return true;
  }

  /////////////////////////// Low-Level Methods ////////////////////////////
  //
  // Most clients will not need the following methods.  They can be slightly
  // more efficient but are harder to use.

  // Returns a pointer to the iterator used internally by this class, in order
  // to avoid the need for clients to create their own iterator.  Clients are
  // allowed to reposition this iterator arbitrarily between method calls.
  Iterator mutableIter() {
    return _it;
  }

  // Low-level helper method that returns true if the given S2ClippedShape
  // referred to by an S2ShapeIndex::Iterator contains the point "p".
  bool shapeContains(in Iterator it, in S2ClippedShape clipped, in S2Point p) const {
    bool inside = clipped.containsCenter();
    int num_edges = clipped.numEdges();
    if (num_edges > 0) {
      // Points and polylines can be ignored unless the vertex model is CLOSED.
      const S2Shape shape = _index.shape(clipped.shapeId());
      if (!shape.hasInterior() && _options.vertexModel() != S2VertexModel.CLOSED) {
        return false;
      }
      // Test containment by drawing a line segment from the cell center to the
      // given point and counting edge crossings.
      auto crosser = new S2CopyingEdgeCrosser(it.center(), p);
      for (int i = 0; i < num_edges; ++i) {
        auto edge = shape.edge(clipped.edge(i));
        int sign = crosser.crossingSign(edge.v0, edge.v1);
        if (sign < 0) continue;
        if (sign == 0) {
          // For the OPEN and CLOSED models, check whether "p" is a vertex.
          if (_options.vertexModel() != S2VertexModel.SEMI_OPEN
              && (edge.v0 == p || edge.v1 == p)) {
            return _options.vertexModel() == S2VertexModel.CLOSED;
          }
          sign = vertexCrossing(crosser.a(), crosser.b(), edge.v0, edge.v1);
        }
        inside ^= cast(bool) sign;
      }
    }
    return inside;
  }

 private:
  IndexT _index;
  Options _options;
  Iterator _it;
}

// Returns an S2ContainsPointQuery for the given S2ShapeIndex.  Note that
// it is efficient to return S2ContainsPointQuery objects by value.
S2ContainsPointQuery!IndexT makeS2ContainsPointQuery(IndexT)(
    IndexT index, S2ContainsPointQueryOptions options = null) {
  return new S2ContainsPointQuery!IndexT(
      index, options !is null ? options : new S2ContainsPointQueryOptions());
}
