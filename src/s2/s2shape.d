// Copyright 2012 Google Inc. All Rights Reserved.
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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2shape;

import s2.s2point;
import s2.s2pointutil : origin;

// The purpose of S2Shape is to represent polygonal geometry in a flexible
// way.  It is organized as a collection of edges that optionally defines an
// interior.  All geometry represented by an S2Shape must have the same
// dimension, which means that an S2Shape can represent either a set of
// points, a set of polylines, or a set of polygons.
//
// S2Shape is defined as an abstract base class in order to give clients
// control over the underlying data representation.  Sometimes an S2Shape does
// not have any data of its own, but instead "wraps" some other class.  There
// are various useful subtypes defined in *_shape.h, and some S2 classes also
// have a nested "Shape" class (e.g., S2Polygon::Shape).  It is easy for
// clients to implement their own subtypes, since the interface is minimal.
//
// S2Shape operations are typically defined on S2ShapeIndex objects rather
// than individual shapes.  An S2ShapeIndex is simply a collection of
// S2Shapes, possibly of different dimensions (e.g. 10 points and 3 polygons),
// organized into a data structure for efficient edge access.
//
// The edges of an S2Shape are identified by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, an S2Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an S2Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BooleanOperation).
abstract class S2Shape {
public:
  // An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
  // allowed, and can be used to represent points.
  struct Edge {
    S2Point v0, v1;

    bool opEquals(in Edge e) const {
      return v0 == e.v0 && v1 == e.v1;
    }

    int opCmp(in Edge e) const {
      if (v0 < e.v0 || (v0 == e.v0 && v1 < e.v1)) {
        return -1;
      } else if (v0 > e.v0 || (v0 == e.v0 && v1 > e.v1)) {
        return 1;
      }
      return 0;
    }
  }

  // A range of edge ids corresponding to a chain of zero or more connected
  // edges, specified as a (start, length) pair.  The chain is defined to
  // consist of edge ids {start, start + 1, ..., start + length - 1}.
  struct Chain {
    int start, length;
  }

  // The position of an edge within a given edge chain, specified as a
  // (chain_id, offset) pair.  Chains are numbered sequentially starting from
  // zero, and offsets are measured from the start of each chain.
  struct ChainPosition {
    int chainId, offset;
  }

  // A ReferencePoint consists of a point P and a boolean indicating whether P
  // is contained by a particular shape.
  struct ReferencePoint {
    S2Point point;
    bool contained;

    // Returns a ReferencePoint with the given "contained" value and a default
    // "point".  It should be used when all points or no points are contained.
    // (Was called "contained")
    static ReferencePoint defaultReferencePoint(bool _contained) {
      return ReferencePoint(origin(), _contained);
    }
  }

  this() {
    _id = -1;
  }

  // Returns the number of edges in this shape.  Edges have ids ranging from 0
  // to num_edges() - 1.
  abstract int numEdges() const;

  // Returns the endpoints of the given edge id.
  //
  // REQUIRES: 0 <= id < num_edges()
  abstract Edge edge(int edge_id) const
  in {
    assert(0 <= edge_id && edge_id < numEdges());
  }

  // Returns the dimension of the geometry represented by this shape.
  //
  //  0 - Point geometry.  Each point is represented as a degenerate edge.
  //
  //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
  //      represent any number of polylines.  Polylines edges may intersect.
  //
  //  2 - Polygon geometry.  Edges should be oriented such that the polygon
  //      interior is always on the left.  In theory the edges may be returned
  //      in any order, but typically the edges are organized as a collection
  //      of edge chains where each chain represents one polygon loop.
  //      Polygons may have degeneracies (e.g., degenerate edges or sibling
  //      pairs consisting of an edge and its corresponding reversed edge).
  //
  // Note that this method allows degenerate geometry of different dimensions
  // to be distinguished, e.g. it allows a point to be distinguished from a
  // polyline or polygon that has been simplified to a single point.
  abstract int dimension() const;

  // Convenience function that returns true if this shape has an interior.
  bool hasInterior() const {
    return dimension() == 2;
  }

  // Returns an arbitrary point P along with a boolean indicating whether P is
  // contained by the shape.  (The boolean value must be false for shapes that
  // do not have an interior.)
  //
  // This ReferencePoint may then be used to compute the containment of other
  // points by counting edge crossings.
  abstract ReferencePoint getReferencePoint() const;

  // Returns the number of contiguous edge chains in the shape.  For example,
  // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
  // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
  // sequentially starting from zero.
  //
  // Note that it is always acceptable to implement this method by returning
  // num_edges() (i.e. every chain consists of a single edge), but this may
  // reduce the efficiency of some algorithms.
  abstract int numChains() const;

  // Returns the range of edge ids corresponding to the given edge chain.  The
  // edge chains must form contiguous, non-overlapping ranges that cover the
  // entire range of edge ids.  This is spelled out more formally below:
  //
  // REQUIRES: 0 <= i < num_chains()
  // REQUIRES: chain(i).length >= 0, for all i
  // REQUIRES: chain(0).start == 0
  // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
  //           for i < num_chains() - 1
  // REQUIRES: chain(i).start + chain(i).length == num_edges(),
  //           for i == num_chains() - 1
  abstract Chain chain(int chain_id) const
  in {
    assert(0 <= chain_id && chain_id < numChains());
  } out (c) {
    assert(c.length >= 0);
    if (chain_id == numChains() - 1) {
      assert(c.start + c.length == numEdges());
    }
  }

  // Returns the edge at offset "offset" within edge chain "chain_id".
  // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
  // but may be more efficient.
  abstract Edge chainEdge(int chain_id, int offset) const;

  // Finds the chain containing the given edge, and returns the position of
  // that edge as a (chain_id, offset) pair.
  //
  // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
  // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
  //
  // where     pos == shape.chain_position(edge_id).
  abstract ChainPosition chainPosition(int edge_id) const;

  // A unique id assigned to this shape by S2ShapeIndex.  Shape ids are
  // assigned sequentially starting from 0 in the order shapes are added.
  int id() const { return _id; }

  // Virtual methods that return pointers of your choice.  These methods are
  // intended to help with the problem of attaching additional data to S2Shape
  // objects.  For example, you could return a pointer to a source object, or
  // a pointer to a bundle of additional data allocated with the S2Shape.
  // Because this method exists in all S2Shapes, you can override it in each
  // type of shape you have and call it without knowing the concrete subtype.
  // For example, if you have polyline and polygon shapes, you can do this:
  //
  //   class MyPolyline : S2Polyline.Shape {
  //   public:
  //     void* mutableUserData() { return &_myData; }
  //   private:
  //     MyData _myData;
  //   }
  //   class MyPolygon : S2Polygon.Shape {
  //   public:
  //     void* mutableUserData() { return &_myData; }
  //   private:
  //     MyData _myData;
  //   };
  //   ...
  //   S2Shape* shape = index.shape(id);
  //   const(MyData)* data = cast(const(MyData)*)(shape.userData());
  //
  // This is not the only way to map from an S2Shape back to your source
  // data.  Other reasonable techniques include:
  //
  //  - Every shape has an id() assigned by S2ShapeIndex.  Ids are assigned
  //    sequentially starting from 0 in the order the shapes are added to the
  //    index.  You can use this id to look up arbitrary data stored in your
  //    own vector.
  //
  //  - If all of your shapes are the same type, then you can create your own
  //    subclass of some existing S2Shape type (such as S2Polyline::Shape) and
  //    add your own methods and fields.  You can access this data by
  //    downcasting the S2Shape pointers returned by S2ShapeIndex methods.
  const(void*) userData() const {
    return null;
  }
  void* mutableUserData() {
    return null;
  }

  // Assigned by MutableS2ShapeIndex when the shape is added.
  package int _id;
}
