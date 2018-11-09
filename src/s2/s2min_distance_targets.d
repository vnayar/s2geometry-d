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
//
// This file defines a collection of classes that are useful for computing
// minimum distances on the sphere.  Their purpose is to allow code to be
// shared among the various query classes that find nearby geometry, such as
// S2ClosestPointQuery and S2ClosestEdgeQuery.

module s2.s2min_distance_targets;

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2contains_point_query;
import s2.s2distance_target;
import s2.s2edge_distances : updateMinDistance, updateEdgePairMinDistance;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2closest_edge_query;

import std.math : sqrt;

/**
 * S2MinDistance is a thin wrapper around S1ChordAngle that is used by classes
 * such as S2ClosestEdgeQuery to compute minimum distances on the sphere (as
 * opposed to maximum distances, ellipsoidal distances, etc).
 *
 * It implements the Distance concept defined by S2DistanceTarget, which is
 * used by query classes such as S2ClosestPointQuery and S2ClosestEdgeQuery.
 * (See s2distance_target.h for details.)
 */
struct S2MinDistance {
public:
  S1ChordAngle s1ChordAngle;

  alias s1ChordAngle this;

  alias Delta = S1ChordAngle;

  static S2MinDistance zero() {
    return S2MinDistance(S1ChordAngle.zero());
  }

  static S2MinDistance infinity() {
    return S2MinDistance(S1ChordAngle.infinity());
  }

  static S2MinDistance negative() {
    return S2MinDistance(S1ChordAngle.negative());
  }

  S1ChordAngle getChordAngleBound() const {
    return plusError(getS1AngleConstructorMaxError());
  }

  // If (dist < *this), updates *this and returns true (used internally).
  bool updateMin(in S2MinDistance dist) {
    if (dist < this) {
      this = dist;
      return true;
    }
    return false;
  }

  int opCmp(in S2MinDistance other) const {
    return s1ChordAngle.opCmp(other.s1ChordAngle);
  }
}

// S2MinDistanceTarget represents a geometric object to which distances are
// measured.  Specifically, it is used to compute minimum distances on the
// sphere (as opposed to maximum distances, ellipsoidal distances, etc).
//
// Subtypes are defined below for measuring the distance to a point, an edge,
// an S2Cell, or an S2ShapeIndex (an arbitrary collection of geometry).
alias S2MinDistanceTarget = S2DistanceTarget!S2MinDistance;

// An S2DistanceTarget subtype for computing the minimum distance to a point.
class S2MinDistancePointTarget : S2MinDistanceTarget {
public:
  this(in S2Point point) {
    _point = point;
  }

  final override
  S2Cap getCapBound() {
    return new S2Cap(_point, S1ChordAngle.zero());
  }

  final override
  bool updateMinDistance(in S2Point p, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(S1ChordAngle(p, _point)));
  }

  final override
  bool updateMinDistance(in S2Point v0, in S2Point v1, ref S2MinDistance min_dist) {
    return .updateMinDistance(_point, v0, v1, min_dist);
  }

  final override
  bool updateMinDistance(in S2Cell cell, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(cell.getDistance(_point)));
  }

  final override
  bool visitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor) {
    return makeS2ContainsPointQuery(index).visitContainingShapes(
        _point, (in S2Shape shape) {
          return visitor(shape, _point);
        });
  }

private:
  S2Point _point;
}

// An S2DistanceTarget subtype for computing the minimum distance to a edge.
class S2MinDistanceEdgeTarget : S2MinDistanceTarget {
public:
  this(in S2Point a, in S2Point b) {
    _a = a;
    _b = b;
  }

  final override
  S2Cap getCapBound() {
    // The following computes a radius equal to half the edge length in an
    // efficient and numerically stable way.
    double d2 = S1ChordAngle(_a, _b).length2();
    double r2 = (0.5 * d2) / (1 + sqrt(1 - 0.25 * d2));
    return new S2Cap((_a + _b).normalize(), S1ChordAngle.fromLength2(r2));
  }

  final override
  bool updateMinDistance(in S2Point p, ref S2MinDistance min_dist) {
    return .updateMinDistance(p, _a, _b, min_dist);
  }

  final override
  bool updateMinDistance(in S2Point v0, in S2Point v1, ref S2MinDistance min_dist) {
    return updateEdgePairMinDistance(_a, _b, v0, v1, min_dist);
  }

  final override
  bool updateMinDistance(in S2Cell cell, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(cell.getDistance(_a, _b)));
  }

  final override
  bool visitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor) {
    // We test the center of the edge in order to ensure that edge targets AB
    // and BA yield identical results (which is not guaranteed by the API but
    // users might expect).  Other options would be to test both endpoints, or
    // return different results for AB and BA in some cases.
    scope target = new S2MinDistancePointTarget((_a + _b).normalize());
    return target.visitContainingShapes(index, visitor);
  }

private:
  S2Point _a;
  S2Point _b;
}

// An S2DistanceTarget subtype for computing the minimum distance to an S2Cell
// (including the interior of the cell).
class S2MinDistanceCellTarget : S2MinDistanceTarget {
public:
  this(in S2Cell cell) {
    _cell = cell;
  }

  final override
  S2Cap getCapBound() {
    return _cell.getCapBound();
  }

  final override
  bool updateMinDistance(in S2Point p, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(_cell.getDistance(p)));
  }

  final override
  bool updateMinDistance(in S2Point v0, in S2Point v1, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(_cell.getDistance(v0, v1)));
  }

  final override
  bool updateMinDistance(in S2Cell cell, ref S2MinDistance min_dist) {
    return min_dist.updateMin(S2MinDistance(_cell.getDistance(cell)));
  }

  final override
  bool visitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor) {
    // The simplest approach is simply to return the polygons that contain the
    // cell center.  Alternatively, if the index cell is smaller than the target
    // cell then we could return all polygons that are present in the
    // S2ShapeIndexCell, but since the index is built conservatively this may
    // include some polygons that don't quite intersect the cell.  So we would
    // either need to recheck for intersection more accurately, or weaken the
    // VisitContainingShapes contract so that it only guarantees approximate
    // intersection, neither of which seems like a good tradeoff.
    scope target = new S2MinDistancePointTarget(_cell.getCenter());
    return target.visitContainingShapes(index, visitor);
  }

private:
  const(S2Cell) _cell;
}

/**
 * An S2DistanceTarget subtype for computing the minimum distance to an
 * S2ShapeIndex (a collection of points, polylines, and/or polygons).
 *
 * Note that ShapeIndexTarget has its own options:
 *
 *   include_interiors()
 *     - specifies that distances are measured to the boundary and interior
 *       of polygons in the S2ShapeIndex.  (If set to false, distance is
 *       measured to the polygon boundary only.)
 *       DEFAULT: true.
 *
 *   brute_force()
 *     - specifies that the distances should be computed by examining every
 *       edge in the S2ShapeIndex (for testing and debugging purposes).
 *       DEFAULT: false.
 *
 * These options are specified independently of the corresponding
 * S2ClosestEdgeQuery options.  For example, if include_interiors is true for
 * an ShapeIndexTarget but false for the S2ClosestEdgeQuery where the target
 * is used, then distances will be measured from the boundary of one
 * S2ShapeIndex to the boundary and interior of the other.
 *
 * Note that when the distance to a ShapeIndexTarget is zero because the
 * target intersects the interior of the query index, you can find a point
 * that achieves this zero distance by calling the VisitContainingShapes()
 * method directly.  For example:
 *
 *   S2ClosestEdgeQuery::ShapeIndexTarget target(&target_index);
 *   target.VisitContainingShapes(
 *       query_index, [](S2Shape* containing_shape,
 *                       const S2Point& target_point) {
 *         ... do something with "target_point" ...
 *         return false;  // Terminate search
 *       }));
 */
class S2MinDistanceShapeIndexTarget : S2MinDistanceTarget {
public:
  this(S2ShapeIndex index) {
    _index = index;
    _query = new S2ClosestEdgeQuery(index);
  }

  // Specifies that distance will be measured to the boundary and interior
  // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
  //
  // DEFAULT: true
  bool includeInteriors() const {
    return _query.options().includeInteriors();
  }

  void setIncludeInteriors(bool include_interiors) {
    _query.mutableOptions().setIncludeInteriors(include_interiors);
  }

  // Specifies that the distances should be computed by examining every edge
  // in the S2ShapeIndex (for testing and debugging purposes).
  //
  // DEFAULT: false
  bool useBruteForce() const {
    return _query.options().useBruteForce();
  }

  void setUseBruteForce(bool use_brute_force) {
    _query.mutableOptions().setUseBruteForce(use_brute_force);
  }

  override
  bool setMaxError(in S1ChordAngle max_error) {
    _query.mutableOptions().setMaxError(max_error.toS1Angle());
    return true;  // Indicates that we may return suboptimal results.
  }

  final override
  S2Cap getCapBound() {
    // TODO: Resume when s2shape_index_region is complete.
    //return makeS2ShapeIndexRegion(_index).getCapBound();
    return null;
  }

  final override
  bool updateMinDistance(in S2Point p, ref S2MinDistance min_dist) {
    _query.mutableOptions().setMaxDistance(min_dist);
    auto target = new S2ClosestEdgeQuery.PointTarget(p);
    S2ClosestEdgeQuery.Result r = _query.findClosestEdge(target);
    if (r.shapeId < 0) return false;
    min_dist = r.distance;
    return true;
  }

  final override
  bool updateMinDistance(in S2Point v0, in S2Point v1, ref S2MinDistance min_dist) {
    _query.mutableOptions().setMaxDistance(min_dist);
    auto target = new S2ClosestEdgeQuery.EdgeTarget(v0, v1);
    S2ClosestEdgeQuery.Result r = _query.findClosestEdge(target);
    if (r.shapeId < 0) return false;
    min_dist = r.distance;
    return true;
  }

  final override
  bool updateMinDistance(in S2Cell cell, ref S2MinDistance min_dist) {
    _query.mutableOptions().setMaxDistance(min_dist);
    auto target = new S2ClosestEdgeQuery.CellTarget(cell);
    S2ClosestEdgeQuery.Result r = _query.findClosestEdge(target);
    if (r.shapeId < 0) return false;
    min_dist = r.distance;
    return true;
  }

  final override
  bool visitContainingShapes(S2ShapeIndex query_index, in ShapeVisitor visitor) {
    // It is sufficient to find the set of chain starts in the target index
    // (i.e., one vertex per connected component of edges) that are contained by
    // the query index, except for one special case to handle full polygons.
    //
    // TODO(ericv): Do this by merge-joining the two S2ShapeIndexes, and share
    // the code with S2BooleanOperation.

    int num_shape_ids = _index.numShapeIds();
    for (int s = 0; s < num_shape_ids; ++s) {
      const(S2Shape) shape = _index.shape(s);
      if (shape is null) continue;
      int num_chains = shape.numChains();
      // Shapes that don't have any edges require a special case (below).
      bool tested_point = false;
      for (int c = 0; c < num_chains; ++c) {
        S2Shape.Chain chain = shape.chain(c);
        if (chain.length == 0) continue;
        tested_point = true;
        S2Point v0 = shape.chainEdge(c, 0).v0;
        auto target = new S2MinDistancePointTarget(v0);
        if (!target.visitContainingShapes(query_index, visitor)) {
          return false;
        }
      }
      if (!tested_point) {
        // Special case to handle full polygons.
        S2Shape.ReferencePoint refPoint = shape.getReferencePoint();
        if (!refPoint.contained) continue;
        auto target = new S2MinDistancePointTarget(refPoint.point);
        if (!target.visitContainingShapes(query_index, visitor)) {
          return false;
        }
      }
    }
    return true;
  }

private:
  const(S2ShapeIndex) _index;
  S2ClosestEdgeQuery _query;
}
