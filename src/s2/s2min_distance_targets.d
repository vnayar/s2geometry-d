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
private:
  S1ChordAngle s1ChordAngle;

public:
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

/+ TODO: Resume when S2ClosestEdgeQuery is implemented.
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

  ~S2MinDistanceShapeIndexTarget() override;

  // Specifies that distance will be measured to the boundary and interior
  // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
  //
  // DEFAULT: true
  bool include_interiors() const;
  void set_include_interiors(bool include_interiors);

  // Specifies that the distances should be computed by examining every edge
  // in the S2ShapeIndex (for testing and debugging purposes).
  //
  // DEFAULT: false
  bool use_brute_force() const;
  void set_use_brute_force(bool use_brute_force);

  bool set_max_error(const S1ChordAngle& max_error) override;
  S2Cap GetCapBound() final;
  bool UpdateMinDistance(const S2Point& p, S2MinDistance* min_dist) final;
  bool UpdateMinDistance(const S2Point& v0, const S2Point& v1,
                         S2MinDistance* min_dist) final;
  bool UpdateMinDistance(const S2Cell& cell,
                         S2MinDistance* min_dist) final;
  bool VisitContainingShapes(const S2ShapeIndex& query_index,
                             const ShapeVisitor& visitor) final;

 private:
  const S2ShapeIndex* index_;
  std::unique_ptr<S2ClosestEdgeQuery> query_;
}
+/
