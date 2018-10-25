// Copyright 2013 Google Inc. All Rights Reserved.
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
// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_edge_query_test.cc.

module s2.s2closest_edge_query_base_test;

import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2closest_edge_query_base;
import s2.s2contains_point_query;
import s2.s2shape;
import s2.s2edge_distances : updateMinDistance;
import s2.s2point;
import s2.s2shape_index;
import s2.s2text_format;

// This is a proof-of-concept prototype of a possible S2FurthestEdgeQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.  (A real implementation would need
// to be more careful about error bounds, it would implement a greater range
// of target types, etc.)
//
// It is based on the principle that for any two geometric objects X and Y on
// the sphere,
//
//     max_dist(X, Y) = Pi - min_dist(-X, Y)
//
// where "-X" denotes the reflection of X through the origin (i.e., to the
// opposite side of the sphere).

// MaxDistance is a class that allows maximum distances to be computed using a
// minimum distance algorithm.  It essentially treats a distance "x" as the
// supplementary distance (Pi - x).
struct MaxDistance {
public:
  // Expose only the methods that are documented as necessary in the API.
  struct Delta {
  public:

    bool opEquals(Delta y) {
      return _a == y._a;
    }

    static Delta zero() {
      Delta r;
      r._a = S1ChordAngle.zero();
      return r;
    }

    // This method is needed to implement Distance::operator-.
    S1ChordAngle opCast(T : S1ChordAngle)() const {
      return _a;
    }

  private:
    S1ChordAngle _a;
  }

  S1ChordAngle getChordAngle() const {
    return _distance;
  }

  static MaxDistance zero() {
    return MaxDistance(S1ChordAngle.straight());
  }

  static MaxDistance infinity() {
    return MaxDistance(S1ChordAngle.negative());
  }

  static MaxDistance negative() {
    return MaxDistance(S1ChordAngle.infinity());
  }

  int opCmp(in MaxDistance y) const {
    return _distance.opCmp(y._distance);
  }

  MaxDistance opBinary(string op)(in Delta delta) const
  if (op == "-") {
    return MaxDistance(_distance + S1ChordAngle(delta._a));
  }

  S1ChordAngle getChordAngleBound() const {
    return S1ChordAngle.straight() - _distance;
  }

  // If (dist < *this), updates *this and returns true (used internally).
  bool updateMin(in MaxDistance dist) {
    if (dist < this) {
      this = dist;
      return true;
    }
    return false;
  }

private:
  S1ChordAngle _distance = S1ChordAngle();
}

alias FurthestEdgeQuery = S2ClosestEdgeQueryBase!MaxDistance;

final class FurthestPointTarget : FurthestEdgeQuery.Target {
public:
  this(in S2Point point) {
    _point = point;
  }

  override
  int maxBruteForceIndexSize() const {
    return 100;
  }

  override
  S2Cap getCapBound() {
    return new S2Cap(-_point, S1ChordAngle.zero());
  }

  override
  bool updateMinDistance(in S2Point p, ref MaxDistance min_dist) {
    return min_dist.updateMin(MaxDistance(S1ChordAngle(p, _point)));
  }

  override
  bool updateMinDistance(in S2Point v0, in S2Point v1, ref MaxDistance min_dist) {
    S1ChordAngle dist180 = S1ChordAngle(min_dist.getChordAngle()).isNegative()
        ? S1ChordAngle.infinity()
        : S1ChordAngle.straight() - S1ChordAngle(min_dist.getChordAngle());
    if (!.updateMinDistance(-_point, v0, v1, dist180)) return false;
    min_dist = MaxDistance(S1ChordAngle.straight() - dist180);
    return true;
  }

  override
  bool updateMinDistance(in S2Cell cell, ref MaxDistance min_dist) {
    return min_dist.updateMin(MaxDistance(S1ChordAngle.straight() - cell.getDistance(-_point)));
  }

  override
  bool visitContainingShapes(S2ShapeIndex index, FurthestEdgeQuery.Target.ShapeVisitor visitor) {
    // For furthest points, we return the polygons whose interior contains the
    // antipode of the target point.  (These are the polygons whose
    // MaxDistance() to the target is MaxDistance::Zero().)
    //
    // For target types consisting of multiple connected components (such as
    // FurthestPointQuery::ShapeIndexTarget), this method should return the
    // polygons containing the antipodal reflection of any connected
    // component.  (It is sufficient to test containment of one vertex per
    // connected component, since the API allows us to also return any polygon
    // whose boundary has MaxDistance::Zero() to the target.)
    return makeS2ContainsPointQuery(index).visitContainingShapes(
        -_point, (in S2Shape shape) {
          return visitor(shape, _point);
        });
  }

private:
  S2Point _point;
}

/+ TODO: Add when makeIndex is implemented.
@("S2ClosestEdgeQueryBase.MaxDistance") unittest {
  auto index = makeIndex("0:0 | 1:0 | 2:0 | 3:0 # #");
  FurthestEdgeQuery query(index.get());
  FurthestPointTarget target(s2textformat::MakePoint("4:0"));
  FurthestEdgeQuery::Options options;
  options.set_max_edges(1);
  auto results = query.FindClosestEdges(&target, options);
  ASSERT_EQ(1, results.size());
  EXPECT_EQ(0, results[0].shape_id);
  EXPECT_EQ(0, results[0].edge_id);
  EXPECT_NEAR(4, S1ChordAngle(results[0].distance).ToAngle().degrees(), 1e-13);
}
+/
