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
// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_point_query_test.cc.

module s2.s2closest_point_query_base_test;

import s2.s2closest_point_query_base;

import s2.logger;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2contains_point_query;
import s2.s2distance_target;
import s2.s2edge_distances;
import s2.s2point;
import s2.s2point_index;
import s2.s2shape_index;
import s2.s2text_format;

import fluent.asserts;


// This is a proof-of-concept prototype of a possible S2FurthestPointQuery
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
  alias Delta = S1ChordAngle;

  this(S1ChordAngle x) {
    _distance = x;
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

  bool opEquals(MaxDistance y) const {
    return _distance == y._distance;
  }

  int opCmp(MaxDistance y) const {
    // Reverses the normal comparison direction.
    if (_distance > y._distance) return -1;
    if (_distance < y._distance) return 1;
    return 0;
  }

  MaxDistance opBinary(string op)(S1ChordAngle delta) {
    return MaxDistance(_distance + delta);
  }

  S1ChordAngle getDistance() {
    return _distance;
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
  S1ChordAngle _distance;
}


alias FurthestPointQueryOptions = S2ClosestPointQueryBaseOptions!MaxDistance;

final class FurthestPointQueryTarget : S2DistanceTarget!MaxDistance {
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
    logger.logFatal("Unimplemented");
    return false;
  }

  override
  bool updateMinDistance(in S2Cell cell, ref MaxDistance min_dist) {
    return min_dist.updateMin(MaxDistance(S1ChordAngle.straight() - cell.getDistance(-_point)));
  }

  override
  bool visitContainingShapes(S2ShapeIndex index, ShapeVisitor visitor) {
    logger.logFatal("Unimplemented");
    return false;
  }

private:
  S2Point _point;
}

template FurthestPointQuery(Data) {
  alias FurthestPointQuery = S2ClosestPointQueryBase!(MaxDistance, Data);
}

@("S2ClosestPointQueryBase.MaxDistance") unittest {
  auto index = new S2PointIndex!int();
  auto points = parsePointsOrDie("0:0, 1:0, 2:0, 3:0");
  for (int i = 0; i < points.length; ++i) {
    index.add(points[i], i);
  }
  auto query = new FurthestPointQuery!int(index);
  auto target = new FurthestPointQueryTarget(makePointOrDie("4:0"));
  auto options = new FurthestPointQueryOptions();
  options.setMaxPoints(1);
  auto results = query.findClosestPoints(target, options);
  Assert.equal(results.length, 1);
  Assert.equal(results[0].point(), points[0]);
  Assert.equal(results[0].data(), 0);
  Assert.approximately(
      S1ChordAngle(results[0].distance().getDistance()).toS1Angle().degrees(), 4.0, 1e-13);
}
