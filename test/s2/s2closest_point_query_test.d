// Copyright 2015 Google Inc. All Rights Reserved.
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

module s2.s2closest_point_query_test;

import s2.s2closest_point_query;

import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2edge_distances;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2loop;
import s2.s2point;
import s2.s2point_index;
import s2.s2pointutil;
import s2.s2testing;
import s2.util.math.matrix3x3;

import fluent.asserts;

import std.algorithm;
import std.conv;
import std.math;

alias TestIndex = S2PointIndex!int;
alias TestQuery = S2ClosestPointQuery!int;

@("S2ClosestPointQuery.NoPoints") unittest {
  auto index = new TestIndex();
  auto query = new TestQuery(index);
  auto target = new S2ClosestPointQueryPointTarget(S2Point(1, 0, 0));
  auto results = query.findClosestPoints(target);
  Assert.equal(results.length, 0);
}

@("S2ClosestPointQuery.ManyDuplicatePoints") unittest {
  enum int kNumPoints = 10000;
  auto kTestPoint = S2Point(1, 0, 0);
  auto index = new TestIndex();
  for (int i = 0; i < kNumPoints; ++i) {
    index.add(kTestPoint, i);
  }
  auto query = new TestQuery(index);
  auto target = new S2ClosestPointQueryPointTarget(kTestPoint);
  const auto results = query.findClosestPoints(target);
  Assert.equal(results.length, kNumPoints);
}

// An abstract class that adds points to an S2PointIndex for benchmarking.
interface PointIndexFactory {
public:
  // Given an index that will be queried using random points from "query_cap",
  // adds approximately "num_points" points to "index".  (Typically the
  // indexed points will occupy some fraction of this cap.)
  void addPoints(in S2Cap query_cap, int num_points, TestIndex index);
}

// Generates points that are regularly spaced along a circle.  The circle is
// centered within the query cap and occupies 25% of its area, so that random
// query points have a 25% chance of being inside the circle.
//
// Points along a circle are nearly the worst case for distance calculations,
// since many points are nearly equidistant from any query point that is not
// immediately adjacent to the circle.
class CirclePointIndexFactory : PointIndexFactory {
  override
  void addPoints(in S2Cap query_cap, int num_points, TestIndex index) const {
    S2Point[] points = S2Testing.makeRegularPoints(
        query_cap.center(), 0.5 * query_cap.getRadius(), num_points);
    for (int i = 0; i < points.length; ++i) {
      index.add(points[i], i);
    }
  }
}

// Generate the vertices of a fractal whose convex hull approximately
// matches the query cap.
class FractalPointIndexFactory : PointIndexFactory {
  override
  void addPoints(in S2Cap query_cap, int num_points, TestIndex index) const {
    auto fractal = new S2Testing.Fractal();
    fractal.setLevelForApproxMaxEdges(num_points);
    fractal.setFractalDimension(1.5);
    auto loop = fractal.makeLoop(S2Testing.getRandomFrameAt(query_cap.center()), query_cap.getRadius());
    for (int i = 0; i < loop.numVertices(); ++i) {
      index.add(loop.vertex(i), i);
    }
  }
}

// Generate vertices on a square grid that includes the entire query cap.
class GridPointIndexFactory : PointIndexFactory {
  override
  void addPoints(in S2Cap query_cap, int num_points, TestIndex index) const {
    int sqrt_num_points = cast(int) ceil(sqrt(cast(double) num_points));
    Matrix3x3_d frame = S2Testing.getRandomFrameAt(query_cap.center());
    double radius = query_cap.getRadius().radians();
    double spacing = 2 * radius / sqrt_num_points;
    for (int i = 0; i < sqrt_num_points; ++i) {
      for (int j = 0; j < sqrt_num_points; ++j) {
        auto point = S2Point(tan((i + 0.5) * spacing - radius),
                      tan((j + 0.5) * spacing - radius), 1.0);
        index.add(fromFrame(frame, point.normalize()), i * sqrt_num_points + j);
      }
    }
  }
}

// The approximate radius of S2Cap from which query points are chosen.
enum S1Angle kRadius = S2Testing.kmToAngle(10);

// An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
enum double kChordAngleError = 1e-15;

// TODO(user): Remove Result and use TestQuery::Result.
alias Result = Pair!(S1Angle, int);

private S1Angle getDistance(TestQuery.Target target, in S2Point point) {
  TestQuery.Distance distance = TestQuery.Distance.infinity();
  target.updateMinDistance(point, distance);
  return distance.toS1Angle();
}

// Use "query" to find the closest point(s) to the given target, and extract
// the query results into the given vector.  Also verify that the results
// satisfy the search criteria.
private void getClosestPoints(TestQuery.Target target, TestQuery query, ref Result[] results) {
  const query_results = query.findClosestPoints(target);
  Assert.notGreaterThan(query_results.length, query.options().maxPoints());
  if (!query.options().region() && query.options().maxDistance() == S1ChordAngle.infinity()) {
    // We can predict exactly how many points should be returned.
    Assert.equal(
        query_results.length, min(query.options().maxPoints(), query.index().numPoints()));
  }
  foreach (const result; query_results) {
    // Check that query->distance() is approximately equal to the
    // S1Angle(point, target) distance.  They may be slightly different
    // because query->distance() is computed using S1ChordAngle.  Note that
    // the error gets considerably larger (1e-7) as the angle approaches Pi.
    Assert.approximately(getDistance(target, result.point()).radians(),
                result.distance().radians(), kChordAngleError);
    // Check that the point satisfies the region() condition.
    if (query.options().region()) {
      Assert.equal(query.mutableOptions().region().contains(result.point()), true);
    }
    // Check that it satisfies the max_distance() condition.
    Assert.notGreaterThan(result.distance(), query.options().maxDistance());
    results ~= Result(result.distance().toS1Angle(), result.data());
  }
}

private void testFindClosestPoints(TestQuery.Target target, TestQuery query) {
  Result[] expected;
  Result[] actual;
  query.mutableOptions().setUseBruteForce(true);
  getClosestPoints(target, query, expected);
  query.mutableOptions().setUseBruteForce(false);
  getClosestPoints(target, query, actual);
  Assert.equal(
      checkDistanceResults(
          expected, actual, query.options().maxPoints(),
          query.options().maxDistance().toS1Angle(), S1Angle.zero()),
      true,
      "max_points=" ~ query.options().maxPoints().to!string
      ~ ", max_distance=" ~ query.options().maxDistance().to!string);
}

// The running time of this test is proportional to
//    num_indexes * num_points * num_queries.
private void testWithIndexFactory(PointIndexFactory factory,
                                 int num_indexes, int num_points,
                                 int num_queries) {
  auto index = new TestIndex();
  for (int i_index = 0; i_index < num_indexes; ++i_index) {
    // Generate a point set and index it.
    S2Cap query_cap = new S2Cap(S2Testing.randomPoint(), kRadius);
    index.clear();
    factory.addPoints(query_cap, num_points, index);
    for (int i_query = 0; i_query < num_queries; ++i_query) {
      // Use a new query each time to avoid resetting default parameters.
      auto options = new TestQuery.Options();
      options.setMaxPoints(1 + S2Testing.rnd.uniform(100));
      if (S2Testing.rnd.oneIn(2)) {
        options.setMaxDistance(S2Testing.rnd.randDouble() * kRadius);
      }
      S2LatLngRect rect = S2LatLngRect.fromCenterSize(
          S2LatLng(S2Testing.samplePoint(query_cap)),
          S2LatLng(S2Testing.rnd.randDouble() * kRadius,
                   S2Testing.rnd.randDouble() * kRadius));
      if (S2Testing.rnd.oneIn(5)) {
        options.setRegion(rect);
      }

      auto query = new TestQuery(index, options);
      if (S2Testing.rnd.oneIn(2)) {
        auto target = new S2ClosestPointQueryPointTarget(S2Testing.samplePoint(query_cap));
        testFindClosestPoints(target, query);
      } else {
        S2Point a = S2Testing.samplePoint(query_cap);
        S2Point b = S2Testing.samplePoint(
            new S2Cap(a, pow(1e-4, S2Testing.rnd.randDouble()) * kRadius));
        auto target = new S2ClosestPointQueryEdgeTarget(a, b);
        testFindClosestPoints(target, query);
      }
    }
  }
}

enum int kNumIndexes = 10;
enum int kNumVertices = 1000;
enum int kNumQueries = 50;

// TODO: Resume here.
@("S2ClosestPointQueryTest.CirclePoints") unittest {
  testWithIndexFactory(new CirclePointIndexFactory(), kNumIndexes, kNumVertices, kNumQueries);
}

@("S2ClosestPointQueryTest.FractalPoints") unittest {
  testWithIndexFactory(new FractalPointIndexFactory(), kNumIndexes, kNumVertices, kNumQueries);
}

@("S2ClosestPointQueryTest.GridPoints") unittest {
  testWithIndexFactory(new GridPointIndexFactory(), kNumIndexes, kNumVertices, kNumQueries);
}
