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

module s2.s2closest_edge_query_test;

import s2.s2closest_edge_query;

//import s2.s2polygon;
import s2.mutable_s2shape_index;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2edge_distances;
import s2.s2edge_vector_shape;
import s2.s2loop;
import s2.s2metrics;
import s2.s2point;
import s2.s2point_vector_shape;
import s2.s2predicates;
import s2.s2testing;
import s2.s2text_format;

import fluent.asserts;

import std.math;
import std.stdio;

@("S2ClosestEdgeQuery.NoEdges") unittest {
  auto index = new MutableS2ShapeIndex();
  auto query = new S2ClosestEdgeQuery(index);
  auto target = new S2ClosestEdgeQuery.PointTarget(S2Point(1, 0, 0));
  const auto edge = query.findClosestEdge(target);
  Assert.equal(edge.distance.s1ChordAngle, S1ChordAngle.infinity());
  Assert.equal(edge.edgeId, -1);
  Assert.equal(edge.shapeId, -1);
  Assert.equal(query.getDistance(target), S1ChordAngle.infinity());
}

@("S2ClosestEdgeQuery.OptionsNotModified") unittest {
  // Tests that FindClosestEdge(), GetDistance(), and IsDistanceLess() do not
  // modify query.options(), even though all of these methods have their own
  // specific options requirements.
  auto options = new S2ClosestEdgeQuery.Options();
  options.setMaxEdges(3);
  options.setMaxDistance(S1ChordAngle.fromDegrees(3));
  options.setMaxError(S1ChordAngle.fromDegrees(0.001).toS1Angle());
  auto index = makeIndexOrDie("1:1 | 1:2 | 1:3 # #");
  auto query = new S2ClosestEdgeQuery(index, options);
  auto target = new S2ClosestEdgeQuery.PointTarget(makePointOrDie("2:2"));
  Assert.equal(query.findClosestEdge(target).edgeId, 1);
  Assert.approximately(query.getDistance(target).degrees(), 1e-15, 1.0);
  Assert.equal(query.isDistanceLess(target, S1ChordAngle.fromDegrees(1.5)), true);

  // Verify that none of the options above were modified.
  Assert.equal(options.maxEdges(), query.options().maxEdges());
  Assert.equal(options.maxDistance(), query.options().maxDistance());
  Assert.equal(options.maxError(), query.options().maxError());
}

@("S2ClosestEdgeQuery.DistanceEqualToLimit") unittest {
  // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
  // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
  // the distance to the target exactly equals the chosen limit.
  S2Point p0 = makePointOrDie("23:12");
  S2Point p1 = makePointOrDie("47:11");
  S2Point[] index_points = [p0];
  auto index = new MutableS2ShapeIndex();
  index.add(new S2PointVectorShape(index_points));
  auto query = new S2ClosestEdgeQuery(index);

  // Start with two identical points and a zero distance.
  auto target0 = new S2ClosestEdgeQuery.PointTarget(p0);
  S1ChordAngle dist0 = S1ChordAngle.zero();
  Assert.equal(query.isDistanceLess(target0, dist0), false);
  Assert.equal(query.isDistanceLessOrEqual(target0, dist0), true);
  Assert.equal(query.isConservativeDistanceLessOrEqual(target0, dist0), true);

  // Now try two points separated by a non-zero distance.
  auto target1 = new S2ClosestEdgeQuery.PointTarget(p1);
  auto dist1 = S1ChordAngle(p0, p1);
  Assert.equal(query.isDistanceLess(target1, dist1), false);
  Assert.equal(query.isDistanceLessOrEqual(target1, dist1), true);
  Assert.equal(query.isConservativeDistanceLessOrEqual(target1, dist1), true);
}

@("S2ClosestEdgeQuery.TrueDistanceLessThanS1ChordAngleDistance") unittest {
  // Tests that IsConservativeDistanceLessOrEqual returns points where the
  // true distance is slightly less than the one computed by S1ChordAngle.
  //
  // The points below had the worst error from among 100,000 random pairs.
  auto p0 = S2Point(0.78516762584829192, -0.50200400690845970, -0.36263449417782678);
  auto p1 = S2Point(0.78563011732429433, -0.50187655940493503, -0.36180828883938054);

  // The S1ChordAngle distance is ~4 ulps greater than the true distance.
  auto dist1 = S1ChordAngle(p0, p1);
  auto limit = dist1.predecessor().predecessor().predecessor().predecessor();
  Assert.lessThan(compareDistance(p0, p1, limit), 0);

  // Verify that IsConservativeDistanceLessOrEqual() still returns "p1".
  S2Point[] index_points = [p0];
  auto index = new MutableS2ShapeIndex();
  index.add(new S2PointVectorShape(index_points));
  auto query = new S2ClosestEdgeQuery(index);
  auto target1 = new S2ClosestEdgeQuery.PointTarget(p1);
  Assert.equal(query.isDistanceLess(target1, limit), false);
  Assert.equal(query.isDistanceLessOrEqual(target1, limit), false);
  Assert.equal(query.isConservativeDistanceLessOrEqual(target1, limit), true);
}

@("S2ClosestEdgeQuery.TestReuseOfQuery") unittest {
  // Tests that between queries, the internal mechanism for de-duplicating
  // results is re-set. See b/71646017.
  auto index = makeIndexOrDie("2:2 # #");
  auto query = new S2ClosestEdgeQuery(index);
  query.mutableOptions().setMaxError(S1Angle.fromDegrees(1));
  auto target_index = makeIndexOrDie("## 0:0, 0:5, 5:5, 5:0");
  auto target = new S2ClosestEdgeQuery.ShapeIndexTarget(target_index);
  auto results1 = query.findClosestEdges(target);
  auto results2 = query.findClosestEdges(target);
  Assert.equal(results1.length, results2.length);
}

@("S2ClosestEdgeQuery.TargetPointInsideIndexedPolygon") unittest {
  // Tests a target point in the interior of an indexed polygon.
  // (The index also includes a polyline loop with no interior.)
  auto index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  auto options = new S2ClosestEdgeQuery.Options();
  options.setIncludeInteriors(true);
  options.setMaxDistance(S1Angle.fromDegrees(1));
  auto query = new S2ClosestEdgeQuery(index, options);
  auto target = new S2ClosestEdgeQuery.PointTarget(makePointOrDie("2:12"));
  auto results = query.findClosestEdges(target);
  Assert.equal(results.length, 1);
  Assert.equal(results[0].distance.s1ChordAngle, S1ChordAngle.zero());
  Assert.equal(results[0].shapeId, 1);
  Assert.equal(results[0].edgeId, -1);
}

@("S2ClosestEdgeQuery.TargetPointOutsideIndexedPolygon") unittest {
  // Tests a target point in the interior of a polyline loop with no
  // interior.  (The index also includes a nearby polygon.)
  auto index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  auto options = new S2ClosestEdgeQuery.Options();
  options.setIncludeInteriors(true);
  options.setMaxDistance(S1Angle.fromDegrees(1));
  auto query = new S2ClosestEdgeQuery(index, options);
  auto target = new S2ClosestEdgeQuery.PointTarget(makePointOrDie("2:2"));
  auto results = query.findClosestEdges(target);
  Assert.equal(results.length, 0);
}

@("S2ClosestEdgeQuery.TargetPolygonContainingIndexedPoints") unittest {
  // Two points are contained within a polyline loop (no interior) and two
  // points are contained within a polygon.
  auto index = makeIndexOrDie("2:2 | 3:3 | 1:11 | 3:13 # #");
  auto query = new S2ClosestEdgeQuery(index);
  query.mutableOptions().setMaxDistance(S1Angle.fromDegrees(1));
  auto target_index = makeIndexOrDie(
      "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  auto target = new S2ClosestEdgeQuery.ShapeIndexTarget(target_index);
  target.setIncludeInteriors(true);
  // TODO: Resume here.
  auto results = query.findClosestEdges(target);
  Assert.equal(results.length, 2);
  Assert.equal(results[0].distance.s1ChordAngle, S1ChordAngle.zero());
  Assert.equal(results[0].shapeId, 0);
  Assert.equal(results[0].edgeId, 2);  // 1:11
  Assert.equal(results[1].distance.s1ChordAngle, S1ChordAngle.zero());
  Assert.equal(results[1].shapeId, 0);
  Assert.equal(results[1].edgeId, 3);  // 3:13
}

@("S2ClosestEdgeQuery.EmptyPolygonTarget") unittest {
  // Verifies that distances are measured correctly to empty polygon targets.
  auto empty_polygon_index = makeIndexOrDie("# # empty");
  auto point_index = makeIndexOrDie("1:1 # #");
  auto full_polygon_index = makeIndexOrDie("# # full");
  auto target = new S2ClosestEdgeQuery.ShapeIndexTarget(empty_polygon_index);
  target.setIncludeInteriors(true);

  auto empty_query = new S2ClosestEdgeQuery(empty_polygon_index);
  empty_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(empty_query.getDistance(target), S1ChordAngle.infinity());

  auto point_query = new S2ClosestEdgeQuery(point_index);
  point_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(point_query.getDistance(target), S1ChordAngle.infinity());

  auto full_query = new S2ClosestEdgeQuery(full_polygon_index);
  full_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(full_query.getDistance(target), S1ChordAngle.infinity());
}

@("S2ClosestEdgeQuery.FullLaxPolygonTarget") unittest {
  // Verifies that distances are measured correctly to full LaxPolygon targets.
  auto empty_polygon_index = makeIndexOrDie("# # empty");
  auto point_index = makeIndexOrDie("1:1 # #");
  auto full_polygon_index = makeIndexOrDie("# # full");
  auto target = new S2ClosestEdgeQuery.ShapeIndexTarget(full_polygon_index);
  target.setIncludeInteriors(true);

  auto empty_query = new S2ClosestEdgeQuery(empty_polygon_index);
  empty_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(empty_query.getDistance(target), S1ChordAngle.infinity());

  auto point_query = new S2ClosestEdgeQuery(point_index);
  point_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(point_query.getDistance(target), S1ChordAngle.zero());

  auto full_query = new S2ClosestEdgeQuery(full_polygon_index);
  full_query.mutableOptions().setIncludeInteriors(true);
  Assert.equal(full_query.getDistance(target), S1ChordAngle.zero());
}

/+ TODO: Add when S2Polygon is added.
@("S2ClosestEdgeQuery.FullS2PolygonTarget") unittest {
  // Verifies that distances are measured correctly to full S2Polygon targets
  // (which use a different representation of "full" than LaxPolygon does).
  auto empty_polygon_index = makeIndexOrDie("# # empty");
  auto point_index = makeIndexOrDie("1:1 # #");
  auto full_polygon_index = makeIndexOrDie("# #");
  full_polygon_index.add(new S2Polygon.OwningShape(
      s2textformat::MakePolygonOrDie("full")));

  S2ClosestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);

  S2ClosestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), empty_query.GetDistance(&target));

  S2ClosestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), point_query.GetDistance(&target));

  S2ClosestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), full_query.GetDistance(&target));
}
+/

@("S2ClosestEdgeQuery.IsConservativeDistanceLessOrEqual") unittest {
  // Test
  int num_tested = 0;
  int num_conservative_needed = 0;
  auto rnd = S2Testing.rnd;
  foreach (int iter; 0..1000) {
    rnd.reset(iter + 1);  // Easier to reproduce a specific case.
    S2Point x = S2Testing.randomPoint();
    S2Point dir = S2Testing.randomPoint();
    S1Angle r = S1Angle.fromRadians(M_PI * pow(1e-30, rnd.randDouble()));
    S2Point y = interpolateAtDistance(r, x, dir);
    auto limit = S1ChordAngle(r);
    if (compareDistance(x, y, limit) <= 0) {
      auto index = new MutableS2ShapeIndex();
      index.add(new S2PointVectorShape([x]));
      auto query = new S2ClosestEdgeQuery(index);
      auto target = new S2ClosestEdgeQuery.PointTarget(y);
      Assert.equal(query.isConservativeDistanceLessOrEqual(target, limit), true);
      ++num_tested;
      if (!query.isDistanceLess(target, limit)) ++num_conservative_needed;
    }
  }
  // Verify that in most test cases, the distance between the target points
  // was close to the desired value.  Also verify that at least in some test
  // cases, the conservative distance test was actually necessary.
  Assert.notLessThan(num_tested, 300);
  Assert.notGreaterThan(num_tested, 700);
  Assert.notLessThan(num_conservative_needed, 25);
}

// An abstract class that adds edges to a MutableS2ShapeIndex for benchmarking.
abstract class ShapeIndexFactory {
public:
  // Requests that approximately "num_edges" edges located within the given
  // S2Cap bound should be added to "index".
  void addEdges(in S2Cap index_cap, int num_edges, MutableS2ShapeIndex index) const;
}

// Generates a regular loop that approximately fills the given S2Cap.
//
// Regular loops are nearly the worst case for distance calculations, since
// many edges are nearly equidistant from any query point that is not
// immediately adjacent to the loop.
class RegularLoopShapeIndexFactory : ShapeIndexFactory {
public:
  override
  void addEdges(in S2Cap index_cap, int num_edges, MutableS2ShapeIndex index) const {
    index.add(new S2Loop.Shape(
        S2Loop.makeRegularLoop(index_cap.center(), index_cap.getRadius(), num_edges)));
  }
}

// TODO: Resume here.
// Generates a fractal loop that approximately fills the given S2Cap.
// class FractalLoopShapeIndexFactory : ShapeIndexFactory {
// public:
//   override
//   void addEdges(in S2Cap index_cap, int num_edges, MutableS2ShapeIndex index) const {
//     S2Testing.Fractal fractal;
//     fractal.SetLevelForApproxMaxEdges(num_edges);
//     index->Add(make_unique<S2Loop::OwningShape>(
//         fractal.MakeLoop(S2Testing::GetRandomFrameAt(index_cap.center()),
//                          index_cap.GetRadius())));
//   }
// };

/+

// Generates a cloud of points that approximately fills the given S2Cap.
class PointCloudShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(const S2Cap& index_cap, int num_edges,
                MutableS2ShapeIndex* index) const override {
    vector<S2Point> points;
    for (int i = 0; i < num_edges; ++i) {
      points.push_back(S2Testing::SamplePoint(index_cap));
    }
    index->Add(make_unique<S2PointVectorShape>(std::move(points)));
  }
};

// The approximate radius of S2Cap from which query edges are chosen.
static const S1Angle kRadius = S2Testing::KmToAngle(10);

// An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
static double kChordAngleError = 1e-15;

using Result = pair<S1Angle, s2shapeutil::ShapeEdgeId>;

// Converts to the format required by CheckDistanceResults() in s2testing.h
vector<Result> ConvertResults(const vector<S2ClosestEdgeQuery::Result>& edges) {
  vector<Result> results;
  for (const auto& edge : edges) {
    results.push_back(
        make_pair(edge.distance.ToAngle(),
                  s2shapeutil::ShapeEdgeId(edge.shape_id, edge.edge_id)));
  }
  return results;
}

// Use "query" to find the closest edge(s) to the given target, then convert
// the query results into two parallel vectors, one for distances and one for
// (shape_id, edge_id) pairs.  Also verify that the results satisfy the search
// criteria.
static void GetClosestEdges(S2ClosestEdgeQuery::Target* target,
                            S2ClosestEdgeQuery *query,
                            vector<S2ClosestEdgeQuery::Result>* edges) {
  query->FindClosestEdges(target, edges);
  EXPECT_LE(edges->size(), query->options().max_edges());
  if (query->options().max_distance() ==
      S2ClosestEdgeQuery::Distance::Infinity()) {
    int min_expected = min(query->options().max_edges(),
                           s2shapeutil::CountEdges(query->index()));
    if (!query->options().include_interiors()) {
      // We can predict exactly how many edges should be returned.
      EXPECT_EQ(min_expected, edges->size());
    } else {
      // All edges should be returned, and possibly some shape interiors.
      EXPECT_LE(min_expected, edges->size());
    }
  }
  for (const auto& edge : *edges) {
    // Check that the edge satisfies the max_distance() condition.
    EXPECT_LE(edge.distance, query->options().max_distance());
  }
}

static S2ClosestEdgeQuery::Result TestFindClosestEdges(
    S2ClosestEdgeQuery::Target* target, S2ClosestEdgeQuery *query) {
  vector<S2ClosestEdgeQuery::Result> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestEdges(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetClosestEdges(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(ConvertResults(expected),
                                   ConvertResults(actual),
                                   query->options().max_edges(),
                                   query->options().max_distance().ToAngle(),
                                   query->options().max_error().ToAngle()))
      << "max_edges=" << query->options().max_edges()
      << ", max_distance=" << query->options().max_distance()
      << ", max_error=" << query->options().max_error();

  if (expected.empty()) return S2ClosestEdgeQuery::Result();

  // Note that when options.max_error() > 0, expected[0].distance may not be
  // the minimum distance.  It is never larger by more than max_error(), but
  // the actual value also depends on max_edges().
  //
  // Here we verify that GetDistance() and IsDistanceLess() return results
  // that are consistent with the max_error() setting.
  S1ChordAngle max_error = query->options().max_error();
  S1ChordAngle min_distance = expected[0].distance;
  EXPECT_LE(query->GetDistance(target), min_distance + max_error);

  // Test IsDistanceLess().
  EXPECT_FALSE(query->IsDistanceLess(target, min_distance - max_error));
  EXPECT_TRUE(query->IsDistanceLess(target, min_distance.Successor()));

  // Return the closest edge result so that we can also test Project.
  return expected[0];
}

// The running time of this test is proportional to
//    (num_indexes + num_queries) * num_edges.
// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(const ShapeIndexFactory& factory,
                                 int num_indexes, int num_edges,
                                 int num_queries) {
  // Build a set of MutableS2ShapeIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<MutableS2ShapeIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    S2Testing::rnd.Reset(FLAGS_s2_random_seed + i);
    index_caps.push_back(S2Cap(S2Testing::RandomPoint(), kRadius));
    indexes.emplace_back(new MutableS2ShapeIndex);
    factory.AddEdges(index_caps.back(), num_edges, indexes.back().get());
  }
  for (int i = 0; i < num_queries; ++i) {
    S2Testing::rnd.Reset(FLAGS_s2_random_seed + i);
    int i_index = S2Testing::rnd.Uniform(num_indexes);
    const S2Cap& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2Cap query_cap(index_cap.center(), query_radius);
    S2ClosestEdgeQuery query(indexes[i_index].get());

    // Occasionally we don't set any limit on the number of result edges.
    // (This may return all edges if we also don't set a distance limit.)
    if (!S2Testing::rnd.OneIn(5)) {
      query.mutable_options()->set_max_edges(1 + S2Testing::rnd.Uniform(10));
    }
    // We set a distance limit 2/3 of the time.
    if (!S2Testing::rnd.OneIn(3)) {
      query.mutable_options()->set_max_distance(
          S2Testing::rnd.RandDouble() * query_radius);
    }
    if (S2Testing::rnd.OneIn(2)) {
      // Choose a maximum error whose logarithm is uniformly distributed over
      // a reasonable range, except that it is sometimes zero.
      query.mutable_options()->set_max_error(S1Angle::Radians(
          pow(1e-4, S2Testing::rnd.RandDouble()) * query_radius.radians()));
    }
    query.mutable_options()->set_include_interiors(S2Testing::rnd.OneIn(2));
    int target_type = S2Testing::rnd.Uniform(4);
    if (target_type == 0) {
      // Find the edges closest to a given point.
      S2Point point = S2Testing::SamplePoint(query_cap);
      S2ClosestEdgeQuery::PointTarget target(point);
      auto closest = TestFindClosestEdges(&target, &query);
      if (!closest.distance.is_infinity()) {
        // Also test the Project method.
        EXPECT_NEAR(
            closest.distance.ToAngle().radians(),
            S1Angle(point, query.Project(point, closest)).radians(),
            kChordAngleError);
      }
    } else if (target_type == 1) {
      // Find the edges closest to a given edge.
      S2Point a = S2Testing::SamplePoint(query_cap);
      S2Point b = S2Testing::SamplePoint(
          S2Cap(a, pow(1e-4, S2Testing::rnd.RandDouble()) * query_radius));
      S2ClosestEdgeQuery::EdgeTarget target(a, b);
      TestFindClosestEdges(&target, &query);
    } else if (target_type == 2) {
      // Find the edges closest to a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = min_level + S2Testing::rnd.Uniform(
          S2CellId::kMaxLevel - min_level + 1);
      S2Point a = S2Testing::SamplePoint(query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      S2ClosestEdgeQuery::CellTarget target(cell);
      TestFindClosestEdges(&target, &query);
    } else {
      DCHECK_EQ(3, target_type);
      // Use another one of the pre-built indexes as the target.
      int j_index = S2Testing::rnd.Uniform(num_indexes);
      S2ClosestEdgeQuery::ShapeIndexTarget target(indexes[j_index].get());
      target.set_include_interiors(S2Testing::rnd.OneIn(2));
      TestFindClosestEdges(&target, &query);
    }
  }
}

static const int kNumIndexes = 50;
static const int kNumEdges = 100;
static const int kNumQueries = 200;

TEST(S2ClosestEdgeQuery, CircleEdges) {
  TestWithIndexFactory(RegularLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, FractalEdges) {
  TestWithIndexFactory(FractalLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, PointCloudEdges) {
  TestWithIndexFactory(PointCloudShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, ConservativeCellDistanceIsUsed) {
  // These specific test cases happen to fail if max_error() is not properly
  // taken into account when measuring distances to S2ShapeIndex cells.
  for (int seed : {42, 681, 894, 1018, 1750, 1759, 2401}) {
    FLAGS_s2_random_seed = seed;  // Automatically restored.
    TestWithIndexFactory(FractalLoopShapeIndexFactory(), 5, 100, 10);
  }
}

+/
