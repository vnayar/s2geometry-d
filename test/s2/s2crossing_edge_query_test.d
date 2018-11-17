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

module s2.s2crossing_edge_query_test;

import s2.s2crossing_edge_query;

import s2.logger;
import s2.mutable_s2shape_index;
import s2.r2point;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2coords : FaceUVtoXYZ, GetNorm, GetUVWAxis, MAX_CELL_LEVEL;
import s2.s2crossing_edge_query : CrossingType;
import s2.s2edge_crossings : crossingSign;
import s2.s2edge_distances : getDistance, interpolateAtDistance;
import s2.s2edge_vector_shape;
import s2.s2metrics : MAX_DIAG;
import s2.s2point;
import s2.s2polyline;
import s2.s2shape;
import s2.s2testing;
import s2.s2text_format;

import std.algorithm : count, isSorted;
import std.conv : to;
import std.range : assumeSorted, back, empty, front;
import std.math : nextafter, pow;

import fluent.asserts;

struct TestEdge {
  S2Point v0, v1;
}

S2Point perturbAtDistance(S1Angle distance, in S2Point a0, in S2Point b0) {
  S2Point x = interpolateAtDistance(distance, a0, b0);
  if (S2Testing.rnd.oneIn(2)) {
    foreach (int i; 0 .. 3) {
      x[i] = nextafter(x[i], S2Testing.rnd.oneIn(2) ? 1 : -1);
    }
    x = x.normalize();
  }
  return x;
}

// Generate sub-edges of some given edge (a0,b0).  The length of the sub-edges
// is distributed exponentially over a large range, and the endpoints may be
// slightly perturbed to one side of (a0,b0) or the other.
void getPerturbedSubEdges(S2Point a0, S2Point b0, int count, TestEdge[] edges) {
  edges.length = 0;
  a0 = a0.normalize();
  b0 = b0.normalize();
  auto length0 = S1Angle(a0, b0);
  for (int i = 0; i < count; ++i) {
    S1Angle length = length0 * pow(1e-15, S2Testing.rnd.randDouble());
    S1Angle offset = (length0 - length) * S2Testing.rnd.randDouble();
    edges ~= TestEdge(
        perturbAtDistance(offset, a0, b0), perturbAtDistance(offset + length, a0, b0));
  }
}

// Generate edges whose center is randomly chosen from the given S2Cap, and
// whose length is randomly chosen up to "max_length".
void getCapEdges(in S2Cap center_cap, S1Angle max_length, int count, TestEdge[] edges) {
  edges.length = 0;
  for (int i = 0; i < count; ++i) {
    S2Point center = S2Testing.samplePoint(center_cap);
    auto edge_cap = new S2Cap(center, 0.5 * max_length);
    S2Point p1 = S2Testing.samplePoint(edge_cap);
    // Compute p1 reflected through "center", and normalize for good measure.
    S2Point p2 = (2 * p1.dotProd(center) * center - p1).normalize();
    edges ~= TestEdge(p1, p2);
  }
}

void checkAllCrossings(in TestEdge[] edges) {
  auto shape = new S2EdgeVectorShape();  // raw pointer since "shape" used below
  foreach (const(TestEdge) edge; edges) {
    shape.add(edge.v0, edge.v1);
  }
  // Force more subdivision than usual to make the test more challenging.
  MutableS2ShapeIndex.Options options;
  options.maxEdgesPerCell(1);
  auto index = new MutableS2ShapeIndex(options);
  index.add(shape);
  // To check that candidates are being filtered reasonably, we count the
  // total number of candidates that the total number of edge pairs that
  // either intersect or are very close to intersecting.
  int num_candidates = 0, num_nearby_pairs = 0;
  int i = 0;
  foreach (const TestEdge edge; edges) {
    logger.logTrace("Iteration ", i++);
    const(S2Point) a = edge.v0;
    const(S2Point) b = edge.v1;
    int[] candidates;
    auto query = new S2CrossingEdgeQuery(index);
    query.getCandidates(a, b, shape, candidates);

    // Verify that the second version of GetCandidates returns the same result.
    S2CrossingEdgeQuery.EdgeMap edge_map;
    query.getCandidates(a, b, edge_map);
    Assert.equal(edge_map.length, 1);
    Assert.equal(edge_map[shape], candidates);
    Assert.equal(candidates.empty(), false);

    // Now check the actual candidates.
    Assert.equal(isSorted(candidates), true);
    Assert.notLessThan(candidates.front(), 0);
    Assert.lessThan(candidates.back(), shape.numEdges());
    num_candidates += candidates.length;
    string missing_candidates;
    int[] expected_crossings;
    int[] expected_interior_crossings;
    for (int j = 0; j < shape.numEdges(); ++j) {
      auto edge2 = shape.edge(j);
      const S2Point c = edge2.v0;
      const S2Point d = edge2.v1;
      int sign = crossingSign(a, b, c, d);
      if (sign >= 0) {
        expected_crossings ~= j;
        if (sign > 0) {
          expected_interior_crossings ~= j;
        }
        ++num_nearby_pairs;
        if (assumeSorted(candidates).equalRange(j).empty) {
          missing_candidates ~= " " ~ to!string(j);
        }
      } else {
        const double kMaxDist = MAX_DIAG.getValue(MAX_CELL_LEVEL);
        if (getDistance(a, c, d).radians() < kMaxDist ||
            getDistance(b, c, d).radians() < kMaxDist ||
            getDistance(c, a, b).radians() < kMaxDist ||
            getDistance(d, a, b).radians() < kMaxDist) {
          ++num_nearby_pairs;
        }
      }
    }
    Assert.equal(missing_candidates.empty(), true, missing_candidates);

    // Test that GetCrossings() returns only the actual crossing edges.
    int[] actual_crossings;
    query.getCrossings(a, b, shape, CrossingType.ALL, actual_crossings);
    Assert.equal(expected_crossings, actual_crossings);

    // Verify that the second version of GetCrossings returns the same result.
    query.getCrossings(a, b, CrossingType.ALL, edge_map);
    Assert.equal(edge_map.length, 1);
    Assert.equal(edge_map[shape], expected_crossings);

    // Verify that CrossingType::INTERIOR returns only the interior crossings.
    query.getCrossings(a, b, shape, CrossingType.INTERIOR, actual_crossings);
    Assert.equal(expected_interior_crossings, actual_crossings);
  }
  // There is nothing magical about this particular ratio; this check exists
  // to catch changes that dramatically increase the number of candidates.
  Assert.notGreaterThan(num_candidates, 3 * num_nearby_pairs);
}

// Test edges that lie in the plane of one of the S2 cube edges.  Such edges
// may lie on the boundary between two cube faces, or pass through a cube
// vertex, or follow a 45 diagonal across a cube face toward its center.
//
// This test is sufficient to demonstrate that padding the cell boundaries is
// necessary for correctness.  (It fails if MutableS2ShapeIndex::kCellPadding
// is set to zero.)
@("GetCrossingCandidates.PerturbedCubeEdges") unittest {
  auto rnd = S2Testing.rnd;
  TestEdge[] edges;
  for (int iter = 0; iter < 10; ++iter) {
    int face = rnd.uniform(6);
    double scale = pow(1e-15, rnd.randDouble());
    auto uv = R2Point(2 * rnd.uniform(2) - 1, 2 * rnd.uniform(2) - 1);  // vertex
    S2Point a0 = FaceUVtoXYZ(face, scale * uv);
    S2Point b0 = a0 - 2 * GetNorm(face);
    // TODO(ericv): This test is currently slow because *every* crossing test
    // needs to invoke s2pred::ExpensiveSign().
    getPerturbedSubEdges(a0, b0, 30, edges);
    checkAllCrossings(edges);
  }
}

// Test edges that lie in the plane of one of the S2 cube face axes.  These
// edges are special because one coordinate is zero, and they lie on the
// boundaries between the immediate child cells of the cube face.
@("GetCrossingCandidates.PerturbedCubeFaceAxes") unittest {
  auto rnd = S2Testing.rnd;
  TestEdge[] edges;
  for (int iter = 0; iter < 5; ++iter) {
    int face = rnd.uniform(6);
    double scale = pow(1e-15, rnd.randDouble());
    S2Point axis = GetUVWAxis(face, rnd.uniform(2));
    S2Point a0 = scale * axis + GetNorm(face);
    S2Point b0 = scale * axis - GetNorm(face);
    getPerturbedSubEdges(a0, b0, 30, edges);
    checkAllCrossings(edges);
  }
}

@("GetCrossingCandidates.CapEdgesNearCubeVertex") unittest {
  // Test a random collection of edges near the S2 cube vertex where the
  // Hilbert curve starts and ends.
  TestEdge[] edges;
  getCapEdges(
      new S2Cap(S2Point(-1, -1, 1).normalize(),
      S1Angle.fromRadians(1e-3)),
      S1Angle.fromRadians(1e-4), 1000, edges);
  checkAllCrossings(edges);
}

@("GetCrossingCandidates.DegenerateEdgeOnCellVertexIsItsOwnCandidate") unittest {
  for (int i = 0; i < 100; ++i) {
    TestEdge[] edges;
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    edges ~= TestEdge(cell.getVertex(0), cell.getVertex(0));
    checkAllCrossings(edges);
  }
}

@("GetCrossingCandidates.CollinearEdgesOnCellBoundaries") unittest {
  const int kNumEdgeIntervals = 8;  // 9*8/2 = 36 edges
  for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
    auto cell = new S2Cell(S2Testing.getRandomCellId(level));
    int k = S2Testing.rnd.uniform(4);
    S2Point p1 = cell.getVertexRaw(k);
    S2Point p2 = cell.getVertexRaw(k + 1);
    S2Point delta = (p2 - p1) / kNumEdgeIntervals;
    TestEdge[] edges;
    for (int i = 0; i <= kNumEdgeIntervals; ++i) {
      for (int j = 0; j < i; ++j) {
        edges ~= TestEdge((p1 + i * delta).normalize(), (p1 + j * delta).normalize());
      }
    }
    checkAllCrossings(edges);
  }
}

// This is the example from the header file, with a few extras.
void checkPolylineCrossings(S2ShapeIndex index, in S2Point a0, in S2Point a1) {
  auto query = new S2CrossingEdgeQuery(index);
  S2CrossingEdgeQuery.EdgeMap edge_map;
  if (!query.getCrossings(a0, a1, CrossingType.ALL, edge_map)) return;
  foreach (S2Shape s, int[] edges; edge_map) {
    S2Polyline.Shape shape = cast(S2Polyline.Shape) s;
    const S2Polyline polyline = shape.polyline();
    // Shapes with no crossings should be filtered out by this method.
    Assert.equal(edges.empty(), false);
    foreach (int edge; edges) {
      const S2Point b0 = polyline.vertex(edge);
      const S2Point b1 = polyline.vertex(edge + 1);
      Assert.notLessThan(crossingSign(a0, a1, b0, b1), 0);
    }
  }
  // Also test that no edges are missing.
  for (int i = 0; i < index.numShapeIds(); ++i) {
    auto shape = cast(S2Polyline.Shape) index.shape(i);
    if (shape !in edge_map) continue;
    const int[] edges = edge_map[shape];
    const S2Polyline polyline = shape.polyline();
    for (int e = 0; e < polyline.numVertices() - 1; ++e) {
      if (crossingSign(a0, a1, polyline.vertex(e), polyline.vertex(e + 1)) >= 0) {
        Assert.equal(count(edges, e), 1);
      }
    }
  }
}

@("GetCrossings.PolylineCrossings") unittest {
  auto index = new MutableS2ShapeIndex();
  // Three zig-zag lines near the equator.
  index.add(new S2Polyline.Shape(makePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")));
  index.add(new S2Polyline.Shape(makePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")));
  index.add(new S2Polyline.Shape(makePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")));
  checkPolylineCrossings(index, makePointOrDie("1:0"), makePointOrDie("1:4"));
  checkPolylineCrossings(index, makePointOrDie("5:5"), makePointOrDie("6:6"));
}

/+ TODO: Add when S@Loop is implemented.
@("GetCrossings.ShapeIdsAreCorrect") unittest {
  // This tests that when some index cells contain only one shape, the
  // intersecting edges are returned with the correct shape id.
  auto index = new MutableS2ShapeIndex();
  index.add(S2Testing.makeRegularPoints(makePointOrDie("0:0"), S1Angle.fromDegrees(5), 100))));
  index.add(S2Testing.makeRegularPoints(makePointOrDie("0:20"), S1Angle.fromDegrees(5), 100))));
  checkPolylineCrossings(index, makePointOrDie("1:-10"), makePointOrDie("1:30"));
}
+/
