// Copyright 2005 Google Inc. All Rights Reserved.
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
// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2cell_test;

import s2.logger;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2coords;
import s2.s2edge_crossings : crossingSign;
import s2.s2edge_distances : interpolate, updateMinDistance, updateMaxDistance;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2latlng_rect_bounder;
import s2.s2loop;
import s2.s2metrics;
import s2.s2point;
import s2.s2pointutil;
import s2.s2testing;
import s2.s2text_format;
import s2.util.coding.coder;
import s2.util.math.s2const;

import algorithm = std.algorithm;
import fluent.asserts;
import math = std.math;
import std.exception;
import std.stdio;

@("S2Cell.TestFaces") unittest {
  int[S2Point] edge_counts;
  int[S2Point] vertex_counts;
  foreach (face; 0 .. 6) {
    S2CellId id = S2CellId.fromFace(face);
    S2Cell cell = new S2Cell(id);
    Assert.equal(id, cell.id());
    Assert.equal(face, cell.face());
    Assert.equal(0, cell.level());
    // Top-level faces have alternating orientations to get RHS coordinates.
    Assert.equal(face & SWAP_MASK, cell.orientation());
    Assert.equal(cell.isLeaf(), false);
    foreach (k; 0 .. 4) {
      edge_counts[cell.getEdgeRaw(k)] += 1;
      vertex_counts[cell.getVertexRaw(k)] += 1;
      Assert.approximately(cell.getVertexRaw(k).dotProd(cell.getEdgeRaw(k)), 0.0, DOUBLE_ERR);
      Assert.approximately(cell.getVertexRaw(k + 1).dotProd(cell.getEdgeRaw(k)), 0.0, DOUBLE_ERR);
      Assert.approximately(
          cell.getVertexRaw(k).crossProd(cell.getVertexRaw(k + 1)).normalize()
              .dotProd(cell.getEdge(k)),
          1.0, DOUBLE_ERR);
    }
  }
  // Check that edges have multiplicity 2 and vertices have multiplicity 3.
  foreach (count; edge_counts.byValue()) {
    Assert.equal(2, count);
  }
  foreach (count; vertex_counts.byValue()) {
    Assert.equal(3, count);
  }
}

struct LevelStats {
  double count = 0;
  double min_area = 100, max_area = 0, avg_area = 0;
  double min_width = 100, max_width = 0, avg_width = 0;
  double min_edge = 100, max_edge = 0, avg_edge = 0, max_edge_aspect = 0;
  double min_diag = 100, max_diag = 0, avg_diag = 0, max_diag_aspect = 0;
  double min_angle_span = 100, max_angle_span = 0, avg_angle_span = 0;
  double min_approx_ratio = 100, max_approx_ratio = 0;
}
private LevelStats[S2CellId.MAX_LEVEL + 1] level_stats;

private void gatherStats(in S2Cell cell) {
  LevelStats* s = &level_stats[cell.level()];
  double exact_area = cell.exactArea();
  double approx_area = cell.approxArea();
  double min_edge = 100, max_edge = 0, avg_edge = 0;
  double min_diag = 100, max_diag = 0;
  double min_width = 100, max_width = 0;
  double min_angle_span = 100, max_angle_span = 0;
  for (int i = 0; i < 4; ++i) {
    double edge = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 1));
    min_edge = algorithm.min(edge, min_edge);
    max_edge = algorithm.max(edge, max_edge);
    avg_edge += 0.25 * edge;
    S2Point mid = cell.getVertexRaw(i) + cell.getVertexRaw(i + 1);
    double width = M_PI_2 - mid.angle(cell.getEdgeRaw(i + 2));
    min_width = algorithm.min(width, min_width);
    max_width = algorithm.max(width, max_width);
    if (i < 2) {
      double diag = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 2));
      min_diag = algorithm.min(diag, min_diag);
      max_diag = algorithm.max(diag, max_diag);
      double angle_span = cell.getEdgeRaw(i).angle(-cell.getEdgeRaw(i + 2));
      min_angle_span = algorithm.min(angle_span, min_angle_span);
      max_angle_span = algorithm.max(angle_span, max_angle_span);
    }
  }
  s.count += 1;
  s.min_area = algorithm.min(exact_area, s.min_area);
  s.max_area = algorithm.max(exact_area, s.max_area);
  s.avg_area += exact_area;
  s.min_width = algorithm.min(min_width, s.min_width);
  s.max_width = algorithm.max(max_width, s.max_width);
  s.avg_width += 0.5 * (min_width + max_width);
  s.min_edge = algorithm.min(min_edge, s.min_edge);
  s.max_edge = algorithm.max(max_edge, s.max_edge);
  s.avg_edge += avg_edge;
  s.max_edge_aspect = algorithm.max(max_edge / min_edge, s.max_edge_aspect);
  s.min_diag = algorithm.min(min_diag, s.min_diag);
  s.max_diag = algorithm.max(max_diag, s.max_diag);
  s.avg_diag += 0.5 * (min_diag + max_diag);
  s.max_diag_aspect = algorithm.max(max_diag / min_diag, s.max_diag_aspect);
  s.min_angle_span = algorithm.min(min_angle_span, s.min_angle_span);
  s.max_angle_span = algorithm.max(max_angle_span, s.max_angle_span);
  s.avg_angle_span += 0.5 * (min_angle_span + max_angle_span);
  double approx_ratio = approx_area / exact_area;
  s.min_approx_ratio = algorithm.min(approx_ratio, s.min_approx_ratio);
  s.max_approx_ratio = algorithm.max(approx_ratio, s.max_approx_ratio);
}

private void testSubdivide(in S2Cell cell) {
  gatherStats(cell);
  if (cell.isLeaf()) return;

  S2Cell[4] children;
  enforce(cell.subdivide(children));
  S2CellId child_id = cell.id().childBegin();
  double exact_area = 0;
  double approx_area = 0;
  double average_area = 0;
  for (int i = 0; i < 4; ++i, child_id = child_id.next()) {
    exact_area += children[i].exactArea();
    approx_area += children[i].approxArea();
    average_area += children[i].averageArea();

    // Check that the child geometry is consistent with its cell ID.
    Assert.equal(child_id, children[i].id());
    Assert.equal(approxEquals(children[i].getCenter(), child_id.toS2Point()), true);
    S2Cell direct = new S2Cell(child_id);
    Assert.equal(direct.face(), children[i].face());
    Assert.equal(direct.level(), children[i].level());
    Assert.equal(direct.orientation(), children[i].orientation());
    Assert.equal(direct.getCenterRaw(), children[i].getCenterRaw());
    for (int k = 0; k < 4; ++k) {
      Assert.equal(direct.getVertexRaw(k), children[i].getVertexRaw(k));
      Assert.equal(direct.getEdgeRaw(k), children[i].getEdgeRaw(k));
    }

    // Test Contains() and MayIntersect().
    Assert.equal(cell.contains(children[i]), true);
    Assert.equal(cell.mayIntersect(children[i]), true);
    Assert.equal(children[i].contains(cell), false);
    Assert.equal(cell.contains(children[i].getCenterRaw()), true);
    for (int j = 0; j < 4; ++j) {
      Assert.equal(cell.contains(children[i].getVertexRaw(j)), true);
      if (j != i) {
        Assert.equal(children[i].contains(children[j].getCenterRaw()), false);
        Assert.equal(children[i].mayIntersect(children[j]), false);
      }
    }

    // Test GetCapBound and GetRectBound.
    S2Cap parent_cap = cell.getCapBound();
    S2LatLngRect parent_rect = cell.getRectBound();
    if (cell.contains(S2Point(0, 0, 1)) || cell.contains(S2Point(0, 0, -1))) {
      Assert.equal(parent_rect.lng().isFull(), true);
    }
    S2Cap child_cap = children[i].getCapBound();
    S2LatLngRect child_rect = children[i].getRectBound();
    Assert.equal(child_cap.contains(children[i].getCenter()), true);
    Assert.equal(child_rect.contains(children[i].getCenterRaw()), true);
    Assert.equal(parent_cap.contains(children[i].getCenter()), true);
    Assert.equal(parent_rect.contains(children[i].getCenterRaw()), true);
    for (int j = 0; j < 4; ++j) {
      Assert.equal(child_cap.contains(children[i].getVertex(j)), true);
      Assert.equal(child_rect.contains(children[i].getVertex(j)), true);
      Assert.equal(child_rect.contains(children[i].getVertexRaw(j)), true);
      Assert.equal(parent_cap.contains(children[i].getVertex(j)), true);
      Assert.equal(parent_rect.contains(children[i].getVertex(j)), true);
      Assert.equal(parent_rect.contains(children[i].getVertexRaw(j)), true);
      if (j != i) {
        // The bounding caps and rectangles should be tight enough so that
        // they exclude at least two vertices of each adjacent cell.
        int cap_count = 0;
        int rect_count = 0;
        for (int k = 0; k < 4; ++k) {
          if (child_cap.contains(children[j].getVertex(k)))
            ++cap_count;
          if (child_rect.contains(children[j].getVertexRaw(k)))
            ++rect_count;
        }
        Assert.notGreaterThan(cap_count, 2);
        if (child_rect.latLo().radians() > -M_PI_2 &&
            child_rect.latHi().radians() < M_PI_2) {
          // Bounding rectangles may be too large at the poles because the
          // pole itself has an arbitrary fixed longitude.
          Assert.notGreaterThan(rect_count, 2);
        }
      }
    }

    // Check all children for the first few levels, and then sample randomly.
    // We also always subdivide the cells containing a few chosen points so
    // that we have a better chance of sampling the minimum and maximum metric
    // values.  kMaxSizeUV is the absolute value of the u- and v-coordinate
    // where the cell size at a given level is maximal.
    const double kMaxSizeUV = 0.3964182625366691;
    const R2Point[] special_uv = [
      R2Point(double.epsilon, double.epsilon),  // Face center
      R2Point(double.epsilon, 1),            // Edge midpoint
      R2Point(1, 1),                      // Face corner
      R2Point(kMaxSizeUV, kMaxSizeUV),    // Largest cell area
      R2Point(double.epsilon, kMaxSizeUV),   // Longest edge/diagonal
    ];
    bool force_subdivide = false;
    foreach (uv; special_uv) {
      if (children[i].getBoundUV().contains(uv))
        force_subdivide = true;
    }
    if (force_subdivide || cell.level() < 5 || S2Testing.rnd.oneIn(5)) {
      testSubdivide(children[i]);
    }
  }

  // Check sum of child areas equals parent area.
  //
  // For ExactArea(), the best relative error we can expect is about 1e-6
  // because the precision of the unit vector coordinates is only about 1e-15
  // and the edge length of a leaf cell is about 1e-9.
  //
  // For ApproxArea(), the areas are accurate to within a few percent.
  //
  // For AverageArea(), the areas themselves are not very accurate, but
  // the average area of a parent is exactly 4 times the area of a child.

  Assert.notGreaterThan(
      math.fabs(math.log(exact_area / cell.exactArea())), math.fabs(math.log(1 + 1e-6)));
  Assert.notGreaterThan(
      math.fabs(math.log(approx_area / cell.approxArea())), math.fabs(math.log(1.03)));
  Assert.notGreaterThan(
      math.fabs(math.log(average_area / cell.averageArea())), math.fabs(math.log(1 + 1e-15)));
}

static void checkMinMaxAvg(int DimV)(
    string label, int level, double count, double abs_error,
    double min_value, double max_value, double avg_value,
    in Metric!DimV min_metric,
    in Metric!DimV max_metric,
    in Metric!DimV avg_metric) {

  // All metrics are minimums, maximums, or averages of differential
  // quantities, and therefore will not be exact for cells at any finite
  // level.  The differential minimum is always a lower bound, and the maximum
  // is always an upper bound, but these minimums and maximums may not be
  // achieved for two different reasons.  First, the cells at each level are
  // sampled and we may miss the most extreme examples.  Second, the actual
  // metric for a cell is obtained by integrating the differential quantity,
  // which is not constant across the cell.  Therefore cells at low levels
  // (bigger cells) have smaller variations.
  //
  // The "tolerance" below is an attempt to model both of these effects.
  // At low levels, error is dominated by the variation of differential
  // quantities across the cells, while at high levels error is dominated by
  // the effects of random sampling.
  double tolerance = (max_metric.getValue(level) - min_metric.getValue(level)) /
      math.sqrt(algorithm.min(count, 0.5 * cast(double)(1 << level)));
  if (tolerance == 0) tolerance = abs_error;

  double min_error = min_value - min_metric.getValue(level);
  double max_error = max_metric.getValue(level) - max_value;
  double avg_error = math.fabs(avg_metric.getValue(level) - avg_value);
  writefln("%-10s (%6.0f samples, tolerance %8.3g) - min %9.4g (%9.3g : %9.3g) "
      ~ "max %9.4g (%9.3g : %9.3g), avg %9.4g (%9.3g : %9.3g)",
      label, count, tolerance,
      min_value, min_error / min_value, min_error / tolerance,
      max_value, max_error / max_value, max_error / tolerance,
      avg_value, avg_error / avg_value, avg_error / tolerance);

  Assert.notGreaterThan(min_metric.getValue(level), min_value + abs_error);
  Assert.greaterThan(min_metric.getValue(level), min_value - tolerance);
  Assert.notGreaterThan(max_metric.getValue(level), max_value + tolerance);
  Assert.greaterThan(max_metric.getValue(level), max_value - abs_error);
  Assert.approximately(avg_metric.getValue(level), avg_value, 10 * tolerance);
}

@("S2Cell.TestSubdivide") unittest {
  // Only test a sample of faces to reduce the runtime.
  testSubdivide(S2Cell.fromFace(0));
  testSubdivide(S2Cell.fromFace(3));
  testSubdivide(S2Cell.fromFace(5));

  // This table is useful in evaluating the quality of the various S2
  // projections.
  //
  // The maximum edge *ratio* is the ratio of the longest edge of any cell to
  // the shortest edge of any cell at the same level (and similarly for the
  // maximum diagonal ratio).
  //
  // The maximum edge *aspect* is the maximum ratio of the longest edge of a
  // cell to the shortest edge of that same cell (and similarly for the
  // maximum diagonal aspect).
  writeln(
      "Ratio:  (Max value for any cell) / (Min value for any cell)\n"
      ~ "Aspect: (Max value / min value) for any cell\n"
      ~ "                   Edge          Diag       Approx Area/    Avg Area/\n"
      ~ "         Area     Length        Length       Exact Area    Exact Area\n"
      ~ "Level   Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n"
      ~ "--------------------------------------------------------------------");
  for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
    LevelStats* s = &level_stats[i];
    if (s.count > 0) {
      s.avg_area /= s.count;
      s.avg_width /= s.count;
      s.avg_edge /= s.count;
      s.avg_diag /= s.count;
      s.avg_angle_span /= s.count;
    }
    writefln("%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f",
           i, s.max_area / s.min_area,
           s.max_edge / s.min_edge, s.max_edge_aspect,
           s.max_diag / s.min_diag, s.max_diag_aspect,
           s.min_approx_ratio, s.max_approx_ratio,
           S2Cell.averageArea(i) / s.max_area,
           S2Cell.averageArea(i) / s.min_area);
  }

  // Now check the validity of the S2 length and area metrics.
  for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
    const(LevelStats)* s = &level_stats[i];
    if (s.count == 0) continue;

    writefln("Level %2d - metric value (error/actual : error/tolerance)", i);

    // The various length calculations are only accurate to 1e-15 or so,
    // so we need to allow for this amount of discrepancy with the theoretical
    // minimums and maximums.  The area calculation is accurate to about 1e-15
    // times the cell width.
    checkMinMaxAvg("area", i, s.count, 1e-15 * s.min_width,
                   s.min_area, s.max_area, s.avg_area,
                   MIN_AREA, MAX_AREA, AVG_AREA);
    checkMinMaxAvg("width", i, s.count, 1e-15,
                   s.min_width, s.max_width, s.avg_width,
                   MIN_WIDTH, MAX_WIDTH, AVG_WIDTH);
    checkMinMaxAvg("edge", i, s.count, 1e-15,
                   s.min_edge, s.max_edge, s.avg_edge,
                   MIN_EDGE, MAX_EDGE, AVG_EDGE);
    checkMinMaxAvg("diagonal", i, s.count, 1e-15,
                   s.min_diag, s.max_diag, s.avg_diag,
                   MIN_DIAG, MAX_DIAG, AVG_DIAG);
    checkMinMaxAvg("angle span", i, s.count, 1e-15,
                   s.min_angle_span, s.max_angle_span, s.avg_angle_span,
                   MIN_ANGLE_SPAN, MAX_ANGLE_SPAN, AVG_ANGLE_SPAN);

    // The aspect ratio calculations are ratios of lengths and are therefore
    // less accurate at higher subdivision levels.
    Assert.notGreaterThan(s.max_edge_aspect, MAX_EDGE_ASPECT + 1e-15 * (1 << i));
    Assert.notGreaterThan(s.max_diag_aspect, MAX_DIAG_ASPECT + 1e-15 * (1 << i));
  }
}

private const int kMaxLevel = 11;

private void expandChildren1(in S2Cell cell) {
  S2Cell[4] children;
  enforce(cell.subdivide(children));
  if (children[0].level() < kMaxLevel) {
    foreach (pos; 0 .. 4) {
      expandChildren1(children[pos]);
    }
  }
}

static void expandChildren2(in S2Cell cell) {
  S2CellId id = cell.id().childBegin();
  for (int pos = 0; pos < 4; ++pos, id = id.next()) {
    S2Cell child = new S2Cell(id);
    if (child.level() < kMaxLevel) expandChildren2(child);
  }
}

/+
// TODO: Wait until proper CPU time measuring is implemented.
TEST(S2Cell, TestPerformance) {
  double subdivide_start = S2Testing::GetCpuTime();
  ExpandChildren1(S2Cell::FromFace(0));
  double subdivide_time = S2Testing::GetCpuTime() - subdivide_start;
  fprintf(stderr, "Subdivide: %.3f seconds\n", subdivide_time);

  double constructor_start = S2Testing::GetCpuTime();
  ExpandChildren2(S2Cell::FromFace(0));
  double constructor_time = S2Testing::GetCpuTime() - constructor_start;
  fprintf(stderr, "Constructor: %.3f seconds\n", constructor_time);
}
+/

@("S2Cell.CellVsLoopRectBound")
unittest {
  // This test verifies that the S2Cell and S2Loop bounds contain each other
  // to within their maximum errors.
  //
  // The S2Cell and S2Loop calculations for the latitude of a vertex can differ
  // by up to 2 * double.epsilon, therefore the S2Cell bound should never exceed
  // the S2Loop bound by more than this (the reverse is not true, because the
  // S2Loop code sometimes thinks that the maximum occurs along an edge).
  // Similarly, the longitude bounds can differ by up to 4 * double.epsilon since
  // the S2Cell bound has an error of 2 * double.epsilon and then expands by this
  // amount, while the S2Loop bound does no expansion at all.

  // Possible additional S2Cell error compared to S2Loop error:
  S2LatLng kCellError = S2LatLng.fromRadians(2 * double.epsilon, 4 * double.epsilon);
  // Possible additional S2Loop error compared to S2Cell error:
  S2LatLng kLoopError = S2LatLngRectBounder.maxErrorForTests();

  for (int iter = 0; iter < 1000; ++iter) {
    S2Cell cell = new S2Cell(S2Testing.getRandomCellId());
    S2Loop loop = new S2Loop(cell);
    S2LatLngRect cell_bound = cell.getRectBound();
    S2LatLngRect loop_bound = loop.getRectBound();
    Assert.equal(true, loop_bound.expanded(kCellError).contains(cell_bound));
    Assert.equal(true, cell_bound.expanded(kLoopError).contains(loop_bound));
  }
}

@("S2Cell.RectBoundIsLargeEnough") unittest {
  // Construct many points that are nearly on an S2Cell edge, and verify that
  // whenever the cell contains a point P then its bound contains S2LatLng(P).

  for (int iter = 0; iter < 1000; /* advanced in loop below */) {
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    int i = S2Testing.rnd.uniform(4);
    S2Point v1 = cell.getVertex(i);
    S2Point v2 = S2Testing.samplePoint(
        new S2Cap(cell.getVertex(i + 1), S1Angle.fromRadians(1e-15)));
    S2Point p = interpolate(S2Testing.rnd.randDouble(), v1, v2);
    if (new S2Loop(cell).contains(p)) {
      Assert.equal(true, cell.getRectBound().contains(S2LatLng(p)));
      ++iter;
    }
  }
}

@("S2Cell.ConsistentWithS2CellIdFromPoint") unittest {
  // Construct many points that are nearly on an S2Cell edge, and verify that
  // S2Cell(S2CellId(p)).contains(p) is always true.

  for (int iter = 0; iter < 1000; ++iter) {
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    int i = S2Testing.rnd.uniform(4);
    S2Point v1 = cell.getVertex(i);
    S2Point v2 = S2Testing.samplePoint(
        new S2Cap(cell.getVertex(i + 1), S1Angle.fromRadians(1e-15)));
    S2Point p = interpolate(S2Testing.rnd.randDouble(), v1, v2);
    Assert.equal(new S2Cell(S2CellId(p)).contains(p), true);
  }
}

@("S2CellId.AmbiguouscontainsPoint") unittest {
  // This tests a case where S2CellId returns the "wrong" cell for a point
  // that is very close to the cell edge. (ConsistentWithS2CellIdFromPoint
  // generates more examples like this.)
  //
  // The S2Point below should have x = 0, but conversion from latlng to
  // (x,y,z) gives x = 6.1e-17.  When xyz is converted to uv, this gives u =
  // -6.1e-17.  However when converting to st, which is centered at 0.5 rather
  // than 0, the low precision bits of u are lost and we wind up with s = 0.5.
  // S2CellId(const S2Point&) then chooses an arbitrary neighboring cell.
  //
  // This tests that S2Cell::contains() expands the cell bounds sufficiently
  // so that the returned cell is still considered to contain "p".
  S2Point p = S2LatLng.fromDegrees(-2, 90).toS2Point();
  S2CellId cell_id = new S2CellId(p).parent(1);
  auto cell = new S2Cell(cell_id);
  Assert.equal(cell.contains(p), true);
}

S1ChordAngle getDistanceToPointBruteForce(in S2Cell cell, in S2Point target) {
  S1ChordAngle min_distance = S1ChordAngle.infinity();
  for (int i = 0; i < 4; ++i) {
    updateMinDistance(target, cell.getVertex(i), cell.getVertex(i + 1), min_distance);
  }
  return min_distance;
}

S1ChordAngle getMaxDistanceToPointBruteForce(in S2Cell cell, in S2Point target) {
  if (cell.contains(-target)) {
    return S1ChordAngle.straight();
  }
  S1ChordAngle max_distance = S1ChordAngle.negative();
  for (int i = 0; i < 4; ++i) {
    updateMaxDistance(target, cell.getVertex(i), cell.getVertex(i + 1), max_distance);
  }
  return max_distance;
}

@("S2Cell.GetDistanceToPoint") unittest {
  S2Testing.rnd.reset(s2RandomSeed);
  for (int iter = 0; iter < 1000; ++iter) {
    logger.logTrace("Iteration ", iter);
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    S2Point target = S2Testing.randomPoint();
    S1Angle expected_to_boundary = getDistanceToPointBruteForce(cell, target).toS1Angle();
    S1Angle expected_to_interior = cell.contains(target) ? S1Angle.zero() : expected_to_boundary;
    S1Angle expected_max = getMaxDistanceToPointBruteForce(cell, target).toS1Angle();
    S1Angle actual_to_boundary = cell.getBoundaryDistance(target).toS1Angle();
    S1Angle actual_to_interior = cell.getDistance(target).toS1Angle();
    S1Angle actual_max = cell.getMaxDistance(target).toS1Angle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    Assert.approximately(actual_to_boundary.radians(), expected_to_boundary.radians(), 1e-12);
    Assert.approximately(actual_to_interior.radians(), expected_to_interior.radians(), 1e-12);
    Assert.approximately(actual_max.radians(), expected_max.radians(), 1e-12);
    if (expected_to_boundary.radians() <= M_PI / 3) {
      Assert.approximately(actual_to_boundary.radians(), expected_to_boundary.radians(), 1e-15);
      Assert.approximately(actual_to_interior.radians(), expected_to_interior.radians(), 1e-15);
    }
    if (expected_max.radians() <= M_PI / 3) {
      Assert.approximately(actual_max.radians(), expected_max.radians(), 1e-15);
    }
  }
}

static void chooseEdgeNearCell(in S2Cell cell, out S2Point a, out S2Point b) {
  S2Cap cap = cell.getCapBound();
  if (S2Testing.rnd.oneIn(5)) {
    // Choose a point anywhere on the sphere.
    a = S2Testing.randomPoint();
  } else {
    // Choose a point inside or somewhere near the cell.
    a = S2Testing.samplePoint(new S2Cap(cap.center(), 1.5 * cap.getRadius()));
  }
  // Now choose a maximum edge length ranging from very short to very long
  // relative to the cell size, and choose the other endpoint.
  double max_length = algorithm.min(
      100.0 * math.pow(1e-4, S2Testing.rnd.randDouble()) * cap.getRadius().radians(), M_PI_2);
  b = S2Testing.samplePoint(new S2Cap(a, S1Angle.fromRadians(max_length)));

  if (S2Testing.rnd.oneIn(20)) {
    // Occasionally replace edge with antipodal edge.
    a = -a;
    b = -b;
  }
}

static S1ChordAngle getDistanceToEdgeBruteForce(in S2Cell cell, in S2Point a, in S2Point b) {
  if (cell.contains(a) || cell.contains(b)) {
    return S1ChordAngle.zero();
  }

  S1ChordAngle min_dist = S1ChordAngle.infinity();
  for (int i = 0; i < 4; ++i) {
    S2Point v0 = cell.getVertex(i);
    S2Point v1 = cell.getVertex(i + 1);
    // If the edge crosses through the cell, max distance is 0.
    if (crossingSign(a, b, v0, v1) >= 0) {
      return S1ChordAngle.zero();
    }
    updateMinDistance(a, v0, v1, min_dist);
    updateMinDistance(b, v0, v1, min_dist);
    updateMinDistance(v0, a, b, min_dist);
  }
  return min_dist;
}

static S1ChordAngle getMaxDistanceToEdgeBruteForce(in S2Cell cell, in S2Point a, in S2Point b) {
  // If any antipodal endpoint is within the cell, the max distance is Pi.
  if (cell.contains(-a) || cell.contains(-b)) {
    return S1ChordAngle.straight();
  }

  S1ChordAngle max_dist = S1ChordAngle.negative();
  for (int i = 0; i < 4; ++i) {
    S2Point v0 = cell.getVertex(i);
    S2Point v1 = cell.getVertex(i + 1);
    // If the antipodal edge crosses through the cell, max distance is Pi.
    if (crossingSign(-a, -b, v0, v1) >= 0) {
      return S1ChordAngle.straight();
    }
    updateMaxDistance(a, v0, v1, max_dist);
    updateMaxDistance(b, v0, v1, max_dist);
    updateMaxDistance(v0, a, b, max_dist);
  }
  return max_dist;
}

@("S2Cell.GetDistanceToEdge") unittest {
  S2Testing.rnd.reset(s2RandomSeed);
  for (int iter = 0; iter < 1000; ++iter) {
    logger.logTrace("Iteration ", iter);
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    S2Point a, b;
    chooseEdgeNearCell(cell, a, b);
    S1Angle expected_min = getDistanceToEdgeBruteForce(cell, a, b).toS1Angle();
    S1Angle expected_max = getMaxDistanceToEdgeBruteForce(cell, a, b).toS1Angle();
    S1Angle actual_min = cell.getDistance(a, b).toS1Angle();
    S1Angle actual_max = cell.getMaxDistance(a, b).toS1Angle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    if (expected_min.radians() > M_PI/2) {
      // Max error for S1ChordAngle as it approaches Pi is about 2e-8.
      Assert.approximately(actual_min.radians(), expected_min.radians(), 3e-8);
    } else if (expected_min.radians() <= M_PI / 3) {
      Assert.approximately(actual_min.radians(), expected_min.radians(), 1e-15);
    } else {
      Assert.approximately(actual_min.radians(), expected_min.radians(), 1e-12);
    }

    Assert.approximately(actual_max.radians(), expected_max.radians(), 1e-12);
    if (expected_max.radians() <= M_PI / 3) {
      Assert.approximately(actual_max.radians(), expected_max.radians(), 1e-15);
    }
  }
}

@("S2Cell.GetMaxDistanceToEdge") unittest {
  // Test an edge for which its antipode crosses the cell. Validates both the
  // standard and brute force implementations for this case.
  auto cell = S2Cell.fromFacePosLevel(0, 0, 20);
  S2Point a = -interpolate(2.0, cell.getCenter(), cell.getVertex(0));
  S2Point b = -interpolate(2.0, cell.getCenter(), cell.getVertex(2));

  S1ChordAngle actual = cell.getMaxDistance(a, b);
  S1ChordAngle expected = getMaxDistanceToEdgeBruteForce(cell, a, b);

  Assert.approximately(expected.radians(), S1ChordAngle.straight().radians(), 1e-15);
  Assert.approximately(actual.radians(), S1ChordAngle.straight().radians(), 1e-15);
}

@("S2Cell.GetMaxDistanceToCellAntipodal") unittest {
  S2Point p = makePointOrDie("0:0");
  auto cell = new S2Cell(p);
  auto antipodal_cell = new S2Cell(-p);
  S1ChordAngle dist = cell.getMaxDistance(antipodal_cell);
  Assert.equal(S1ChordAngle.straight(), dist);
}

@("S2Cell.GetMaxDistanceToCell") unittest {
  for (int i = 0; i < 1000; i++) {
    auto cell = new S2Cell(S2Testing.getRandomCellId());
    auto test_cell = new S2Cell(S2Testing.getRandomCellId());
    auto antipodal_leaf_id = S2CellId(-test_cell.getCenter());
    auto antipodal_test_cell = new S2Cell(antipodal_leaf_id.parent(test_cell.level()));

    auto dist_from_min = S1ChordAngle.straight() - cell.getDistance(antipodal_test_cell);
    auto dist_from_max = cell.getMaxDistance(test_cell);
    Assert.approximately(dist_from_min.radians(), dist_from_max.radians(), 1e-8);
  }
}

@("S2Cell.EncodeDecode") unittest {
  auto orig_cell = new S2Cell(S2LatLng.fromDegrees(40.7406264, -74.0029963));
  auto encoder = makeEncoder();
  orig_cell.encode(encoder);

  auto decoded_cell = new S2Cell(S2LatLng.fromDegrees(51.494987, -0.146585));
  auto decoder = makeDecoder(encoder.buffer().data());
  decoded_cell.decode(decoder);

  Assert.equal(orig_cell, decoded_cell);
  Assert.equal(orig_cell.face(), decoded_cell.face());
  Assert.equal(orig_cell.level(), decoded_cell.level());
  Assert.equal(orig_cell.orientation(), decoded_cell.orientation());
  Assert.equal(orig_cell.id(), decoded_cell.id());
  Assert.equal(orig_cell.getBoundUV(), decoded_cell.getBoundUV());
}
