// Copyright 2016 Google Inc. All Rights Reserved.
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
// The bulk of this file consists of "tests" that attempt to construct worst
// cases for the various constants used in S2CellIdSnapFunction and
// IntLatLngSnapFunction implementations.  For all of these constants I have
// done hand analysis of the planar configurations, but sometimes the
// spherical case is slightly better or worse because of the spherical
// distortion.

module s2.builder.util.snap_functions_test;

import s2.builder.util.snap_functions;

import s2.r2point;
import s2.s1angle;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2edge_distances;
import s2.s2latlng;
import s2.s2measures;
import s2.s2metrics;
import s2.s2testing;
import s2.s2text_format;
import s2.s2point;

import std.algorithm;
import std.array;
import std.math;
import std.stdio;

import fluent.asserts;

@("S2CellIdSnapFunction.LevelToFromSnapRadius") unittest {
  for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
    S1Angle radius = S2CellIdSnapFunction.minSnapRadiusForLevel(level);
    Assert.equal(level, S2CellIdSnapFunction.levelForMaxSnapRadius(radius));
    Assert.equal(min(level + 1, S2CellId.MAX_LEVEL),
              S2CellIdSnapFunction.levelForMaxSnapRadius(0.999 * radius));
  }
  Assert.equal(0, S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.fromRadians(5)));
  Assert.equal(S2CellId.MAX_LEVEL,
      S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.fromRadians(1e-30)));
}

@("S2CellIdSnapFunction.SnapPoint") unittest {
  for (int iter = 0; iter < 1000; ++iter) {
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      // This checks that points are snapped to the correct level, since
      // S2CellId centers at different levels are always different.
      auto f = new S2CellIdSnapFunction(level);
      S2Point p = S2Testing.getRandomCellId(level).toS2Point();
      Assert.equal(p, f.snapPoint(p));
    }
  }
}

@("IntLatLngSnapFunction.ExponentToFromSnapRadius") unittest {
  for (int exponent = IntLatLngSnapFunction.MIN_EXPONENT;
       exponent <= IntLatLngSnapFunction.MAX_EXPONENT; ++exponent) {
    S1Angle radius = IntLatLngSnapFunction.minSnapRadiusForExponent(exponent);
    Assert.equal(exponent,
        IntLatLngSnapFunction.exponentForMaxSnapRadius(radius));
    Assert.equal(min(exponent + 1, IntLatLngSnapFunction.MAX_EXPONENT),
        IntLatLngSnapFunction.exponentForMaxSnapRadius(0.999 * radius));
  }
  Assert.equal(IntLatLngSnapFunction.MIN_EXPONENT,
      IntLatLngSnapFunction.exponentForMaxSnapRadius(
          S1Angle.fromRadians(5)));
  Assert.equal(IntLatLngSnapFunction.MAX_EXPONENT,
      IntLatLngSnapFunction.exponentForMaxSnapRadius(
          S1Angle.fromRadians(1e-30)));
}

@("IntLatLngSnapFunction.SnapPoint") unittest {
  for (int iter = 0; iter < 1000; ++iter) {
    // Test that IntLatLngSnapFunction does not modify points that were
    // generated using the S2LatLng::From{E5,E6,E7} methods.  This ensures
    // that both functions are using bitwise-compatible conversion methods.
    S2Point p = S2Testing.randomPoint();
    auto ll = S2LatLng(p);
    S2Point p5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5()).toS2Point();
    Assert.equal(p5, new IntLatLngSnapFunction(5).snapPoint(p5));
    S2Point p6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6()).toS2Point();
    Assert.equal(p6, new IntLatLngSnapFunction(6).snapPoint(p6));
    S2Point p7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7()).toS2Point();
    Assert.equal(p7, new IntLatLngSnapFunction(7).snapPoint(p7));

    // Make sure that we're not snapping using some lower exponent.
    S2Point p7not6 = S2LatLng.fromE7(10 * ll.lat().e6() + 1,
        10 * ll.lng().e6() + 1).toS2Point();
    Assert.notEqual(p7not6, new IntLatLngSnapFunction(6).snapPoint(p7not6));
  }
}

private S2CellId kSearchRootId = S2CellId.fromFace(0);
private S2CellId kSearchFocusId = S2CellId.fromFace(0).child(3);

private S1Angle getMaxVertexDistance(in S2Point p, S2CellId id) {
  auto cell = new S2Cell(id);
  return max(max(S1Angle(p, cell.getVertex(0)), S1Angle(p, cell.getVertex(1))),
      max(S1Angle(p, cell.getVertex(2)), S1Angle(p, cell.getVertex(3))));
}

struct Score {
  double first;
  S2CellId second;

  int opCmp(in Score o) {
    if (first < o.first) return -1;
    else if (first > o.first) return 1;
    else if (second < o.second) return -1;
    else if (second > o.second) return 1;
    return 0;
  }
}

/// Helper function that computes the vertex separation between "id0" and its
/// neighbors.
private void updateS2CellIdMinVertexSeparation(S2CellId id0, ref Score[] scores) {
  S2Point site0 = id0.toS2Point();
  S2CellId[] nbrs;
  id0.appendAllNeighbors(id0.level(), nbrs);
  foreach (S2CellId id1; nbrs) {
    S2Point site1 = id1.toS2Point();
    auto vertex_sep = S1Angle(site0, site1);
    S1Angle max_snap_radius = getMaxVertexDistance(site0, id1);
    assert(max_snap_radius >= S2CellIdSnapFunction.minSnapRadiusForLevel(id0.level()));
    double r = (vertex_sep / max_snap_radius).radians();
    scores ~= Score(r, id0);
  }
}

private double getS2CellIdMinVertexSeparation(int level, ref bool[S2CellId] best_cells) {
  // The worst-case separation ratios always occur when the snap_radius is not
  // much larger than the minimum, since this allows the site spacing to be
  // reduced by as large a fraction as possible.
  //
  // For the minimum vertex separation ratio, we choose a site and one of its
  // 8-way neighbors, then look at the ratio of the distance to the center of
  // that neighbor to the distance to the furthest corner of that neighbor
  // (which is the largest possible snap radius for this configuration).
  Score[] scores;
  if (level == 0) {
    updateS2CellIdMinVertexSeparation(kSearchRootId, scores);
  } else {
    foreach (S2CellId parent; best_cells.keys) {
      for (S2CellId id0 = parent.childBegin(); id0 != parent.childEnd(); id0 = id0.next()) {
        updateS2CellIdMinVertexSeparation(id0, scores);
      }
    }
  }
  // Now sort the entries, print out the "num_to_print" best ones, and keep
  // the best "num_to_keep" of them to seed the next round.
  scores = sort(scores).uniq().array;
  best_cells.clear();
  int num_to_keep = 300;
  int num_to_print = 1;
  foreach (entry; scores) {
    S2CellId id = entry.second;
    if (--num_to_print >= 0) {
      R2Point uv = id.getCenterUV();
      writefln("Level %2d: min_vertex_sep_ratio = %.15f u=%.6f v=%.6f %s\n",
             level, entry.first, uv[0], uv[1], id.toToken());
    }
    if (kSearchFocusId.contains(id) || id.contains(kSearchFocusId)) {
      if (id !in best_cells) {
        best_cells[id] = true;
        if (--num_to_keep <= 0) break;
      }
    }
  }
  return scores[0].first;
}

@("S2CellIdSnapFunction.MinVertexSeparationSnapRadiusRatio") unittest {
  // The purpose of this "test" is to compute a lower bound to the fraction
  // (min_vertex_separation() / snap_radius()).  Essentially this involves
  // searching for two adjacent cells A and B such when one of the corner
  // vertices of B is snapped to the center of B, the distance to the center
  // of A decreases as much as possible.  In other words, we want the ratio
  //
  //   distance(center(A), center(B)) / distance(center(A), vertex(B))
  //
  // to be as small as possible.  We do this by considering one cell level at
  // a time, and remembering the cells that had the lowest ratios.  When we
  // proceed from one level to the next, we consider all the children of those
  // cells and keep the best ones.
  //
  // The reason we can restrict the search to children of cells at the
  // previous level is that the ratio above is essentially a function of the
  // local distortions created by projecting the S2 cube space onto the
  // sphere.  These distortions change smoothly over the sphere, so by keeping
  // a fairly large number of candidates ("num_to_keep"), we are essentially
  // keeping all the neighbors of the optimal cell as well.
  double best_score = 1e10;
  bool[S2CellId] best_cells;
  for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
    double score = getS2CellIdMinVertexSeparation(level, best_cells);
    best_score = min(best_score, score);
  }
  writefln("min_vertex_sep / snap_radius ratio: %.15f\n", best_score);
}

private S1Angle getCircumRadius(in S2Point a, in S2Point b, in S2Point c) {
  // We return this value is the circumradius is very large.
  S1Angle kTooBig = S1Angle.fromRadians(M_PI);
  double turn_angle = turnAngle(a, b, c);
  if (abs(remainder(turn_angle, M_PI)) < 1e-2) return kTooBig;

  real a2 = (b - c).norm2();
  real b2 = (c - a).norm2();
  real c2 = (a - b).norm2();
  if (a2 > 2 || b2 > 2 || c2 > 2) return kTooBig;
  real ma = a2 * (b2 + c2 - a2);
  real mb = b2 * (c2 + a2 - b2);
  real mc = c2 * (a2 + b2 - c2);
  S2Point p = (ma * a + mb * b + mc * c) / (ma + mb + mc);
  return S1Angle(p, a);
}

private S2CellId[] getNeighbors(S2CellId id) {
  const int kNumLayers = 2;
  S2CellId[] nbrs;
  nbrs ~= id;
  for (int layer = 0; layer < kNumLayers; ++layer) {
    S2CellId[] new_nbrs;
    foreach (S2CellId nbr; nbrs) {
      nbr.appendAllNeighbors(id.level(), new_nbrs);
    }
    nbrs ~= new_nbrs;
    nbrs = nbrs.remove!(a => a == id)().sort().uniq().array();
  }
  return nbrs;
}

// S2CellIdMinEdgeSeparationFunction defines an objective function that will
// be optimized by GetS2CellIdMinEdgeSeparation() by finding worst-case
// configurations of S2CellIds.  We use this to find the worst cases under
// various conditions (e.g., when the minimum snap radius at a given level is
// being used).  The objective function is called for a specific configuration
// of vertices that are snapped at the given S2CellId level.  "edge_sep" is
// the edge-vertex distance that is achieved by this configuration, and
// "min_snap_radius" and "max_snap_radius" are the minimum and maximum snap
// radii for which this configuration is valid (i.e., where the desired
// snapping will take place).
alias S2CellIdMinEdgeSeparationFunction =
    double function(int level, S1Angle edge_sep, S1Angle min_snap_radius, S1Angle max_snap_radius);

/+ TODO: Resume here.

// Returns the minimum value of the given objective function over sets of
// nearby vertices that are designed to minimize the edge-vertex separation
// when an edge is snapped.
static double GetS2CellIdMinEdgeSeparation(
    const char* label, S2CellIdMinEdgeSeparationFunction objective,
    int level, set<S2CellId>* best_cells) {
  // To find minimum edge separations, we choose a cell ("id0") and two nearby
  // cells ("id1" and "id2"), where "nearby" is defined by GetNeighbors().
  // Let "site0", "site1", and "site2" be the centers of these cells.  The
  // idea is to consider an input edge E that intersects the Voronoi regions
  // of "site1" and "site2" (and therefore snaps to an edge E' between these
  // sites) but does not not intersect the Voronoi region of "site0" (and
  // therefore can't be snapped to site0).  The goal is to search for snapped
  // edges E' that approach site0 as closely as possible.
  //
  // To do this, we first compute the circumradius of the three cell centers
  // ("site0", "site1", and "site2"); this is the minimum snap radius in order
  // for it to be possible to construct an edge E that snaps to "site1" and
  // "site2" but not to "site0".  We also compute the distance from "site0" to
  // the snapped edge.  Next we find the corner vertex of "id1" and "id2" that
  // is furthest from "site0"; the smaller of these two distances is the
  // maximum snap radius such that "site1" and "site2" can be chosen as
  // sites after choosing "site0".  If the maximum is less than the minimum,
  // then this configuration is rejected; otherwise we evaluate the given
  // objective function and keep the configurations that result in the
  // smallest values.
  //
  // The optimization process works by keeping track of the set of S2CellIds
  // that yielded the best results at the previous level, and exploring all
  // the nearby neighbor combinations of the children of those cells at the
  // next level.  In order to get better coverage, we keep track of the best
  // score and configuration (i.e. the two neighboring cells "id1" and "id2")
  // for each initial cell "id0".
  map<S2CellId, double> best_scores;
  map<S2CellId, pair<S2CellId, S2CellId>> best_configs;
  for (S2CellId parent : *best_cells) {
    for (S2CellId id0 = parent.child_begin(level);
         id0 != parent.child_end(level); id0 = id0.next()) {
      S2Point site0 = id0.ToPoint();
      vector<S2CellId> nbrs = GetNeighbors(id0);
      for (S2CellId id1 : nbrs) {
        S2Point site1 = id1.ToPoint();
        S1Angle max_v1 = GetMaxVertexDistance(site0, id1);
        for (S2CellId id2 : nbrs) {
          if (id2 <= id1) continue;
          S2Point site2 = id2.ToPoint();
          S1Angle min_snap_radius = GetCircumRadius(site0, site1, site2);
          if (min_snap_radius > S2Builder::SnapFunction::kMaxSnapRadius()) {
            continue;
          }
          // Note that it is only the original points *before* snapping that
          // need to be at least "snap_radius" away from "site0".  The points
          // after snapping ("site1" and "site2") may be closer.
          S1Angle max_v2 = GetMaxVertexDistance(site0, id2);
          S1Angle max_snap_radius = min(max_v1, max_v2);
          if (min_snap_radius > max_snap_radius) continue;
          DCHECK_GE(max_snap_radius,
                    S2CellIdSnapFunction::MinSnapRadiusForLevel(level));

          // This is a valid configuration, so evaluate it.
          S1Angle edge_sep = S2::GetDistance(site0, site1, site2);
          double score = objective(level, edge_sep,
                                   min_snap_radius, max_snap_radius);
          double& best_score = best_scores[id0];
          if (best_score == 0 || best_score > score) {
            best_score = score;
            best_configs[id0] = make_pair(id1, id2);
          }
        }
      }
    }
  }
  // Now sort the entries, print out the "num_to_print" best ones, and
  // generate a set of candidates for the next round by generating all the
  // 8-way neighbors of the best candidates, and keeping up to"num_to_keep" of
  // them.  The results vary slightly according to how many candidates we
  // keep, but the variations are much smaller than the conservative
  // assumptions made by the S2CellIdSnapFunction implementation.
  int num_to_keep = google::DEBUG_MODE ? 20 : 100;
  int num_to_print = 3;
  vector<pair<double, S2CellId>> sorted;
  for (const auto& entry : best_scores) {
    sorted.push_back(make_pair(entry.second, entry.first));
  }
  std::sort(sorted.begin(), sorted.end());
  best_cells->clear();
  printf("Level %d:\n", level);
  for (const auto& entry : sorted) {
    S2CellId id = entry.second;
    if (--num_to_print >= 0) {
      R2Point uv = id.GetCenterUV();
      const pair<S2CellId, S2CellId>& nbrs = best_configs.find(id)->second;
      printf("  %s = %.15f u=%7.4f v=%7.4f %s %s %s\n",
             label, entry.first, uv[0], uv[1], id.ToToken().c_str(),
             nbrs.first.ToToken().c_str(), nbrs.second.ToToken().c_str());
    }
    vector<S2CellId> nbrs(1, id);
    id.AppendAllNeighbors(id.level(), &nbrs);
    for (S2CellId nbr : nbrs) {
      // The S2Cell hierarchy has many regions that are symmetrical.  We can
      // eliminate most of the "duplicates" by restricting the search to cells
      // in kS2CellIdFocus.
      if (kSearchFocusId.contains(nbr) || nbr.contains(kSearchFocusId)) {
        if (best_cells->insert(nbr).second && --num_to_keep <= 0) {
          return sorted[0].first;
        }
      }
    }
  }
  return sorted[0].first;
}

static double GetS2CellIdMinEdgeSeparation(
    const char* label, S2CellIdMinEdgeSeparationFunction objective) {
  double best_score = 1e10;
  set<S2CellId> best_cells;
  best_cells.insert(kSearchRootId);
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    double score = GetS2CellIdMinEdgeSeparation(label, objective, level,
                                                &best_cells);
    best_score = min(best_score, score);
  }
  return best_score;
}

TEST(S2CellIdSnapFunction, MinEdgeVertexSeparationForLevel) {
  // Computes the minimum edge separation (as a fraction of kMinDiag) for any
  // snap radius at each level.
  double score = GetS2CellIdMinEdgeSeparation("min_sep_for_level",
                                              [](int level, S1Angle edge_sep,
                                                 S1Angle min_snap_radius,
                                                 S1Angle max_snap_radius) {
    return edge_sep.radians() / S2::kMinDiag.GetValue(level);
  });
  printf("min_edge_vertex_sep / kMinDiag ratio: %.15f\n", score);
}
TEST(S2CellIdSnapFunction, MinEdgeVertexSeparationAtMinSnapRadius) {
  // Computes the minimum edge separation (as a fraction of kMinDiag) for the
  // special case where the minimum snap radius is being used.
  double score = GetS2CellIdMinEdgeSeparation("min_sep_at_min_radius",
                                              [](int level, S1Angle edge_sep,
                                                 S1Angle min_snap_radius,
                                                 S1Angle max_snap_radius) {
    double min_radius_at_level = S2::kMaxDiag.GetValue(level) / 2;
    return (min_snap_radius.radians() <= (1 + 1e-10) * min_radius_at_level) ?
        (edge_sep.radians() / S2::kMinDiag.GetValue(level)) : 100.0;
  });
  printf("min_edge_vertex_sep / kMinDiag at MinSnapRadiusForLevel: %.15f\n",
         score);
}

TEST(S2CellIdSnapFunction, MinEdgeVertexSeparationSnapRadiusRatio) {
  // Computes the minimum edge separation expressed as a fraction of the
  // maximum snap radius that could yield that edge separation.
  double score = GetS2CellIdMinEdgeSeparation("min_sep_snap_radius_ratio",
                                              [](int level, S1Angle edge_sep,
                                                 S1Angle min_snap_radius,
                                                 S1Angle max_snap_radius) {
    return edge_sep.radians() / max_snap_radius.radians();
  });
  printf("min_edge_vertex_sep / snap_radius ratio: %.15f\n", score);
}

// A scaled S2LatLng with integer coordinates, similar to E7 coordinates,
// except that the scale is variable (see LatLngConfig below).
using IntLatLng = Vector2<int64>;

static bool IsValid(const IntLatLng& ll, int64 scale) {
  // A coordinate value of "scale" corresponds to 180 degrees.
  return (std::abs(ll[0]) <= scale / 2 && std::abs(ll[1]) <= scale);
}

static bool HasValidVertices(const IntLatLng& ll, int64 scale) {
  // Like IsValid, but excludes latitudes of 90 and longitudes of 180.
  // A coordinate value of "scale" corresponds to 180 degrees.
  return (std::abs(ll[0]) < scale / 2 && std::abs(ll[1]) < scale);
}

static IntLatLng Rescale(const IntLatLng&ll, double scale_factor) {
  return IntLatLng(MathUtil::FastInt64Round(scale_factor * ll[0]),
                   MathUtil::FastInt64Round(scale_factor * ll[1]));
}

static S2Point ToPoint(const IntLatLng& ll, int64 scale) {
  return S2LatLng::FromRadians(ll[0] * (M_PI / scale),
                               ll[1] * (M_PI / scale)).ToPoint();
}

static S2Point GetVertex(const IntLatLng& ll, int64 scale, int i) {
  // Return the points in CCW order starting from the lower left.
  int dlat = (i == 0 || i == 3) ? -1 : 1;
  int dlng = (i == 0 || i == 1) ? -1 : 1;
  return ToPoint(2 * ll + IntLatLng(dlat, dlng), 2 * scale);
}

static S1Angle GetMaxVertexDistance(const S2Point& p,
                                    const IntLatLng& ll, int64 scale) {
  return max(max(S1Angle(p, GetVertex(ll, scale, 0)),
                 S1Angle(p, GetVertex(ll, scale, 1))),
             max(S1Angle(p, GetVertex(ll, scale, 2)),
                 S1Angle(p, GetVertex(ll, scale, 3))));
}

static double GetLatLngMinVertexSeparation(int64 old_scale, int64 scale,
                                           set<IntLatLng>* best_configs) {
  // The worst-case separation ratios always occur when the snap_radius is not
  // much larger than the minimum, since this allows the site spacing to be
  // reduced by as large a fraction as possible.
  //
  // For the minimum vertex separation ratio, we choose a site and one of its
  // 8-way neighbors, then look at the ratio of the distance to the center of
  // that neighbor to the distance to the furthest corner of that neighbor
  // (which is the largest possible snap radius for this configuration).
  S1Angle min_snap_radius_at_scale = S1Angle::Radians(M_SQRT1_2 * M_PI / scale);
  vector<pair<double, IntLatLng>> scores;
  double scale_factor = static_cast<double>(scale) / old_scale;
  for (const IntLatLng& parent : *best_configs) {
    IntLatLng new_parent = Rescale(parent, scale_factor);
    for (int dlat0 = -7; dlat0 <= 7; ++dlat0) {
      IntLatLng ll0 = new_parent + IntLatLng(dlat0, 0);
      if (!IsValid(ll0, scale) || ll0[0] < 0) continue;
      S2Point site0 = ToPoint(ll0, scale);
      for (int dlat1 = 0; dlat1 <= 2; ++dlat1) {
        for (int dlng1 = 0; dlng1 <= 5; ++dlng1) {
          IntLatLng ll1 = ll0 + IntLatLng(dlat1, dlng1);
          if (ll1 == ll0 || !HasValidVertices(ll1, scale)) continue;
          S1Angle max_snap_radius = GetMaxVertexDistance(site0, ll1, scale);
          if (max_snap_radius < min_snap_radius_at_scale) continue;
          S2Point site1 = ToPoint(ll1, scale);
          S1Angle vertex_sep(site0, site1);
          double r = vertex_sep / max_snap_radius;
          scores.push_back(make_pair(r, ll0));
        }
      }
    }
  }
  // Now sort the entries, print out the "num_to_print" best ones, and keep
  // the best "num_to_keep" of them to seed the next round.
  std::sort(scores.begin(), scores.end());
  scores.erase(std::unique(scores.begin(), scores.end()), scores.end());
  best_configs->clear();
  int num_to_keep = 100;
  int num_to_print = 1;
  for (const auto& entry : scores) {
    if (--num_to_print >= 0) {
      printf("Scale %14" PRId64 ": min_vertex_sep_ratio = %.15f, %s\n",
             static_cast<int64_t>(scale), entry.first,
             s2textformat::ToString(ToPoint(entry.second, scale)).c_str());
    }
    if (best_configs->insert(entry.second).second && --num_to_keep <= 0) break;
  }
  return scores[0].first;
}

TEST(IntLatLngSnapFunction, MinVertexSeparationSnapRadiusRatio) {
  double best_score = 1e10;
  set<IntLatLng> best_configs;
  int64 scale = 18;
  for (int lat0 = 0; lat0 <= 9; ++lat0) {
    best_configs.insert(IntLatLng(lat0, 0));
  }
  for (int exp = 0; exp <= 10; ++exp, scale *= 10) {
    double score = GetLatLngMinVertexSeparation(scale, 10 * scale,
                                                &best_configs);
    best_score = min(best_score, score);
  }
  printf("min_vertex_sep / snap_radius ratio: %.15f\n", best_score);
}

// A triple of scaled S2LatLng coordinates.  The coordinates are multiplied by
// (M_PI / scale) to convert them to radians.
struct LatLngConfig {
  int64 scale;
  IntLatLng ll0, ll1, ll2;

  LatLngConfig(int64 _scale, const IntLatLng& _ll0,
               const IntLatLng& _ll1, const IntLatLng& _ll2) :
      scale(_scale), ll0(_ll0), ll1(_ll1), ll2(_ll2) {
  }
  bool operator<(const LatLngConfig& other) const {
    DCHECK_EQ(scale, other.scale);
    return (make_pair(ll0, make_pair(ll1, ll2)) <
            make_pair(other.ll0, make_pair(other.ll1, other.ll2)));
  }
  bool operator==(const LatLngConfig& other) const {
    return (ll0 == other.ll0 && ll1 == other.ll1 && ll2 == other.ll2);
  }
};

typedef double LatLngMinEdgeSeparationFunction(int64 scale, S1Angle edge_sep,
                                               S1Angle max_snap_radius);

static double GetLatLngMinEdgeSeparation(
    const char* label, LatLngMinEdgeSeparationFunction objective,
    int64 scale, vector<LatLngConfig>* best_configs) {
  S1Angle min_snap_radius_at_scale = S1Angle::Radians(M_SQRT1_2 * M_PI / scale);
  vector<pair<double, LatLngConfig>> scores;
  for (LatLngConfig parent : *best_configs) {
    // To reduce duplicates, we require that site0 always has longitude 0.
    DCHECK_EQ(0, parent.ll0[1]);
    double scale_factor = static_cast<double>(scale) / parent.scale;
    parent.ll0 = Rescale(parent.ll0, scale_factor);
    parent.ll1 = Rescale(parent.ll1, scale_factor);
    parent.ll2 = Rescale(parent.ll2, scale_factor);
    for (int dlat0 = -1; dlat0 <= 1; ++dlat0) {
      IntLatLng ll0 = parent.ll0 + IntLatLng(dlat0, 0);
      // To reduce duplicates, we require that site0.latitude >= 0.
      if (!IsValid(ll0, scale) || ll0[0] < 0) continue;
      S2Point site0 = ToPoint(ll0, scale);
      for (int dlat1 = -1; dlat1 <= 1; ++dlat1) {
        for (int dlng1 = -2; dlng1 <= 2; ++dlng1) {
          IntLatLng ll1 = parent.ll1 + IntLatLng(dlat0 + dlat1, dlng1);
          if (ll1 == ll0 || !HasValidVertices(ll1, scale)) continue;
          // Only consider neighbors within 2 latitude units of site0.
          if (std::abs(ll1[0] - ll0[0]) > 2) continue;

          S2Point site1 = ToPoint(ll1, scale);
          S1Angle max_v1 = GetMaxVertexDistance(site0, ll1, scale);
          for (int dlat2 = -1; dlat2 <= 1; ++dlat2) {
            for (int dlng2 = -2; dlng2 <= 2; ++dlng2) {
              IntLatLng ll2 = parent.ll2 + IntLatLng(dlat0 + dlat2, dlng2);
              if (!HasValidVertices(ll2, scale)) continue;
              // Only consider neighbors within 2 latitude units of site0.
              if (std::abs(ll2[0] - ll0[0]) > 2) continue;
              // To reduce duplicates, we require ll1 < ll2 lexicographically
              // and site2.longitude >= 0.  (It's *not* okay to
              // require site1.longitude >= 0, because then some configurations
              // with site1.latitude == site2.latitude would be missed.)
              if (ll2 <= ll1 || ll2[1] < 0) continue;

              S2Point site2 = ToPoint(ll2, scale);
              S1Angle min_snap_radius = GetCircumRadius(site0, site1, site2);
              if (min_snap_radius > S2Builder::SnapFunction::kMaxSnapRadius()) {
                continue;
              }
              // Only the original points *before* snapping that need to be at
              // least "snap_radius" away from "site0".  The points after
              // snapping ("site1" and "site2") may be closer.
              S1Angle max_v2 = GetMaxVertexDistance(site0, ll2, scale);
              S1Angle max_snap_radius = min(max_v1, max_v2);
              if (min_snap_radius > max_snap_radius) continue;
              if (max_snap_radius < min_snap_radius_at_scale) continue;

              // This is a valid configuration, so evaluate it.
              S1Angle edge_sep = S2::GetDistance(site0, site1, site2);
              double score = objective(scale, edge_sep, max_snap_radius);
              LatLngConfig config(scale, ll0, ll1, ll2);
              scores.push_back(make_pair(score, config));
            }
          }
        }
      }
    }
  }
  // Now sort the entries, print out the "num_to_print" best ones, and keep
  // the best "num_to_keep" of them to seed the next round.
  std::sort(scores.begin(), scores.end());
  scores.erase(std::unique(scores.begin(), scores.end()), scores.end());
  best_configs->clear();
  int num_to_keep = google::DEBUG_MODE ? 50 : 200;
  int num_to_print = 3;
  printf("Scale %" PRId64 ":\n", static_cast<int64_t>(scale));
  for (const auto& entry : scores) {
    const LatLngConfig& config = entry.second;
    int64 scale = config.scale;
    if (--num_to_print >= 0) {
      printf("  %s = %.15f %s %s %s\n",
             label, entry.first,
             s2textformat::ToString(ToPoint(config.ll0, scale)).c_str(),
             s2textformat::ToString(ToPoint(config.ll1, scale)).c_str(),
             s2textformat::ToString(ToPoint(config.ll2, scale)).c_str());
    }
    // Optional: filter the candidates to concentrate on a specific region
    // (e.g., the north pole).
    best_configs->push_back(config);
    if (--num_to_keep <= 0) break;
  }
  return scores[0].first;
}

static double GetLatLngMinEdgeSeparation(
    const char* label, LatLngMinEdgeSeparationFunction objective) {
  double best_score = 1e10;
  vector<LatLngConfig> best_configs;
  int64 scale = 6;  // Initially points are 30 degrees apart.
  int max_lng = scale;
  int max_lat = scale / 2;
  for (int lat0 = 0; lat0 <= max_lat; ++lat0) {
    for (int lat1 = lat0 - 2; lat1 <= min(max_lat, lat0 + 2); ++lat1) {
      for (int lng1 = 0; lng1 <= max_lng; ++lng1) {
        for (int lat2 = lat1; lat2 <= min(max_lat, lat0 + 2); ++lat2) {
          for (int lng2 = 0; lng2 <= max_lng; ++lng2) {
            IntLatLng ll0(lat0, 0);
            IntLatLng ll1(lat1, lng1);
            IntLatLng ll2(lat2, lng2);
            if (ll2 <= ll1) continue;
            best_configs.push_back(LatLngConfig(scale, ll0, ll1, ll2));
          }
        }
      }
    }
  }
  LOG(INFO) << "Starting with " << best_configs.size() << " configurations";
  int64 target_scale = 180;
  for (int exp = 0; exp <= 10; ++exp, target_scale *= 10) {
    while (scale < target_scale) {
      scale = min(static_cast<int64>(1.8 * scale), target_scale);
      double score = GetLatLngMinEdgeSeparation(label, objective, scale,
                                                &best_configs);
      if (scale == target_scale) {
        best_score = min(best_score, score);
      }
    }
  }
  return best_score;
}

TEST(IntLatLngSnapFunction, MinEdgeVertexSeparationForLevel) {
  // Computes the minimum edge separation (as a fraction of kMinDiag) for any
  // snap radius at each level.
  double score = GetLatLngMinEdgeSeparation("min_sep_for_level",
                                            [](int64 scale, S1Angle edge_sep,
                                               S1Angle max_snap_radius) {
    double e_unit = M_PI / scale;
    return edge_sep.radians() / e_unit;
  });
  printf("min_edge_vertex_sep / e_unit ratio: %.15f\n", score);
}

TEST(IntLatLngSnapFunction, MinEdgeVertexSeparationSnapRadiusRatio) {
  // Computes the minimum edge separation expressed as a fraction of the
  // maximum snap radius that could yield that edge separation.
  double score = GetLatLngMinEdgeSeparation("min_sep_snap_radius_ratio",
                                              [](int64 scale, S1Angle edge_sep,
                                                 S1Angle max_snap_radius) {
    return edge_sep.radians() / max_snap_radius.radians();
  });
  printf("min_edge_vertex_sep / snap_radius ratio: %.15f\n", score);
}
+/
