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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2metrics;

// The following are various constants that describe the shapes and sizes of
// S2Cells (see s2coords.h and s2cell_id.h).  They are useful for deciding
// which cell level to use in order to satisfy a given condition (e.g. that
// cell vertices must be no further than "x" apart).  All of the raw constants
// are differential quantities; you can use the GetValue(level) method to
// compute the corresponding length or area on the unit sphere for cells at a
// given level.  The minimum and maximum bounds are valid for cells at all
// levels, but they may be somewhat conservative for very large cells
// (e.g. face cells).

import s2coords = s2.s2coords;
import algorithm = std.algorithm;
import math = std.math;
import s2.util.math.s2const;

// Defines a cell metric of the given dimension (1 == length, 2 == area).
struct Metric(int DimV) {
private:
  immutable double _deriv;

public:
  this(double deriv) {
    _deriv = deriv;
  }

  // The "deriv" value of a metric is a derivative, and must be multiplied by
  // a length or area in (s,t)-space to get a useful value.
  @property
  double deriv() const {
    return _deriv;
  }

  // Return the value of a metric for cells at the given level. The value is
  // either a length or an area on the unit sphere, depending on the
  // particular metric.
  double getValue(int level) const {
    return math.ldexp(_deriv, -DimV * level);
  }

  // Return the level at which the metric has approximately the given
  // value.  For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the
  // level at which the average cell edge length is approximately 0.1.
  // The return value is always a valid level.
  int getClosestLevel(double value) const {
    return getLevelForMaxValue((DimV == 1 ? math.SQRT2 : 2) * value);
  }

  // Return the minimum level such that the metric is at most the given
  // value, or S2CellId::kMaxLevel if there is no such level.  For example,
  // S2::kMaxDiag.GetLevelForMaxValue(0.1) returns the minimum level such
  // that all cell diagonal lengths are 0.1 or smaller.  The return value
  // is always a valid level.
  int getLevelForMaxValue(double value) const
  out (level) {
    assert(level == s2coords.MAX_CELL_LEVEL || getValue(level) <= value);
    assert(level == 0 || getValue(level - 1) > value);
  } do {
    if (value <= 0) {
      return s2coords.MAX_CELL_LEVEL;
    }

    // This code is equivalent to computing a floating-point "level"
    // value and rounding up.  ilogb() returns the exponent corresponding to a
    // fraction in the range [1,2).
    int level = math.ilogb(value / _deriv);
    level = algorithm.max(0, algorithm.min(s2coords.MAX_CELL_LEVEL, -(level >> (DimV - 1))));
    return level;
  }

  // Return the maximum level such that the metric is at least the given
  // value, or zero if there is no such level.  For example,
  // S2::kMinWidth.GetLevelForMinValue(0.1) returns the maximum level such
  // that all cells have a minimum width of 0.1 or larger.  The return value
  // is always a valid level.
  int getLevelForMinValue(double value) const
  out (level) {
    assert(level == 0 || getValue(level) >= value);
    assert(level == s2coords.MAX_CELL_LEVEL || getValue(level + 1) < value);
  } do {
    if (value <= 0) {
      return s2coords.MAX_CELL_LEVEL;
    }

    // This code is equivalent to computing a floating-point "level"
    // value and rounding down.
    int level = math.ilogb(_deriv / value);
    level = algorithm.max(0, algorithm.min(s2coords.MAX_CELL_LEVEL, level >> (DimV - 1)));
    return level;
  }
}

alias LengthMetric = Metric!1;
alias AreaMetric = Metric!2;

version = S2_QUADRATIC_PROJECTION;

// Each cell is bounded by four planes passing through its four edges and
// the center of the sphere.  These metrics relate to the angle between each
// pair of opposite bounding planes, or equivalently, between the planes
// corresponding to two different s-values or two different t-values.  For
// example, the maximum angle between opposite bounding planes for a cell at
// level k is kMaxAngleSpan.GetValue(k), and the average angle span for all
// cells at level k is approximately kAvgAngleSpan.GetValue(k).
version (S2_LINEAR_PROJECTION) {
  immutable LengthMetric MIN_ANGLE_SPAN = LengthMetric(1.0);           // 1.000
  immutable LengthMetric MAX_ANGLE_SPAN = LengthMetric(2.0);           // 2.000
}
version (S2_TAN_PROJECTION) {
  immutable LengthMetric MIN_ANGLE_SPAN = LengthMetric(M_PI / 2);   // 1.571
  immutable LengthMetric MAX_ANGLE_SPAN = LengthMetric(M_PI / 2);   // 1.571
}
version (S2_QUADRATIC_PROJECTION) {
  immutable LengthMetric MIN_ANGLE_SPAN = LengthMetric(4.0 / 3);       // 1.333
  immutable LengthMetric MAX_ANGLE_SPAN = LengthMetric(1.704897179199218452); // 1.705
}

// This is true for all projections.
immutable LengthMetric AVG_ANGLE_SPAN = LengthMetric(M_PI / 2);     // 1.571

// The width of geometric figure is defined as the distance between two
// parallel bounding lines in a given direction.  For cells, the minimum
// width is always attained between two opposite edges, and the maximum
// width is attained between two opposite vertices.  However, for our
// purposes we redefine the width of a cell as the perpendicular distance
// between a pair of opposite edges.  A cell therefore has two widths, one
// in each direction.  The minimum width according to this definition agrees
// with the classic geometric one, but the maximum width is different.  (The
// maximum geometric width corresponds to kMaxDiag defined below.)
//
// For a cell at level k, the distance between opposite edges is at least
// kMinWidth.GetValue(k) and at most kMaxWidth.GetValue(k).  The average
// width in both directions for all cells at level k is approximately
// kAvgWidth.GetValue(k).
//
// The width is useful for bounding the minimum or maximum distance from a
// point on one edge of a cell to the closest point on the opposite edge.
// For example, this is useful when "growing" regions by a fixed distance.
//
// Note that because S2Cells are not usually rectangles, the minimum width of
// a cell is generally smaller than its minimum edge length.  (The interior
// angles of an S2Cell range from 60 to 120 degrees.)
version (S2_LINEAR_PROJECTION) {
  immutable LengthMetric MIN_WIDTH = LengthMetric(math.sqrt(2.0 / 3));             // 0.816
  immutable LengthMetric AVG_WIDTH = LengthMetric(1.411459345844456965);           // 1.411
}

version (S2_TAN_PROJECTION) {
  immutable LengthMetric MIN_WIDTH = LengthMetric(M_PI / (2 * math.sqrt(2.0))); // 1.111
  immutable LengthMetric AVG_WIDTH = LengthMetric(1.437318638925160885);           // 1.437
}

version (S2_QUADRATIC_PROJECTION) {
  immutable LengthMetric MIN_WIDTH = LengthMetric(2 * math.sqrt(2.0) / 3);         // 0.943
  immutable LengthMetric AVG_WIDTH = LengthMetric(1.434523672886099389);           // 1.435
}

// This is true for all projections.
immutable LengthMetric MAX_WIDTH = LengthMetric(MAX_ANGLE_SPAN.deriv());

// The minimum edge length of any cell at level k is at least
// kMinEdge.GetValue(k), and the maximum is at most kMaxEdge.GetValue(k).
// The average edge length is approximately kAvgEdge.GetValue(k).
//
// The edge length metrics can also be used to bound the minimum, maximum,
// or average distance from the center of one cell to the center of one of
// its edge neighbors.  In particular, it can be used to bound the distance
// between adjacent cell centers along the space-filling Hilbert curve for
// cells at any given level.
version (S2_LINEAR_PROJECTION) {
  immutable LengthMetric MIN_EDGE = LengthMetric(2 * math.sqrt(2.0) / 3);        // 0.943
  immutable LengthMetric AVG_EDGE = LengthMetric(1.440034192955603643);          // 1.440
}
version (S2_TAN_PROJECTION) {
  immutable LengthMetric MIN_EDGE = LengthMetric(M_PI / (2 * math.sqrt(2.0)));  // 1.111
  immutable LengthMetric AVG_EDGE = LengthMetric(1.461667032546739266);         // 1.462
}
version (S2_QUADRATIC_PROJECTION) {
  immutable LengthMetric MIN_EDGE = LengthMetric(2 * math.sqrt(2.0) / 3);       // 0.943
  immutable LengthMetric AVG_EDGE = LengthMetric(1.459213746386106062);         // 1.459
}
// This is true for all projections.
immutable LengthMetric MAX_EDGE = LengthMetric(MAX_ANGLE_SPAN.deriv());

// The minimum diagonal length of any cell at level k is at least
// kMinDiag.GetValue(k), and the maximum is at most kMaxDiag.GetValue(k).
// The average diagonal length is approximately kAvgDiag.GetValue(k).
//
// The maximum diagonal also happens to be the maximum diameter of any cell,
// and also the maximum geometric width (see the discussion above).  So for
// example, the distance from an arbitrary point to the closest cell center
// at a given level is at most half the maximum diagonal length.
version (S2_LINEAR_PROJECTION) {
  immutable LengthMetric MIN_DIAG = LengthMetric(2 * math.sqrt(2.0) / 3);        // 0.943
  immutable LengthMetric MAX_DIAG = LengthMetric(2 * math.sqrt(2.0));            // 2.828
  immutable LengthMetric AVG_DIAG = LengthMetric(2.031817866418812674);          // 2.032
}
version (S2_TAN_PROJECTION) {
  immutable LengthMetric MIN_DIAG = LengthMetric(M_PI * math.sqrt(2.0) / 3);     // 1.481
  immutable LengthMetric MAX_DIAG = LengthMetric(M_PI * math.sqrt(2.0 / 3));     // 2.565
  immutable LengthMetric AVG_DIAG = LengthMetric(2.063623197195635753);          // 2.064
}
version (S2_QUADRATIC_PROJECTION) {
  immutable LengthMetric MIN_DIAG = LengthMetric(8 * math.sqrt(2.0) / 9);        // 1.257
  immutable LengthMetric MAX_DIAG = LengthMetric(2.438654594434021032);          // 2.439
  immutable LengthMetric AVG_DIAG = LengthMetric(2.060422738998471683);          // 2.060
}

// The minimum area of any cell at level k is at least kMinArea.GetValue(k),
// and the maximum is at most kMaxArea.GetValue(k).  The average area of all
// cells at level k is exactly kAvgArea.GetValue(k).
version (S2_LINEAR_PROJECTION) {
  immutable AreaMetric MIN_AREA = AreaMetric(4 / (3 * math.sqrt(3)));        // 0.770
  immutable AreaMetric MAX_AREA = AreaMetric(4);                             // 4.000
}
version (S2_TAN_PROJECTION) {
  immutable AreaMetric MIN_AREA = AreaMetric((M_PI*M_PI) / (4*math.sqrt(2.0))); // 1.745
  immutable AreaMetric MAX_AREA = AreaMetric(M_PI * M_PI / 4);               // 2.467
}
version (S2_QUADRATIC_PROJECTION) {
  immutable AreaMetric MIN_AREA = AreaMetric(8 * math.sqrt(2.0) / 9);        // 1.257
  immutable AreaMetric MAX_AREA = AreaMetric(2.635799256963161491);          // 2.636
}
immutable AreaMetric AVG_AREA = AreaMetric(4 * M_PI / 6);                    // 2.094

// This is the maximum edge aspect ratio over all cells at any level, where
// the edge aspect ratio of a cell is defined as the ratio of its longest
// edge length to its shortest edge length.
version (S2_LINEAR_PROJECTION) {
  immutable double MAX_EDGE_ASPECT = math.sqrt(2);                               // 1.414
}
version (S2_TAN_PROJECTION) {
  immutable double MAX_EDGE_ASPECT = math.sqrt(2);                               // 1.414
}
version (S2_QUADRATIC_PROJECTION) {
  immutable double MAX_EDGE_ASPECT = 1.442615274452682920;                       // 1.443
}

// This is the maximum diagonal aspect ratio over all cells at any level,
// where the diagonal aspect ratio of a cell is defined as the ratio of its
// longest diagonal length to its shortest diagonal length.
immutable double MAX_DIAG_ASPECT = math.sqrt(3.0);                               // 1.732
