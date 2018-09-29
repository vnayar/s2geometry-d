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

module s2.s2padded_cell_test;

import s2.r1interval;
import s2.r2point;
import s2.r2rect;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2padded_cell;
import s2.s2testing;

import std.math : pow;
import std.algorithm : min;

import fluent.asserts;

void compareS2CellToPadded(in S2Cell cell, S2PaddedCell pcell, double padding) {
  Assert.equal(cell.id(), pcell.id());
  Assert.equal(cell.level(), pcell.level());
  Assert.equal(padding, pcell.padding());
  Assert.equal(cell.getBoundUV().expanded(padding), pcell.bound());
  R2Point center_uv = cell.id().getCenterUV();
  Assert.equal(R2Rect.fromPoint(center_uv).expanded(padding), pcell.middle());
  Assert.equal(cell.getCenter(), pcell.getCenter());
}

@("S2PaddedCell.S2CellMethods")
unittest {
  // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
  enum int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2CellId id = S2Testing.getRandomCellId();
    double padding = pow(1e-15, S2Testing.rnd.randDouble());
    auto cell = new S2Cell(id);
    auto pcell = new S2PaddedCell(id, padding);
    compareS2CellToPadded(cell, pcell, padding);
    if (!id.isLeaf()) {
      S2Cell[4] children;
      Assert.equal(cell.subdivide(children), true);
      for (int pos = 0; pos < 4; ++pos) {
        int i, j;
        pcell.getChildIJ(pos, i, j);
        compareS2CellToPadded(children[pos], new S2PaddedCell(pcell, i, j), padding);
      }
    }
  }
}

@("S2PaddedCell.GetEntryExitVertices")
unittest {
  const int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2CellId id = S2Testing.getRandomCellId();
    // Check that entry/exit vertices do not depend on padding.
    Assert.equal(
        new S2PaddedCell(id, 0).getEntryVertex(), new S2PaddedCell(id, 0.5).getEntryVertex());
    Assert.equal(
        new S2PaddedCell(id, 0).getExitVertex(), new S2PaddedCell(id, 0.5).getExitVertex());

    // Check that the exit vertex of one cell is the same as the entry vertex
    // of the immediately following cell.  (This also tests wrapping from the
    // end to the start of the S2CellId curve with high probability.)
    Assert.equal(
        new S2PaddedCell(id, 0).getExitVertex(),
        new S2PaddedCell(id.nextWrap(), 0).getEntryVertex());

    // Check that the entry vertex of a cell is the same as the entry vertex
    // of its first child, and similarly for the exit vertex.
    if (!id.isLeaf()) {
      Assert.equal(
          new S2PaddedCell(id, 0).getEntryVertex(),
          new S2PaddedCell(id.child(0), 0).getEntryVertex());
      Assert.equal(
          new S2PaddedCell(id, 0).getExitVertex(),
          new S2PaddedCell(id.child(3), 0).getExitVertex());
    }
  }
}

static double sampleInterval(in R1Interval x) {
  return S2Testing.rnd.uniformDouble(x.lo(), x.hi());
}

@("S2PaddedCell.ShrinkToFit")
unittest {
  const int kIters = 1000;
  auto rnd = S2Testing.rnd;
  for (int iter = 0; iter < kIters; ++iter) {
    // Start with the desired result and work backwards.
    S2CellId result = S2Testing.getRandomCellId();
    R2Rect result_uv = result.getBoundUV();
    R2Point size_uv = result_uv.getSize();

    // Find the biggest rectangle that fits in "result" after padding.
    // (These calculations ignore numerical errors.)
    double max_padding = 0.5 * min(size_uv[0], size_uv[1]);
    double padding = max_padding * rnd.randDouble();
    R2Rect max_rect = result_uv.expanded(-padding);

    // Start with a random subset of the maximum rectangle.
    auto a = R2Point(sampleInterval(max_rect[0]), sampleInterval(max_rect[1]));
    auto b = R2Point(sampleInterval(max_rect[0]), sampleInterval(max_rect[1]));
    if (!result.isLeaf()) {
      // If the result is not a leaf cell, we must ensure that no child of
      // "result" also satisfies the conditions of ShrinkToFit().  We do this
      // by ensuring that "rect" intersects at least two children of "result"
      // (after padding).
      size_t axis = rnd.uniform(2);
      double center = result.getCenterUV()[axis];

      // Find the range of coordinates that are shared between child cells
      // along that axis.
      auto sharedInterval = R1Interval(center - padding, center + padding);
      double midInterval = sampleInterval(sharedInterval.intersection(max_rect[axis]));
      a[axis] = sampleInterval(R1Interval(max_rect[axis].lo(), midInterval));
      b[axis] = sampleInterval(R1Interval(midInterval, max_rect[axis].hi()));
    }
    R2Rect rect = R2Rect.fromPointPair(a, b);

    // Choose an arbitrary ancestor as the S2PaddedCell.
    S2CellId initial_id = result.parent(rnd.uniform(result.level() + 1));
    Assert.equal(result, new S2PaddedCell(initial_id, padding).shrinkToFit(rect));
  }
}
