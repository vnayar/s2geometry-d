// Copyright 2018 Google Inc. All Rights Reserved.
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
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.s2shape_index_region_test;

import s2.mutable_s2shape_index;
import s2.r2rect;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2coords;
import s2.s2edge_clipping;
import s2.s2latlng_rect;
import s2.s2lax_loop_shape;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index_region;

import std.algorithm;

import fluent.asserts;

S2CellId makeCellId(string str) {
  return S2CellId.fromDebugString(str);
}

// Pad by at least twice the maximum error for reliable results.
private enum double kPadding = 2 * (FACE_CLIP_ERROR_UV_COORD + INTERSECTS_RECT_ERROR_UV_DIST);

S2Shape newPaddedCell(S2CellId id, double padding_uv) {
  int[2] ij;
  int orientation;
  int face = id.toFaceIJOrientation(ij[0], ij[1], orientation);
  R2Rect uv = S2CellId.IJLevelToBoundUV(ij, id.level()).expanded(padding_uv);
  S2Point[4] vertices;
  for (int i = 0; i < 4; ++i) {
    vertices[i] = FaceUVtoXYZ(face, uv.getVertex(i)).normalize();
  }
  return new S2LaxLoopShape(vertices);
}

@("S2ShapeIndexRegion.GetCapBound") unittest {
  auto id = S2CellId.fromDebugString("3/0123012301230123012301230123");

  // Add a polygon that is slightly smaller than the cell being tested.
  auto index = new MutableS2ShapeIndex();
  index.add(newPaddedCell(id, -kPadding));
  S2Cap cell_bound = new S2Cell(id).getCapBound();
  S2Cap index_bound = makeS2ShapeIndexRegion(index).getCapBound();
  Assert.equal(index_bound.contains(cell_bound), true);

  // Note that S2CellUnion::GetCapBound returns a slightly larger bound than
  // S2Cell::GetBound even when the cell union consists of a single S2CellId.
  Assert.notGreaterThan(index_bound.getRadius(), 1.00001 * cell_bound.getRadius());
}

@("S2ShapeIndexRegion.GetRectBound") unittest {
  auto id = S2CellId.fromDebugString("3/0123012301230123012301230123");

  // Add a polygon that is slightly smaller than the cell being tested.
  auto index = new MutableS2ShapeIndex();
  index.add(newPaddedCell(id, -kPadding));
  S2LatLngRect cell_bound = new S2Cell(id).getRectBound();
  S2LatLngRect index_bound = makeS2ShapeIndexRegion(index).getRectBound();
  Assert.equal(index_bound, cell_bound);
}

@("S2ShapeIndexRegion.GetCellUnionBoundMultipleFaces") unittest {
  S2CellId[] ids = [ makeCellId("3/00123"), makeCellId("2/11200013") ];
  auto index = new MutableS2ShapeIndex();
  foreach (id; ids) index.add(newPaddedCell(id, -kPadding));
  S2CellId[] covering;
  makeS2ShapeIndexRegion(index).getCellUnionBound(covering);
  sort(ids);
  Assert.equal(covering, ids);
}

@("S2ShapeIndexRegion.GetCellUnionBoundOneFace") unittest {
  // This tests consists of 3 pairs of S2CellIds.  Each pair is located within
  // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
  // We expect GetCellUnionBound to compute the smallest cell that bounds the
  // pair on each face.
  S2CellId[] input = [
    makeCellId("5/010"), makeCellId("5/0211030"),
    makeCellId("5/110230123"), makeCellId("5/11023021133"),
    makeCellId("5/311020003003030303"), makeCellId("5/311020023")
  ];
  S2CellId[] expected = [
    makeCellId("5/0"), makeCellId("5/110230"), makeCellId("5/3110200")
  ];
  auto index = new MutableS2ShapeIndex();
  foreach (id; input) {
    // Add each shape 3 times to ensure that the S2ShapeIndex subdivides.
    for (int copy = 0; copy < 3; ++copy) {
      index.add(newPaddedCell(id, -kPadding));
    }
  }
  S2CellId[] actual;
  makeS2ShapeIndexRegion(index).getCellUnionBound(actual);
  Assert.equal(actual, expected);
}

@("S2ShapeIndexRegion.ContainsCellMultipleShapes") unittest {
  auto id = S2CellId.fromDebugString("3/0123012301230123012301230123");

  // Add a polygon that is slightly smaller than the cell being tested.
  auto index = new MutableS2ShapeIndex();
  index.add(newPaddedCell(id, -kPadding));
  Assert.equal(makeS2ShapeIndexRegion(index).contains(new S2Cell(id)), false);

  // Add a second polygon that is slightly larger than the cell being tested.
  // Note that Contains() should return true if *any* shape contains the cell.
  index.add(newPaddedCell(id, kPadding));
  Assert.equal(makeS2ShapeIndexRegion(index).contains(new S2Cell(id)), true);

  // Verify that all children of the cell are also contained.
  for (S2CellId child = id.childBegin(); child != id.childEnd(); child = child.next()) {
    Assert.equal(makeS2ShapeIndexRegion(index).contains(new S2Cell(child)), true);
  }
}

@("S2ShapeIndexRegion.IntersectsShrunkenCell") unittest {
  auto target = S2CellId.fromDebugString("3/0123012301230123012301230123");

  // Add a polygon that is slightly smaller than the cell being tested.
  auto index = new MutableS2ShapeIndex();
  index.add(newPaddedCell(target, -kPadding));
  auto region = makeS2ShapeIndexRegion(index);

  // Check that the index intersects the cell itself, but not any of the
  // neighboring cells.
  Assert.equal(region.mayIntersect(new S2Cell(target)), true);
  S2CellId[] nbrs;
  target.appendAllNeighbors(target.level(), nbrs);
  foreach (S2CellId id; nbrs) {
    Assert.equal(region.mayIntersect(new S2Cell(id)), false);
  }
}

@("S2ShapeIndexRegion.IntersectsExactCell") unittest {
  auto target = S2CellId.fromDebugString("3/0123012301230123012301230123");

  // Adds a polygon that exactly follows a cell boundary.
  auto index = new MutableS2ShapeIndex();
  index.add(newPaddedCell(target, 0.0));
  auto region = makeS2ShapeIndexRegion(index);

  // Check that the index intersects the cell and all of its neighbors.
  S2CellId[] ids = [ target ];
  target.appendAllNeighbors(target.level(), ids);
  foreach (S2CellId id; ids) {
    Assert.equal(region.mayIntersect(new S2Cell(id)), true);
  }
}
