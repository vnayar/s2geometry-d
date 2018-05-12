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

module s2.s2cell_id_test;

import algorithm = std.algorithm;
import array = std.array;
import fluent.asserts;
import math = std.math;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell_id;
import s2.s2coords;
import s2.s2latlng;
import s2.s2point;
import s2.s2testing;
import s2.util.math.matrix3x3;
import s2coords = s2.s2coords;
import s2metrics = s2.s2metrics;

private S2CellId getCellId(double lat_degrees, double lng_degrees) {
  auto id = S2CellId(S2LatLng.fromDegrees(lat_degrees, lng_degrees));
  return id;
}

@("DefaultConstructor")
unittest {
  S2CellId id;
  Assert.equal(id.id(), 0);
  Assert.equal(id.isValid(), false);
}

@("S2CellIdHash")
unittest {
  Assert.equal(getCellId(0, 90).toHash(), getCellId(0, 90).toHash());
}

@("FaceDefinitions")
unittest {
  Assert.equal(getCellId(0, 0).face(), 0);
  Assert.equal(getCellId(0, 90).face(), 1);
  Assert.equal(getCellId(90, 0).face(), 2);
  Assert.equal(getCellId(0, 180).face(), 3);
  Assert.equal(getCellId(0, -90).face(), 4);
  Assert.equal(getCellId(-90, 0).face(), 5);
}

@("FromFace")
unittest {
  for (int face = 0; face < 6; ++face) {
    Assert.equal(S2CellId.fromFace(face), S2CellId.fromFacePosLevel(face, 0, 0));
  }
}

@("ParentChildRelationships")
unittest {
  S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, S2CellId.MAX_LEVEL - 4);
  Assert.equal(id.isValid(), true);
  Assert.equal(id.face(), 3);
  Assert.equal(id.pos(), 0x12345700);
  Assert.equal(id.level(), S2CellId.MAX_LEVEL - 4);
  Assert.equal(id.isLeaf(), false);

  Assert.equal(id.childBegin(id.level() + 2).pos(), 0x12345610);
  Assert.equal(id.childBegin().pos(), 0x12345640);
  Assert.equal(id.parent().pos(), 0x12345400);
  Assert.equal(id.parent(id.level() - 2).pos(), 0x12345000);

  // Check ordering of children relative to parents.
  Assert.lessThan(id.childBegin(), id);
  Assert.greaterThan(id.childEnd(), id);
  Assert.equal(id.childEnd(), id.childBegin().next().next().next().next());
  Assert.equal(id.rangeMin(), id.childBegin(S2CellId.MAX_LEVEL));
  Assert.equal(id.rangeMax().next(), id.childEnd(S2CellId.MAX_LEVEL));

  // Check that cells are represented by the position of their center
  // along the Hilbert curve.
  Assert.equal(2 * id.id(), id.rangeMin().id() + id.rangeMax().id());
}

@("SentinelRangeMinMax")
unittest {
  Assert.equal(S2CellId.sentinel(), S2CellId.sentinel().rangeMin());
  Assert.equal(S2CellId.sentinel(), S2CellId.sentinel().rangeMax());
}

@("CenterSiTi")
unittest {
  S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, S2CellId.MAX_LEVEL);
  // Check that the (si, ti) coordinates of the center end in a
  // 1 followed by (30 - level) 0s.
  int si, ti;

  // Leaf level, 30.
  id.getCenterSiTi(si, ti);
  Assert.equal(1 << 0, si & 1);
  Assert.equal(1 << 0, ti & 1);

  // Level 29.
  id.parent(S2CellId.MAX_LEVEL - 1).getCenterSiTi(si, ti);
  Assert.equal(1 << 1, si & 3);
  Assert.equal(1 << 1, ti & 3);

  // Level 28.
  id.parent(S2CellId.MAX_LEVEL - 2).getCenterSiTi(si, ti);
  Assert.equal(1 << 2, si & 7);
  Assert.equal(1 << 2, ti & 7);

  // Level 20.
  id.parent(S2CellId.MAX_LEVEL - 10).getCenterSiTi(si, ti);
  Assert.equal(1 << 10, si & ((1 << 11) - 1));
  Assert.equal(1 << 10, ti & ((1 << 11) - 1));

  // Level 10.
  id.parent(S2CellId.MAX_LEVEL - 20).getCenterSiTi(si, ti);
  Assert.equal(1 << 20, si & ((1 << 21) - 1));
  Assert.equal(1 << 20, ti & ((1 << 21) - 1));

  // Level 0.
  id.parent(0).getCenterSiTi(si, ti);
  Assert.equal(1 << 30, si & ((1U << 31) - 1));
  Assert.equal(1 << 30, ti & ((1U << 31) - 1));
}

@("Wrapping")
unittest {
  // Check wrapping from beginning of Hilbert curve to end and vice versa.
  Assert.equal(S2CellId.end(0).prev(), S2CellId.begin(0).prevWrap());

  Assert.equal(
      S2CellId.fromFacePosLevel(5, ~0uL >> S2CellId.FACE_BITS, S2CellId.MAX_LEVEL),
      S2CellId.begin(S2CellId.MAX_LEVEL).prevWrap());
  Assert.equal(
      S2CellId.fromFacePosLevel(5, ~0uL >> S2CellId.FACE_BITS, S2CellId.MAX_LEVEL),
      S2CellId.begin(S2CellId.MAX_LEVEL).advanceWrap(-1));

  Assert.equal(S2CellId.begin(4), S2CellId.end(4).prev().nextWrap());
  Assert.equal(S2CellId.begin(4), S2CellId.end(4).advance(-1).advanceWrap(1));

  Assert.equal(
      S2CellId.fromFacePosLevel(0, 0, S2CellId.MAX_LEVEL),
      S2CellId.end(S2CellId.MAX_LEVEL).prev().nextWrap());
  Assert.equal(
      S2CellId.fromFacePosLevel(0, 0, S2CellId.MAX_LEVEL),
      S2CellId.end(S2CellId.MAX_LEVEL).advance(-1).advanceWrap(1));
}

@("Advance")
unittest {
  S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, S2CellId.MAX_LEVEL - 4);

  // Check basic properties of advance().
  Assert.equal(S2CellId.end(0), S2CellId.begin(0).advance(7));
  Assert.equal(S2CellId.end(0), S2CellId.begin(0).advance(12));
  Assert.equal(S2CellId.begin(0), S2CellId.end(0).advance(-7));
  Assert.equal(S2CellId.begin(0), S2CellId.end(0).advance(-12000000));

  int num_level_5_cells = 6 << (2 * 5);
  Assert.equal(
      S2CellId.end(5).advance(500 - num_level_5_cells), S2CellId.begin(5).advance(500));
  Assert.equal(
      id.next().childBegin(S2CellId.MAX_LEVEL),
      id.childBegin(S2CellId.MAX_LEVEL).advance(256));
  Assert.equal(
      S2CellId.fromFacePosLevel(5, 0, S2CellId.MAX_LEVEL),
      S2CellId.fromFacePosLevel(1, 0, S2CellId.MAX_LEVEL)
          .advance(4L << (2 * S2CellId.MAX_LEVEL)));

  // Check basic properties of advance_wrap().
  Assert.equal(S2CellId.fromFace(1), S2CellId.begin(0).advanceWrap(7));
  Assert.equal(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(12));
  Assert.equal(S2CellId.fromFace(4), S2CellId.fromFace(5).advanceWrap(-7));
  Assert.equal(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(-12000000));
  Assert.equal(S2CellId.begin(5).advanceWrap(6644), S2CellId.begin(5).advanceWrap(-11788));
  Assert.equal(
      id.next().childBegin(S2CellId.MAX_LEVEL),
      id.childBegin(S2CellId.MAX_LEVEL).advanceWrap(256));
  Assert.equal(
      S2CellId.fromFacePosLevel(1, 0, S2CellId.MAX_LEVEL),
      S2CellId.fromFacePosLevel(5, 0, S2CellId.MAX_LEVEL)
          .advanceWrap(2L << (2 * S2CellId.MAX_LEVEL)));
}

@("DistanceFromBegin")
unittest {
  Assert.equal(6L, S2CellId.end(0).distanceFromBegin());
  Assert.equal(
      6 * (1L << (2 * S2CellId.MAX_LEVEL)),
      S2CellId.end(S2CellId.MAX_LEVEL).distanceFromBegin());

  Assert.equal(0L, S2CellId.begin(0).distanceFromBegin());
  Assert.equal(0L, S2CellId.begin(S2CellId.MAX_LEVEL).distanceFromBegin());

  S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, S2CellId.MAX_LEVEL - 4);
  Assert.equal(id, S2CellId.begin(id.level()).advance(id.distanceFromBegin()));
}

@("MaximumTile")
unittest {
  // This method is tested more thoroughly in s2cell_union_test.cc.
  for (int iter = 0; iter < 1000; ++iter) {
    S2CellId id = S2Testing.getRandomCellId(10);

    // Check that "limit" is returned for tiles at or beyond "limit".
    Assert.equal(id, id.maximumTile(id));
    Assert.equal(id, id.child(0).maximumTile(id));
    Assert.equal(id, id.child(1).maximumTile(id));
    Assert.equal(id, id.next().maximumTile(id));
    Assert.equal(id.child(0), id.maximumTile(id.child(0)));

    // Check that the tile size is increased when possible.
    Assert.equal(id, id.child(0).maximumTile(id.next()));
    Assert.equal(id, id.child(0).maximumTile(id.next().child(0)));
    Assert.equal(id, id.child(0).maximumTile(id.next().child(1).child(0)));
    Assert.equal(id, id.child(0).child(0).maximumTile(id.next()));
    Assert.equal(id, id.child(0).child(0).child(0).maximumTile(id.next()));

    // Check that the tile size is decreased when necessary.
    Assert.equal(id.child(0), id.maximumTile(id.child(0).next()));
    Assert.equal(id.child(0), id.maximumTile(id.child(0).next().child(0)));
    Assert.equal(id.child(0), id.maximumTile(id.child(0).next().child(1)));
    Assert.equal(id.child(0).child(0), id.maximumTile(id.child(0).child(0).next()));
    Assert.equal(
        id.child(0).child(0).child(0),
        id.maximumTile(id.child(0).child(0).child(0).next()));

    // Check that the tile size is otherwise unchanged.
    Assert.equal(id, id.maximumTile(id.next()));
    Assert.equal(id, id.maximumTile(id.next().child(0)));
    Assert.equal(id, id.maximumTile(id.next().child(1).child(0)));
  }
}

@("GetCommonAncestorLevel")
unittest {
  // Two identical cell ids.
  Assert.equal(0, S2CellId.fromFace(0).getCommonAncestorLevel(S2CellId.fromFace(0)));
  Assert.equal(30, S2CellId.fromFace(0).childBegin(30)
      .getCommonAncestorLevel(S2CellId.fromFace(0).childBegin(30)));

  // One cell id is a descendant of the other.
  Assert.equal(0, S2CellId.fromFace(0).childBegin(30)
      .getCommonAncestorLevel(S2CellId.fromFace(0)));
  Assert.equal(0, S2CellId.fromFace(5)
      .getCommonAncestorLevel(S2CellId.fromFace(5).childEnd(30).prev()));

  // Two cells that have no common ancestor.
  Assert.equal(-1, S2CellId.fromFace(0)
      .getCommonAncestorLevel(S2CellId.fromFace(5)));
  Assert.equal(-1, S2CellId.fromFace(2).childBegin(30)
      .getCommonAncestorLevel(S2CellId.fromFace(3).childEnd(20)));

  // Two cells that have a common ancestor distinct from both of them.
  Assert.equal(8, S2CellId.fromFace(5).childBegin(9).next().childBegin(15)
      .getCommonAncestorLevel(S2CellId.fromFace(5).childBegin(9).childBegin(20)));
  Assert.equal(1, S2CellId.fromFace(0).childBegin(2).childBegin(30)
      .getCommonAncestorLevel(S2CellId.fromFace(0).childBegin(2).next().childBegin(5)));
}

@("Inverses")
unittest {
  // Check the conversion of random leaf cells to S2LatLngs and back.
  foreach (int i; 0 .. 200_000) {
    S2CellId id = S2Testing.getRandomCellId(S2CellId.MAX_LEVEL);
    Assert.equal(id.isLeaf(), true);
    Assert.equal(S2CellId.MAX_LEVEL, id.level());
    S2LatLng center = id.toLatLng();
    Assert.equal(id.id(), S2CellId(center).id());
  }
}

@("Tokens")
unittest {
  // Test random cell ids at all levels.
  foreach (int i; 0 .. 10_000) {
    S2CellId id = S2Testing.getRandomCellId();
    string token = id.toToken();
    Assert.notGreaterThan(token.length, 16);
    Assert.equal(id, S2CellId.fromToken(token));
  }
  // Check that invalid cell ids can be encoded, and round-trip is
  // the identity operation.
  string token = S2CellId.none().toToken();
  Assert.equal(S2CellId.none(), S2CellId.fromToken(token));

  // Sentinel is invalid.
  token = S2CellId.sentinel().toToken();
  Assert.equal(S2CellId.fromToken(token), S2CellId.sentinel());

  // Check an invalid face.
  token = S2CellId.fromFace(7).toToken();
  Assert.equal(S2CellId.fromToken(token), S2CellId.fromFace(7));

  // Check that supplying tokens with non-alphanumeric characters
  // returns S2CellId::None().
  Assert.equal(S2CellId.none(), S2CellId.fromToken("876b e99"));
  Assert.equal(S2CellId.none(), S2CellId.fromToken("876bee99\n"));
  Assert.equal(S2CellId.none(), S2CellId.fromToken("876[ee99"));
  Assert.equal(S2CellId.none(), S2CellId.fromToken(" 876bee99"));
}

/+
TEST(S2CellId, EncodeDecode) {
  S2CellId id(0x7837423);
  Encoder encoder;
  id.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2CellId decoded_id;
  EXPECT_TRUE(decoded_id.Decode(&decoder));
  EXPECT_EQ(id, decoded_id);
}

TEST(S2CellId, EncodeDecodeNoneCell) {
  S2CellId none_id = S2CellId::None();
  Encoder encoder;
  none_id.Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2CellId decoded_id;
  EXPECT_TRUE(decoded_id.Decode(&decoder));
  EXPECT_EQ(none_id, decoded_id);
}

TEST(S2CellId, DecodeFailsWithTruncatedBuffer) {
  S2CellId id(0x7837423);
  Encoder encoder;
  id.Encode(&encoder);

  // Truncate encoded buffer.
  Decoder decoder(encoder.base(), encoder.length() - 2);
  S2CellId decoded_id;
  EXPECT_FALSE(decoded_id.Decode(&decoder));
}
+/


private immutable int MAX_EXPAND_LEVEL = 3;

private void expandCell(in S2CellId parent, S2CellId[] cells, S2CellId[S2CellId] parent_map) {
  cells ~= parent;
  if (parent.level() == MAX_EXPAND_LEVEL) return;
  int i, j, orientation;
  int face = parent.toFaceIJOrientation(i, j, orientation);
  Assert.equal(parent.face(), face);

  S2CellId child = parent.childBegin();
  for (int pos = 0; child != parent.childEnd(); child = child.next(), ++pos) {
    parent_map[child] = parent;
    // Do some basic checks on the children.
    Assert.equal(child, parent.child(pos));
    Assert.equal(pos, child.childPosition());
    // Test child_position(level) on all the child's ancestors.
    for (S2CellId ancestor = child; ancestor.level() >= 1; ancestor = parent_map[ancestor]) {
      Assert.equal(child.childPosition(ancestor.level()), ancestor.childPosition());
    }
    Assert.equal(pos, child.childPosition(child.level()));
    Assert.equal(parent.level() + 1, child.level());
    Assert.equal(child.isLeaf(), false);
    int child_orientation;
    Assert.equal(face, child.toFaceIJOrientation(i, j, child_orientation));
    Assert.equal(orientation ^ s2coords.POS_TO_ORIENTATION[pos], child_orientation);
    expandCell(child, cells, parent_map);
  }
}

@("Containment")
unittest {
  // Test contains() and intersects().
  S2CellId[S2CellId] parent_map;
  S2CellId[] cells;
  for (int face = 0; face < 6; ++face) {
    expandCell(S2CellId.fromFace(face), cells, parent_map);
  }
  foreach (S2CellId end_id; cells) {
    foreach (S2CellId begin_id; cells) {
      bool contained = true;
      for (S2CellId id = begin_id; id != end_id; id = parent_map[id]) {
        if (id !in parent_map) {
          contained = false;
          break;
        }
      }
      Assert.equal(contained, end_id.contains(begin_id));
      Assert.equal(contained,
                begin_id >= end_id.rangeMin() &&
                begin_id <= end_id.rangeMax());
      Assert.equal(end_id.intersects(begin_id),
                end_id.contains(begin_id) || begin_id.contains(end_id));
    }
  }
}

immutable int MAX_WALK_LEVEL = 8;

@("Continuity")
unittest {
  // Make sure that sequentially increasing cell ids form a continuous
  // path over the surface of the sphere, i.e. there are no
  // discontinuous jumps from one region to another.

  double max_dist = s2metrics.MAX_EDGE.getValue(MAX_WALK_LEVEL);
  S2CellId end = S2CellId.end(MAX_WALK_LEVEL);
  S2CellId id = S2CellId.begin(MAX_WALK_LEVEL);
  for (; id != end; id = id.next()) {
    Assert.notGreaterThan(id.toS2PointRaw().angle(id.nextWrap().toS2PointRaw()), max_dist);
    Assert.equal(id.nextWrap(), id.advanceWrap(1));
    Assert.equal(id, id.nextWrap().advanceWrap(-1));

    // Check that the ToPointRaw() returns the center of each cell
    // in (s,t) coordinates.
    double u, v;
    s2coords.XYZtoFaceUV(id.toS2PointRaw(), u, v);
    immutable double kCellSize = 1.0 / (1 << MAX_WALK_LEVEL);
    Assert.approximately(math.remainder(s2coords.UVtoST(u), 0.5 * kCellSize), 0.0, 1e-15);
    Assert.approximately(math.remainder(s2coords.UVtoST(v), 0.5 * kCellSize), 0.0, 1e-15);
  }
}

@("Coverage")
unittest {
  // Make sure that random points on the sphere can be represented to the
  // expected level of accuracy, which in the worst case is sqrt(2/3) times
  // the maximum arc length between the points on the sphere associated with
  // adjacent values of "i" or "j".  (It is sqrt(2/3) rather than 1/2 because
  // the cells at the corners of each face are stretched -- they have 60 and
  // 120 degree angles.)

  double max_dist = 0.5 * s2metrics.MAX_DIAG.getValue(S2CellId.MAX_LEVEL);
  for (int i = 0; i < 1_000_000; ++i) {
    S2Point p = S2Testing.randomPoint();
    S2Point q = S2CellId(p).toS2PointRaw();
    Assert.notGreaterThan(p.angle(q), max_dist);
  }
}

private void testAllNeighbors(S2CellId id, int level)
in {
  assert(level >= id.level());
  assert(level < S2CellId.MAX_LEVEL);
} body {

  // We compute AppendAllNeighbors, and then add in all the children of "id"
  // at the given level.  We then compare this against the result of finding
  // all the vertex neighbors of all the vertices of children of "id" at the
  // given level.  These should give the same result.
  S2CellId[] all;
  S2CellId[] expected;
  id.appendAllNeighbors(level, all);
  S2CellId end = id.childEnd(level + 1);
  for (S2CellId c = id.childBegin(level + 1); c != end; c = c.next()) {
    all ~= c.parent();
    c.appendVertexNeighbors(level, expected);
  }
  // Sort the results and eliminate duplicates.
  all = array.array(algorithm.uniq(algorithm.sort(all)));
  expected = array.array(algorithm.uniq(algorithm.sort(expected)));
  Assert.equal(expected, all);
}

@("Neighbors")
unittest {
  // Check the edge neighbors of face 1.
  immutable int[] out_faces = [ 5, 3, 2, 0 ];
  S2CellId[4] face_nbrs;
  S2CellId.fromFace(1).getEdgeNeighbors(face_nbrs);
  for (int i = 0; i < 4; ++i) {
    Assert.equal(face_nbrs[i].isFace(), true);
    Assert.equal(out_faces[i], face_nbrs[i].face());
  }

  // Check the edge neighbors of the corner cells at all levels.  This case is
  // trickier because it requires projecting onto adjacent faces.
  immutable int kMaxIJ = S2CellId.MAX_SIZE - 1;
  for (int level = 1; level <= S2CellId.MAX_LEVEL; ++level) {
    S2CellId id = S2CellId.fromFaceIJ(1, 0, 0).parent(level);
    S2CellId[4] nbrs;
    id.getEdgeNeighbors(nbrs);
    // These neighbors were determined manually using the face and axis
    // relationships defined in s2coords.cc.
    int size_ij = S2CellId.getSizeIJ(level);
    Assert.equal(S2CellId.fromFaceIJ(5, kMaxIJ, kMaxIJ).parent(level), nbrs[0]);
    Assert.equal(S2CellId.fromFaceIJ(1, size_ij, 0).parent(level), nbrs[1]);
    Assert.equal(S2CellId.fromFaceIJ(1, 0, size_ij).parent(level), nbrs[2]);
    Assert.equal(S2CellId.fromFaceIJ(0, kMaxIJ, 0).parent(level), nbrs[3]);
  }

  // Check the vertex neighbors of the center of face 2 at level 5.
  S2CellId[] nbrs;
  S2CellId(S2Point(0, 0, 1)).appendVertexNeighbors(5, nbrs);
  algorithm.sort(nbrs);
  for (int i = 0; i < 4; ++i) {
    Assert.equal(
        S2CellId.fromFaceIJ(2, (1 << 29) - (i < 2), (1 << 29) - (i == 0 || i == 3)).parent(5),
        nbrs[i]);
  }
  nbrs.length = 0;

  // Check the vertex neighbors of the corner of faces 0, 4, and 5.
  S2CellId id = S2CellId.fromFacePosLevel(0, 0, S2CellId.MAX_LEVEL);
  id.appendVertexNeighbors(0, nbrs);
  algorithm.sort(nbrs);
  Assert.equal(3uL, nbrs.length);
  Assert.equal(S2CellId.fromFace(0), nbrs[0]);
  Assert.equal(S2CellId.fromFace(4), nbrs[1]);
  Assert.equal(S2CellId.fromFace(5), nbrs[2]);

  // Check that AppendAllNeighbors produces results that are consistent
  // with AppendVertexNeighbors for a bunch of random cells.
  for (int i = 0; i < 1000; ++i) {
    id = S2Testing.getRandomCellId();
    if (id.isLeaf()) id = id.parent();

    // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell ids,
    // so it's not reasonable to use large values of "diff".
    int max_diff = algorithm.min(6, S2CellId.MAX_LEVEL - id.level() - 1);
    int level = id.level() + S2Testing.rnd.uniform(max_diff);
    testAllNeighbors(id, level);
  }
}

// Returns a random point on the boundary of the given rectangle.
private R2Point sampleBoundary(in R2Rect rect) {
  R2Point uv;
  int d = S2Testing.rnd.uniform(2);
  uv[d] = S2Testing.rnd.uniformDouble(rect[d][0], rect[d][1]);
  uv[1-d] = S2Testing.rnd.oneIn(2) ? rect[1-d][0] : rect[1-d][1];
  return uv;
}

// Returns the closest point to "uv" on the boundary of "rect".
private R2Point projectToBoundary(in R2Point uv, in R2Rect rect) {
  double du0 = math.abs(uv[0] - rect[0][0]);
  double du1 = math.abs(uv[0] - rect[0][1]);
  double dv0 = math.abs(uv[1] - rect[1][0]);
  double dv1 = math.abs(uv[1] - rect[1][1]);
  double dmin = algorithm.min(algorithm.min(du0, du1), algorithm.min(dv0, dv1));
  if (du0 == dmin) return R2Point(rect[0][0], rect[1].project(uv[1]));
  if (du1 == dmin) return R2Point(rect[0][1], rect[1].project(uv[1]));
  if (dv0 == dmin) return R2Point(rect[0].project(uv[0]), rect[1][0]);
  Assert.equal(dmin, dv1, "Bug in ProjectToBoundary");
  return R2Point(rect[0].project(uv[0]), rect[1][1]);
}

private void testExpandedByDistanceUV(in S2CellId id, in S1Angle distance) {
  R2Rect bound = id.getBoundUV();
  R2Rect expanded = S2CellId.expandedByDistanceUV(bound, distance);
  for (int iter = 0; iter < 100; ++iter) {
    // Choose a point on the boundary of the rectangle.
    int face = S2Testing.rnd.uniform(6);
    R2Point center_uv = sampleBoundary(bound);
    S2Point center = s2coords.FaceUVtoXYZ(face, center_uv).normalize();

    // Now sample a point from a disc of radius (2 * distance).
    S2Point p = S2Testing.samplePoint(new S2Cap(center, 2 * distance.abs()));

    // Find the closest point on the boundary to the sampled point.
    R2Point uv;
    if (!s2coords.FaceXYZtoUV(face, p, uv)) {
      continue;
    }
    R2Point closest_uv = projectToBoundary(uv, bound);
    S2Point closest = s2coords.FaceUVtoXYZ(face, closest_uv).normalize();
    S1Angle actual_dist = S1Angle(p, closest);

    if (distance >= S1Angle.zero()) {
      // "expanded" should contain all points in the original bound, and also
      // all points within "distance" of the boundary.
      if (bound.contains(uv) || actual_dist < distance) {
        Assert.equal(expanded.contains(uv), true);
      }
    } else {
      // "expanded" should not contain any points within "distance" of the
      // original boundary.
      if (actual_dist < -distance) {
        Assert.equal(expanded.contains(uv), false);
      }
    }
  }
}

@("ExpandedByDistanceUV")
unittest {
  double max_dist_degrees = 10;
  for (int iter = 0; iter < 100; ++iter) {
    S2CellId id = S2Testing.getRandomCellId();
    double dist_degrees = S2Testing.rnd.uniformDouble(-max_dist_degrees, max_dist_degrees);
    testExpandedByDistanceUV(id, S1Angle.fromDegrees(dist_degrees));
  }
}

@("toString")
unittest {
  Assert.equal("3/", S2CellId.fromFace(3).toString());
  Assert.equal("4/000000000000000000000000000000", S2CellId.fromFace(4).rangeMin().toString());
  Assert.equal("Invalid: 0000000000000000", S2CellId.none().toString());
}

@("fromDebugString")
unittest {
  Assert.equal(S2CellId.fromFace(3), S2CellId.fromDebugString("3/"));
  Assert.equal(S2CellId.fromFace(0).child(2).child(1),
      S2CellId.fromDebugString("0/21"));
  Assert.equal(S2CellId.fromFace(4).rangeMin(),
      S2CellId.fromDebugString("4/000000000000000000000000000000"));
  Assert.equal(S2CellId.none(),
      S2CellId.fromDebugString("4/0000000000000000000000000000000"));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString(""));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString("7/"));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString(" /"));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString("3:0"));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString("3/ 12"));
  Assert.equal(S2CellId.none(), S2CellId.fromDebugString("3/1241"));
}

/+
////
// TODO: A method to measure CPU usage is needed before converting this.
////

TEST(S2CellId, ToPointBenchmark) {
  // This "test" is really a benchmark, so skip it unless we're optimized.
  if (google::DEBUG_MODE) return;

  // Test speed of conversions from points to leaf cells.
  double control_start = S2Testing::GetCpuTime();
  S2CellId begin = S2CellId::Begin(S2CellId::kMaxLevel);
  S2CellId end = S2CellId::End(S2CellId::kMaxLevel);
  uint64 delta = (end.id() - begin.id()) / FLAGS_iters;
  delta &= ~static_cast<uint64>(1);  // Make sure all ids are leaf cells.

  S2CellId id = begin;
  double sum = 0;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += static_cast<double>(id.id());
    id = S2CellId(id.id() + delta);
  }
  double control_time = S2Testing::GetCpuTime() - control_start;
  printf("\tControl:    %8.3f usecs\n", 1e6 * control_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.

  double test_start = S2Testing::GetCpuTime();
  sum = 0;
  id = begin;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += id.ToPointRaw()[0];
    id = S2CellId(id.id() + delta);
  }
  double test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tToPointRaw: %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.

  test_start = S2Testing::GetCpuTime();
  sum = 0;
  id = begin;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += id.ToPoint()[0];
    id = S2CellId(id.id() + delta);
  }
  test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tToPoint:    %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.
}

TEST(S2CellId, FromPointBenchmark) {
  // This "test" is really a benchmark, so skip it unless we're optimized.
  if (google::DEBUG_MODE) return;

  // The sample points follow a spiral curve that completes one revolution
  // around the z-axis every 1/dt samples.  The z-coordinate increases
  // from -4 to +4 over FLAGS_iters samples.

  S2Point start(1, 0, -4);
  double dz = (-2 * start.z()) / FLAGS_iters;
  double dt = 1.37482937133e-4;

  // Test speed of conversions from leaf cells to points.
  double control_start = S2Testing::GetCpuTime();
  uint64 isum = 0;
  S2Point p = start;
  for (int i = FLAGS_iters; i > 0; --i) {
    // Cheap rotation around the z-axis (spirals inward slightly
    // each revolution).
    p += S2Point(-dt * p.y(), dt * p.x(), dz);
    isum += MathUtil::FastIntRound(p[0] + p[1] + p[2]);
  }
  double control_time = S2Testing::GetCpuTime() - control_start;
  printf("\tControl:    %8.3f usecs\n", 1e6 * control_time / FLAGS_iters);
  EXPECT_NE(isum, 0);  // Don't let the loop get optimized away.

  double test_start = S2Testing::GetCpuTime();
  isum = 0;
  p = start;
  for (int i = FLAGS_iters; i > 0; --i) {
    p += S2Point(-dt * p.y(), dt * p.x(), dz);
    isum += S2CellId(p).id();
  }
  double test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tFromPoint:  %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(isum, 0);  // Don't let the loop get optimized away.
}
+/
