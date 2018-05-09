module s2.s2coords_test;

// Original author: ericv@google.com (Eric Veach)

import fluent.asserts;
import s2.s2coords;
import s2.s2cellid;
import s2.s2testing;
import s2.s2point;
import math = std.math;

private int swapAxes(int ij) {
  return ((ij >> 1) & 1) + ((ij & 1) << 1);
}

private int invertBits(int ij) {
  return ij ^ 3;
}

enum DOUBLE_ERR = 0.0001;

@("S2.TraversalOrder")
unittest {
  foreach (r; 0 .. 4) {
    foreach (i; 0 .. 4) {
      // Check consistency with respect to swapping axes.
      Assert.equal(IJ_TO_POS[r][i], IJ_TO_POS[r ^ SWAP_MASK][swapAxes(i)]);
      Assert.equal(POS_TO_IJ[r][i], swapAxes(POS_TO_IJ[r ^ SWAP_MASK][i]));

      // Check consistency with respect to reversing axis directions.
      Assert.equal(IJ_TO_POS[r][i], IJ_TO_POS[r ^ INVERT_MASK][invertBits(i)]);
      Assert.equal(POS_TO_IJ[r][i], invertBits(POS_TO_IJ[r ^ INVERT_MASK][i]));

      // Check that the two tables are inverses of each other.
      Assert.equal(IJ_TO_POS[r][POS_TO_IJ[r][i]], i);
      Assert.equal(POS_TO_IJ[r][IJ_TO_POS[r][i]], i);
    }
  }
}

@("S2.ST_UV_Converstions")
unittest {
  // Check boundary conditions.
  for (double s = 0; s <= 1; s += 0.5) {
    double u = STtoUV(s);
    Assert.equal(u, 2 * s - 1);
  }
  for (double u = -1; u <= 1; ++u) {
    double s = UVtoST(u);
    Assert.equal(s, 0.5 * (u + 1));
  }
  // Check that UVtoST and STtoUV are inverses.
  for (double x = 0; x <= 1; x += 0.0001) {
    Assert.approximately(UVtoST(STtoUV(x)), x, 1e-15);
    Assert.approximately(STtoUV(UVtoST(2*x-1)), 2*x-1, 1e-15);
  }
}

@("S2.FaceUVtoXYZ")
unittest {
  // Check that each face appears exactly once.
  S2Point sum;
  foreach (face; 0 .. 6) {
    S2Point center = FaceUVtoXYZ(face, 0, 0);
    Assert.equal(GetNorm(face), center);
    Assert.equal(math.fabs(center[center.largestAbsComponent()]), 1);
    sum += center.fabs();
  }
  Assert.equal(sum, S2Point(2, 2, 2));

  // Check that each face has a right-handed coordinate system.
  foreach (face; 0 .. 6) {
    Assert.equal(GetUAxis(face).crossProd(GetVAxis(face)).dotProd(FaceUVtoXYZ(face, 0, 0)), 1);
  }

  // Check that the Hilbert curves on each face combine to form a
  // continuous curve over the entire cube.
  foreach (face; 0 .. 6) {
    // The Hilbert curve on each face starts at (-1,-1) and terminates
    // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
    int sign = (face & SWAP_MASK) ? -1 : 1;
    Assert.equal(FaceUVtoXYZ(face, sign, -sign), FaceUVtoXYZ((face + 1) % 6, -1, -1));
  }
}

@("S2.FaceXYZtoUVW")
unittest {
  foreach (face; 0 .. 6) {
    Assert.equal(S2Point( 0,  0,  0), FaceXYZtoUVW(face, S2Point(0, 0, 0)));
    Assert.equal(S2Point( 1,  0,  0), FaceXYZtoUVW(face,  GetUAxis(face)));
    Assert.equal(S2Point(-1,  0,  0), FaceXYZtoUVW(face, -GetUAxis(face)));
    Assert.equal(S2Point( 0,  1,  0), FaceXYZtoUVW(face,  GetVAxis(face)));
    Assert.equal(S2Point( 0, -1,  0), FaceXYZtoUVW(face, -GetVAxis(face)));
    Assert.equal(S2Point( 0,  0,  1), FaceXYZtoUVW(face,  GetNorm(face)));
    Assert.equal(S2Point( 0,  0, -1), FaceXYZtoUVW(face, -GetNorm(face)));
  }
}

@("S2.XYZToFaceSiTi")
unittest {
  // Check the conversion of random cells to center points and back.
  for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
    for (int i = 0; i < 1000; ++i) {
      S2CellId id = S2Testing.getRandomCellId(level);
      int face;
      uint si, ti;
      int actual_level = XYZtoFaceSiTi(id.toPoint(), face, si, ti);
      Assert.equal(level, actual_level);
      S2CellId actual_id = S2CellId.fromFaceIJ(face, si / 2, ti / 2).parent(level);
      Assert.equal(id, actual_id);

      // Now test a point near the cell center but not equal to it.
      S2Point p_moved = id.toPoint() + S2Point(1e-13, 1e-13, 1e-13);
      int face_moved;
      uint si_moved, ti_moved;
      actual_level = XYZtoFaceSiTi(p_moved, face_moved, si_moved, ti_moved);
      Assert.equal(actual_level, -1);
      Assert.equal(face_moved, face);
      Assert.equal(si_moved, si);
      Assert.equal(ti_moved, ti);

      // Finally, test some random (si,ti) values that may be at different
      // levels, or not at a valid level at all (for example, si == 0).
      int face_random = S2Testing.rnd.uniform(S2CellId.NUM_FACES);
      uint si_random, ti_random;
      uint mask = -1 << (S2CellId.MAX_LEVEL - level);
      do {
        si_random = S2Testing.rnd.rand32() & mask;
        ti_random = S2Testing.rnd.rand32() & mask;
      } while (si_random > MAX_SI_TI || ti_random > MAX_SI_TI);
      S2Point p_random = FaceSiTitoXYZ(face_random, si_random, ti_random);
      actual_level = XYZtoFaceSiTi(p_random, face, si, ti);
      if (face != face_random) {
        // The chosen point is on the edge of a top-level face cell.
        Assert.equal(actual_level, -1);
        Assert.equal(si == 0 || si == MAX_SI_TI || ti == 0 || ti == MAX_SI_TI, true);
      } else {
        Assert.equal(si_random, si);
        Assert.equal(ti_random, ti);
        if (actual_level >= 0) {
          Assert.equal(
              p_random,
              S2CellId.fromFaceIJ(face, si / 2, ti / 2).parent(actual_level).toPoint());
        }
      }
    }
  }
}

@("S2.UVNorms")
unittest {
  // Check that GetUNorm and GetVNorm compute right-handed normals for
  // an edge in the increasing U or V direction.
  foreach (face; 0 .. 6) {
    for (double x = -1; x <= 1; x += 1/1024.0) {
      Assert.approximately(
          FaceUVtoXYZ(face, x, -1)
              .crossProd(FaceUVtoXYZ(face, x, 1))
              .angle(GetUNorm(face, x)),
          0, DOUBLE_ERR);
      Assert.approximately(
          FaceUVtoXYZ(face, -1, x)
              .crossProd(FaceUVtoXYZ(face, 1, x))
              .angle(GetVNorm(face, x)),
          0, DOUBLE_ERR);
    }
  }
}

@("S2.UVWAxis")
unittest {
  foreach (face; 0 .. 6) {
    // Check that axes are consistent with FaceUVtoXYZ.
    Assert.equal(FaceUVtoXYZ(face, 1, 0) - FaceUVtoXYZ(face, 0, 0), GetUAxis(face));
    Assert.equal(FaceUVtoXYZ(face, 0, 1) - FaceUVtoXYZ(face, 0, 0), GetVAxis(face));
    Assert.equal(FaceUVtoXYZ(face, 0, 0), GetNorm(face));

    // Check that every face coordinate frame is right-handed.
    Assert.equal(GetUAxis(face).crossProd(GetVAxis(face)).dotProd(GetNorm(face)), 1);

    // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
    Assert.equal(GetUAxis(face), GetUVWAxis(face, 0));
    Assert.equal(GetVAxis(face), GetUVWAxis(face, 1));
    Assert.equal(GetNorm(face), GetUVWAxis(face, 2));
  }
}

@("S2.UVWFace")
unittest {
  // Check that GetUVWFace is consistent with GetUVWAxis.
  foreach (face; 0 .. 6) {
    foreach (axis; 0 .. 3) {
      Assert.equal(GetFace(-GetUVWAxis(face, axis)), GetUVWFace(face, axis, 0));
      Assert.equal(GetFace(GetUVWAxis(face, axis)), GetUVWFace(face, axis, 1));
    }
  }
}
