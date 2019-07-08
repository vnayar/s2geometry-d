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

module s2.s2cap_test;

import fluent.asserts;
import math = std.math;
import s2.r1interval;
import s2.s1angle;
import s2.s1chord_angle;
import s2.s1interval;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2coords;
import s2.s2latlng;
import s2.s2latlng_rect;
import s2.s2point;
import s2.s2testing;
import s2.util.coding.coder;
import s2.util.math.vector;
import s2metrics = s2.s2metrics;

static S2Point getLatLngPoint(double lat_degrees, double lng_degrees) {
  return S2LatLng.fromDegrees(lat_degrees, lng_degrees).toS2Point();
}

enum {
  // About 9 times the double-precision roundoff relative error.
  double EPSILON = 1e-15,
  double DOUBLE_ERR = 0.0001
}

@("S2Cap.Basic")
unittest {
  // Test basic properties of empty and full caps.
  S2Cap empty = S2Cap.empty();
  S2Cap full = S2Cap.full();
  Assert.equal(empty.isValid(), true);
  Assert.equal(empty.isEmpty(), true);
  Assert.equal(empty.complement().isFull(), true);
  Assert.equal(full.isValid(), true);
  Assert.equal(full.isFull(), true);
  Assert.equal(full.complement().isEmpty(), true);
  Assert.equal(full.height(), 2);
  Assert.approximately(full.getRadius().degrees(), 180.0, DOUBLE_ERR);

  // Test the S1Angle constructor using out-of-range arguments.
  Assert.equal(new S2Cap(S2Point(1, 0, 0), S1Angle.fromRadians(-20)).isEmpty(), true);
  Assert.equal(new S2Cap(S2Point(1, 0, 0), S1Angle.fromRadians(5)).isFull(), true);
  Assert.equal(new S2Cap(S2Point(1, 0, 0), S1Angle.infinity()).isFull(), true);

  // Check that the default S2Cap is identical to Empty().
  S2Cap default_empty = new S2Cap();
  Assert.equal(default_empty.isValid(), true);
  Assert.equal(default_empty.isEmpty(), true);
  Assert.equal(empty.center(), default_empty.center());
  Assert.equal(empty.height(), default_empty.height());

  // Containment and intersection of empty and full caps.
  Assert.equal(empty.contains(empty), true);
  Assert.equal(full.contains(empty), true);
  Assert.equal(full.contains(full), true);
  Assert.equal(empty.interiorIntersects(empty), false);
  Assert.equal(full.interiorIntersects(full), true);
  Assert.equal(full.interiorIntersects(empty), false);

  // Singleton cap containing the x-axis.
  S2Cap xaxis = S2Cap.fromPoint(S2Point(1, 0, 0));
  Assert.equal(xaxis.contains(S2Point(1, 0, 0)), true);
  Assert.equal(xaxis.contains(S2Point(1, 1e-20, 0)), false);
  Assert.equal(xaxis.getRadius().radians(), 0.0);

  // Singleton cap containing the y-axis.
  S2Cap yaxis = S2Cap.fromPoint(S2Point(0, 1, 0));
  Assert.equal(yaxis.contains(xaxis.center()), false);
  Assert.equal(xaxis.height(), 0.0);

  // Check that the complement of a singleton cap is the full cap.
  S2Cap xcomp = xaxis.complement();
  Assert.equal(xcomp.isValid(), true);
  Assert.equal(xcomp.isFull(), true);
  Assert.equal(xcomp.contains(xaxis.center()), true);

  // Check that the complement of the complement is *not* the original.
  Assert.equal(xcomp.complement().isValid(), true);
  Assert.equal(xcomp.complement().isEmpty(), true);
  Assert.equal(xcomp.complement().contains(xaxis.center()), false);

  // Check that very small caps can be represented accurately.
  // Here "kTinyRad" is small enough that unit vectors perturbed by this
  // amount along a tangent do not need to be renormalized.
  static const double kTinyRad = 1e-10;
  S2Cap tiny = new S2Cap(S2Point(1, 2, 3).normalize(), S1Angle.fromRadians(kTinyRad));
  Vector3_d tangent = tiny.center().crossProd(S2Point(3, 2, 1)).normalize();
  Assert.equal(tiny.contains(tiny.center() + 0.99 * kTinyRad * tangent), true);
  Assert.equal(tiny.contains(tiny.center() + 1.01 * kTinyRad * tangent), false);

  // Basic tests on a hemispherical cap.
  S2Cap hemi = S2Cap.fromCenterHeight(S2Point(1, 0, 1).normalize(), 1);
  Assert.equal(-hemi.center(), S2Point(hemi.complement().center()));
  Assert.equal(hemi.complement().height(), 1);
  Assert.equal(hemi.contains(S2Point(1, 0, 0)), true);
  Assert.equal(hemi.complement().contains(S2Point(1, 0, 0)), false);
  Assert.equal(hemi.contains(S2Point(1, 0, -(1 - EPSILON)).normalize()), true);
  Assert.equal(hemi.interiorContains(S2Point(1, 0, -(1 + EPSILON)).normalize()), false);

  // A concave cap.  Note that the error bounds for point containment tests
  // increase with the cap angle, so we need to use a larger error bound
  // here.  (It would be painful to do this everywhere, but this at least
  // gives an example of how to compute the maximum error.)
  S2Point center = getLatLngPoint(80, 10);
  S1ChordAngle radius = S1ChordAngle(S1Angle.fromDegrees(150));
  double max_error = radius.getS2PointConstructorMaxError()
      + radius.getS1AngleConstructorMaxError() + 3 * double.epsilon;  // GetLatLngPoint() error
  S2Cap concave = new S2Cap(center, radius);
  S2Cap concave_min = new S2Cap(center, radius.plusError(-max_error));
  S2Cap concave_max = new S2Cap(center, radius.plusError(max_error));
  Assert.equal(concave_max.contains(getLatLngPoint(-70, 10)), true);
  Assert.equal(concave_min.contains(getLatLngPoint(-70, 10)), false);
  Assert.equal(concave_max.contains(getLatLngPoint(-50, -170)), true);
  Assert.equal(concave_min.contains(getLatLngPoint(-50, -170)), false);

  // Cap containment tests.
  Assert.equal(empty.contains(xaxis), false);
  Assert.equal(empty.interiorIntersects(xaxis), false);
  Assert.equal(full.contains(xaxis), true);
  Assert.equal(full.interiorIntersects(xaxis), true);
  Assert.equal(xaxis.contains(full), false);
  Assert.equal(xaxis.interiorIntersects(full), false);
  Assert.equal(xaxis.contains(xaxis), true);
  Assert.equal(xaxis.interiorIntersects(xaxis), false);
  Assert.equal(xaxis.contains(empty), true);
  Assert.equal(xaxis.interiorIntersects(empty), false);
  Assert.equal(hemi.contains(tiny), true);
  Assert.equal(
      hemi.contains(new S2Cap(S2Point(1, 0, 0), S1Angle.fromRadians(math.PI_4 - EPSILON))), true);
  Assert.equal(
      hemi.contains(new S2Cap(S2Point(1, 0, 0), S1Angle.fromRadians(math.PI_4 + EPSILON))),
      false);
  Assert.equal(concave.contains(hemi), true);
  Assert.equal(concave.interiorIntersects(hemi.complement()), true);
  Assert.equal(concave.contains(S2Cap.fromCenterHeight(-concave.center(), 0.1)), false);
}

@("S2Cap.GetRectBound")
unittest {
  // Empty and full caps.
  Assert.equal(S2Cap.empty().getRectBound().isEmpty(), true);
  Assert.equal(S2Cap.full().getRectBound().isFull(), true);

  static const double kDegreeEps = 1e-13;
  // Maximum allowable error for latitudes and longitudes measured in
  // degrees.  (EXPECT_DOUBLE_EQ isn't sufficient.)

  // Cap that includes the south pole.
  S2LatLngRect rect = new S2Cap(getLatLngPoint(-45, 57), S1Angle.fromDegrees(50)).getRectBound();
  Assert.approximately(rect.latLo().degrees(), -90.0, kDegreeEps);
  Assert.approximately(rect.latHi().degrees(), 5.0, kDegreeEps);
  Assert.equal(rect.lng().isFull(), true);

  // Cap that is tangent to the north pole.
  rect = new S2Cap(
      S2Point(1, 0, 1).normalize(), S1Angle.fromRadians(math.PI_4 + 1e-15)).getRectBound();
  Assert.approximately(rect.lat().lo(), 0, EPSILON);
  Assert.approximately(rect.lat().hi(), math.PI_2, EPSILON);
  Assert.equal(rect.lng().isFull(), true);

  rect = new S2Cap(S2Point(1, 0, 1).normalize(), S1Angle.fromDegrees(45 + 2e-14)).getRectBound();
  Assert.approximately(rect.latLo().degrees(), 0, kDegreeEps);
  Assert.approximately(rect.latHi().degrees(), 90, kDegreeEps);
  Assert.equal(true, rect.lng().isFull());

  // The eastern hemisphere.
  rect = new S2Cap(S2Point(0, 1, 0),
               S1Angle.fromRadians(math.PI_2 + 2e-16)).getRectBound();
  Assert.approximately(rect.latLo().degrees(), -90, kDegreeEps);
  Assert.approximately(rect.latHi().degrees(), 90, kDegreeEps);
  Assert.equal(rect.lng().isFull(), true);

  // A cap centered on the equator.
  rect = new S2Cap(getLatLngPoint(0, 50), S1Angle.fromDegrees(20)).getRectBound();
  Assert.approximately(rect.latLo().degrees(), -20, kDegreeEps);
  Assert.approximately(rect.latHi().degrees(), 20, kDegreeEps);
  Assert.approximately(rect.lngLo().degrees(), 30, kDegreeEps);
  Assert.approximately(rect.lngHi().degrees(), 70, kDegreeEps);

  // A cap centered on the north pole.
  rect = new S2Cap(getLatLngPoint(90, 123), S1Angle.fromDegrees(10)).getRectBound();
  Assert.approximately(rect.latLo().degrees(), 80, kDegreeEps);
  Assert.approximately(rect.latHi().degrees(), 90, kDegreeEps);
  Assert.equal(rect.lng().isFull(), true);
}

@("S2Cap.S2CellMethods")
unittest {
  // For each cube face, we construct some cells on
  // that face and some caps whose positions are relative to that face,
  // and then check for the expected intersection/containment results.

  // The distance from the center of a face to one of its vertices.
  const double kFaceRadius = math.atan(math.sqrt(2.0));

  for (int face = 0; face < 6; ++face) {
    // The cell consisting of the entire face.
    S2Cell root_cell = S2Cell.fromFace(face);

    // A leaf cell at the midpoint of the v=1 edge.
    S2Cell edge_cell = new S2Cell(FaceUVtoXYZ(face, 0, 1 - EPSILON));

    // A leaf cell at the u=1, v=1 corner.
    S2Cell corner_cell = new S2Cell(FaceUVtoXYZ(face, 1 - EPSILON, 1 - EPSILON));

    // Quick check for full and empty caps.
    Assert.equal(S2Cap.full().contains(root_cell), true);
    Assert.equal(S2Cap.empty().mayIntersect(root_cell), false);

    // Check intersections with the bounding caps of the leaf cells that are
    // adjacent to 'corner_cell' along the Hilbert curve.  Because this corner
    // is at (u=1,v=1), the curve stays locally within the same cube face.
    S2CellId first = corner_cell.id().advance(-3);
    S2CellId last = corner_cell.id().advance(4);
    for (S2CellId id = first; id < last; id = id.next()) {
      S2Cell cell = new S2Cell(id);
      Assert.equal(id == corner_cell.id(), cell.getCapBound().contains(corner_cell));
      Assert.equal(
          id.parent().contains(corner_cell.id()), cell.getCapBound().mayIntersect(corner_cell));
    }

    int anti_face = (face + 3) % 6;  // Opposite face.
    for (int cap_face = 0; cap_face < 6; ++cap_face) {
      // A cap that barely contains all of 'cap_face'.
      S2Point center = GetNorm(cap_face);
      S2Cap covering = new S2Cap(center, S1Angle.fromRadians(kFaceRadius + EPSILON));
      Assert.equal(cap_face == face, covering.contains(root_cell));
      Assert.equal(cap_face != anti_face, covering.mayIntersect(root_cell));
      Assert.equal(center.dotProd(edge_cell.getCenter()) > 0.1,
                covering.contains(edge_cell));
      Assert.equal(covering.mayIntersect(edge_cell), covering.contains(edge_cell));
      Assert.equal(cap_face == face, covering.contains(corner_cell));
      Assert.equal(center.dotProd(corner_cell.getCenter()) > 0,
                covering.mayIntersect(corner_cell));

      // A cap that barely intersects the edges of 'cap_face'.
      S2Cap bulging = new S2Cap(center, S1Angle.fromRadians(math.PI_4 + EPSILON));
      Assert.equal(false, bulging.contains(root_cell));
      Assert.equal(cap_face != anti_face, bulging.mayIntersect(root_cell));
      Assert.equal(cap_face == face, bulging.contains(edge_cell));
      Assert.equal(center.dotProd(edge_cell.getCenter()) > 0.1,
                bulging.mayIntersect(edge_cell));
      Assert.equal(false, bulging.contains(corner_cell));
      Assert.equal(false, bulging.mayIntersect(corner_cell));

      // A singleton cap.
      S2Cap singleton = new S2Cap(center, S1Angle.zero());
      Assert.equal(cap_face == face, singleton.mayIntersect(root_cell));
      Assert.equal(false, singleton.mayIntersect(edge_cell));
      Assert.equal(false, singleton.mayIntersect(corner_cell));
    }
  }
}

@("S2Cap.GetCellUnionBoundLevel1Radius")
unittest {
  // Check that a cap whose radius is approximately the width of a level 1
  // S2Cell can be covered by only 3 faces.
  S2Cap cap = new S2Cap(
      S2Point(1, 1, 1).normalize(), S1Angle.fromRadians(s2metrics.MIN_WIDTH.getValue(1)));
  S2CellId[] covering;
  cap.getCellUnionBound(covering);
  Assert.equal(covering.length, 3);
}

@("S2Cap.Expanded")
unittest {
  Assert.equal(S2Cap.empty().expanded(S1Angle.fromRadians(2)).isEmpty(), true);
  Assert.equal(S2Cap.full().expanded(S1Angle.fromRadians(2)).isFull(), true);
  S2Cap cap50 = new S2Cap(S2Point(1, 0, 0), S1Angle.fromDegrees(50));
  S2Cap cap51 = new S2Cap(S2Point(1, 0, 0), S1Angle.fromDegrees(51));
  Assert.equal(cap50.expanded(S1Angle.fromRadians(0)).approxEquals(cap50), true);
  Assert.equal(cap50.expanded(S1Angle.fromDegrees(1)).approxEquals(cap51), true);
  Assert.equal(cap50.expanded(S1Angle.fromDegrees(129.99)).isFull(), false);
  Assert.equal(cap50.expanded(S1Angle.fromDegrees(130.01)).isFull(), true);
}

@("S2Cap.GetCentroid")
unittest {
  // Empty and full caps.
  Assert.equal(S2Point(), S2Cap.empty().getCentroid());
  Assert.notGreaterThan(S2Cap.full().getCentroid().norm(), 1e-15);

  // Random caps.
  for (int i = 0; i < 100; ++i) {
    S2Point center = S2Testing.randomPoint();
    double height = S2Testing.rnd.uniformDouble(0.0, 2.0);
    S2Cap cap = S2Cap.fromCenterHeight(center, height);
    S2Point centroid = cap.getCentroid();
    S2Point expected = center * (1.0 - height / 2.0) * cap.getArea();
    Assert.notGreaterThan((expected - centroid).norm(), 1e-15);
  }
}

@("S2Cap.Unite")
unittest {
  // Two caps which have the same center but one has a larger radius.
  S2Cap a = new S2Cap(getLatLngPoint(50.0, 10.0), S1Angle.fromDegrees(0.2));
  S2Cap b = new S2Cap(getLatLngPoint(50.0, 10.0), S1Angle.fromDegrees(0.3));
  Assert.equal(b.contains(a), true);
  Assert.equal(a.unite(b), b);

  // Two caps where one is the full cap.
  Assert.equal(a.unite(S2Cap.full()).isFull(), true);

  // Two caps where one is the empty cap.
  Assert.equal(a.unite(S2Cap.empty()), a);

  // Two caps which have different centers, one entirely encompasses the other.
  S2Cap c = new S2Cap(getLatLngPoint(51.0, 11.0), S1Angle.fromDegrees(1.5));
  Assert.equal(c.contains(a), true);
  Assert.equal(a.unite(c).center(), c.center());
  Assert.equal(a.unite(c).getRadius(), c.getRadius());

  // Two entirely disjoint caps.
  S2Cap d = new S2Cap(getLatLngPoint(51.0, 11.0), S1Angle.fromDegrees(0.1));
  Assert.equal(d.contains(a), false);
  Assert.equal(d.intersects(a), false);
  Assert.equal(a.unite(d).approxEquals(d.unite(a)), true);
  Assert.approximately(S2LatLng(a.unite(d).center()).lat().degrees(), 50.4588, 0.001);
  Assert.approximately(S2LatLng(a.unite(d).center()).lng().degrees(), 10.4525, 0.001);
  Assert.approximately(a.unite(d).getRadius().degrees(), 0.7425, 0.001);

  // Two partially overlapping caps.
  S2Cap e = new S2Cap(getLatLngPoint(50.3, 10.3), S1Angle.fromDegrees(0.2));
  Assert.equal(e.contains(a), false);
  Assert.equal(e.intersects(a), true);
  Assert.equal(a.unite(e).approxEquals(e.unite(a)), true);
  Assert.approximately(S2LatLng(a.unite(e).center()).lat().degrees(), 50.1500, 0.001);
  Assert.approximately(S2LatLng(a.unite(e).center()).lng().degrees(), 10.1495, 0.001);
  Assert.approximately(a.unite(e).getRadius().degrees(), 0.3781, 0.001);

  // Two very large caps, whose radius sums to in excess of 180 degrees, and
  // whose centers are not antipodal.
  S2Cap f = new S2Cap(S2Point(0, 0, 1).normalize(), S1Angle.fromDegrees(150));
  S2Cap g = new S2Cap(S2Point(0, 1, 0).normalize(), S1Angle.fromDegrees(150));
  Assert.equal(f.unite(g).isFull(), true);

  // Two non-overlapping hemisphere caps with antipodal centers.
  S2Cap hemi = S2Cap.fromCenterHeight(S2Point(0, 0, 1).normalize(), 1);
  Assert.equal(hemi.unite(hemi.complement()).isFull(), true);
}

@("S2Cap.EncodeDecode")
unittest {
  S2Cap cap = S2Cap.fromCenterHeight(S2Point(3, 2, 1).normalize(), 1);
  auto encoder = makeEncoder();
  cap.encode(encoder);
  auto decoder = makeDecoder(encoder.buffer().data());
  auto decoded_cap = new S2Cap();
  Assert.equal(decoded_cap.decode(decoder), true);
  Assert.equal(decoded_cap, cap);
}
