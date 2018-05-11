module s2.s2latlngrectbounder;

// Original author: ericv@google.com (Eric Veach)

import algorithm = std.algorithm;
import math = std.math;
import s2.r1interval;
import s2.s1chordangle;
import s2.s1interval;
import s2.s2latlng;
import s2.s2latlngrect;
import s2.s2point;
import s2.s2pointutil;
import s2.util.math.vector;

/**
 * This class computes a bounding rectangle that contains all edges defined
 * by a vertex chain v0, v1, v2, ...  All vertices must be unit length.
 * Note that the bounding rectangle of an edge can be larger than the
 * bounding rectangle of its endpoints, e.g. consider an edge that passes
 * through the north pole.
 *
 * The bounds are calculated conservatively to account for numerical errors
 * when S2Points are converted to S2LatLngs.  More precisely, this class
 * guarantees the following.  Let L be a closed edge chain (loop) such that
 * the interior of the loop does not contain either pole.  Now if P is any
 * point such that L.Contains(P), then RectBound(L).Contains(S2LatLng(P)).
 */
class S2LatLngRectBounder {
private:
  // Common back end for AddPoint() and AddLatLng().  b and b_latlng
  // must refer to the same vertex.
  void addInternal(in S2Point b, in S2LatLng b_latlng)
  in {
    // Simple consistency check to verify that b and b_latlng are alternate
    // representations of the same vertex.
    assert(approxEquals(b, b_latlng.toS2Point()));
  } body {
    if (_bound.isEmpty()) {
      _bound.addPoint(b_latlng);
    } else {
      // First compute the cross product N = A x B robustly.  This is the normal
      // to the great circle through A and B.  We don't use S2::RobustCrossProd()
      // since that method returns an arbitrary vector orthogonal to A if the two
      // vectors are proportional, and we want the zero vector in that case.
      Vector3_d n = (_a - b).crossProd(_a + b);  // N = 2 * (A x B)

      // The relative error in N gets large as its norm gets very small (i.e.,
      // when the two points are nearly identical or antipodal).  We handle this
      // by choosing a maximum allowable error, and if the error is greater than
      // this we fall back to a different technique.  Since it turns out that
      // the other sources of error in converting the normal to a maximum
      // latitude add up to at most 1.16 * DBL_EPSILON (see below), and it is
      // desirable to have the total error be a multiple of DBL_EPSILON, we have
      // chosen to limit the maximum error in the normal to 3.84 * DBL_EPSILON.
      // It is possible to show that the error is less than this when
      //
      //   n.Norm() >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * DBL_EPSILON
      //            = 1.91346e-15 (about 8.618 * DBL_EPSILON)
      double n_norm = n.norm();
      if (n_norm < 1.91346e-15) {
        // A and B are either nearly identical or nearly antipodal (to within
        // 4.309 * DBL_EPSILON, or about 6 nanometers on the earth's surface).
        if (_a.dotProd(b) < 0) {
          // The two points are nearly antipodal.  The easiest solution is to
          // assume that the edge between A and B could go in any direction
          // around the sphere.
          _bound = S2LatLngRect.full();
        } else {
          // The two points are nearly identical (to within 4.309 * DBL_EPSILON).
          // In this case we can just use the bounding rectangle of the points,
          // since after the expansion done by GetBound() this rectangle is
          // guaranteed to include the (lat,lng) values of all points along AB.
          _bound = _bound.unite(S2LatLngRect.fromPointPair(_aLatLng, b_latlng));
        }
      } else {
        // Compute the longitude range spanned by AB.
        S1Interval lng_ab = S1Interval.fromPointPair(
            _aLatLng.lng().radians(), b_latlng.lng().radians());
        if (lng_ab.getLength() >= math.PI - 2 * double.epsilon) {
          // The points lie on nearly opposite lines of longitude to within the
          // maximum error of the calculation.  (Note that this test relies on
          // the fact that M_PI is slightly less than the true value of Pi, and
          // that representable values near M_PI are 2 * DBL_EPSILON apart.)
          // The easiest solution is to assume that AB could go on either side
          // of the pole.
          lng_ab = S1Interval.full();
        }

        // Next we compute the latitude range spanned by the edge AB.  We start
        // with the range spanning the two endpoints of the edge:
        R1Interval lat_ab = R1Interval.fromPointPair(
            _aLatLng.lat().radians(), b_latlng.lat().radians());

        // This is the desired range unless the edge AB crosses the plane
        // through N and the Z-axis (which is where the great circle through A
        // and B attains its minimum and maximum latitudes).  To test whether AB
        // crosses this plane, we compute a vector M perpendicular to this
        // plane and then project A and B onto it.
        Vector3_d m = n.crossProd(S2Point(0, 0, 1));
        double m_a = m.dotProd(_a);
        double m_b = m.dotProd(b);

        // We want to test the signs of "m_a" and "m_b", so we need to bound
        // the error in these calculations.  It is possible to show that the
        // total error is bounded by
        //
        //  (1 + sqrt(3)) * DBL_EPSILON * n_norm + 8 * sqrt(3) * (DBL_EPSILON**2)
        //    = 6.06638e-16 * n_norm + 6.83174e-31

        double m_error = 6.06638e-16 * n_norm + 6.83174e-31;
        if (m_a * m_b < 0 || math.fabs(m_a) <= m_error || math.fabs(m_b) <= m_error) {
          // Minimum/maximum latitude *may* occur in the edge interior.
          //
          // The maximum latitude is 90 degrees minus the latitude of N.  We
          // compute this directly using atan2 in order to get maximum accuracy
          // near the poles.
          //
          // Our goal is compute a bound that contains the computed latitudes of
          // all S2Points P that pass the point-in-polygon containment test.
          // There are three sources of error we need to consider:
          //  - the directional error in N (at most 3.84 * DBL_EPSILON)
          //  - converting N to a maximum latitude
          //  - computing the latitude of the test point P
          // The latter two sources of error are at most 0.955 * DBL_EPSILON
          // individually, but it is possible to show by a more complex analysis
          // that together they can add up to at most 1.16 * DBL_EPSILON, for a
          // total error of 5 * DBL_EPSILON.
          //
          // We add 3 * DBL_EPSILON to the bound here, and GetBound() will pad
          // the bound by another 2 * DBL_EPSILON.
          double max_lat = algorithm.min(
              math.atan2(math.sqrt(n[0]*n[0] + n[1]*n[1]), math.fabs(n[2])) + 3 * double.epsilon,
              math.PI_2);

          // In order to get tight bounds when the two points are close together,
          // we also bound the min/max latitude relative to the latitudes of the
          // endpoints A and B.  First we compute the distance between A and B,
          // and then we compute the maximum change in latitude between any two
          // points along the great circle that are separated by this distance.
          // This gives us a latitude change "budget".  Some of this budget must
          // be spent getting from A to B; the remainder bounds the round-trip
          // distance (in latitude) from A or B to the min or max latitude
          // attained along the edge AB.
          double lat_budget = 2 * math.asin(0.5 * (_a - b).norm() * math.sin(max_lat));
          double max_delta = 0.5 * (lat_budget - lat_ab.getLength()) + double.epsilon;

          // Test whether AB passes through the point of maximum latitude or
          // minimum latitude.  If the dot product(s) are small enough then the
          // result may be ambiguous.
          if (m_a <= m_error && m_b >= -m_error) {
            lat_ab.setHi(algorithm.min(max_lat, lat_ab.hi() + max_delta));
          }
          if (m_b <= m_error && m_a >= -m_error) {
            lat_ab.setLo(algorithm.max(-max_lat, lat_ab.lo() - max_delta));
          }
        }
        _bound = _bound.unite(new S2LatLngRect(lat_ab, lng_ab));
      }
    }
    _a = b;
    _aLatLng = b_latlng;
  }

  S2Point _a;             // The previous vertex in the chain.
  S2LatLng _aLatLng;     // The corresponding latitude-longitude.
  S2LatLngRect _bound;    // The current bounding rectangle.

public:
  this() {
    _bound = S2LatLngRect.empty();
  }

  /**
   * This method is called to add a vertex to the chain when the vertex is
   * represented as an S2Point.  Requires that 'b' has unit length.  Repeated
   * vertices are ignored.
   */
  void addPoint(in S2Point b)
  in {
    assert(isUnitLength(b));
  } body {
    addInternal(b, S2LatLng(b));
  }

  /**
   * This method is called to add a vertex to the chain when the vertex is
   * represented as an S2LatLng.  Repeated vertices are ignored.
   */
  void addLatLng(in S2LatLng b_latlng) {
    addInternal(b_latlng.toS2Point(), b_latlng);
  }

  /**
   * Returns the bounding rectangle of the edge chain that connects the
   * vertices defined so far.  This bound satisfies the guarantee made
   * above, i.e. if the edge chain defines a loop, then the bound contains
   * the S2LatLng coordinates of all S2Points contained by the loop.
   */
  S2LatLngRect getBound() const {
    // To save time, we ignore numerical errors in the computed S2LatLngs while
    // accumulating the bounds and then account for them here.
    //
    // S2LatLng(S2Point) has a maximum error of 0.955 * DBL_EPSILON in latitude.
    // In the worst case, we might have rounded "inwards" when computing the
    // bound and "outwards" when computing the latitude of a contained point P,
    // therefore we expand the latitude bounds by 2 * DBL_EPSILON in each
    // direction.  (A more complex analysis shows that 1.5 * DBL_EPSILON is
    // enough, but the expansion amount should be a multiple of DBL_EPSILON in
    // order to avoid rounding errors during the expansion itself.)
    //
    // S2LatLng(S2Point) has a maximum error of DBL_EPSILON in longitude, which
    // is simply the maximum rounding error for results in the range [-Pi, Pi].
    // This is true because the Gnu implementation of atan2() comes from the IBM
    // Accurate Mathematical Library, which implements correct rounding for this
    // instrinsic (i.e., it returns the infinite precision result rounded to the
    // nearest representable value, with ties rounded to even values).  This
    // implies that we don't need to expand the longitude bounds at all, since
    // we only guarantee that the bound contains the *rounded* latitudes of
    // contained points.  The *true* latitudes of contained points may lie up to
    // DBL_EPSILON outside of the returned bound.

    const S2LatLng kExpansion = S2LatLng.fromRadians(2 * double.epsilon, 0);
    return _bound.expanded(kExpansion).polarClosure();
  }

  /**
   * Expands a bound returned by GetBound() so that it is guaranteed to
   * contain the bounds of any subregion whose bounds are computed using
   * this class.  For example, consider a loop L that defines a square.
   * GetBound() ensures that if a point P is contained by this square, then
   * S2LatLng(P) is contained by the bound.  But now consider a diamond
   * shaped loop S contained by L.  It is possible that GetBound() returns a
   * *larger* bound for S than it does for L, due to rounding errors.  This
   * method expands the bound for L so that it is guaranteed to contain the
   * bounds of any subregion S.
   *
   * More precisely, if L is a loop that does not contain either pole, and S
   * is a loop such that L.Contains(S), then
   *
   *   ExpandForSubregions(RectBound(L)).Contains(RectBound(S)).
   */
  static S2LatLngRect expandForSubregions(in S2LatLngRect bound) {
    // Empty bounds don't need expansion.
    if (bound.isEmpty()) return new S2LatLngRect(bound);

    // First we need to check whether the bound B contains any nearly-antipodal
    // points (to within 4.309 * DBL_EPSILON).  If so then we need to return
    // S2LatLngRect::Full(), since the subregion might have an edge between two
    // such points, and AddPoint() returns Full() for such edges.  Note that
    // this can happen even if B is not Full(); for example, consider a loop
    // that defines a 10km strip straddling the equator extending from
    // longitudes -100 to +100 degrees.
    //
    // It is easy to check whether B contains any antipodal points, but checking
    // for nearly-antipodal points is trickier.  Essentially we consider the
    // original bound B and its reflection through the origin B', and then test
    // whether the minimum distance between B and B' is less than 4.309 *
    // DBL_EPSILON.

    // "lng_gap" is a lower bound on the longitudinal distance between B and its
    // reflection B'.  (2.5 * DBL_EPSILON is the maximum combined error of the
    // endpoint longitude calculations and the GetLength() call.)
    double lng_gap = algorithm.max(0.0, math.PI - bound.lng().getLength() - 2.5 * double.epsilon);

    // "min_abs_lat" is the minimum distance from B to the equator (if zero or
    // negative, then B straddles the equator).
    double min_abs_lat = algorithm.max(bound.lat().lo(), -bound.lat().hi());

    // "lat_gap1" and "lat_gap2" measure the minimum distance from B to the
    // south and north poles respectively.
    double lat_gap1 = math.PI_2 + bound.lat().lo();
    double lat_gap2 = math.PI_2 - bound.lat().hi();

    if (min_abs_lat >= 0) {
      // The bound B does not straddle the equator.  In this case the minimum
      // distance is between one endpoint of the latitude edge in B closest to
      // the equator and the other endpoint of that edge in B'.  The latitude
      // distance between these two points is 2*min_abs_lat, and the longitude
      // distance is lng_gap.  We could compute the distance exactly using the
      // Haversine formula, but then we would need to bound the errors in that
      // calculation.  Since we only need accuracy when the distance is very
      // small (close to 4.309 * DBL_EPSILON), we substitute the Euclidean
      // distance instead.  This gives us a right triangle XYZ with two edges of
      // length x = 2*min_abs_lat and y ~= lng_gap.  The desired distance is the
      // length of the third edge "z", and we have
      //
      //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
      //
      // Therefore the region may contain nearly antipodal points only if
      //
      //  2*min_abs_lat + lng_gap  <  sqrt(2) * 4.309 * DBL_EPSILON
      //                           ~= 1.354e-15
      //
      // Note that because the given bound B is conservative, "min_abs_lat" and
      // "lng_gap" are both lower bounds on their true values so we do not need
      // to make any adjustments for their errors.
      if (2 * min_abs_lat + lng_gap < 1.354e-15) {
        return S2LatLngRect.full();
      }
    } else if (lng_gap >= math.PI_2) {
      // B spans at most Pi/2 in longitude.  The minimum distance is always
      // between one corner of B and the diagonally opposite corner of B'.  We
      // use the same distance approximation that we used above; in this case
      // we have an obtuse triangle XYZ with two edges of length x = lat_gap1
      // and y = lat_gap2, and angle Z >= Pi/2 between them.  We then have
      //
      //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
      //
      // Unlike the case above, "lat_gap1" and "lat_gap2" are not lower bounds
      // (because of the extra addition operation, and because M_PI_2 is not
      // exactly equal to Pi/2); they can exceed their true values by up to
      // 0.75 * DBL_EPSILON.  Putting this all together, the region may
      // contain nearly antipodal points only if
      //
      //   lat_gap1 + lat_gap2  <  (sqrt(2) * 4.309 + 1.5) * DBL_EPSILON
      //                        ~= 1.687e-15
      if (lat_gap1 + lat_gap2 < 1.687e-15) {
        return S2LatLngRect.full();
      }
    } else {
      // Otherwise we know that (1) the bound straddles the equator and (2) its
      // width in longitude is at least Pi/2.  In this case the minimum
      // distance can occur either between a corner of B and the diagonally
      // opposite corner of B' (as in the case above), or between a corner of B
      // and the opposite longitudinal edge reflected in B'.  It is sufficient
      // to only consider the corner-edge case, since this distance is also a
      // lower bound on the corner-corner distance when that case applies.

      // Consider the spherical triangle XYZ where X is a corner of B with
      // minimum absolute latitude, Y is the closest pole to X, and Z is the
      // point closest to X on the opposite longitudinal edge of B'.  This is a
      // right triangle (Z = Pi/2), and from the spherical law of sines we have
      //
      //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
      //     sin(max_lat_gap) / 1  =  sin(d_min) / sin(lng_gap)
      //     sin(d_min)  =  sin(max_lat_gap) * sin(lng_gap)
      //
      // where "max_lat_gap" = max(lat_gap1, lat_gap2) and "d_min" is the
      // desired minimum distance.  Now using the facts that sin(t) >= (2/Pi)*t
      // for 0 <= t <= Pi/2, that we only need an accurate approximation when
      // at least one of "max_lat_gap" or "lng_gap" is extremely small (in
      // which case sin(t) ~= t), and recalling that "max_lat_gap" has an error
      // of up to 0.75 * DBL_EPSILON, we want to test whether
      //
      //   max_lat_gap * lng_gap  <  (4.309 + 0.75) * (Pi/2) * DBL_EPSILON
      //                          ~= 1.765e-15
      if (algorithm.max(lat_gap1, lat_gap2) * lng_gap < 1.765e-15) {
        return S2LatLngRect.full();
      }
    }
    // Next we need to check whether the subregion might contain any edges that
    // span (M_PI - 2 * DBL_EPSILON) radians or more in longitude, since AddPoint
    // sets the longitude bound to Full() in that case.  This corresponds to
    // testing whether (lng_gap <= 0) in "lng_expansion" below.

    // Otherwise, the maximum latitude error in AddPoint is 4.8 * DBL_EPSILON.
    // In the worst case, the errors when computing the latitude bound for a
    // subregion could go in the opposite direction as the errors when computing
    // the bound for the original region, so we need to double this value.
    // (More analysis shows that it's okay to round down to a multiple of
    // DBL_EPSILON.)
    //
    // For longitude, we rely on the fact that atan2 is correctly rounded and
    // therefore no additional bounds expansion is necessary.

    double lat_expansion = 9 * double.epsilon;
    double lng_expansion = (lng_gap <= 0) ? math.PI : 0;
    return bound.expanded(S2LatLng.fromRadians(lat_expansion, lng_expansion)).polarClosure();
  }

  // Returns the maximum error in GetBound() provided that the result does
  // not include either pole.  It is only to be used for testing purposes
  // (e.g., by passing it to S2LatLngRect::ApproxEquals).
  static S2LatLng maxErrorForTests() {
    // The maximum error in the latitude calculation is
    //    3.84 * DBL_EPSILON   for the RobustCrossProd calculation
    //    0.96 * DBL_EPSILON   for the Latitude() calculation
    //    5    * DBL_EPSILON   added by AddPoint/GetBound to compensate for error
    //    ------------------
    //    9.80 * DBL_EPSILON   maximum error in result
    //
    // The maximum error in the longitude calculation is DBL_EPSILON.  GetBound
    // does not do any expansion because this isn't necessary in order to
    // bound the *rounded* longitudes of contained points.
    return S2LatLng.fromRadians(10 * double.epsilon, 1 * double.epsilon);
  }
}
