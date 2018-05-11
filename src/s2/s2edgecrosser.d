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

module s2.s2edgecrosser;

import math = std.math;
import s2.logger;
import s2.s2edgecrossings;
import s2.s2point;
import s2.s2pointutil;
import s2.s2predicates;
import s2.util.math.vector;

/**
 * This class allows edges to be efficiently tested for intersection with a
 * given fixed edge AB.  It is especially efficient when testing for
 * intersection with an edge chain connecting vertices v0, v1, v2, ...
 *
 * Example usage:
 *
 *   void CountIntersections(const S2Point& a, const S2Point& b,
 *                           const vector<pair<S2Point, S2Point>>& edges) {
 *     int count = 0;
 *     S2EdgeCrosser crosser(&a, &b);
 *     for (const auto& edge : edges) {
 *       if (crosser.CrossingSign(&edge.first, &edge.second) >= 0) {
 *         ++count;
 *       }
 *     }
 *     return count;
 *   }
 *
 * This class expects that the client already has all the necessary vertices
 * stored in memory, so that this class can refer to them with pointers and
 * does not need to make its own copies.  If this is not the case (e.g., you
 * want to pass temporary objects as vertices), see S2CopyingEdgeCrosser.
 */
class S2EdgeCrosser {
public:
  // Default constructor; must be followed by a call to Init().
  this() {}

  // Convenience constructor that calls Init() with the given fixed edge AB.
  // The arguments "a" and "b" must point to values that persist for the
  // lifetime of the S2EdgeCrosser object (or until the next Init() call).
  this(ref in S2Point a, ref in S2Point b) {
    if (!isUnitLength(a)) logger.logWarn("Invalid");
    if (!isUnitLength(b)) logger.logWarn("Invalid");

    _a = &a;
    _b = &b;
    _aCrossB = _a.crossProd(*_b);
    _haveTangents = false;
    _c = null;
  }

  // Accessors for the constructor arguments.
  const (S2Point*) a() {
    return _a;
  }

  const (S2Point*) b() {
    return _b;
  }

  // Initialize the S2EdgeCrosser with the given fixed edge AB.  The arguments
  // "a" and "b" must point to values that persist for the lifetime of the
  // S2EdgeCrosser object (or until the next Init() call).
  void init(ref in S2Point a, ref in S2Point b) {
    _a = &a;
    _b = &b;
    _aCrossB = a.crossProd(*_b);
    _haveTangents = false;
    _c = null;
  }

  /**
   * This function determines whether the edge AB intersects the edge CD.
   * Returns +1 if AB crosses CD at a point that is interior to both edges.
   * Returns  0 if any two vertices from different edges are the same.
   * Returns -1 otherwise.
   *
   * Note that if an edge is degenerate (A == B or C == D), the return value
   * is 0 if two vertices from different edges are the same and -1 otherwise.
   *
   * Properties of CrossingSign:
   *
   *  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
   *  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
   *  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
   *  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
   *
   * This function implements an exact, consistent perturbation model such
   * that no three points are ever considered to be collinear.  This means
   * that even if you have 4 points A, B, C, D that lie exactly in a line
   * (say, around the equator), C and D will be treated as being slightly to
   * one side or the other of AB.  This is done in a way such that the
   * results are always consistent (see s2pred::Sign).
   *
   * Note that if you want to check an edge against a chain of other edges,
   * it is slightly more efficient to use the single-argument version of
   * CrossingSign below.
   *
   * The arguments must point to values that persist until the next call.
   */
  int crossingSign(ref in S2Point c, ref in S2Point d) {
    if (_c != &c) {
      restartAt(c);
    }
    return crossingSign(d);
  }

  /**
   * This method extends the concept of a "crossing" to the case where AB
   * and CD have a vertex in common.  The two edges may or may not cross,
   * according to the rules defined in VertexCrossing() below.  The rules
   * are designed so that point containment tests can be implemented simply
   * by counting edge crossings.  Similarly, determining whether one edge
   * chain crosses another edge chain can be implemented by counting.
   *
   * Returns true if CrossingSign(c, d) > 0, or AB and CD share a vertex
   * and VertexCrossing(a, b, c, d) returns true.
   *
   * The arguments must point to values that persist until the next call.
   */
  bool edgeOrVertexCrossing(ref in S2Point c, ref in S2Point d) {
    if (_c != &c) {
      restartAt(c);
    }
    return edgeOrVertexCrossing(d);
  }

  ///////////////////////// Edge Chain Methods ///////////////////////////
  //
  // You don't need to use these unless you're trying to squeeze out every
  // last drop of performance.  Essentially all you are saving is a test
  // whether the first vertex of the current edge is the same as the second
  // vertex of the previous edge.  Example usage:
  //
  //   vector<S2Point> chain;
  //   crosser.RestartAt(&chain[0]);
  //   for (int i = 1; i < chain.size(); ++i) {
  //     if (crosser.EdgeOrVertexCrossing(&chain[i])) { ++count; }
  //   }

  /**
   * Convenience constructor that uses AB as the fixed edge, and C as the
   * first vertex of the vertex chain (equivalent to calling RestartAt(c)).
   *
   * The arguments must point to values that persist until the next call.
   */
  this(ref in S2Point a, ref in S2Point b, ref in S2Point c) {
    if (!isUnitLength(a)) logger.logWarn("Invalid");
    if (!isUnitLength(b)) logger.logWarn("Invalid");

    _a = &a;
    _b = &b;
    _aCrossB = _a.crossProd(*_b);
    _haveTangents = false;
    restartAt(c);
  }


  /**
   * Call this method when your chain 'jumps' to a new place.
   * The argument must point to a value that persists until the next call.
   */
  void restartAt(ref in S2Point c) {
    if (!isUnitLength(c)) logger.logWarn("Invalid");
    _c = &c;
    _acb = -triageSign(*_a, *_b, *_c, _aCrossB);
  }


  /**
   * Like CrossingSign above, but uses the last vertex passed to one of
   * the crossing methods (or RestartAt) as the first vertex of the current
   * edge.
   *
   * The argument must point to a value that persists until the next call.
   */
  int crossingSign(ref in S2Point d) {
    if (!isUnitLength(d)) logger.logWarn("Invalid");
    // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
    // all be oriented the same way (CW or CCW).  We keep the orientation of ACB
    // as part of our state.  When each new point D arrives, we compute the
    // orientation of BDA and check whether it matches ACB.  This checks whether
    // the points C and D are on opposite sides of the great circle through AB.

    // Recall that TriageSign is invariant with respect to rotating its
    // arguments, i.e. ABD has the same orientation as BDA.
    int bda = triageSign(*_a, *_b, d, _aCrossB);
    if (_acb == -bda && bda != 0) {
      // The most common case -- triangles have opposite orientations.  Save the
      // current vertex D as the next vertex C, and also save the orientation of
      // the new triangle ACB (which is opposite to the current triangle BDA).
      _c = &d;
      _acb = -bda;
      return -1;
    }
    _bda = bda;
    return crossingSignInternal(d);
  }

  /**
   * Like EdgeOrVertexCrossing above, but uses the last vertex passed to one
   * of the crossing methods (or RestartAt) as the first vertex of the
   * current edge.
   *
   * The argument must point to a value that persists until the next call.
   */
  bool edgeOrVertexCrossing(ref in S2Point d) {
    // We need to copy c_ since it is clobbered by CrossingSign().
    const S2Point c = *_c;
    int crossing = crossingSign(d);
    if (crossing < 0) {
      return false;
    }
    if (crossing > 0) {
      return true;
    }
    return vertexCrossing(*_a, *_b, c, d);
  }

  /**
   * Returns the last vertex of the current edge chain being tested, i.e. the
   * C vertex that will be used to construct the edge CD when one of the
   * methods above is called.
   */
  @property
  const(S2Point)* c() {
    return _c;
  }

private:
  /** These functions handle the "slow path" of CrossingSign(). */
  int crossingSignInternal(ref in S2Point d) {
    // Compute the actual result, and then save the current vertex D as the next
    // vertex C, and save the orientation of the next triangle ACB (which is
    // opposite to the current triangle BDA).
    int result = crossingSignInternal2(d);
    _c = &d;
    _acb = -_bda;
    return result;
  }

  int crossingSignInternal2(ref in S2Point d) {
    // At this point, a very common situation is that A,B,C,D are four points on
    // a line such that AB does not overlap CD.  (For example, this happens when
    // a line or curve is sampled finely, or when geometry is constructed by
    // computing the union of S2CellIds.)  Most of the time, we can determine
    // that AB and CD do not intersect by computing the two outward-facing
    // tangents at A and B (parallel to AB) and testing whether AB and CD are on
    // opposite sides of the plane perpendicular to one of these tangents.  This
    // is moderately expensive but still much cheaper than s2pred::ExpensiveSign.
    if (!_haveTangents) {
      S2Point norm = robustCrossProd(*_a, *_b).normalize();
      _aTangent = _a.crossProd(norm);
      _bTangent = norm.crossProd(*_b);
      _haveTangents = true;
    }
    // The error in robustCrossProd() is insignificant.  The maximum error in
    // the call to crossProd() (i.e., the maximum norm of the error vector) is
    // (0.5 + 1/sqrt(3)) * double.epsilon.  The maximum error in each call to
    // DotProd() below is double.epsilon.  (There is also a small relative error
    // term that is insignificant because we are comparing the result against a
    // constant that is very close to zero.)
    static const double kError = (1.5 + 1 / math.sqrt(3.0)) * double.epsilon;
    if ((_c.dotProd(_aTangent) > kError && d.dotProd(_aTangent) > kError) ||
        (_c.dotProd(_bTangent) > kError && d.dotProd(_bTangent) > kError)) {
      return -1;
    }

    // Otherwise, eliminate the cases where two vertices from different edges
    // are equal.  (These cases could be handled in the code below, but we would
    // rather avoid calling ExpensiveSign whenever possible.)
    if (*_a == *_c || *_a == d || *_b == *_c || *_b == d) {
      return 0;
    }

    // Eliminate cases where an input edge is degenerate.  (Note that in most
    // cases, if CD is degenerate then this method is not even called because
    // _acb and bda have different signs.)
    if (*_a == *_b || *_c == d) {
      return -1;
    }

    // Otherwise it's time to break out the big guns.
    if (_acb == 0) {
      _acb = -expensiveSign(*_a, *_b, *_c);
    }
    assert(_acb != 0);
    if (_bda == 0) {
      _bda = expensiveSign(*_a, *_b, d);
    }
    assert(_bda != 0);
    if (_bda != _acb) {
      return -1;
    }

    Vector3_d c_cross_d = _c.crossProd(d);
    int cbd = -sign(*_c, d, *_b, c_cross_d);
    assert(cbd != 0);
    if (cbd != _acb) {
      return -1;
    }
    int dac = sign(*_c, d, *_a, c_cross_d);
    assert(dac != 0);
    return (dac != _acb) ? -1 : 1;
  }

  /** Used internally by S2CopyingEdgeCrosser.  Updates "_c" only. */
  void setC(ref S2Point c) {
    _c = &c;
  }

  // The fields below are constant after the call to Init().
  const(S2Point)* _a;
  const(S2Point)* _b;
  Vector3_d _aCrossB;

  // To reduce the number of calls to s2pred::ExpensiveSign(), we compute an
  // outward-facing tangent at A and B if necessary.  If the plane
  // perpendicular to one of these tangents separates AB from CD (i.e., one
  // edge on each side) then there is no intersection.
  bool _haveTangents;  // True if the tangents have been computed.
  S2Point _aTangent;   // Outward-facing tangent at A.
  S2Point _bTangent;   // Outward-facing tangent at B.

  // The fields below are updated for each vertex in the chain.
  const(S2Point)* _c;      // Previous vertex in the vertex chain.
  int _acb;                // The orientation of triangle ACB.

  // The field below is a temporary used by CrossingSignInternal().
  int _bda;                // The orientation of triangle BDA.
}

/**
 * S2CopyingEdgeCrosser is exactly like S2EdgeCrosser, except that it makes its
 * own copy of all arguments so that they do not need to persist between
 * calls.  This is less efficient, but makes it possible to use points that
 * are generated on demand and cannot conveniently be stored by the client.
 */
class S2CopyingEdgeCrosser {
public:
  // These methods are all exactly like S2EdgeCrosser, except that the
  // arguments can be temporaries.
  this() {}
  this(in S2Point a, in S2Point b) {
    _a = a;
    _b = b;
    _c = S2Point();
    _crosser = new S2EdgeCrosser(_a, _b);
  }

  S2Point a() {
    return _a;
  }

  S2Point b() {
    return _b;
  }

  S2Point c() {
    return _c;
  }

  void init(in S2Point a, in S2Point b) {
    _a = a;
    _b = b;
    _c = S2Point();
    _crosser.init(_a, _b);
  }

  int crossingSign(ref in S2Point c, ref in S2Point d) {
    if (c != _c || _crosser._c == null) {
      restartAt(c);
    }
    return crossingSign(d);
  }

  bool edgeOrVertexCrossing(in S2Point c, in S2Point d) {
    if (c != _c || _crosser._c == null) {
      restartAt(c);
    }
    return edgeOrVertexCrossing(d);
  }

  this(in S2Point a, in S2Point b, in S2Point c) {
    _a = a;
    _b = b;
    _c = c;
    _crosser = new S2EdgeCrosser(_a, _b, _c);
  }

  void restartAt(in S2Point c) {
    _c = c;
    _crosser.restartAt(_c);
  }

  int crossingSign(ref in S2Point d) {
    int result = _crosser.crossingSign(d);
    _c = d;
    _crosser.setC(_c);
    return result;
  }

  bool edgeOrVertexCrossing(in S2Point d) {
    bool result = _crosser.edgeOrVertexCrossing(d);
    _c = d;
    _crosser.setC(_c);
    return result;
  }

private:
  S2Point _a;
  S2Point _b;
  S2Point _c;
  // TODO(ericv): It would be more efficient to implement S2CopyingEdgeCrosser
  // directly rather than as a wrapper around S2EdgeCrosser.
  S2EdgeCrosser _crosser;
}
