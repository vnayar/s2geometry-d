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

module s2.s2lax_polygon_shape;

import s2.s2loop;
import s2.s2point;
import s2.s2polygon;
import s2.s2shape;
import s2.shapeutil.get_reference_point : getReferencePoint;

import std.algorithm : copy;
import std.range : assumeSorted, enumerate, SortedRange;

alias ChainPosition = S2Shape.ChainPosition;

import std.stdio;

/**
 * S2LaxPolygonShape represents a region defined by a collection of zero or
 * more closed loops.  The interior is the region to the left of all loops.
 * This is similar to S2Polygon::Shape except that this class supports
 * polygons with degeneracies.  Degeneracies are of two types: degenerate
 * edges (from a vertex to itself) and sibling edge pairs (consisting of two
 * oppositely oriented edges).  Degeneracies can represent either "shells" or
 * "holes" depending on the loop they are contained by.  For example, a
 * degenerate edge or sibling pair contained by a "shell" would be interpreted
 * as a degenerate hole.  Such edges form part of the boundary of the polygon.
 *
 * Loops with fewer than three vertices are interpreted as follows:
 *  - A loop with two vertices defines two edges (in opposite directions).
 *  - A loop with one vertex defines a single degenerate edge.
 *  - A loop with no vertices is interpreted as the "full loop" containing
 *    all points on the sphere.  If this loop is present, then all other loops
 *    must form degeneracies (i.e., degenerate edges or sibling pairs).  For
 *    example, two loops {} and {X} would be interpreted as the full polygon
 *    with a degenerate single-point hole at X.
 *
 * S2LaxPolygonShape does not have any error checking, and it is perfectly
 * fine to create S2LaxPolygonShape objects that do not meet the requirements
 * below (e.g., in order to analyze or fix those problems).  However,
 * S2LaxPolygonShapes must satisfy some additional conditions in order to
 * perform certain operations:
 *
 *  - In order to be valid for point containment tests, the polygon must
 *    satisfy the "interior is on the left" rule.  This means that there must
 *    not be any crossing edges, and if there are duplicate edges then all but
 *    at most one of thm must belong to a sibling pair (i.e., the number of
 *    edges in opposite directions must differ by at most one).
 *
 *  - To be valid for boolean operations (S2BooleanOperation), degenerate
 *    edges and sibling pairs cannot coincide with any other edges.  For
 *    example, the following situations are not allowed:
 *
 *      {AA, AA}      // degenerate edge coincides with another edge
 *      {AA, AB}      // degenerate edge coincides with another edge
 *      {AB, BA, AB}  // sibling pair coincides with another edge
 *
 * Note that S2LaxPolygonShape is must faster to initialize and is more
 * compact than S2Polygon, but unlike S2Polygon it does not have any built-in
 * operations.  Instead you should use S2ShapeIndex operations
 * (S2BooleanOperation, S2ClosestEdgeQuery, etc).
 */
class S2LaxPolygonShape : S2Shape {
public:
  // Constructs an empty polygon.
  this() {
    _numLoops = 0;
    _numVertices = 0;
  }

  // Constructs an S2LaxPolygonShape from the given vertex loops.
  alias Loop = S2Point[];

  this(in Loop[] loops) {
    initialize(loops);
  }

  /**
   * Constructs an S2LaxPolygonShape from an S2Polygon, by copying its data.
   * Full and empty S2Polygons are supported.
   */
  this(in S2Polygon polygon) {
    initialize(polygon);
  }

  /// Initializes an S2LaxPolygonShape from the given vertex loops.
  void initialize(in S2Point[][] loops) {
    writeln("S2LaxPolygonShape 0: loops.length=", loops.length);
    _numLoops = cast(int) loops.length;
    if (_numLoops == 0) {
      _numVertices = 0;
      _vertices = null;
    } else if (_numLoops == 1) {
      _numVertices = cast(int) loops[0].length;
      _vertices = loops[0].dup;
    } else {
      writeln("S2LaxPolygonShape 1: _numLoops=", _numLoops);
      _cumulativeVertices = new int[_numLoops + 1];
      int num_vertices = 0;
      for (int i = 0; i < _numLoops; ++i) {
        _cumulativeVertices[i] = num_vertices;
        num_vertices += loops[i].length;
      }
      _cumulativeVertices[_numLoops] = num_vertices;
      _vertices = new S2Point[num_vertices];
      writeln("S2LaxPolygonShape 2: _cumulativeVertices=", _cumulativeVertices);
      for (int i = 0; i < _numLoops; ++i) {
        writeln("S2LaxPolygonShape 3: _vertices=", _vertices);
        copy(loops[i], _vertices[_cumulativeVertices[i] .. $]);
      }
    }
  }

  /**
   * Initializes an S2LaxPolygonShape from an S2Polygon, by copying its data.
   * Full and empty S2Polygons are supported.
   */
  void initialize(in S2Polygon polygon) {
    writeln("S2LaxPolygonShape.initialize >");
    const(S2Point[])[] spans;
    for (int i = 0; i < polygon.numLoops(); ++i) {
      writeln("S2LaxPolygonShape.initialize 1:");
      const S2Loop loop = polygon.loop(i);
      if (loop.isFull()) {
        spans ~= new S2Point[0];  // Empty span.
      } else {
        spans ~= loop.vertices();
      }
    }
    initialize(spans);
  }

  // Returns the number of loops.
  int numLoops() const {
    return _numLoops;
  }

  // Returns the total number of vertices in all loops.
  int numVertices() const {
    if (numLoops() <= 1) {
      return _numVertices;
    } else {
      return _cumulativeVertices[numLoops()];
    }
  }

  // Returns the number of vertices in the given loop.
  int numLoopVertices(int i) const
  in {
    assert(i < numLoops());
  } body {
    if (numLoops() == 1) {
      return _numVertices;
    } else {
      return _cumulativeVertices[i + 1] - _cumulativeVertices[i];
    }
  }

  // Returns the vertex from loop "i" at index "j".
  S2Point loopVertex(int i, int j) const
  in {
    assert(i < numLoops());
    assert(j < numLoopVertices(i));
  } body {
    if (numLoops() == 1) {
      return _vertices[j];
    } else {
      return _vertices[_cumulativeVertices[i] + j];
    }
  }

  const(S2Point[]) loopVertices(int i) const
  in {
    assert(i < numLoops());
  } body {
    if (numLoops() == 1) {
      return _vertices;
    } else {
      size_t indexStart = _cumulativeVertices[i];
      return _vertices[indexStart .. indexStart + numLoopVertices(i)];
    }
  }

  // S2Shape interface:
  final override
  int numEdges() const {
    return numVertices();
  }

  final override
  Edge edge(int e0) const
  in {
    assert(e0 < numEdges());
  } body {
    int e1 = e0 + 1;
    if (numLoops() == 1) {
      if (e1 == _numVertices) { e1 = 0; }
    } else {
      // Find the index of the first vertex of the loop following this one.
      const int kMaxLinearSearchLoops = 12;  // From benchmarks.
      //int* next = _cumulativeVertices + 1;
      size_t nextIndex = 1;
      if (numLoops() <= kMaxLinearSearchLoops) {
        while (_cumulativeVertices[nextIndex] <= e0) ++nextIndex;
      } else {
        //next = std::lower_bound(next, next + num_loops(), e1);
        nextIndex += _cumulativeVertices[nextIndex .. nextIndex + numLoops()]
            .assumeSorted
            .lowerBound(e1)
            .length;
      }
      // Wrap around to the first vertex of the loop if necessary.
      if (e1 == _cumulativeVertices[nextIndex]) { e1 = _cumulativeVertices[nextIndex - 1]; }
    }
    writeln("S2LaxPolygonShape 0: e0=", e0, ", e1=", e1);
    return Edge(_vertices[e0], _vertices[e1]);
  }

  final override
  int dimension() const {
    return 2;
  }

  final override
  ReferencePoint getReferencePoint() const {
    return .getReferencePoint(this);
  }

  final override
  int numChains() const {
    return numLoops();
  }

  final override
  Chain chain(int i) const
  in {
    assert(i < numLoops());
  } body {
    if (numLoops() == 1) {
      return Chain(0, _numVertices);
    } else {
      int start = _cumulativeVertices[i];
      return Chain(start, _cumulativeVertices[i + 1] - start);
    }
  }

  final override
  Edge chainEdge(int i, int j) const
  in {
    assert(i < numLoops());
    assert(j < numLoopVertices(i));
  } body {
    int n = numLoopVertices(i);
    int k = (j + 1 == n) ? 0 : j + 1;
    if (numLoops() == 1) {
      return Edge(_vertices[j], _vertices[k]);
    } else {
      int base = _cumulativeVertices[i];
      return Edge(_vertices[base + j], _vertices[base + k]);
    }
  }

  final override
  ChainPosition chainPosition(int e) const
  in {
    assert(e < numEdges());
  } body {
    const int kMaxLinearSearchLoops = 12;  // From benchmarks.
    if (numLoops() == 1) {
      return ChainPosition(0, e);
    } else {
      // Find the index of the first vertex of the loop following this one.
      //int* nextIndex = _cumulativeVertices + 1;
      int nextIndex = 1;
      if (numLoops() <= kMaxLinearSearchLoops) {
        while (_cumulativeVertices[nextIndex] <= e) ++nextIndex;
      } else {
        nextIndex = cast(int) _cumulativeVertices[nextIndex .. $]
            .assumeSorted!("a <= b")
            .lowerBound(e)
            .length + 1;
      }
      return ChainPosition(nextIndex - 1, e - _cumulativeVertices[nextIndex - 1]);
    }
  }

private:

  int _numLoops;
  S2Point[] _vertices;
  // If num_loops_ <= 1, this union stores the number of vertices.
  // Otherwise it points to an array of size (num_loops + 1) where element "i"
  // is the total number of vertices in loops 0..i-1.
  union {
    int _numVertices;
    int[] _cumulativeVertices;  // Don't use unique_ptr in unions.
  }
}
