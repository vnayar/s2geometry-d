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

module s2.shapeutil.get_reference_point;

import s2.s2contains_vertex_query;
import s2.s2point;
import s2.s2shape;
import s2.s2shape_index;

import std.algorithm : sort;
import std.exception : enforce;

import std.stdio;

// This is a helper function for implementing S2Shape::GetReferencePoint().
//
// Given a shape consisting of closed polygonal loops, the interior of the
// shape is defined as the region to the left of all edges (which must be
// oriented consistently).  This function then chooses an arbitrary point and
// returns true if that point is contained by the shape.
//
// Unlike S2Loop and S2Polygon, this method allows duplicate vertices and
// edges, which requires some extra care with definitions.  The rule that we
// apply is that an edge and its reverse edge "cancel" each other: the result
// is the same as if that edge pair were not present.  Therefore shapes that
// consist only of degenerate loop(s) are either empty or full; by convention,
// the shape is considered full if and only if it contains an empty loop (see
// S2LaxPolygonShape for details).
//
// Determining whether a loop on the sphere contains a point is harder than
// the corresponding problem in 2D plane geometry.  It cannot be implemented
// just by counting edge crossings because there is no such thing as a "point
// at infinity" that is guaranteed to be outside the loop.
S2Shape.ReferencePoint getReferencePoint(in S2Shape shape)
in {
  assert(shape.hasInterior());
} body {
  writeln("getReferencePoint >");
  scope(exit) writeln("getReferencePoint <");
  if (shape.numEdges() == 0) {
    // A shape with no edges is defined to be "full" if and only if it
    // contains an empty loop.
    return S2Shape.ReferencePoint(shape.numChains() > 0);
  }
  // Define a "matched" edge as one that can be paired with a corresponding
  // reversed edge.  Define a vertex as "balanced" if all of its edges are
  // matched. In order to determine containment, we must find an unbalanced
  // vertex.  Often every vertex is unbalanced, so we start by trying an
  // arbitrary vertex.
  auto edge = shape.edge(0);
  S2Shape.ReferencePoint result;
  if (getReferencePointAtVertex(shape, edge.v0, result)) {
    return result;
  }
  // That didn't work, so now we do some extra work to find an unbalanced
  // vertex (if any).  Essentially we gather a list of edges and a list of
  // reversed edges, and then sort them.  The first edge that appears in one
  // list but not the other is guaranteed to be unmatched.
  int n = shape.numEdges();
  auto edges = new S2Shape.Edge[n];
  auto rev_edges = new S2Shape.Edge[n];
  for (int i = 0; i < n; ++i) {
    auto edge2 = shape.edge(i);
    edges[i] = edge2;
    rev_edges[i] = S2Shape.Edge(edge2.v1, edge2.v0);
  }
  sort(edges);
  sort(rev_edges);
  for (int i = 0; i < n; ++i) {
    if (edges[i] < rev_edges[i]) {  // edges[i] is unmatched
      enforce(getReferencePointAtVertex(shape, edges[i].v0, result));
      return result;
    }
    if (rev_edges[i] < edges[i]) {  // rev_edges[i] is unmatched
      enforce(getReferencePointAtVertex(shape, rev_edges[i].v0, result));
      return result;
    }
  }
  // All vertices are balanced, so this polygon is either empty or full.  By
  // convention it is defined to be "full" if it contains any empty loop.
  for (int i = 0; i < shape.numChains(); ++i) {
    if (shape.chain(i).length == 0) return S2Shape.ReferencePoint(true);
  }
  return S2Shape.ReferencePoint(false);
}

// This is a helper function for GetReferencePoint() below.
//
// If the given vertex "vtest" is unbalanced (see definition below), sets
// "result" to a ReferencePoint indicating whther "vtest" is contained and
// returns true.  Otherwise returns false.
private bool getReferencePointAtVertex(
    in S2Shape shape, in S2Point vtest, ref S2Shape.ReferencePoint result) {
  // Let P be an unbalanced vertex.  Vertex P is defined to be inside the
  // region if the region contains a particular direction vector starting from
  // P, namely the direction S2::Ortho(P).  This can be calculated using
  // S2ContainsVertexQuery.
  auto contains_query = new S2ContainsVertexQuery(vtest);
  int n = shape.numEdges();
  for (int e = 0; e < n; ++e) {
    auto edge = shape.edge(e);
    if (edge.v0 == vtest) contains_query.addEdge(edge.v1, 1);
    if (edge.v1 == vtest) contains_query.addEdge(edge.v0, -1);
  }
  int contains_sign = contains_query.containsSign();
  if (contains_sign == 0) {
    return false;  // There are no unmatched edges incident to this vertex.
  }
  result.point = vtest;
  result.contained = contains_sign > 0;
  return true;
}

