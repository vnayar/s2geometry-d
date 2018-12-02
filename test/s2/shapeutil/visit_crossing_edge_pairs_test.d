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

// Original author: ericv@google.com (Eric Veach)
// Converted to D:  madric@gmail.com (Vijay Nayar)

module s2.shapeutil.visit_crossing_edge_pairs_test;

import s2.mutable_s2shape_index;
import s2.s2crossing_edge_query : CrossingType;
import s2.s2edge_crossings : crossingSign;
import s2.s2edge_vector_shape;
import s2.s2latlng;
import s2.s2shape_index;
import s2.shapeutil.contains_brute_force;
import s2.shapeutil.edge_iterator;
import s2.shapeutil.shape_edge;
import s2.shapeutil.shape_edge_id;
import s2.shapeutil.visit_crossing_edge_pairs;

import std.algorithm : count, sort, uniq;
import std.array;
import std.range;

import fluent.asserts;

// A set of edge pairs within an S2ShapeIndex.
struct EdgePair {
  ShapeEdgeId first;
  ShapeEdgeId second;

  int opCmp(EdgePair v) const {
    int cmp = first.opCmp(v.first);
    return cmp != 0 ? cmp : second.opCmp(v.second);
  }
}

alias EdgePairVector = EdgePair[];

EdgePairVector getCrossings(S2ShapeIndex index, CrossingType type) {
  EdgePairVector edge_pairs;
  visitCrossingEdgePairs(
      index, type, (in ShapeEdge a, in ShapeEdge b, bool) {
        edge_pairs ~= EdgePair(a.id(), b.id());
        return true;  // Continue visiting.
      });
  if (edge_pairs.length > 1) {
    edge_pairs = edge_pairs.sort().uniq().array();
  }
  return edge_pairs;
}

EdgePairVector getCrossingEdgePairsBruteForce(S2ShapeIndex index, CrossingType type) {
  EdgePairVector result;
  int min_sign = (type == CrossingType.ALL) ? 0 : 1;
  for (auto a_iter = EdgeIterator(index); !a_iter.done(); a_iter.next()) {
    auto a = a_iter.edge();
    EdgeIterator b_iter = a_iter;
    for (b_iter.next(); !b_iter.done(); b_iter.next()) {
      auto b = b_iter.edge();
      if (crossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
        result ~= EdgePair(a_iter.shapeEdgeId(), b_iter.shapeEdgeId());
      }
    }
  }
  return result;
}

void checkGetCrossingEdgePairs(S2ShapeIndex index, CrossingType type) {
  import std.stdio;
  import std.conv : to;
  EdgePairVector expected = getCrossingEdgePairsBruteForce(index, type);
  EdgePairVector actual = getCrossings(index, type);
  if (actual != expected) {
    Assert.equal(actual.length, expected.length, "Unexpected edge pairs; see details below."
        ~ "\nExpected number of edge pairs: " ~ to!string(expected.length)
        ~ "\nActual number of edge pairs: " ~ to!string(actual.length));
    foreach (const edge_pair; expected) {
      if (count(actual, edge_pair) != 1) {
        writeln("Missing value: ", edge_pair);
      }
    }
    foreach (const edge_pair; actual) {
      if (count(expected, edge_pair) != 1) {
        writeln("Extra value: ", edge_pair);
      }
    }
  }
}

@("GetCrossingEdgePairs.NoIntersections") unittest {
  auto index = new MutableS2ShapeIndex();
  checkGetCrossingEdgePairs(index, CrossingType.ALL);
  checkGetCrossingEdgePairs(index, CrossingType.INTERIOR);
}

@("GetCrossingEdgePairs.EdgeGrid") unittest {
  const int kGridSize = 10;  // (kGridSize + 1) * (kGridSize + 1) crossings
  auto index = new MutableS2ShapeIndex();
  auto shape = new S2EdgeVectorShape();
  for (int i = 0; i <= kGridSize; ++i) {
    shape.add(S2LatLng.fromDegrees(0, i).toS2Point(),
        S2LatLng.fromDegrees(kGridSize, i).toS2Point());
    shape.add(S2LatLng.fromDegrees(i, 0).toS2Point(),
        S2LatLng.fromDegrees(i, kGridSize).toS2Point());
  }
  index.add(shape);
  checkGetCrossingEdgePairs(index, CrossingType.ALL);
  checkGetCrossingEdgePairs(index, CrossingType.INTERIOR);
}
