// Copyright 2012 Google Inc. All Rights Reserved.
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

module s2.mutable_s2shape_index_test;

import s2.logger;
import s2.mutable_s2shape_index;
import s2.r2point;
import s2.r2rect;
import s2.s1angle;
import s2.s2cap;
import s2.s2cell;
import s2.s2cell_id;
import s2.s2cell_union;
import s2.s2debug;
import s2.s2edge_clipping : clipToPaddedFace, intersectsRect, INTERSECTS_RECT_ERROR_UV_DIST;
import s2.s2edge_crosser;
import s2.s2edge_vector_shape;
import s2.s2error;
import s2.s2loop;
import s2.s2point;
import s2.s2pointutil;
import s2.s2polygon;
import s2.s2polyline;
import s2.s2shape;
import s2.s2shape_index;
import s2.s2testing;
import s2.s2text_format;
import s2.shapeutil.contains_brute_force : containsBruteForce;
import s2.shapeutil.visit_crossing_edge_pairs;

import fluent.asserts;

import std.range;
import std.algorithm;
import std.array;
import core.thread : Thread;
import core.sync.condition : Condition;
import core.sync.mutex : Mutex;

// Verify that that every cell of the index contains the correct edges, and
// that no cells are missing from the index.  The running time of this
// function is quadratic in the number of edges.
void quadraticValidate(MutableS2ShapeIndex index) {
  import s2.s2text_format : toString;
  static int callNum = 1;

  // Iterate through a sequence of nonoverlapping cell ids that cover the
  // sphere and include as a subset all the cell ids used in the index.  For
  // each cell id, verify that the expected set of edges is present.

  // "min_cellid" is the first S2CellId that has not been validated yet.
  S2CellId min_cellid = S2CellId.begin(S2CellId.MAX_LEVEL);
  for (auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN); ;
       it.next()) {
    // Generate a list of S2CellIds ("skipped cells") that cover the gap
    // between the last cell we validated and the next cell in the index.
    auto skipped = new S2CellUnion();
    if (!it.done()) {
      S2CellId cellid = it.id();
      Assert.notLessThan(cellid, min_cellid);
      skipped.initFromBeginEnd(min_cellid, cellid.rangeMin());
      min_cellid = cellid.rangeMax().next();
    } else {
      // Validate the empty cells beyond the last cell in the index.
      skipped.initFromBeginEnd(min_cellid, S2CellId.end(S2CellId.MAX_LEVEL));
    }
    // Iterate through all the shapes, simultaneously validating the current
    // index cell and all the skipped cells.
    int short_edges = 0;  // number of edges counted toward subdivision
    for (int id = 0; id < index.numShapeIds(); ++id) {
      const(S2Shape) shape = index.shape(id);
      const(S2ClippedShape)* clipped = null;

      if (!it.done()) {
        clipped = it.cell().findClipped(id);
      }

      // First check that contains_center() is set correctly.
      foreach (S2CellId skipped_id; skipped.cellIds()) {
        validateInterior(shape, skipped_id, false);
      }
      if (!it.done()) {
        bool contains_center = clipped && clipped.containsCenter();
        validateInterior(shape, it.id(), contains_center);
      }
      // If this shape has been released, it should not be present at all.
      if (shape is null) {
        Assert.equal(clipped, null);
        continue;
      }
      // Otherwise check that the appropriate edges are present.
      for (int e = 0; e < shape.numEdges(); ++e) {
        auto edge = shape.edge(e);
        for (int j = 0; j < skipped.numCells(); ++j) {
          validateEdge(edge.v0, edge.v1, skipped.cellId(j), false);
        }
        if (!it.done()) {
          bool has_edge = clipped && clipped.containsEdge(e);
          validateEdge(edge.v0, edge.v1, it.id(), has_edge);
          int max_level = index.getEdgeMaxLevel(edge);
          if (has_edge && it.id().level() < max_level) {
            ++short_edges;
          }
        }
      }
    }
    Assert.notGreaterThan(short_edges, index.options().maxEdgesPerCell());
    if (it.done()) break;
  }
}

// Given an edge and a cell id, determine whether or not the edge should be
// present in that cell and verify that this matches "index_has_edge".
void validateEdge(in S2Point a, in S2Point b, S2CellId id, bool index_has_edge) {
  // Expand or shrink the padding slightly to account for errors in the
  // function we use to test for intersection (IntersectsRect).
  double padding = MutableS2ShapeIndex.CELL_PADDING;
  padding += (index_has_edge ? 1 : -1) * INTERSECTS_RECT_ERROR_UV_DIST;
  R2Rect bound = id.getBoundUV().expanded(padding);
  R2Point a_uv, b_uv;
  // Assert.equal(clipToPaddedFace(a, b, id.face(), padding, a_uv, b_uv)
  //     && intersectsRect(a_uv, b_uv, bound),
  //     index_has_edge);
  assert(
      (clipToPaddedFace(a, b, id.face(), padding, a_uv, b_uv) && intersectsRect(a_uv, b_uv, bound))
      == index_has_edge);
}

// Given a shape and a cell id, determine whether or not the shape contains
// the cell center and verify that this matches "index_contains_center".
void validateInterior(in S2Shape shape, S2CellId id, bool index_contains_center) {
  static int callNum = 1;

  if (shape is null) {
    assert(index_contains_center == false);
  } else {
    assert(containsBruteForce(shape, id.toS2Point()) == index_contains_center);
  }
}

void checkIteratorMethods(MutableS2ShapeIndex index) {
  auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
  Assert.equal(it.prev(), false);
  it.finish();
  Assert.equal(it.done(), true);
  S2CellId[] ids;
  auto it2 = new MutableS2ShapeIndex.Iterator(index);
  S2CellId min_cellid = S2CellId.begin(S2CellId.MAX_LEVEL);
  for (it.begin(); !it.done(); it.next()) {
    S2CellId cellid = it.id();
    auto skipped = S2CellUnion.fromBeginEnd(min_cellid, cellid.rangeMin());
    foreach (S2CellId skipped_id; skipped.cellIds) {
      Assert.equal(it2.locate(skipped_id.toS2Point()), false);
      Assert.equal(S2ShapeIndex.CellRelation.DISJOINT, it2.locate(skipped_id));
      it2.begin();
      it2.seek(skipped_id);
      Assert.equal(cellid, it2.id());
    }
    if (!ids.empty()) {
      it2.copy(it);
      Assert.equal(it2.prev(), true);
      Assert.equal(ids.back(), it2.id());
      it2.next();
      Assert.equal(cellid, it2.id());
      it2.seek(ids.back());
      Assert.equal(ids.back(), it2.id());
    }
    it2.begin();
    Assert.equal(cellid.toS2Point(), it.center());
    Assert.equal(it2.locate(it.center()), true);
    Assert.equal(cellid, it2.id());
    it2.begin();
    Assert.equal(S2ShapeIndex.CellRelation.INDEXED, it2.locate(cellid));
    Assert.equal(cellid, it2.id());
    if (!cellid.isFace()) {
      it2.begin();
      Assert.equal(S2ShapeIndex.CellRelation.SUBDIVIDED, it2.locate(cellid.parent()));
      Assert.notGreaterThan(it2.id(), cellid);
      Assert.notLessThan(it2.id(), cellid.parent().rangeMin());
    }
    if (!cellid.isLeaf()) {
      for (int i = 0; i < 4; ++i) {
        it2.begin();
        Assert.equal(S2ShapeIndex.CellRelation.INDEXED, it2.locate(cellid.child(i)));
        Assert.equal(cellid, it2.id());
      }
    }
    ids ~= cellid;
    min_cellid = cellid.rangeMax().next();
  }
}

@("MutableS2ShapeIndexTest.SpaceUsed") unittest {
  auto index = new MutableS2ShapeIndex();
  index.add(new S2EdgeVectorShape(S2Point(1, 0, 0), S2Point(0, 1, 0)));
  Assert.equal(index.isFresh(), false);
  size_t size_before = index.spaceUsed();
  Assert.equal(index.isFresh(), false);

  quadraticValidate(index);
  size_t size_after = index.spaceUsed();

  Assert.equal(index.isFresh(), true);

  // TODO: After size is implemented.
  //EXPECT_TRUE(size_after > size_before);
}

@("MutableS2ShapeIndexTest.NoEdges") unittest {
  auto index = new MutableS2ShapeIndex();
  auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
  Assert.equal(it.done(), true);
  checkIteratorMethods(index);
}

@("MutableS2ShapeIndexTest.OneEdge") unittest {
  auto index = new MutableS2ShapeIndex();
  Assert.equal(
      index.add(new S2EdgeVectorShape(S2Point(1, 0, 0), S2Point(0, 1, 0))), 0);
  quadraticValidate(index);
  checkIteratorMethods(index);
}

@("MutableS2ShapeIndexTest.ShrinkToFitOptimization") unittest {
  // This used to trigger a bug in the ShrinkToFit optimization.  The loop
  // below contains almost all of face 0 except for a small region in the
  // 0/00000 subcell.  That subcell is the only one that contains any edges.
  // This caused the index to be built only in that subcell.  However, all the
  // other cells on that face should also have index entries, in order to
  // indicate that they are contained by the loop.
  auto index = new MutableS2ShapeIndex();
  S2Loop loop = S2Loop.makeRegularLoop(
      S2Point(1, 0.5, 0.5).normalize(), S1Angle.fromDegrees(89.0), 100);
  index.add(new S2Loop.Shape(loop));
  quadraticValidate(index);
}

@("MutableS2ShapeIndexTest.LoopsSpanningThreeFaces") unittest {
  auto index = new MutableS2ShapeIndex();
  auto polygon = new S2Polygon();
  const int kNumEdges = 100;  // Validation is quadratic
  // Construct two loops consisting of kNumEdges vertices each, centered
  // around the cube vertex at the start of the Hilbert curve.
  S2Testing.concentricLoopsPolygon(
      S2Point(1, -1, -1).normalize(), 2, kNumEdges, polygon);
  S2Loop[] loops = polygon.release();
  foreach (loop; loops) {
    index.add(new S2Loop.Shape(loop));
  }
  quadraticValidate(index);
  checkIteratorMethods(index);
}

@("MutableS2ShapeIndexTest.ManyIdenticalEdges") unittest {
  auto index = new MutableS2ShapeIndex();
  const int kNumEdges = 100;  // Validation is quadratic
  auto a = S2Point(0.99, 0.99, 1).normalize();
  auto b = S2Point(-0.99, -0.99, 1).normalize();
  for (int i = 0; i < kNumEdges; ++i) {
    Assert.equal(i, index.add(new S2EdgeVectorShape(a, b)));
  }
  quadraticValidate(index);
  checkIteratorMethods(index);
  // Since all edges span the diagonal of a face, no subdivision should
  // have occurred (with the default index options).
  for (auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
       !it.done(); it.next()) {
    Assert.equal(0, it.id().level());
  }
}

@("MutableS2ShapeIndexTest.DegenerateEdge") unittest {
  auto index = new MutableS2ShapeIndex();
  // This test verifies that degenerate edges are supported.  The following
  // point is a cube face vertex, and so it should be indexed in 3 cells.
  S2Point a = S2Point(1, 1, 1).normalize();
  auto shape = new S2EdgeVectorShape();
  shape.add(a, a);
  index.add(shape);
  quadraticValidate(index);
  // Check that exactly 3 index cells contain the degenerate edge.
  int count = 0;
  for (auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
       !it.done(); it.next(), ++count) {
    Assert.equal(it.id().isLeaf(), true);
    Assert.equal(it.cell().numClipped(), 1);
    Assert.equal(it.cell().clipped(0).numEdges(), 1);
  }
  Assert.equal(count, 3);
}

@("MutableS2ShapeIndexTest.ManyTinyEdges") unittest {
  auto index = new MutableS2ShapeIndex();
  // This test adds many edges to a single leaf cell, to check that
  // subdivision stops when no further subdivision is possible.
  enum int kNumEdges = 100;  // Validation is quadratic
  // Construct two points in the same leaf cell.
  S2Point a = S2CellId(S2Point(1, 0, 0)).toS2Point();
  S2Point b = (a + S2Point(0, 1e-12, 0)).normalize();
  auto shape = new S2EdgeVectorShape();
  for (int i = 0; i < kNumEdges; ++i) {
    shape.add(a, b);
  }
  index.add(shape);
  quadraticValidate(index);
  // Check that there is exactly one index cell and that it is a leaf cell.
  auto it = new MutableS2ShapeIndex.Iterator(index, S2ShapeIndex.InitialPosition.BEGIN);
  Assert.equal(!it.done(), true);
  Assert.equal(it.id().isLeaf(), true);
  it.next();
  Assert.equal(it.done(), true);
}

@("MutableS2ShapeIndexTest.SimpleUpdates") unittest {
  auto index = new MutableS2ShapeIndex();
  // Add 5 loops one at a time, then release them one at a time,
  // validating the index at each step.
  auto polygon = new S2Polygon();
  S2Testing.concentricLoopsPolygon(S2Point(1, 0, 0), 5, 20, polygon);
  for (int i = 0; i < polygon.numLoops(); ++i) {
    index.add(new S2Loop.Shape(polygon.loop(i)));
    quadraticValidate(index);
  }
  for (int id = 0; id < polygon.numLoops(); ++id) {
    index.release(id);
    quadraticValidate(index);
  }
}

@("MutableS2ShapeIndexTest.RandomUpdates") unittest {
  auto index = new MutableS2ShapeIndex();

  // Allow the seed to be varied from the command line.
  S2Testing.rnd.reset(s2RandomSeed);

  // A few polylines.
  index.add(new S2Polyline.Shape(
      makePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6")));
  index.add(new S2Polyline.Shape(
      makePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6")));
  index.add(new S2Polyline.Shape(
      makePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6")));

  // A loop that used to trigger an indexing bug.
  index.add(new S2Loop.Shape(S2Loop.makeRegularLoop(
      S2Point(1, 0.5, 0.5).normalize(), S1Angle.fromDegrees(89.0), 20)));

  // Five concentric loops.
  auto polygon5 = new S2Polygon();
  S2Testing.concentricLoopsPolygon(S2Point(1, -1, -1).normalize(), 5, 20, polygon5);
  for (int i = 0; i < polygon5.numLoops(); ++i) {
    index.add(new S2Loop.Shape(polygon5.loop(i)));
  }

  // Two clockwise loops around S2Cell cube vertices.
  index.add(new S2Loop.Shape(S2Loop.makeRegularLoop(
      S2Point(-1, 1, 1).normalize(), S1Angle.fromRadians(M_PI - 0.001), 10)));
  index.add(new S2Loop.Shape(S2Loop.makeRegularLoop(
      S2Point(-1, -1, -1).normalize(), S1Angle.fromRadians(M_PI - 0.001), 10)));

  // A shape with no edges and no interior.
  index.add(new S2Loop.Shape(new S2Loop(S2Loop.empty())));

  // A shape with no edges that covers the entire sphere.
  index.add(new S2Loop.Shape(
      new S2Loop(S2Loop.full())));

  S2Shape[] released;
  int[] added = iota(0, index.numShapeIds()).array;
  quadraticValidate(index);
  for (int iter = 0; iter < 100; ++iter) {
    logger.logDebug("Iteration: ", iter);
    // Choose some shapes to add and release.
    int num_updates = 1 + S2Testing.rnd.skewed(5);
    for (int n = 0; n < num_updates; ++n) {
      if (S2Testing.rnd.oneIn(2) && !added.empty()) {
        int i = S2Testing.rnd.uniform(cast(int) added.length);
        logger.logDebug("  Released shape ", added[i], " (", index.shape(added[i]), ")");
        released ~= index.release(added[i]);
        added = added.remove(i);
      } else if (!released.empty()) {
        int i = S2Testing.rnd.uniform(cast(int) released.length);
        S2Shape shape = released[i];
        index.add(released[i]);  // Changes shape->id().
        released = released.remove(i);
        added ~= shape.id();
        logger.logDebug("  Added shape ", shape.id(), " (", shape, ")");
      }
    }
    quadraticValidate(index);
  }
}

// Return true if any loop crosses any other loop (including vertex crossings
// and duplicate edges), or any loop has a self-intersection (including
// duplicate vertices).
private bool hasSelfIntersection(MutableS2ShapeIndex index) {
  S2Error error;
  if (findSelfIntersection(index, error)) {
    logger.logDebug(error);
    return true;
  }
  return false;
}

// This function recursively verifies that HasCrossing returns the given
// result for all possible cyclic permutations of the loop vertices for the
// given set of loops.
void checkHasCrossingPermutations(S2Loop[] loops, int i, bool has_crossing) {
  static int callNum = 1;
  if (i == loops.length) {
    auto index = new MutableS2ShapeIndex();
    auto polygon = new S2Polygon(loops);
    index.add(new S2Polygon.Shape(polygon));
    Assert.equal(has_crossing, hasSelfIntersection(index));
    loops = polygon.release();
  } else {
    S2Loop orig_loop = loops[i];
    for (int j = 0; j < orig_loop.numVertices(); ++j) {
      S2Point[] vertices;
      for (int k = 0; k < orig_loop.numVertices(); ++k) {
        vertices ~= orig_loop.vertex(j + k);
      }
      loops[i] = new S2Loop(vertices);
      checkHasCrossingPermutations(loops, i+1, has_crossing);
    }
    loops[i] = orig_loop;
  }
}

// Given a string reprsenting a polygon, and a boolean indicating whether this
// polygon has any self-intersections or loop crossings, verify that all
// hasSelfIntersection returns the expected result for all possible cyclic
// permutations of the loop vertices.
void checkHasCrossing(string polygon_str, bool has_crossing) {
  flagsS2Debug = false; // Allow invalid polygons (restore on exit).
  scope(exit) flagsS2Debug = true;

  auto polygon = makePolygonOrDie(polygon_str);
  S2Loop[] loops = polygon.release();
  checkHasCrossingPermutations(loops, 0, has_crossing);
}

@("MutableS2ShapeIndexTest.HasCrossing") unittest {
  // Coordinates are (lat,lng), which can be visualized as (y,x).
  checkHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
  checkHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true);  // duplicate vertex
  checkHasCrossing("0:0, 0:1, 1:0, 1:1", true);  // edge crossing
  checkHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true);  // duplicate edge
  checkHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true);  // reversed edge
  checkHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0", true);  // vertex crossing
}

// A test that repeatedly updates "index_" in one thread and attempts to
// concurrently read the index_ from several other threads.  When all threads
// have finished reading, the first thread makes another update.
//
// Note that we only test concurrent read access, since MutableS2ShapeIndex
// requires all updates to be single-threaded and not concurrent with any
// reads.
class LazyUpdatesTest {
 public:
  this() {
    _index = new MutableS2ShapeIndex();
    _lock = new Mutex();
    _updateReady = new Condition(_lock);
    _allReadersDone = new Condition(_lock);

    _numUpdates = 0;
    _numReadersLeft = 0;
  }

  // The function executed by each reader thread.
  void readerThread() {
    _lock.lock();
    for (int last_update = 0; ; last_update = _numUpdates) {
      while (_numUpdates == last_update) {
        _updateReady.wait();
      }
      if (_numUpdates < 0) break;

      // The index is built on demand the first time we attempt to use it.
      // We intentionally release the lock so that many threads have a chance
      // to access the MutableS2ShapeIndex in parallel.
      _lock.unlock();
      for (auto it = new MutableS2ShapeIndex.Iterator(_index, S2ShapeIndex.InitialPosition.BEGIN);
           !it.done(); it.next()) {
        continue;
      }
      _lock.lock();
      if (--_numReadersLeft == 0) {
        _allReadersDone.notify();
      }
    }
    _lock.unlock();
  }

  class ReaderThreadPool {
  public:
    this(int num_threads) {
      _threads = new Thread[](num_threads);
      for (int i = 0; i < _threads.length; ++i) {
        _threads[i] = new Thread(&LazyUpdatesTest.readerThread);
      }
    }
    ~this() {
      foreach (thread; _threads) {
        thread.join();
      }
    }

   private:
    Thread[] _threads;
  }

  MutableS2ShapeIndex _index;
  // The following fields are guarded by lock_.
  Mutex _lock;
  int _numUpdates;
  int _numReadersLeft;

  // Signalled when a new update is ready to be processed.
  Condition _updateReady;
  // Signalled when all readers have processed the latest update.
  Condition _allReadersDone;
}

/+ TODO: Implement when S2Loop is implemented.
@("LazyUpdatesTest.ConstMethodsThreadSafe") unittest {
  auto test = new LazyUpdatesTest();
  // Ensure that lazy updates are thread-safe.  In other words, make sure that
  // nothing bad happens when multiple threads call "const" methods that
  // cause pending updates to be applied.

  // The number of readers should be large enough so that it is likely that
  // several readers will be running at once (with a multiple-core CPU).
  const int kNumReaders = 8;
  auto pool = test.new ReaderThreadPool(kNumReaders);
  test._lock.lock();
  const int kIters = 100;
  for (int iter = 0; iter < kIters; ++iter) {
    // Loop invariant: lock_ is held and num_readers_left_ == 0.
    Assert.equal(test._numReadersLeft, 0);
    // Since there are no readers, it is safe to modify the index.
    test._index.clear();
    int num_vertices = 4 * S2Testing.rnd.skewed(10);  // Up to 4K vertices
    S2Loop loop(S2Loop::MakeRegularLoop(
        S2Testing::RandomPoint(), S2Testing::KmToAngle(5), num_vertices));
    index_.Add(make_unique<S2Loop::Shape>(loop.get()));
    num_readers_left_ = kNumReaders;
    ++num_updates_;
    update_ready_.SignalAll();
    while (num_readers_left_ > 0) {
      all_readers_done_.Wait(&lock_);
    }
  }
  // Signal the readers to exit.
  num_updates_ = -1;
  update_ready_.SignalAll();
  lock_.Unlock();
  // ReaderThreadPool destructor waits for all threads to complete.
}
+/

@("MutableS2ShapeIndex.MixedGeometry") unittest {
  // This test used to trigger a bug where the presence of a shape with an
  // interior could cause shapes that don't have an interior to suddenly
  // acquire one.  This would cause extra S2ShapeIndex cells to be created
  // that are outside the bounds of the given geometry.
  S2Polyline[] polylines;
  polylines ~= makePolylineOrDie("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6");
  polylines ~= makePolylineOrDie("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6");
  polylines ~= makePolylineOrDie("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6");
  auto index = new MutableS2ShapeIndex();
  foreach (polyline; polylines) {
    index.add(new S2Polyline.Shape(polyline));
  }
  auto loop = new S2Loop(new S2Cell(S2CellId.begin(S2CellId.MAX_LEVEL)));
  index.add(new S2Loop.Shape(loop));
  auto it = new MutableS2ShapeIndex.Iterator(index);
  // No geometry intersects face 1, so there should be no index cells there.
  Assert.equal(it.locate(S2CellId.fromFace(1)), S2ShapeIndex.CellRelation.DISJOINT);
}

@("S2Shape.user_data") unittest {
  struct MyData {
    int x, y;
  }
  class MyEdgeVectorShape : S2EdgeVectorShape {
  public:
    this(in MyData data) {
      super();
      _data = data;
    }

    override
    const(void*) userData() const { return &_data; }

    override
    void* mutableUserData() { return &_data; }

  private:
    MyData _data;
  }
  auto shape = new MyEdgeVectorShape(MyData(3, 5));
  MyData* data = cast(MyData*)(shape.mutableUserData());
  Assert.equal(data.x, 3);
  data.y = 10;
  Assert.equal((cast(const(MyData*)) shape.userData()).y, 10);
}
