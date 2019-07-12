// Copyright 2011 Google Inc. All Rights Reserved.
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

//
// Given a sequence of S2Points assumed to be the center of level-k cells,
// compresses it into a stream using the following method:
// - decompose the points into (face, si, ti) tuples (see s2coords.h)
// - run-length encode the faces, combining face number and count into a
//     varint32.  See the Faces class in s2point_compression.cc.
// - right shift the (si, ti) to remove the part that's constant for all cells
//     of level-k.  The result is called the (pi, qi) space.
// - 2nd derivative encode the pi and qi sequences (linear prediction)
// - zig-zag encode all derivative values but the first, which cannot be
//     negative
// - interleave the zig-zag encoded values
// - encode the first interleaved value in a fixed length encoding
//     (varint would make this value larger)
// - encode the remaining interleaved values as varint64s, as the
//     derivative encoding should make the values small.
// In addition, provides a lossless method to compress a sequence of points even
// if some points are not the center of level-k cells. These points are stored
// exactly, using 3 double precision values, after the above encoded string,
// together with their index in the sequence (this leads to some redundancy - it
// is expected that only a small fraction of the points are not cell centers).
//
// Require that the encoder was constructed with the no-arg constructor, as
// Ensure() will be called to allocate space.

//
// To encode leaf cells, this requires 8 bytes for the first vertex plus
// an average of 3.8 bytes for each additional vertex, when computed on
// Google's geographic repository.

module s2.s2point_compression;

import s2.s2point;
import s2.util.coding.coder;

/**
 * The XYZ and face,si,ti coordinates of an S2Point and, if this point is equal
 * to the center of an S2Cell, the level of this cell (-1 otherwise).
 */
struct S2XYZFaceSiTi {
  S2Point xyz;
  int face;
  uint si;
  uint ti;
  int cellLevel;
}

/+ TODO

// Encode the points in the encoder, using an optimized compressed format for
// points at the center of a cell at 'level', plus 3 double values for the
// others.
void encodeS2PointsCompressed(ORangeT)(
    in S2XYZFaceSiTi[] points, int level, Encoder!ORangeT encoder) {
  int[2][] vertices_pi_qi = new int[2][points.length];
  int[] off_center;
  Faces faces;
  for (int i = 0; i < points.length; ++i) {
    faces.addFace(points[i].face);
    vertices_pi_qi[i][0] = SiTitoPiQi(points[i].si, level);
    vertices_pi_qi[i][1] = SiTitoPiQi(points[i].ti, level);
    if (points[i].cellLevel != level) {
      off_center ~= i;
    }
  }
  faces.encode(encoder);
  encodePointsCompressed(vertices_pi_qi, level, encoder);
  int num_off_center = off_center.length;
  encoder.ensure(Encoder::kVarintMax32 +
                  (Encoder::kVarintMax32 + sizeof(S2Point)) * num_off_center);
  encoder->put_varint32(num_off_center);
  DCHECK_GE(encoder->avail(), 0);
  for (int index : off_center) {
    encoder->put_varint32(index);
    encoder->putn(&points[index].xyz, sizeof(points[index].xyz));
    DCHECK_GE(encoder->avail(), 0);
  }
}

// Decode points encoded with S2EncodePointsCompressed. Requires that the
// level is the level that was used in S2EncodePointsCompressed. Ensures
// that the decoded points equal the encoded points. Returns true on success.
bool decodeS2PointsCompressed(Decoder* decoder, int level,
                              absl::Span<S2Point> points);

+/
