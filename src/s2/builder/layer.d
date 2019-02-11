// Copyright 2016 Google Inc. All Rights Reserved.
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

module s2.builder.layer;

import s2.builder.graph;
import s2.s2builder;
import s2.s2error;

/**
 * This class is not needed by ordinary S2Builder clients.  It is only
 * necessary if you wish to implement a new S2Builder::Layer subtype.
 */
abstract class Layer {
public:
  // Convenience declarations for layer subtypes.
  alias EdgeType = S2Builder.EdgeType;
  alias GraphOptions = .GraphOptions;
  alias Graph = .Graph;
  alias Label = S2Builder.Label;
  alias LabelSetId = S2Builder.LabelSetId;

  /// Defines options for building the edge graph that is passed to Build().
  abstract GraphOptions graphOptions() const;

  // Assembles a graph of snapped edges into the geometry type implemented by
  // this layer.  If an error is encountered, sets "error" appropriately.
  //
  // Note that when there are multiple layers, the Graph objects passed to all
  // layers are guaranteed to be valid until the last Build() method returns.
  // This makes it easier to write algorithms that gather the output graphs
  // from several layers and process them all at once (such as
  // s2builderutil::ClosedSetNormalizer).
  abstract void build(in Graph g, ref S2Error error);
}
