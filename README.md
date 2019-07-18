# S2 Geometry Library

## Overview

This project is a migration of Google's S2 Geometric Library from C++ to the D programming language.

Converted from [commit cdf1a05](https://github.com/google/s2geometry/commit/cdf1a05acf462ad60e1110a354c86e0b606c8931) on February 2nd 2018.

The S2 Geometric Library is a package for manipulating geometric shapes. Unlike many geometry
libraries, S2 is primarily designed to work with _spherical geometry_, i.e., shapes drawn on a
sphere rather than on a planar 2D map. This makes it especially suitable for working with geographic
data.

If you want to learn more about the library, start by reading the
[overview](http://s2geometry.io/about/overview) and [quick start
document](http://s2geometry.io/devguide/cpp/quickstart), then read the
introduction to the [basic types](http://s2geometry.io/devguide/basic_types).

S2 documentation can be found on [s2geometry.io](http://s2geometry.io).

## Requirements for End Users

* A [D compiler](https://dlang.org/).
* The [dub](https://code.dlang.org/download) build tool.

## Build and Install

You must [clone the git repository](https://help.github.com/articles/cloning-a-repository/).

```sh
cd [parent of directory where you want to put S2]
git clone https://github.com/vnayar/s2geometry-d.git
cd s2geometry-d
```
### Testing

The unittests may be run using the following command:

```sh
dub test
```

To run a single test, use its name along with the parameters "-s" and "-d", to run in
single-threaded mode and to enable debug output respectively.

```sh
dub test -- -s -d s2.s2predicates_test
```

### Building

From the appropriate directory depending on how you got the source, run the command:

```sh
dub
```

### Generating Documentation

Documentation is generated from the [ddoc](https://dlang.org/spec/ddoc.html) formatted comments in
the code. The best looking output comes from using the tool `ddox` and can be generated with the
command:

```sh
dub build --build=ddox
```

A local server can be run to allow you to preview the documentation using the command:

```sh
dub run --build=ddox
```

## Disclaimer

This is not an official Google product.
