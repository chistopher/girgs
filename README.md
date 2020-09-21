
[![Build Status](https://travis-ci.org/chistopher/girgs.svg)](https://travis-ci.org/chistopher/girgs)

This project contains an implementation of a linear-time sampling algorithm 
for geometric inhomogeneous random graphs (GIRG)
with a special case implementation for hyperbolic random graphs (HRG).
Details on the algorithms as well as the models can be found in our corresponding paper:
Efficiently Generating Geometric Inhomogeneous and Hyperbolic Random Graphs ([conference version](http://dx.doi.org/10.4230/LIPIcs.ESA.2019.21), [full version](https://arxiv.org/abs/1905.06706)).
To cite this work, please use the bibtex entry for our conference version.

We provide a C++ library and a command line client for each of the two generators.

# Installation

## Ubuntu 

The most recent version of the generators can be installed from our personal package archive ([PPA link](https://launchpad.net/~chistopher/+archive/ubuntu/ppa)).
```
sudo add-apt-repository ppa:chistopher/ppa
sudo apt-get update
sudo apt-get install libgirgs-all 
```

The `libgirgs-all` package depends on all sub-packages: 
- the runtime `libgirgs`
- the cli's `libgirgs-cli`
- the headers and cmake files `libgirgs-dev`
- and the debug runtime `libgirgs-dbg`

If you installed `libgirgs-cli` you can verify the install by running `genhrg --version`.
The result should look something like this.
```
HyperGIRGs command line interface.

girgs v1.0.0 (2dd23aad456e)
Generator for Geometric Inhomogeneous Random Graphs
Hasso Plattner Institute
https://github.com/chistopher/girgs/
christopher.weyand@hpi.de
```
 
## Windows

A Win64 installer for each release can be downloaded from the release page.
It features similar components as the ubuntu PPA.
If a more recent version is desired, the project can be build and installed from source. 

## macOS

A packaged build is not available. Feel free to open a pull request.
However, building from source should work just as detailed below.
We will add osx to the travis builds as soon as their images have CMake >=3.12.

## NetworKit

An integration into the [NetworKit](https://networkit.github.io) framework is planned.

## Build from Source

To build the project from source you need
- CMake 3.9
- C++11
- OpenMP
- OPTIONAL: CPU with BMI2 instruction set

The optional development components use
- [Google Test](https://github.com/google/googletest)
- [Google Benchmark](https://github.com/google/benchmark)
- [Doxygen](https://github.com/google/benchmark), which uses [Graphviz](https://www.graphviz.org/) and latex
- [Boost](https://www.boost.org/)

The simplest way to build the project is
```
mkdir build
cd build
cmake ..
make
```

Otherwise, our `configure` script provides reasonable defaults
for release, debug, and packaging builds.
For example, to install the project from source use
```
./configure
./configure pack
cmake --build ./build/
cmake --build ./build/ --target install
```
Specifying the `pack` target instead of `install` generates a debian package that can be installed via aptitude or dpkg.
 
# CLI

The project contains two CLI's: `gengirg` and `genhrg`.
You can use `--version` to determine the version of your respective generator.

The GIRG generator features the following input parameters.
```
./gengirg --help
usage: ./gengirg
		[-n anInt]          // number of nodes                          default 10000 
        	[-d anInt]          // dimension of geometry    range [1,5]     default 1
		[-ple aFloat]       // power law exponent       range (2,3]     default 2.5
		[-alpha aFloat]     // model parameter          range (1,inf]   default infinity
		[-deg aFloat]       // average degree           range [1,n)     default 10
		[-wseed anInt]      // weight seed                              default 12
		[-pseed anInt]      // position seed                            default 130
		[-sseed anInt]      // sampling seed                            default 1400
		[-threads anInt]    // number of threads to use                 default 1
		[-file aString]     // file name for output (w/o ext)           default "graph"
		[-dot 0|1]          // write result as dot (.dot)               default 0
		[-edge 0|1]         // write result as edgelist (.txt)          default 0
```

The HRG generator features the following input parameters.
```
./genhrg --help
usage: ./genhrg
		[-n anInt]          // number of nodes                          default 10000
		[-alpha aFloat]     // model parameter          range [0.5,1]   default 0.75
		[-t aFloat]         // temperature parameter    range [0,1)     default 0
		[-deg aFloat]       // average degree           range [1,n)     default 10
		[-rseed anInt]      // radii seed                               default 12
		[-aseed anInt]      // angle seed                               default 130
		[-sseed anInt]      // sampling seed                            default 1400
		[-threads anInt]    // number of threads to use                 default 1
		[-nkr 0|1]          // use NetworKit R estimation               default 0
		[-file aString]     // file name for output (w/o ext)           default "graph"
		[-edge 0|1]         // write result as edgelist (.txt)          default 0
		[-coord 0|1]        // write hyp. coordinates (.hyp)            default 0
```

# C++ Library

The library is based in the [cmake-init](https://github.com/cginternals/cmake-init) project template.
Detailed information on the build system design can be found there.


For the standard use case, we provide a buffered edge list generation supporting parallel execution in the headers `hypergirgs/Generator.h` and `girgs/Generator.h`.
```cpp
#include <hypergirgs/Generator.h>

auto R = hypergirgs::calculateRadius(n, alpha, T, deg);
auto radii = hypergirgs::sampleRadii(n, alpha, R, rseed);
auto angles = hypergirgs::sampleAngles(n, aseed);
// generates edge list as vector<pair<int,int>>
auto hrg_edges = hypergirgs::generateEdges(radii, angles, T, R, sseed);
```

Internally, the algorithm is templated with a callback that is called for each emitted edge.
Using lambdas, a custom callback can be used as follows.
```cpp
#include <hypergirgs/HyperbolicTree.h>

auto callback = [] (int a, int b, int tid) { ... }; // tid is thread id
auto generator = hypergirgs::makeHyperbolicTree(radii, angles, T, R, callback);
generator.generate(seed);
```

For details we refer to our example applications in `source/examples/` or the CLI's in `source/cli/`.

