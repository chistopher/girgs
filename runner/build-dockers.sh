#!/usr/bin/env bash

# the common image installs C++, Python and tools used by all other containers
docker build -t girgs_common common

docker build -t girgs_networkit networkit

git clone .. girgs/source
docker build --no-cache  -t girgs_girgs girgs
rm -rf girgs/source

git clone .. hypergirgs/source
docker build --no-cache  -t girgs_hypergirgs hypergirgs
rm -rf hypergirgs/source

docker build --no-cache  -t girgs_hypergen  hypergen

docker build --no-cache  -t girgs_embedder  embedder

docker build --no-cache  -t girgs_nkbin     nkbin

