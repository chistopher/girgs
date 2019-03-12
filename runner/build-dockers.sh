#!/usr/bin/env bash

# the common image installs C++, Python and tools used by all other containers
docker build -t girgs_common common

docker build -t girgs_networkit networkit

cd girgs; git clone ../.. source
docker build -t girgs_girgs

cd hypergirgs; git clone ../.. source
docker build -t girgs_hypergirgs hypergirgs

docker build -t girgs_hypergen  hypergen
