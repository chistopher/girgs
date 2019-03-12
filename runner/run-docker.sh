#!/usr/bin/env bash
DATADIR=$(pwd)/data
mkdir -p $DATADIR

DOCKER_CMD="docker run --rm -it -u $(id -u ${USER}):$(id -g ${USER})"

$DOCKER_CMD -v $DATADIR:/data -t girgs_networkit
$DOCKER_CMD -v $DATADIR:/data -t girgs_girgs
$DOCKER_CMD -v $DATADIR:/data -t girgs_hypergirgs
