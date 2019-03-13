#!/usr/bin/env bash
DATADIR=$(pwd)/data
mkdir -p $DATADIR

DOCKER_CMD="docker run --rm -it -u $(id -u ${USER}):$(id -g ${USER}) --hostname $(hostname)"

$DOCKER_CMD -v $DATADIR:/data -t girgs_hypergen
$DOCKER_CMD -v $DATADIR:/data -t girgs_hypergirgs
$DOCKER_CMD -v $DATADIR:/data -t girgs_networkit
$DOCKER_CMD -v $DATADIR:/data -t girgs_girgs
