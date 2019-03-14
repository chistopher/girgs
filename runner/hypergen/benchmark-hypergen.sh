#!/bin/bash
cd /data
LD_PRELOAD=libtbbmalloc_proxy.so.2 /runner/build_gcc/benchmark_hypergen 2> hypergen_bench.csv