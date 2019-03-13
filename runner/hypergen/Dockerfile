FROM girgs_common

RUN apt-get -y install scons

RUN mkdir /runner /data
RUN cd /runner; git clone https://github.com/manpen/hypergen.git . && git checkout bench_hypergirgs
RUN cd /runner/; ./compile.sh

CMD cd /data; LD_PRELOAD=libtbbmalloc_proxy.so.2 /runner/build_gcc/benchmark_hypergen 2> hypergen_bench.csv