FROM girgs_common

RUN apt-get -y install libgflags-dev libgoogle-glog-dev libgsl-dev

RUN mkdir /runner /data
RUN cd /runner; \
   git clone https://github.com/manpen/hyperbolic-embedder.git . && \
   git checkout benchmark_hyper
RUN cd /runner; mkdir build; cd build; \
   cmake -DCMAKE_BUILD_TYPE="Release" -GNinja ..; \
   ninja

VOLUME /data
CMD cd /data; /runner/build/benchmark_embedder 2> embedder_bench.csv
