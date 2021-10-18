FROM girgs_common

RUN mkdir /runner /data
COPY source /source
VOLUME /data

RUN cd /runner; \
 cmake -DOPTION_BUILD_EXAMPLES=ON -DOPTION_USE_BMI2=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native" -GNinja /source && \
 ninja bench_hyper

CMD cd /data; /runner/bench_hyper 2> hypergirgs_bench.csv
