FROM girgs_common

RUN mkdir /runner /data
RUN cd /runner; git clone https://github.com/manpen/hypergen.git . && git checkout benchmark_girg
RUN cd /runner/; ./build-benchmark.sh

CMD cd /data; /runner/run-benchmark.sh