FROM girgs_common

RUN mkdir /runner /data
RUN cd /runner; git clone https://github.com/chistopher/networkit-general-polylog.git . && git checkout benchmark_hyper
RUN cd /runner; ./build-benchmark.sh

VOLUME /data
CMD cd /data; . /runner/pyenv/bin/activate; python /runner/benchmark.py
