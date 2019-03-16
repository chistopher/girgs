FROM girgs_common

RUN mkdir /runner /data
COPY source /source
VOLUME /data

RUN cd /runner; \
 cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native" -GNinja /source && \
 ninja

CMD cd /data; . /runner/pyenv/bin/activate; python /runner/benchmark.py
