FROM girgs_common

RUN apt-get -y install scons

RUN mkdir /runner /data
RUN cd /runner; git clone https://github.com/manpen/hypergen.git . && git checkout bench_hypergirgs
RUN cd /runner/; ./compile.sh

COPY *.sh /usr/bin/
