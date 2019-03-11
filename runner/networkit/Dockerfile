FROM ubuntu:18.04

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
   git g++-7 ninja-build cmake \
   python3-pip python3-venv

RUN mkdir /runner /data
RUN cd /runner; git clone https://github.com/manpen/networkit.git . && git checkout benchmark_hyper
RUN cd /runner; ./build-benchmark.sh

VOLUME /data
CMD cd /data; . /runner/pyenv/bin/activate; python /runner/benchmark.py