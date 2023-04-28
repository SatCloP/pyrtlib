FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

ADD . /home/dev/pyrtlib

WORKDIR /home/dev/pyrtlib
RUN python3 setup.py install

ENTRYPOINT ["/bin/bash"]