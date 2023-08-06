FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev wget \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN python3 -m pip install cftime

# TODO: uncomment when repository will be public
WORKDIR /tmp
RUN wget https://github.com/SatCloP/pyrtlib/archive/refs/heads/main.zip && \
    unzip -a pyrtlib-main.zip && \
    cd pyrtlib-main
# ADD . /home/dev/pyrtlib
# WORKDIR /home/dev/pyrtlib
# RUN python3 setup.py install

RUN python3 setup.py install

ENTRYPOINT ["/bin/bash"]