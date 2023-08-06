FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev unzip wget \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN python3 -m pip install cftime

# WORKDIR /tmp
RUN wget https://github.com/SatCloP/pyrtlib/archive/refs/heads/main.zip && \
    unzip -a main.zip && \
    cd pyrtlib-main && \
    python3 setup.py install

ENTRYPOINT ["/bin/bash"]