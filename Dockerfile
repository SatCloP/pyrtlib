FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3-pip python3-dev unzip wget build-essential \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && cd /usr/local/bin 

RUN python3 -m pip install --upgrade pip setuptools --break-system-packages
RUN python3 -m pip install cftime --break-system-packages

# WORKDIR /tmp
RUN wget https://github.com/SatCloP/pyrtlib/archive/refs/heads/main.zip && \
    unzip -a main.zip && \
    cd pyrtlib-main && \
    python3 setup.py install

ENTRYPOINT ["/bin/bash"]