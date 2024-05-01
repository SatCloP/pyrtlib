FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3-pip python3-cftime python3-dev unzip git wget build-essential \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/SatCloP/pyrtlib.git && \
    cd pyrtlib && \
    python3 setup.py install

ENTRYPOINT ["/bin/bash"]