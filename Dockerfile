FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends python3-pip python3-cftime python3-dev python3-venv unzip git wget build-essential \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"

RUN git clone https://github.com/SatCloP/pyrtlib.git && \
    cd pyrtlib && \
    pip install .

ENTRYPOINT ["/bin/bash"]