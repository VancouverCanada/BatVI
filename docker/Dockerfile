# syntax=docker/dockerfile:1.6
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV BATVI_HOME=/opt/batvi
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    bash \
    bedtools \
    build-essential \
    bwa \
    ca-certificates \
    default-jdk \
    default-jre-headless \
    gawk \
    libtool \
    make \
    ncbi-blast+ \
    perl \
    picard-tools \
    pkg-config \
    samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR ${BATVI_HOME}
COPY . ${BATVI_HOME}

RUN chmod +x \
    ${BATVI_HOME}/*.sh \
    ${BATVI_HOME}/bin/*.sh \
    ${BATVI_HOME}/msapipeline/*.sh \
    ${BATVI_HOME}/picard_samtofastq.sh \
    ${BATVI_HOME}/docker/entrypoint.sh \
    && bash ${BATVI_HOME}/build.sh

ENV PATH="${BATVI_HOME}/bin:${PATH}"

ENTRYPOINT ["/opt/batvi/docker/entrypoint.sh"]
CMD ["shell"]
