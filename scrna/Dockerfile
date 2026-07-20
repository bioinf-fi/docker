# syntax=docker/dockerfile:1.7

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# The bioinformatics packages come from Ubuntu 24.04 so the same Dockerfile
# builds natively on both linux/amd64 and linux/arm64.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        bzip2 \
        ca-certificates \
        coreutils \
        curl \
        file \
        fonts-dejavu-core \
        gzip \
        hisat2 \
        less \
        nano \
        pigz \
        procps \
        rna-star \
        samtools \
        sra-toolkit \
        subread \
        time \
        unzip \
        wget \
        x11-utils \
        xauth \
        xpdf \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

COPY scripts/show-tools /usr/local/bin/show-tools
COPY scripts/download-sra /usr/local/bin/download-sra
COPY scripts/download-paper /usr/local/bin/download-paper

RUN chmod 0755 \
        /usr/local/bin/show-tools \
        /usr/local/bin/download-sra \
        /usr/local/bin/download-paper \
    && install -d -m 1777 /tmp/workshop-home \
    && install -d -m 0777 /data \
    && command -v STAR \
    && command -v hisat2 \
    && command -v hisat2-build \
    && command -v samtools \
    && command -v featureCounts \
    && command -v prefetch \
    && command -v fasterq-dump \
    && command -v xpdf

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    HOME=/tmp/workshop-home

WORKDIR /data

LABEL org.opencontainers.image.title="scRNA-seq workshop" \
      org.opencontainers.image.description="Interactive STAR, HISAT2, SAMtools, featureCounts and SRA Toolkit workshop image with X11 PDF viewing"

CMD ["/bin/bash"]
