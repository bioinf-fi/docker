# syntax=docker/dockerfile:1.7

FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# The bioinformatics packages come from Ubuntu 24.04 so the same Dockerfile
# builds natively on both linux/amd64 and linux/arm64.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ssh \
        openssh-client \
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
        python3 \
        python3-pip \
        python3-venv \
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

# Isolated Python env for Jupyter and the workshop notebook kernel.
COPY requirements-notebook.txt /tmp/requirements-notebook.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements-notebook.txt \
    && /opt/venv/bin/python -m ipykernel install \
        --name=pydeseq2 \
        --display-name="pydeseq2" \
        --prefix=/usr/local \
    && rm -f /tmp/requirements-notebook.txt

COPY scripts/show-tools /usr/local/bin/show-tools
COPY scripts/download-sra /usr/local/bin/download-sra
COPY scripts/download-paper /usr/local/bin/download-paper
COPY scripts/start-notebook /usr/local/bin/start-notebook
COPY notebooks/ /opt/workshop/notebooks/

RUN chmod 0755 \
        /usr/local/bin/show-tools \
        /usr/local/bin/download-sra \
        /usr/local/bin/download-paper \
        /usr/local/bin/start-notebook \
    && install -d -m 1777 /tmp/workshop-home \
    && install -d -m 0777 /data \
    && command -v STAR \
    && command -v hisat2 \
    && command -v hisat2-build \
    && command -v samtools \
    && command -v featureCounts \
    && command -v prefetch \
    && command -v fasterq-dump \
    && command -v xpdf \
    && /opt/venv/bin/jupyter kernelspec list | grep -F pydeseq2 \
    && /opt/venv/bin/python -c "import pydeseq2, pandas, numpy, matplotlib, seaborn, scipy, scanpy, decoupler, pertpy" \
    && test -f /opt/workshop/notebooks/pydeseq2_notebook.ipynb

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    HOME=/tmp/workshop-home \
    PATH="/opt/venv/bin:${PATH}" \
    NOTEBOOK_SRC=/opt/workshop/notebooks \
    NOTEBOOK_DST=/data/notebooks \
    JUPYTER_PORT=8888

WORKDIR /data

EXPOSE 8888

LABEL org.opencontainers.image.title="scRNA-seq workshop" \
      org.opencontainers.image.description="Interactive STAR, HISAT2, SAMtools, featureCounts, SRA Toolkit and Jupyter Lab workshop image"

CMD ["/bin/bash"]
