# =========
# BUILDER
# =========
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
# BUILDER stage
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    curl \
    wget \
    pkg-config \
    autoconf automake libtool \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libdeflate-dev \
    libncurses-dev \
    bison flex \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/build
ENV PATH="/opt/build/bin:${PATH}"
RUN mkdir -p /opt/build/bin

# --- HTSlib (for samtools build) ---
ARG HTSLIB_VERSION=1.20
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} && autoheader && autoconf || true \
    && ./configure --prefix=/opt/build && make -j"$(nproc)" && make install

# --- samtools ---
ARG SAMTOOLS_VERSION=1.20
RUN curl -fsSL -o samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --with-htslib=/opt/build --prefix=/opt/build \
    && make -j"$(nproc)" && make install

# --- minimap2 ---
ARG MINIMAP2_VERSION=2.28
RUN curl -fsSL -o minimap2-${MINIMAP2_VERSION}.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && tar --no-same-owner -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && cd minimap2-${MINIMAP2_VERSION} \
    && arch="$(uname -m)"; \
        if [ "$arch" = "aarch64" ] || [ "$arch" = "arm64" ]; then \
            make -j"$(nproc)" aarch64=1; \
        else \
            make -j"$(nproc)"; \
        fi \
    && install -m 0755 minimap2 /opt/build/bin/minimap2

# --- seqtk ---
ARG SEQTK_VERSION=1.4
RUN git clone --branch v${SEQTK_VERSION} --depth 1 https://github.com/lh3/seqtk.git \
    && cd seqtk && make -j"$(nproc)" \
    && install -m 0755 seqtk /opt/build/bin/seqtk

# --- bioawk ---
ARG BIOAWK_REF=
RUN git clone --depth 1 https://github.com/lh3/bioawk.git \
    && cd bioawk \
    # generate parser files first (Makefile expects ytab.* to exist)
    && bison -y -d awkgram.y \
    && mv y.tab.c ytab.c && mv y.tab.h ytab.h \
    # build (serialize to avoid races)
    && make YACC="bison -y" -j1 \
    && install -m 0755 bioawk /opt/build/bin/bioawk

# deps for rust/cargo builds
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl pkg-config libssl-dev zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# install rustup (minimal profile)
RUN curl -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

# build Modkit from a tagged release (set to the current latest)
ARG MODKIT_TAG=v0.5.0
RUN git clone --depth 1 --branch ${MODKIT_TAG} https://github.com/nanoporetech/modkit.git \
 && cd modkit \
 # install the 'modkit' package from the workspace into /opt/build/bin
 && cargo install --path modkit --root /opt/build \
 && /opt/build/bin/modkit --version

# --- pomfret installation deferred ---
# Will be installed in runtime stage using a different approach
# Issues with current build: multiple definition errors, hardcoded paths in Makefile

# --- modbamtools ---
# Install via pip since git clone is having network issues in build environment
# modbamtools is available as a Python package, so we'll handle it in the pybuilder stage

# =========
# PYTHON WHEEL BUILDER
# =========
FROM ubuntu:24.04 AS pybuilder


ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
  && add-apt-repository ppa:deadsnakes/ppa \
  && apt-get update && apt-get install -y --no-install-recommends \
    python3.9 python3.9-venv python3.9-dev python3-pip \
    build-essential git \
    cython3 \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    zlib1g-dev libbz2-dev liblzma-dev libffi-dev \
  && rm -rf /var/lib/apt/lists/*

# Make a temp venv just for building wheels using Python 3.9 (required for modbamtools)
RUN python3.9 -m venv /opt/pwb && . /opt/pwb/bin/activate \
 && pip install --upgrade pip setuptools wheel

# Build wheels for all requested packages + their deps
# (pip wheel downloads sdists/wheels for deps and builds anything missing)
COPY check_wheels.py /tmp/check_wheels.py
RUN . /opt/pwb/bin/activate \
 && mkdir -p /wheelhouse \
 && pip wheel -w /wheelhouse \
      pysam \
      moddotplot \
      whatshap \
      pybedtools \
      pybigwig \
      ndindex \
      NanoPlot \
 && ( \
      pip wheel -w /wheelhouse methylartist \
      || pip wheel -w /wheelhouse git+https://github.com/adamewing/methylartist \
    ) \
 && ( \
      pip wheel -w /wheelhouse modbamtools \
      || echo "modbamtools not available via pip, will need manual installation" \
    ) \
 && python3 /tmp/check_wheels.py

# =========
# RUNTIME
# =========
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

# Core runtime only (no compilers, no *-dev), incl. ncurses for `samtools tview`
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
  && add-apt-repository ppa:deadsnakes/ppa \
  && apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    locales \
    python3.9 python3.9-venv libpython3.9 \
    r-base \
    libcurl4 libxml2 libssl3 \
    libzstd1 libbz2-1.0 liblzma5 libdeflate0 \
    libncursesw6 \
    fonts-dejavu-core \
    gv \
    default-jre-headless \
    wget \
    unzip \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8
ENV LANG=C.UTF-8

# Copy compiled C/C++ tools from the builder stage (htslib, samtools, minimap2, seqtk, bioawk)
COPY --from=builder /opt/build /usr/local
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
RUN ldconfig
ENV PATH="/usr/local/bin:${PATH}"

# --- Install IGV ---
ARG IGV_VERSION=2.18.2
RUN wget -O /tmp/IGV_Linux_${IGV_VERSION}_WithJava.zip \
    "https://data.broadinstitute.org/igv/projects/downloads/2.18/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
    && unzip /tmp/IGV_Linux_${IGV_VERSION}_WithJava.zip -d /opt/ \
    && ln -s /opt/IGV_Linux_${IGV_VERSION}/igv.sh /usr/local/bin/igv \
    && rm /tmp/IGV_Linux_${IGV_VERSION}_WithJava.zip

# --- Install pomfret (using symlinks to satisfy hardcoded paths) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /usr/bin/htslib \
    && ln -sf /usr/local/include/htslib /usr/bin/htslib/ \
    && ln -sf /usr/local/lib/libhts.* /usr/bin/htslib/ \
    && curl -fsSL -o /tmp/pomfret.tar.gz "https://github.com/nanoporetech/pomfret/archive/refs/heads/main.tar.gz" \
    && cd /tmp && tar -xzf pomfret.tar.gz \
    && cd pomfret-main \
    && make CFLAGS="-g -O2 -Wall -Wno-error -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-label -fcommon" \
    && install -m 0755 pomfret /usr/local/bin/pomfret \
    && cd / && rm -rf /tmp/pomfret* \
    && rm -rf /usr/bin/htslib \
    && apt-get remove -y build-essential curl zlib1g-dev \
    && apt-get autoremove -y

# Note: modbamtools attempted via Python wheels

# ---- Python tools via offline wheels from pybuilder ----
# (pybuilder stage must have: pip wheel -w /wheelhouse pysam methylartist moddotplot whatshap pybedtools pybigwig ndindex NanoPlot)
COPY --from=pybuilder /wheelhouse /wheels

# Create runtime venv and install from local wheels (no internet/compilers here)
# after copying /wheels
COPY verify.py /tmp/verify.py
RUN set -eux; \
  python3.9 -m venv /opt/venv; \
  . /opt/venv/bin/activate; \
  pip install --upgrade pip; \
  ls -1 /wheels | sed 's/^/WHEEL: /'; \
  # install everything from the wheelhouse (names allow pip to resolve deps inside /wheels)
  pip install --no-index --find-links=/wheels \
      cython pysam moddotplot whatshap pybedtools pyBigWig ndindex methylartist NanoPlot; \
  # Try to install modbamtools if wheel exists, otherwise install directly from PyPI
  pip install --no-index --find-links=/wheels modbamtools || pip install modbamtools; \
  # verify core Python libs by import
  python /tmp/verify.py
# clean up wheels after successful install
RUN rm -rf /wheels
ENV PATH="/opt/venv/bin:${PATH}"

# Non-root user for podman/docker
ARG USER=worker
ARG UID=2000
ARG GID=2000
RUN groupadd -g ${GID} ${USER} \
 && useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USER}

# --- Interactive greeting + `tools` helper for bash shells ---
# System-wide function (root + all users)
# Ensure venv/bin is first so CLIs resolve (belt-and-suspenders)
ENV PATH="/opt/venv/bin:/usr/local/bin:${PATH}"

# Replace the tools helper to avoid importing methylartist (use CLI check instead)
COPY genomics-tools.sh /etc/profile.d/genomics-tools.sh

# Make sure the worker shell loads the helper and shows the tools summary at login
RUN echo '. /etc/profile.d/genomics-tools.sh 2>/dev/null' >> /home/worker/.bashrc \
 && echo '[ -z "$NO_GREETING" ] && [[ $- == *i* ]] && tools || true' >> /home/worker/.bashrc

USER ${USER}
WORKDIR /data
ENV HOME=/home/${USER}

# Make `bash` the default so interactive shells get the greeting and function
CMD ["bash"]
