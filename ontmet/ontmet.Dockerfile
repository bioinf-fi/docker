# =========================================
# PARALLEL OPTIMIZED GENOMICS DOCKERFILE
# =========================================
# This version uses parallel stages to maximize build speed

# =========
# STAGE 1: C/C++ TOOLS BUILDER (Parallel Group A)
# =========
FROM ubuntu:24.04 AS c_tools_builder

ENV DEBIAN_FRONTEND=noninteractive
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

# --- HTSlib (foundation for samtools) ---
ARG HTSLIB_VERSION=1.20
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} && autoheader && autoconf || true \
    && ./configure --prefix=/opt/build && make -j"$(nproc)" && make install

# --- samtools (depends on HTSlib) ---
ARG SAMTOOLS_VERSION=1.20
RUN curl -fsSL -o samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --with-htslib=/opt/build --prefix=/opt/build \
    && make -j"$(nproc)" && make install

# --- minimap2 (independent) ---
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

# --- seqtk (independent) ---
ARG SEQTK_VERSION=1.4
RUN git clone --branch v${SEQTK_VERSION} --depth 1 https://github.com/lh3/seqtk.git \
    && cd seqtk && make -j"$(nproc)" \
    && install -m 0755 seqtk /opt/build/bin/seqtk

# --- bioawk (independent) ---
RUN git clone --depth 1 https://github.com/lh3/bioawk.git \
    && cd bioawk \
    && bison -y -d awkgram.y \
    && mv y.tab.c ytab.c && mv y.tab.h ytab.h \
    && make YACC="bison -y" -j1 \
    && install -m 0755 bioawk /opt/build/bin/bioawk

# --- pomfret (will use HTSlib from this stage) ---
RUN mkdir -p /usr/bin/htslib \
    && ln -sf /opt/build/include/htslib /usr/bin/htslib/ \
    && ln -sf /opt/build/lib/libhts.* /usr/bin/htslib/ \
    && curl -fsSL -o /tmp/pomfret.tar.gz "https://github.com/nanoporetech/pomfret/archive/refs/heads/main.tar.gz" \
    && cd /tmp && tar -xzf pomfret.tar.gz \
    && cd pomfret-main \
    && make CFLAGS="-g -O2 -Wall -Wno-error -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-label -fcommon" \
    && install -m 0755 pomfret /opt/build/bin/pomfret \
    && cd / && rm -rf /tmp/pomfret* /usr/bin/htslib

# =========
# STAGE 2: RUST BUILDER (Parallel Group B)
# =========
FROM ubuntu:24.04 AS rust_builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl pkg-config libssl-dev zlib1g-dev build-essential git \
    && rm -rf /var/lib/apt/lists/*

# Install rustup (minimal profile)
RUN curl -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /opt/build
RUN mkdir -p /opt/build/bin

# --- modkit (Rust tool) ---
ARG MODKIT_TAG=v0.5.0
RUN git clone --depth 1 --branch ${MODKIT_TAG} https://github.com/nanoporetech/modkit.git \
    && cd modkit \
    && cargo install --path modkit --root /opt/build \
    && /opt/build/bin/modkit --version

# =========
# STAGE 3: PYTHON WHEEL BUILDER (Parallel Group C)
# =========
FROM ubuntu:24.04 AS python_builder

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

# Create build venv using Python 3.9 (required for modbamtools)
RUN python3.9 -m venv /opt/pwb && . /opt/pwb/bin/activate \
    && pip install --upgrade pip setuptools wheel

# Build wheels for all Python packages
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
# STAGE 4: BINARY DOWNLOADER (Parallel Group D)
# =========
FROM ubuntu:24.04 AS binary_downloader

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/downloads

# --- IGV Download ---
ARG IGV_VERSION=2.18.2
RUN wget -O IGV_Linux_${IGV_VERSION}_WithJava.zip \
    "https://data.broadinstitute.org/igv/projects/downloads/2.18/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
    && unzip IGV_Linux_${IGV_VERSION}_WithJava.zip \
    && rm IGV_Linux_${IGV_VERSION}_WithJava.zip

# =========
# STAGE 5: RUNTIME ASSEMBLY
# =========
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install minimal runtime dependencies
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
    && rm -rf /var/lib/apt/lists/* \
    && locale-gen en_US.UTF-8
ENV LANG=C.UTF-8

# Copy artifacts from all parallel stages
# Stage 1: C/C++ tools (htslib, samtools, minimap2, seqtk, bioawk, pomfret)
COPY --from=c_tools_builder /opt/build /usr/local

# Stage 2: Rust tools (modkit)
COPY --from=rust_builder /opt/build/bin/modkit /usr/local/bin/modkit

# Stage 3: Python wheels
COPY --from=python_builder /wheelhouse /wheels

# Stage 4: Binary downloads (IGV)
COPY --from=binary_downloader /opt/downloads/IGV_Linux_* /opt/

# Configure library paths and create IGV symlink
ENV LD_LIBRARY_PATH="/usr/local/lib"
RUN ldconfig
ENV PATH="/usr/local/bin:${PATH}"

# Link IGV (find the actual path since directory name may vary)
RUN find /opt -name 'igv.sh' -exec ln -sf {} /usr/local/bin/igv \;

# Create runtime Python environment and install wheels
COPY verify.py /tmp/verify.py
RUN set -eux; \
    python3.9 -m venv /opt/venv; \
    . /opt/venv/bin/activate; \
    pip install --upgrade pip; \
    ls -1 /wheels | sed 's/^/WHEEL: /'; \
    pip install --no-index --find-links=/wheels \
        cython pysam moddotplot whatshap pybedtools pyBigWig ndindex methylartist NanoPlot; \
    pip install --no-index --find-links=/wheels modbamtools || pip install modbamtools; \
    python /tmp/verify.py

# Clean up wheels and set Python environment
RUN rm -rf /wheels
ENV PATH="/opt/venv/bin:${PATH}"

# Create non-root user
ARG USER=worker
ARG UID=2000
ARG GID=2000
RUN groupadd -g ${GID} ${USER} \
    && useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USER}

# Configure interactive environment
ENV PATH="/opt/venv/bin:/usr/local/bin:${PATH}"
COPY genomics-tools.sh /etc/profile.d/genomics-tools.sh

RUN echo '. /etc/profile.d/genomics-tools.sh 2>/dev/null' >> /home/worker/.bashrc \
    && echo '[ -z "$NO_GREETING" ] && [[ $- == *i* ]] && tools || true' >> /home/worker/.bashrc

USER ${USER}
WORKDIR /data
ENV HOME=/home/${USER}

CMD ["bash"]
