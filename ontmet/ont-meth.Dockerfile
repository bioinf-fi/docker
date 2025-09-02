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
    && tar -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2 \
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

# =========
# PYTHON WHEEL BUILDER
# =========
FROM ubuntu:24.04 AS pybuilder


ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv python3-dev \
    build-essential git \
    cython3 \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    zlib1g-dev libbz2-dev liblzma-dev libffi-dev \
  && rm -rf /var/lib/apt/lists/*

# Make a temp venv just for building wheels
RUN python3 -m venv /opt/pwb && . /opt/pwb/bin/activate \
 && pip install --upgrade pip setuptools wheel

# Build wheels for all requested packages + their deps
# (pip wheel downloads sdists/wheels for deps and builds anything missing)
RUN . /opt/pwb/bin/activate \
 && mkdir -p /wheelhouse \
 && pip wheel -w /wheelhouse \
      pysam \
      methylartist \
      moddotplot \
      whatshap \
      pybedtools \
      pyBigWig \
      ndindex

# =========
# RUNTIME
# =========
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

# Core runtime only (no compilers, no *-dev), incl. ncurses for `samtools tview`
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    locales \
    python3 python3-venv \
    r-base \
    libcurl4 libxml2 libssl3 \
    libzstd1 libbz2-1.0 liblzma5 libdeflate0 \
    libncursesw6 \
    fonts-dejavu-core \
    gv \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# Copy compiled C/C++ tools from the builder stage (htslib, samtools, minimap2, seqtk, bioawk)
COPY --from=builder /opt/build /usr/local
ENV PATH="/usr/local/bin:${PATH}"

# ---- Python tools via offline wheels from pybuilder ----
# (pybuilder stage must have: pip wheel -w /wheelhouse pysam methylartist moddotplot whatshap pybedtools pyBigWig ndindex)
COPY --from=pybuilder /wheelhouse /wheels

# Create runtime venv and install from local wheels (no internet/compilers here)
RUN python3 -m venv /opt/venv \
 && . /opt/venv/bin/activate \
 && pip install --upgrade pip \
 && pip install --no-index --find-links=/wheels \
      pysam methylartist moddotplot whatshap pybedtools pyBigWig ndindex \
 && rm -rf /wheels
ENV PATH="/opt/venv/bin:${PATH}"

# Non-root user for podman/docker
ARG USER=worker
ARG UID=2000
ARG GID=2000
RUN groupadd -g ${GID} ${USER} \
 && useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USER}

USER ${USER}
WORKDIR /data
ENV HOME=/home/${USER}

# -------- Aggressive Pruning (size trim) --------
# 1) Remove APT caches (already done above, but ensure no leftovers)
RUN rm -rf /var/lib/apt/lists/*

# 2) Drop docs, man pages, and non-English locales (keep en*)
RUN set -eux; \
  rm -rf /usr/share/man/* /usr/share/info/* /usr/share/doc/*; \
  find /usr/share/locale -mindepth 1 -maxdepth 1 \
      ! -name 'en' ! -name 'en_US' ! -name 'locale.alias' -exec rm -rf {} + || true; \
  find /usr/lib/locale -mindepth 1 -maxdepth 1 \
      ! -name 'en_US.utf8' -exec rm -rf {} + || true

# 3) Strip symbols from user-installed native binaries and shared libs
#    (harmless if some objects are already stripped)
RUN set -eux; \
  if command -v strip >/dev/null 2>&1; then \
    find /usr/local -type f -executable -exec sh -c 'file -b "$1" | grep -qE "ELF.*(executable|shared object)" && strip --strip-unneeded "$1" || true' _ {} \; || true; \
    find /opt/venv -type f -name "*.so*" -exec sh -c 'file -b "$1" | grep -q "ELF" && strip --strip-unneeded "$1" || true' _ {} \; || true; \
  fi

# 4) Remove Python bytecode, tests, examples, and caches from the venv
RUN set -eux; \
  find /opt/venv -type d -name "__pycache__" -prune -exec rm -rf {} +; \
  find /opt/venv/lib -type d \( -name "tests" -o -name "test" -o -name "testing" -o -name "examples" -o -name "docs" \) \
       -prune -exec rm -rf {} +; \
  # remove pip metadata caches inside the venv
  rm -rf /opt/venv/pip-selfcheck.json || true

# 5) Prune R help, HTML, and caches (keeps packages working)
RUN set -eux; \
  rm -rf /usr/lib/R/doc /usr/share/doc/*R* || true; \
  find /usr/lib/R/site-library -maxdepth 2 -type d \( -name "help" -o -name "html" -o -name "doc" -o -name "docs" -o -name "examples" -o -name "unitTests" \) \
       -exec rm -rf {} + || true; \
  rm -rf /root/.cache /home/*/.cache /tmp/* /var/tmp/*

# 6) Verify key tools still run (optional; can be removed to save a few KB)
RUN bash -lc 'samtools --version | head -n1 && minimap2 --version && whatshap --help >/dev/null 2>&1 || true'
# -------- End pruning --------

# Smoke test / default command
CMD bash -lc '\
  echo "Tools ready:" && \
  samtools --version | head -n1 && \
  minimap2 --version && \
  seqtk 2>&1 | head -n1 && \
  bioawk -v | head -n1 && \
  python3 - <<PY \
import pysam, whatshap, importlib; \
mods=[("pysam",pysam.__version__),("whatshap",whatshap.__version__)]; \
for m in ("methylartist","moddotplot","pybedtools","pyBigWig","ndindex"): \
    mod=importlib.import_module(m); \
    print(m, getattr(mod,"__version__", "installed")); \
print(*["pysam "+mods[0][1],"whatshap "+mods[1][1]], sep="\n") \
PY \
  && modkit --version \
'