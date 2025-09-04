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
      moddotplot \
      whatshap \
      pybedtools \
      pybigwig \
      ndindex \
 && ( \
      pip wheel -w /wheelhouse methylartist \
      || pip wheel -w /wheelhouse git+https://github.com/adamewing/methylartist \
    ) \
 && python - <<'PY'
import glob, sys
wheels = glob.glob("/wheelhouse/*.whl")
need = {"pysam","moddotplot","whatshap","pybedtools","pybigwig","ndindex","methylartist"}
have = {w.split("/")[-1].split("-")[0].lower() for w in wheels}
missing = sorted(need - have)
if missing:
    print("Missing wheels:", missing)
    sys.exit(1)
print("Wheelhouse OK.")
PY

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
    libpython3.12 \  
    r-base \
    libcurl4 libxml2 libssl3 \
    libzstd1 libbz2-1.0 liblzma5 libdeflate0 \
    libncursesw6 \
    fonts-dejavu-core \
    gv \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8
ENV LANG=C.UTF-8

# Copy compiled C/C++ tools from the builder stage (htslib, samtools, minimap2, seqtk, bioawk)
COPY --from=builder /opt/build /usr/local
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
RUN ldconfig
ENV PATH="/usr/local/bin:${PATH}"

# ---- Python tools via offline wheels from pybuilder ----
# (pybuilder stage must have: pip wheel -w /wheelhouse pysam methylartist moddotplot whatshap pybedtools pybigwig ndindex)
COPY --from=pybuilder /wheelhouse /wheels

# Create runtime venv and install from local wheels (no internet/compilers here)
# after copying /wheels
RUN set -eux; \
  python3 -m venv /opt/venv; \
  . /opt/venv/bin/activate; \
  pip install --upgrade pip; \
  ls -1 /wheels | sed 's/^/WHEEL: /'; \
  # install everything from the wheelhouse (names allow pip to resolve deps inside /wheels)
  pip install --no-index --find-links=/wheels \
      cython pysam moddotplot whatshap pybedtools pyBigWig ndindex methylartist; \
  # verify core Python libs by import
  python - <<'PY'
import importlib, sys
mods=("pysam","whatshap","moddotplot","pybedtools","pyBigWig","ndindex")
bad=[]
for m in mods:
    try: importlib.import_module(m)
    except Exception as e: bad.append((m,e.__class__.__name__,str(e)))
if bad: print("Import failures:",bad); sys.exit(1)
print("Core Python modules import OK.")
PY
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
RUN tee /etc/profile.d/genomics-tools.sh >/dev/null <<'SH'
# Genomics Toolkit helpers

tools() {
  echo "Tools ready:"
  command -v samtools >/dev/null 2>&1 \
    && samtools --version 2>/dev/null | head -n1 \
    || echo "samtools: not found (or missing libhts)"
  # minimap2 (print name + version)
  if command -v minimap2 >/dev/null 2>&1; then
    v="$(minimap2 --version 2>/dev/null || true)"
    [ -n "$v" ] && echo "minimap2 $v" || echo "minimap2: installed"
  else
    echo "minimap2: not found"
  fi
  # seqtk (print name + version)
    if command -v seqtk >/dev/null 2>&1; then
    ver="$(seqtk 2>&1 | awk '/^Version:/{v=$2; if ($3) v=v"-"$3; print v; exit}')"
    if [ -n "$ver" ]; then
    echo "seqtk $ver"
    else
    # fallback: first non-empty line of help
    line="$(seqtk 2>&1 | awk 'NF{print;exit}')"
    [ -n "$line" ] && echo "seqtk $line" || echo "seqtk: installed"
    fi
    else
    echo "seqtk: not found"
    fi
  command -v bioawk >/dev/null 2>&1 \
    && echo "bioawk: installed" \
    || echo "bioawk: not found"

  # Import-check for core Python libs (exclude methylartist)
  /opt/venv/bin/python - <<'PY'
import importlib
mods=("pysam","whatshap","moddotplot","pybedtools","pyBigWig","ndindex")
for m in mods:
    try:
        mod=importlib.import_module(m)
        print(f"{m} {getattr(mod,'__version__','installed')}")
    except Exception as e:
        print(f"{m}: import failed ({e.__class__.__name__}: {e})")
PY

  # Methylartist via CLI
  if command -v methylartist >/dev/null 2>&1; then
    # some versions print help on --version; fall back to a simple presence message
    methylartist --version 2>/dev/null || echo "methylartist: installed (CLI)"
  else
    echo "methylartist: not found"
  fi

  command -v modkit >/dev/null 2>&1 \
    && modkit --version \
    || echo "modkit: not found"
}

# Greeting in interactive shells (disable with NO_GREETING=1)
if [ -n "$BASH_VERSION" ] && [[ $- == *i* ]] && [ -z "$NO_GREETING" ]; then
  echo "Welcome to the ONT MET container. Type 'tools' to see versions."
fi
SH

# Make sure the worker shell loads the helper and shows the tools summary at login
RUN echo '. /etc/profile.d/genomics-tools.sh 2>/dev/null' >> /home/worker/.bashrc \
 && echo '[ -z "$NO_GREETING" ] && [[ $- == *i* ]] && tools || true' >> /home/worker/.bashrc

USER ${USER}
WORKDIR /data
ENV HOME=/home/${USER}

# Make `bash` the default so interactive shells get the greeting and function
CMD ["bash"]