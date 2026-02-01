# Multi-stage Dockerfile for Genome Assembly
# Focus: Long-read genome assembly tools

# =========================================
# Stage 1: Install all bioconda tools
# =========================================
FROM docker.io/condaforge/mambaforge:latest AS builder

WORKDIR /tmp/build

# Configure conda channels (bioconda for bioinformatics tools)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Update conda to ensure latest package metadata
RUN mamba update -n base -y mamba conda && \
    mamba clean -afy

# Install bioinformatics tools available in bioconda
# Following verkko's official installation instructions
# verkko pinned to v2.2.1 for reproducibility
RUN mamba create -n assembly -y -c conda-forge -c bioconda -c defaults \
    verkko=2.2.1 \
    hifiasm \
    seqtk \
    nanoplot \
    quast \
    bioawk \
    samtools \
    mashmap \
    minimap2 \
    igv \
    bandage \
    && mamba clean -afy

# =========================================
# Stage 2: Runtime environment (final image)
# =========================================
FROM docker.io/ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
# perl, bash, bc - needed by verkko
# X11 libraries - needed for GUI applications (IGV, Bandage)
# tree, git, openssh-client - utility tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    libgomp1 \
    wget \
    ca-certificates \
    perl \
    bash \
    python3 \
    bc \
    tree \
    git \
    openssh-client \
    libxrender1 \
    libxtst6 \
    libxi6 \
    libxrandr2 \
    libxcursor1 \
    libxdamage1 \
    libxcomposite1 \
    libxfixes3 \
    && rm -rf /var/lib/apt/lists/*

# Fix: verkko uses bash syntax but may call /bin/sh
# Ubuntu's /bin/sh is dash which doesn't support [[
# Point /bin/sh to bash for compatibility
RUN ln -sf /bin/bash /bin/sh

# Copy conda environment from builder
COPY --from=builder /opt/conda /opt/conda

# Copy tool listing script
COPY show-tools.sh /usr/local/bin/show-tools
RUN chmod +x /usr/local/bin/show-tools

# Set up environment variables
ENV PATH="/opt/conda/envs/assembly/bin:$PATH" \
    CONDA_DEFAULT_ENV=assembly \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Configure bash to show tools on interactive login
RUN echo '[ -z "$TOOLS_SHOWN" ] && [[ $- == *i* ]] && export TOOLS_SHOWN=1 && show-tools' >> /etc/bash.bashrc

# Create a working directory for analyses
WORKDIR /data

# Add labels for metadata
LABEL maintainer="bionf-fi" \
      description="Docker image for genome assembly using long reads" \
      version="1.0"

# Document installed tools and their purposes
RUN echo "=== Installed tools ===" && \
    echo "  - verkko: Telomere-to-telomere assembly pipeline" && \
    echo "  - hifiasm: HiFi read assembler" && \
    echo "  - seqtk: Toolkit for processing sequences in FASTA/Q formats" && \
    echo "  - nanoplot: Quality control for long-read sequencing data" && \
    echo "  - quast: Quality assessment tool for genome assemblies" && \
    echo "  - bioawk: AWK with biological data extensions" && \
    echo "  - samtools: Tools for manipulating SAM/BAM/CRAM files" && \
    echo "  - mashmap: Fast approximate aligner for long sequences" && \
    echo "  - minimap2: Versatile sequence alignment program" && \
    echo "  - igv: Integrative Genomics Viewer (GUI)" && \
    echo "  - bandage: Assembly graph visualization tool (GUI)"

# Default command (interactive shell)
CMD ["/bin/bash"]
