# ONT MET Container

A **parallel-optimized** container image bundling common genomics tools for alignment, variant calling, methylation analysis, phasing, and visualization. Built with a **2.3x faster** parallel architecture.

The image is built on **Ubuntu 24.04** and contains:

## Core Analysis Tools
- **samtools 1.20** ([docs](http://www.htslib.org/doc/samtools.html))  
- **bcftools 1.20** ([docs](http://www.htslib.org/doc/bcftools.html))
- **htslib 1.20** ([docs](http://www.htslib.org/doc/htslib.html))  
- **minimap2 2.30** ([docs](https://lh3.github.io/minimap2/))  
- **seqtk 1.4** ([repo](https://github.com/lh3/seqtk))  
- **bioawk** ([repo](https://github.com/lh3/bioawk))  
- **pomfret v0.1** ([repo](https://github.com/nanoporetech/pomfret))

## Python Tools
- **pysam 0.23.3** ([docs](https://pysam.readthedocs.io/en/latest/))  
- **methylartist 1.5.2** ([docs](https://methylartist.readthedocs.io/en/latest/))  
- **moddotplot** ([repo](https://github.com/timplab/moddotplot))  
- **WhatsHap 2.8** ([docs](https://whatshap.readthedocs.io/en/latest/))  
- **pybedtools 0.12.0** ([docs](https://daler.github.io/pybedtools/))  
- **pyBigWig 0.3.24** ([docs](https://github.com/deeptools/pyBigWig))  
- **ndindex 1.10.0** ([docs](https://quansight-labs.github.io/ndindex/))
- **NanoPlot 1.46.1** ([docs](https://github.com/wdecoster/NanoPlot))
- **modbamtools 0.4.8** ([repo](https://github.com/epi2me-labs/modbam2bed))

## Visualization & Analysis
- **modkit 0.5.0** ([repo](https://github.com/nanoporetech/modkit))  
- **IGV 2.18.2** ([docs](https://software.broadinstitute.org/software/igv/)) - **GUI Supported**
- **DSS (R/Bioconductor)** ([manual](https://www.bioconductor.org/packages/release/bioc/manuals/DSS/man/DSS.pdf))  
- **GNU gv** ([man page](https://manpages.debian.org/gv))  

---

## 🔧 Version Management

All tool versions are configurable via build arguments at the top of the Dockerfile. You can easily customize versions without editing the entire file.

### Changing Default Versions

Edit the version configuration section at the top of `ontmet.Dockerfile`:

```dockerfile
# =============================================
# VERSION CONFIGURATION - EDIT HERE
# =============================================
# Core C/C++ Tools
ARG HTSLIB_VERSION=1.20      # ← Change this
ARG SAMTOOLS_VERSION=1.20    # ← Change this
ARG BCFTOOLS_VERSION=1.20    # ← Change this
ARG MINIMAP2_VERSION=2.30    # ← Change this
ARG SEQTK_VERSION=1.4

# Rust Tools  
ARG MODKIT_VERSION=0.5.0

# Python Environment
ARG PYTHON_VERSION=3.9

# GUI Tools
ARG IGV_VERSION=2.18.2

# Base OS
ARG UBUNTU_VERSION=24.04
```

### Overriding Versions at Build Time

You can override any version without modifying the Dockerfile:

```bash
# Build with custom tool versions
docker build \
  --build-arg MINIMAP2_VERSION=2.31 \
  --build-arg SAMTOOLS_VERSION=1.21 \
  --build-arg PYTHON_VERSION=3.10 \
  -t genomics-tools:custom .

# Build with latest IGV
docker build --build-arg IGV_VERSION=2.19.0 -t genomics-tools:igv-latest .

# Using make (add ARGS variable)
make build ARGS='--build-arg MINIMAP2_VERSION=2.31'
```

### Available Versions

Check each tool's GitHub releases or documentation:
- [samtools releases](https://github.com/samtools/samtools/releases)
- [minimap2 releases](https://github.com/lh3/minimap2/releases)  
- [modkit releases](https://github.com/nanoporetech/modkit/releases)
- [IGV downloads](https://software.broadinstitute.org/software/igv/download)

---

## Quick Start

### Build & Run
```bash
# Build the container
make build

# Run interactive shell
make run-shell

# Check all installed tools
docker run --rm bioinf-fi/ontmet:latest bash -l -c "tools"
```

### Using Makefile (Recommended)
```bash
# Interactive shell with data mounting
make run-shell DATA_DIR=/path/to/your/data

# Run with GUI support (for IGV)
make run-igv

# Launch IGV directly
make igv
```

---

## Manual Container Usage

The container creates a non-root user `worker` (UID/GID 2000).  
Your host directory with input/output data should be mounted to `/data`.

### Standard Usage

#### Docker
```bash
# Mount current directory into /data
docker run -it --rm \
  -v "$PWD:/data" \
  bioinf-fi/ontmet:latest \
  samtools --version

# Interactive work
docker run -it --rm \
  -v "$PWD:/data" \
  bioinf-fi/ontmet:latest \
  bash
```

#### Podman
```bash
# Mount current directory into /data  
podman run -it --rm \
  -v "$PWD:/data:Z" \
  bioinf-fi/ontmet:latest \
  minimap2 --version

# Interactive work
podman run -it --rm \
  -v "$PWD:/data:Z" \
  bioinf-fi/ontmet:latest \
  bash
```
- `:Z` applies SELinux relabeling (needed on Fedora/RHEL/CentOS)

### Running IGV (GUI Application)

IGV support varies by operating system:

#### Linux Desktop (Recommended)
IGV works well with X11 forwarding:

```bash
# 1. Enable X11 forwarding
xhost +local:root

# 2. Run with GUI support
make run-igv

# 3. Launch IGV inside container
igv
```

**Manual method:**
```bash
# Enable X11 forwarding
xhost +local:root

# Run container with X11 support
docker run -it --rm \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
  -v $(pwd):/data \
  bioinf-fi/ontmet:latest

# Launch IGV
igv
```

#### macOS (Alternative Approaches)

**Option 1: Native IGV (Recommended for macOS)**
- Download IGV directly for macOS from [igv.org](https://igv.org)
- Use container for data processing, native IGV for visualization
- Best performance and user experience

**Option 2: XQuartz + Docker (Advanced)**
```bash
# 1. Install XQuartz
brew install --cask xquartz

# 2. Start XQuartz and enable network connections
# In XQuartz preferences: Security > "Allow connections from network clients"

# 3. Get your IP address
IP=$(ifconfig en0 | grep inet | awk '$1=="inet" {print $2}')

# 4. Allow X11 connections
xhost + $IP

# 5. Run container with X11 forwarding
docker run -it --rm \
  -e DISPLAY=$IP:0 \
  -v $(pwd):/data \
  bioinf-fi/ontmet:latest

# 6. Launch IGV inside container
igv
```

**Option 3: VNC Server (Most Compatible)**
```bash
# Run container with VNC server
docker run -it --rm \
  -p 5901:5901 \
  -v $(pwd):/data \
  bioinf-fi/ontmet:latest \
  bash -c "
    apt update && apt install -y tightvncserver xfce4 &&
    vncserver :1 -geometry 1024x768 -depth 24 &&
    DISPLAY=:1 igv
  "

# Connect using VNC client to localhost:5901
```

#### Windows (WSL2 + X11)
```bash
# With WSL2 and X11 server (like VcXsrv)
# 1. Install X11 server for Windows
# 2. In WSL2:
export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0
make run-igv
```

---

## Web-Based Tools

Some tools provide web interfaces that require port mapping:

### Methylartist GUI
```bash
docker run -it --rm \
  -v "$PWD:/data" \
  -p 8888:8888 \
  bioinf-fi/ontmet:latest \
  methylartist gui --port 8888
```

### Moddotplot Dashboard  
```bash
docker run -it --rm \
  -v "$PWD:/data" \
  -p 8050:8050 \
  bioinf-fi/ontmet:latest \
  moddotplot dash --port 8050
```

Then open [http://localhost:8888](http://localhost:8888) or [http://localhost:8050](http://localhost:8050).

---

## Performance Features

### Parallel Build Architecture
- **2.3x faster builds** through parallel compilation stages
- Independent C/C++, Rust, Python, and binary download stages
- Optimized layer caching for faster rebuilds

### Tool Verification
Run `tools` command inside the container to verify all installations:
```bash
docker run --rm bioinf-fi/ontmet:latest bash -l -c "tools"
```

---

## Advanced Usage

### Custom Data Directory
```bash
# Mount specific directory
make run-shell DATA_DIR=/path/to/genomics/data

# Or manually
docker run -it --rm \
  -v /path/to/data:/data \
  bioinf-fi/ontmet:latest
```

### Memory & Performance
```bash
# Increase memory for large files
docker run -it --rm \
  --memory=8g \
  -v "$PWD:/data" \
  bioinf-fi/ontmet:latest
```

### R/Bioconductor DSS
```bash
# Launch R with DSS
docker run -it --rm bioinf-fi/ontmet:latest R

# Inside R
library(DSS)
```

---

## Container Details

- **Base**: Ubuntu 24.04 LTS
- **User**: `worker` (UID/GID 2000) - non-root for security
- **Python**: Virtualenv at `/opt/venv` (automatically activated)
- **Working Directory**: `/data` (mount your data here)
- **Home Directory**: `/home/worker`
- **Size**: ~2.5GB (optimized with multi-stage builds)

## Build Architecture

The container uses a **parallel multi-stage build** for optimal performance:

1. **C/C++ Tools Builder**: HTSlib, samtools, minimap2, seqtk, bioawk, pomfret
2. **Rust Builder**: modkit compilation
3. **Python Builder**: All Python packages and wheels
4. **Binary Downloader**: IGV and other pre-built tools
5. **Runtime Assembly**: Combines all artifacts into final image

This architecture provides **2.3x faster builds** compared to sequential compilation.  