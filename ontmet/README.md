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

## 🏗️ Tools Installation Architecture

This container uses a **simplified single-source approach** for maximum reliability and maintainability:

```
🔧 TOOLS INSTALLATION SOURCES:
├── Python Tools (pip + Python 3.8 venv)
│   ├── pysam, moddotplot, whatshap, pybedtools 
│   ├── pyBigWig, ndindex, methylartist, NanoPlot
│   └── modbamtools (libmodbampy resolved)
├── C/C++ Tools (source compilation)  
│   ├── samtools, bcftools, htslib, minimap2
│   ├── seqtk, bioawk
│   └── pomfret (with HTSlib integration)
├── Rust Tools (source compilation)
│   └── modkit (cargo build)
└── Java Tools (pre-built binaries)
    └── IGV (with Java runtime)
```

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

---

## Quick Start

### Container Engine Selection

This container works with both **Docker** and **Podman**. Use the `ENGINE` variable to choose:

```bash
# Using Docker (default)
make build
make run-shell

# Using Podman
make build ENGINE=podman
make run-shell ENGINE=podman

# Set as default for session
export ENGINE=podman  # All subsequent make commands will use Podman
```

### Build & Run
```bash
# Build the container
make build

# Run interactive shell
make run-shell

# Check all installed tools
make tools
```

### Using Makefile (Recommended)
```bash
# Interactive shell with data mounting
make run-shell DATA_DIR=/path/to/your/data

# Run with GUI support (Linux X11)
make run-igv

# Launch IGV directly (Linux X11)
make igv

# List all available tools and versions
make tools
```

---

## Manual Container Usage

The container creates a non-root user `worker` (UID/GID 2000).  
Your host directory with input/output data should be mounted to `/data`.

### ⚠️ File Permission Fix for Workshops

This container is designed to be **workshop-friendly** with minimal permission issues. For **bulletproof workshop setup**, use this complete command sequence:

```bash
# Complete workshop setup (one-time per directory)
mkdir -p workshop-data
chmod 755 workshop-data
cd workshop-data

# Run container with permission mapping
docker run -it --rm --user $(id -u):$(id -g) -v "$PWD:/data:Z" bioinf-fi/ontmet:latest

# With Podman (requires U flag for write permissions)  
podman run -it --rm --userns=keep-id:uid=2000,gid=2000 -v "$PWD:/data:Z,U" bioinf-fi/ontmet:latest
```

**Why This Works:**
- ✅ **User ID mapping**: `--user $(id -u):$(id -g)` matches host user permissions
- ✅ **Host directory setup**: `chmod 755` ensures the mounted directory is accessible
- ✅ **Permissive umask (000)**: All created files inside container are world-writable
- ✅ **No build required**: Works with pre-built images from any source

**Alternative for Existing Directory:**
If you already have data in a directory and see permission errors:
```bash
chmod 755 .  # Make current directory accessible
docker run -it --rm --user $(id -u):$(id -g) -v "$PWD:/data" bioinf-fi/ontmet:latest
```

### 📊 NanoPlot Static Plot Issues

**If you see Chrome/Kaleido warnings with NanoPlot**, add `--no_static` to disable static plot generation:

```bash
# Instead of:
NanoPlot --fastq file.fastq -o output

# Use:
NanoPlot --fastq file.fastq -o output --no_static
```

This avoids Chrome dependency issues in containers while still generating interactive HTML plots.

### Standard Usage

#### Docker/Podman (interchangeable)
```bash
# Mount current directory into /data
docker run -it --rm \
  --user $(id -u):$(id -g) \
  -v "$PWD:/data" \
  bioinf-fi/ontmet:latest \
  samtools --version

# Interactive work
docker run -it --rm \
  --user $(id -u):$(id -g) \
  -v "$PWD:/data" \
  bioinf-fi/ontmet:latest
```

#### Platform-Specific Notes

**For SELinux systems (Fedora/RHEL/CentOS) - IMPORTANT for Linux workshops:**

If you get "Permission denied" errors on Linux systems with SELinux enabled:

```bash
# Docker on SELinux systems
docker run -it --rm --user $(id -u):$(id -g) -v "$PWD:/data:Z" bioinf-fi/ontmet:latest

# Podman on SELinux systems (requires additional U flag for write permissions)
podman run -it --rm --userns=keep-id:uid=2000,gid=2000 -v "$PWD:/data:Z,U" bioinf-fi/ontmet:latest
```

**Volume Mount Flags Explained:**
- `:Z` - SELinux relabeling for container access (required for both Docker & Podman)
- `,U` - **Podman-specific**: Changes ownership of mounted volume to container user
- **Docker**: Only needs `:Z` 
- **Podman**: Needs `:Z,U` for write permissions
- Not needed on macOS/Windows, but harmless if included

### Running IGV (GUI Application)

#### Linux (Recommended)
IGV works seamlessly with X11 forwarding on Linux:

```bash
# 1. Enable X11 forwarding
xhost +local:root

# 2. Launch IGV directly
make igv

# Or run interactive container with GUI support
make run-igv
```

**Manual Docker command:**
```bash
# Enable X11 forwarding
xhost +local:root

# Run container with X11 support
docker run -it --rm \
  -e DISPLAY=$DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix:rw \
  -v $(pwd):/data \
  bioinf-fi/ontmet:latest

# Launch IGV inside container
igv
```

#### macOS/Windows (Alternative)
For optimal performance on macOS and Windows, use native IGV:

1. **Download native IGV** from [igv.org](https://igv.org)
2. **Use container for data processing** (samtools, minimap2, etc.)
3. **Load results in native IGV** for visualization


**Example workflow:**
```bash
# Process data in container
make run-shell DATA_DIR=/path/to/data
# Inside container: run samtools, minimap2, etc.

# Visualize results in native IGV
# Open IGV natively and load files from /path/to/data
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
