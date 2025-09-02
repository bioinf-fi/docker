# Genomics Toolkit Container

A container image bundling common genomics tools for alignment, variant calling, methylation analysis, phasing, and visualization.

The image is built on **Ubuntu 24.04** and contains:

- **samtools** ([docs](http://www.htslib.org/doc/samtools.html))  
- **htslib** ([docs](http://www.htslib.org/doc/htslib.html))  
- **minimap2** ([docs](https://lh3.github.io/minimap2/))  
- **seqtk** ([repo](https://github.com/lh3/seqtk))  
- **bioawk** ([repo](https://github.com/lh3/bioawk))  
- **pysam** ([docs](https://pysam.readthedocs.io/en/latest/))  
- **methylartist** ([docs](https://methylartist.readthedocs.io/en/latest/))  
- **moddotplot** ([repo](https://github.com/timplab/moddotplot))  
- **modkit** ([repo](https://github.com/nanoporetech/modkit))  
- **WhatsHap** ([docs](https://whatshap.readthedocs.io/en/latest/))  
- **pybedtools** ([docs](https://daler.github.io/pybedtools/))  
- **pyBigWig** ([docs](https://github.com/deeptools/pyBigWig))  
- **ndindex** ([docs](https://quansight-labs.github.io/ndindex/))  
- **DSS (R/Bioconductor)** ([manual](https://www.bioconductor.org/packages/release/bioc/manuals/DSS/man/DSS.pdf))  
- **GNU gv** ([man page](https://manpages.debian.org/gv))  

---

## Running the Image

The container creates a non-root user `worker` (UID/GID 2000).  
Your host directory with input/output data should be mounted to `/data`.

### Podman

    # Mount current directory into /data
    podman run -it --rm \
      -v "$PWD:/data:Z" \
      localhost/genomics-toolkit:latest \
      samtools --version

    # For interactive work:
    podman run -it --rm \
      -v "$PWD:/data:Z" \
      localhost/genomics-toolkit:latest \
      bash

- `:Z` applies SELinux relabeling (needed on Fedora/RHEL/CentOS).

### Docker

    # Mount current directory into /data
    docker run -it --rm \
      -v "$PWD:/data" \
      genomics-toolkit:latest \
      minimap2 --version

    # For interactive work:
    docker run -it --rm \
      -v "$PWD:/data" \
      genomics-toolkit:latest \
      bash

---

## Ports

Most tools are **command-line only** and don’t require network ports.  
If you run Python-based visualization servers (e.g. `methylartist gui`, `moddotplot dash`), you need to map a port:

### Podman

    podman run -it --rm \
      -v "$PWD:/data:Z" \
      -p 8888:8888 \
      localhost/genomics-toolkit:latest \
      methylartist gui --port 8888

### Docker

    docker run -it --rm \
      -v "$PWD:/data" \
      -p 8888:8888 \
      genomics-toolkit:latest \
      methylartist gui --port 8888

Then open [http://localhost:8888](http://localhost:8888).

---

## Notes

- All Python tools are installed in a virtualenv at `/opt/venv` and exposed on `PATH`.  
- R + Bioconductor DSS is available by launching R inside the container:

        podman run -it --rm localhost/genomics-toolkit:latest R

  and then inside R:

        library(DSS)

- User home directory inside the container: `/home/worker`.  
- Default working directory: `/data` (your mounted host folder).  
- Image has been aggressively pruned to reduce size (stripped binaries, removed docs/locales/caches).  