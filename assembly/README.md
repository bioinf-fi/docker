# Genome Assembly Container

Container image for genome assembly with long reads. Works with Docker or Podman.

## Tools Included

**Assembly:** verkko, hifiasm, minimap2, mashmap  
**Quality:** nanoplot, quast  
**Processing:** seqtk, bioawk, samtools  
**Visualization:** IGV, Bandage (GUI)

## Quick Start

```bash
# Build
make build

# Run command-line tools
make run

# Run with GUI support (for IGV/Bandage)
make gui
```

## Setup

### Prerequisites
- Linux desktop with X11
- Docker or Podman installed

### First Time
```bash
cd assembly
make build
```

This takes a few minutes. You only do it once.

## Usage

### Command-Line Tools

```bash
make run
```

Inside the container:
```bash
# Quality control
NanoPlot --fastq reads.fastq.gz -o qc_output

# Assembly
hifiasm -o assembly -t 8 reads.fastq.gz

# Quality assessment
quast assembly.fasta -o results

# View installed tools
show-tools
```

### GUI Applications

```bash
make gui
```

Inside the container:
```bash
# View assembly graphs
Bandage

# View alignments
igv
```

Your current directory is mounted at `/data` in the container.

## Troubleshooting

### GUI Won't Start

Try relaxed X11 mode:
```bash
make gui-simple
```

Or manually:
```bash
xhost +local:
make gui
xhost -local:
```

### Permission Errors

With Docker:
```bash
docker run -it --rm -u $(id -u):$(id -g) -v $(pwd):/data genome-assembly:latest
```

With Podman, this is automatic.

### Check Your Setup

```bash
make test-x11    # Test X11 for GUI
make test        # Test tools installed
make help        # See all commands
```

## Examples

```bash
# Assembly with hifiasm
hifiasm -o asm -t 8 hifi_reads.fastq.gz
awk '/^S/{print ">"$2;print $3}' asm.bp.p_ctg.gfa > assembly.fasta

# Hybrid assembly with verkko
verkko -d output --hifi hifi.fastq.gz --nano ont.fastq.gz

# Align reads
minimap2 -ax map-hifi ref.fasta reads.fastq.gz | samtools sort -o aligned.bam

# View graph
Bandage load asm.bp.p_ctg.gfa
```

## Technical Notes

- Auto-detects Docker or Podman
- SELinux-compatible (Fedora/RHEL)
- X11 forwarding for GUI apps
- Volume mounts at `/data`
- Tools from bioconda (latest versions)

## Make Commands

```bash
make build       # Build image
make run         # Start container (CLI)
make gui         # Start container (GUI)
make gui-simple  # Start with relaxed X11
make test-x11    # Test X11 setup
make test        # Test tool installation
make clean       # Remove image
make help        # Show all commands
```

## Docker vs Podman

Both work. The Makefile detects which you have. Run `make help` to see.

**Podman users on Fedora/RHEL:** Everything is automatic. SELinux is handled.

## Resource Limits

```bash
# Docker
docker run -it --rm --memory=32g --cpus=8 -v $(pwd):/data genome-assembly:latest

# Podman
podman run -it --rm --memory=32g --cpus=8 -v $(pwd):/data:Z genome-assembly:latest
```

## File Organization

All your data files go in the directory where you run `make gui` or `make run`. They appear at `/data` inside the container.

## Common Issues

**"Cannot open display"**  
Use `make gui-simple`

**"Permission denied" on files**  
Podman handles this automatically. Docker users add `-u $(id -u):$(id -g)`

**Wayland instead of X11**  
```bash
export DISPLAY=:0
xhost +local:
make gui
```

**Check which runtime you're using**  
```bash
make help  # First line shows docker or podman
```

## Support

- Tool versions: `show-tools` (inside container)
- IGV docs: https://software.broadinstitute.org/software/igv/
- Bandage help: `Bandage --help`

---

Built with bioconda. Multi-stage build for minimal size. Compatible with Docker and Podman.
