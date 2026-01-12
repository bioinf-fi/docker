# Genome Assembly Docker Image

Docker image for interactive bioinformatics class focused on genome assembly using long reads.

## Included Tools

### Core Assembly Tools
- **verkko** - Telomere-to-telomere assembly pipeline for accurate long-read assembly
- **hifiasm** - Fast haplotype-resolved de novo assembler for PacBio HiFi reads
- **minimap2** - Versatile pairwise aligner for genomic sequences

### Quality Control & Assessment
- **nanoplot** - Visualization and quality control for long-read sequencing data
- **quast** - Quality assessment tool for genome assemblies

### Sequence Processing
- **seqtk** - Fast toolkit for processing sequences in FASTA/FASTQ formats
- **bioawk** - AWK with biological data format extensions
- **samtools** - Tools for manipulating sequence alignment/map formats

### Alignment & Mapping
- **mashmap** - Fast approximate aligner for long DNA sequences

## Building the Image

```bash
# Build from the assembly directory
docker build -f assembly.Dockerfile -t genome-assembly:latest .

# Or build from the repository root
docker build -f assembly/assembly.Dockerfile -t genome-assembly:latest .
```

## Running the Container

### Interactive Shell (Recommended for Class)
```bash
# Run with current directory mounted
# Tool versions are displayed automatically on startup
docker run -it --rm -v $(pwd):/data genome-assembly:latest

# Run with specific data directory
docker run -it --rm -v /path/to/data:/data genome-assembly:latest

# Inside the container, you can run 'show-tools' anytime to see available tools
```

### Running Specific Commands
```bash
# Check tool versions
docker run --rm genome-assembly:latest hifiasm --version

# Run assembly workflow
docker run --rm -v $(pwd):/data genome-assembly:latest \
    hifiasm -o output -t 8 reads.fastq.gz
```

## Usage Examples

### Quality Control with NanoPlot
```bash
docker run -it --rm -v $(pwd):/data genome-assembly:latest
# Inside container:
NanoPlot --fastq reads.fastq.gz -o qc_output
```

### Assembly with Hifiasm
```bash
# HiFi reads assembly
hifiasm -o assembly.asm -t 8 hifi_reads.fastq.gz

# Extract primary assembly
awk '/^S/{print ">"$2;print $3}' assembly.asm.bp.p_ctg.gfa > assembly.fasta
```

### Assembly with Verkko
```bash
# Hybrid assembly (HiFi + ONT)
verkko -d output_dir \
    --hifi hifi_reads.fastq.gz \
    --nano ont_reads.fastq.gz
```

### Quality Assessment with QUAST
```bash
# Evaluate assembly
quast assembly.fasta -o quast_results

# Compare multiple assemblies
quast assembly1.fasta assembly2.fasta -o comparison
```

### Sequence Processing with seqtk
```bash
# Sample 10% of reads
seqtk sample reads.fastq.gz 0.1 > subset.fastq

# Convert FASTQ to FASTA
seqtk seq -a reads.fastq.gz > reads.fasta

# Get sequence statistics
seqtk comp reads.fasta
```

### Alignment with Minimap2
```bash
# Align long reads to reference
minimap2 -ax map-hifi reference.fasta reads.fastq.gz | samtools sort -o aligned.bam

# Align assembly to reference
minimap2 -ax asm5 reference.fasta assembly.fasta > alignment.sam
```

## Resource Recommendations

- **CPU**: Most tools benefit from multiple cores (use `-t` flag)
- **Memory**: 
  - Small genomes (<100 Mb): 8-16 GB
  - Bacterial genomes: 16-32 GB
  - Large genomes (>1 Gb): 64-128 GB or more
- **Storage**: Depends on dataset size, but allow 3-5x the input data size

## Docker Resource Limits

```bash
# Run with memory limit
docker run -it --rm --memory=32g -v $(pwd):/data genome-assembly:latest

# Run with CPU limit
docker run -it --rm --cpus=8 -v $(pwd):/data genome-assembly:latest

# Combined
docker run -it --rm --memory=32g --cpus=8 -v $(pwd):/data genome-assembly:latest
```

## Notes

- GUI tools (Bandage, IGV) are not included in this command-line focused image
- All tools are installed via bioconda/conda for reproducibility
- The image uses a multistage build to minimize final size
- Data directory `/data` is set as the working directory
- NanoPlot is a Python tool with heavy dependencies (numpy, pandas, matplotlib) so it may take a few seconds to start

## Troubleshooting

### Permission Issues
If you encounter permission issues with mounted volumes:
```bash
docker run -it --rm -u $(id -u):$(id -g) -v $(pwd):/data genome-assembly:latest
```

### Check Installed Versions
```bash
docker run --rm genome-assembly:latest bash -c \
    "hifiasm --version && verkko --version && samtools --version"
```
