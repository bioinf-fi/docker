# scRNA-seq workshop container

This image provides an interactive Linux terminal with STAR, HISAT2,
SAMtools, featureCounts, the NCBI SRA Toolkit, `curl`/`wget`, and the Xpdf
graphical PDF viewer. It follows the X11 workflow used by the earlier
`bioinf-fi/docker` assembly workshop.

The current host folder is mounted at `/data`, so downloaded reads, reference
files, indexes, alignments, and counts persist after the container exits.

## Requirements

- Linux desktop with X11 or XWayland
- Docker or Podman
- At least 16 GB RAM for the supplied STAR indexing command
- Substantial free disk space: each SRA object is roughly 3 GB, the plain
  FASTQ is much larger, and `fasterq-dump` needs additional scratch space

On Docker Desktop, raise the virtual machine's memory above the 12 GB passed
to `--limitGenomeGenerateRAM` before building a STAR mouse-genome index.

## Build and start

```bash
make build
make test
make run
```

Inside the container:

```bash
show-tools
```

To use a specific runtime when both are installed:

```bash
make run CONTAINER_RUNTIME=podman
```

## Jupyter notebook

The image embeds `notebooks/pydeseq2_notebook.ipynb` with a preinstalled
`pydeseq2` kernel (pydeseq2, pandas, scanpy, decoupler, pertpy, etc.). After
building:

```bash
make notebook
```

Open http://127.0.0.1:8888/lab. Missing notebooks are copied from the image
into `./notebooks` on the host mount so edits persist (including the embedded
`metadata.tsv`). The notebook also expects count tables at `../hisat2/`
relative to that folder (i.e. `./hisat2` on the host). Override the published
port with `make notebook NOTEBOOK_PORT=8889`.

## X11 and the paper

Check the host, then start the image with X11 forwarding:

```bash
make test-x11
make gui
```

Inside the container:

```bash
download-paper
xpdf /data/paper.pdf &
```

Alternatively, download the paper before starting the GUI container:

```bash
make paper
make gui
# inside the container
xpdf paper.pdf &
```

If authenticated X11 forwarding is unavailable, use `make gui-simple`. This
temporarily runs `xhost +local:` and restores access control when the
container exits.

This setup targets Linux workshop machines, as did the earlier container.
macOS requires an X server such as XQuartz and different `DISPLAY` handling.

## Download the workshop runs

Inside the container, download one run first to confirm disk capacity and
network performance:

```bash
download-sra SRR11449097
```

Then fetch the remaining listed runs if needed:

```bash
download-sra SRR11449098 SRR11449111 SRR11449112
```

The output is placed in `/data/data`, the SRA cache in `/data/sra-cache`, and
scratch files in `/data/tmp`. The defaults can be changed, for example:

```bash
THREADS=8 TMPDIR=/data/fast-scratch download-sra SRR11449097
```

These 10x Genomics v2 SRA records contain three read segments. NCBI reports
8 + 26 + 57 bp for SRR11449097/98 and 8 + 27 + 50 bp for
SRR11449111/12. `download-sra` intentionally uses
`fasterq-dump --concatenate-reads --include-technical`, producing one 91 or
85 bp record per spot. The workshop's HISAT2 option `-5 34` then reproduces
the supplied trimming and leaves 57 or 51 bp, respectively. A normal
`fasterq-dump` call would omit the technical segments and would therefore
not reproduce the shown command.

The example alignment accessions in the supplied notes are different from
the four runs in the DATA section. To reproduce those examples, also run:

```bash
download-sra SRR11449109 SRR11449101
```

## Reference files

Put the instructor-provided mouse reference files here on the host:

```text
reference/GRCm39_genomic.fna
reference/GRCm39_genomic.gtf
```

They will appear at the same paths under `/data` in the container. Reference
data are deliberately not copied into the image.

## STAR index

Inside the container:

```bash
mkdir -p reference/star

STAR \
  --runThreadN 6 \
  --limitGenomeGenerateRAM 12000000000 \
  --runMode genomeGenerate \
  --genomeDir /data/reference/star \
  --genomeFastaFiles /data/reference/GRCm39_genomic.fna \
  --sjdbGTFfile /data/reference/GRCm39_genomic.gtf \
  --sjdbOverhang 50
```

The supplied `--sjdbOverhang 50` corresponds to the 51-base post-trimming
read length in the SRR11449109/SRR11449101 examples. If the workshop is
changed to use only the 57-base reads from SRR11449097/98, use 56 when
building a new STAR index.

## HISAT2, SAMtools, and featureCounts

Build the HISAT2 index:

```bash
mkdir -p reference/hisat2 alignments counts
hisat2-build \
  /data/reference/GRCm39_genomic.fna \
  /data/reference/hisat2/GRCm39
```

Align the examples. The `-x` before the index prefix is required by HISAT2:

```bash
hisat2 \
  -x /data/reference/hisat2/GRCm39 \
  -p 6 \
  -U /data/data/SRR11449109.fastq \
  -5 34 \
  -S /data/alignments/SRR11449109.sam

hisat2 \
  -x /data/reference/hisat2/GRCm39 \
  -p 6 \
  -U /data/data/SRR11449101.fastq \
  -5 34 \
  -S /data/alignments/SRR11449101.sam
```

Convert, sort, index, and count:

```bash
samtools view -@ 6 -bS /data/alignments/SRR11449101.sam \
  | samtools sort -@ 6 -o /data/alignments/SRR11449101_sorted.bam

samtools index /data/alignments/SRR11449101_sorted.bam

featureCounts \
  -T 6 \
  -a /data/reference/GRCm39_genomic.gtf \
  -o /data/counts/read_counts_101.txt \
  /data/alignments/SRR11449101_sorted.bam
```

## Cell Ranger

Cell Ranger is not included. It is not needed for the planned exercise, is a
large download, and 10x Genomics requires users to accept its End User
Software License Agreement before downloading it. If it becomes necessary,
download it separately from 10x Genomics and place the extracted directory
under this mounted workshop folder.

## Links

- [BioProject PRJNA616273](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA616273)
- [Paper (PMID 33526011)](https://pubmed.ncbi.nlm.nih.gov/33526011/)
- [10x Chromium 3′ library structure](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)
- [Cell Ranger installation](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in)
- [NCBI fasterq-dump guide](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
