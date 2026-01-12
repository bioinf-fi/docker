.PHONY: build run test clean help

IMAGE_NAME = genome-assembly
IMAGE_TAG = latest
FULL_IMAGE = $(IMAGE_NAME):$(IMAGE_TAG)

help:
	@echo "Genome Assembly Docker Image - Makefile Commands"
	@echo ""
	@echo "Available targets:"
	@echo "  make build    - Build the Docker image"
	@echo "  make run      - Run interactive container (shows tools on startup)"
	@echo "  make test     - Test that all tools are installed correctly"
	@echo "  make shell    - Alias for 'run'"
	@echo "  make versions - Show versions of all installed tools"
	@echo "  make clean    - Remove the Docker image"
	@echo ""
	@echo "Inside container: run 'show-tools' to display tool versions"
	@echo ""

build:
	@echo "Building $(FULL_IMAGE)..."
	docker build -f assembly.Dockerfile -t $(FULL_IMAGE) .
	@echo ""
	@echo "Build complete! Run 'make test' to verify installation."

run:
	@echo "Starting interactive container..."
	docker run -it --rm -v $$(pwd):/data $(FULL_IMAGE)

shell: run

test:
	@echo "Testing tool installations..."
	@docker run --rm $(FULL_IMAGE) show-tools

versions:
	@docker run --rm $(FULL_IMAGE) show-tools

clean:
	@echo "Removing Docker image $(FULL_IMAGE)..."
	docker rmi $(FULL_IMAGE)
	@echo "Image removed."

# Quick command examples
example-qc:
	@echo "Example: Quality control with NanoPlot"
	@echo "docker run --rm -v \$$(pwd):/data $(FULL_IMAGE) NanoPlot --fastq reads.fastq.gz -o qc_output"

example-assembly:
	@echo "Example: Assembly with hifiasm"
	@echo "docker run --rm -v \$$(pwd):/data $(FULL_IMAGE) hifiasm -o assembly -t 8 reads.fastq.gz"

example-quast:
	@echo "Example: Assembly quality assessment"
	@echo "docker run --rm -v \$$(pwd):/data $(FULL_IMAGE) quast assembly.fasta -o quast_results"
