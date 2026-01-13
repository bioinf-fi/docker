.PHONY: build run test clean help

# Auto-detect container runtime (docker or podman)
CONTAINER_RUNTIME := $(shell command -v podman 2>/dev/null || command -v docker 2>/dev/null)
ifeq ($(CONTAINER_RUNTIME),)
    $(error Neither Docker nor Podman found. Please install one of them.)
endif

IMAGE_NAME = genome-assembly
IMAGE_TAG = latest
FULL_IMAGE = $(IMAGE_NAME):$(IMAGE_TAG)

help:
	@echo "Genome Assembly Container Image - Makefile Commands"
	@echo ""
	@echo "Detected container runtime: $(CONTAINER_RUNTIME)"
	@echo ""
	@echo "Available targets:"
	@echo "  make build      - Build the container image"
	@echo "  make run        - Run interactive container (shows tools on startup)"
	@echo "  make gui        - Run container with GUI support (for IGV, Bandage)"
	@echo "  make gui-simple - Run GUI with relaxed X11 security (if gui fails)"
	@echo "  make test-x11   - Test if your system is ready for GUI apps"
	@echo "  make test       - Test that all tools are installed correctly"
	@echo "  make shell      - Alias for 'run'"
	@echo "  make versions   - Show versions of all installed tools"
	@echo "  make clean      - Remove the container image"
	@echo ""
	@echo "Inside container: run 'show-tools' to display tool versions"
	@echo "For GUI apps: Run 'make test-x11' first, then 'make gui'"
	@echo "If X11 auth fails: Try 'make gui-simple' instead"
	@echo ""

build:
	@echo "Building $(FULL_IMAGE) using $(CONTAINER_RUNTIME)..."
	$(CONTAINER_RUNTIME) build -f assembly.Dockerfile -t $(FULL_IMAGE) .
	@echo ""
	@echo "Build complete! Run 'make test' to verify installation."

run:
	@echo "Starting interactive container..."
	$(CONTAINER_RUNTIME) run -it --rm -v $$(pwd):/data:Z $(FULL_IMAGE)

gui:
	@echo "Starting interactive container with GUI support..."
	@echo "You can now run: igv, Bandage"
	@if [ -f "$$HOME/.Xauthority" ]; then \
		$(CONTAINER_RUNTIME) run -it --rm \
			-v $$(pwd):/data:Z \
			-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
			-v $$HOME/.Xauthority:/root/.Xauthority:ro \
			-e DISPLAY=$(DISPLAY) \
			-e XAUTHORITY=/root/.Xauthority \
			--network=host \
			--security-opt label=disable \
			$(FULL_IMAGE); \
	else \
		echo "Warning: .Xauthority not found, using relaxed X11 access"; \
		xhost +local: 2>/dev/null || true; \
		$(CONTAINER_RUNTIME) run -it --rm \
			-v $$(pwd):/data:Z \
			-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
			-e DISPLAY=$(DISPLAY) \
			--network=host \
			--security-opt label=disable \
			$(FULL_IMAGE); \
	fi

gui-simple:
	@echo "Starting container with GUI support (relaxed X11 security)..."
	@echo "Note: Running 'xhost +local:' to allow container X11 access"
	@xhost +local: 2>/dev/null || echo "Warning: xhost command not found or failed"
	@echo "You can now run: igv, Bandage"
	$(CONTAINER_RUNTIME) run -it --rm \
		-v $$(pwd):/data:Z \
		-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
		-e DISPLAY=$(DISPLAY) \
		--network=host \
		--security-opt label=disable \
		$(FULL_IMAGE)
	@echo "Restoring X11 security..."
	@xhost -local: 2>/dev/null || true

shell: run

test:
	@echo "Testing tool installations..."
	@$(CONTAINER_RUNTIME) run --rm $(FULL_IMAGE) show-tools

test-x11:
	@echo "Testing X11/GUI setup on your system..."
	@./test-gui.sh

versions:
	@$(CONTAINER_RUNTIME) run --rm $(FULL_IMAGE) show-tools

clean:
	@echo "Removing container image $(FULL_IMAGE)..."
	$(CONTAINER_RUNTIME) rmi $(FULL_IMAGE)
	@echo "Image removed."

# Quick command examples
example-qc:
	@echo "Example: Quality control with NanoPlot"
	@echo "$(CONTAINER_RUNTIME) run --rm -v \$$(pwd):/data:Z $(FULL_IMAGE) NanoPlot --fastq reads.fastq.gz -o qc_output"

example-assembly:
	@echo "Example: Assembly with hifiasm"
	@echo "$(CONTAINER_RUNTIME) run --rm -v \$$(pwd):/data:Z $(FULL_IMAGE) hifiasm -o assembly -t 8 reads.fastq.gz"

example-quast:
	@echo "Example: Assembly quality assessment"
	@echo "$(CONTAINER_RUNTIME) run --rm -v \$$(pwd):/data:Z $(FULL_IMAGE) quast assembly.fasta -o quast_results"
