# Makefile for building and exporting genomics container

IMAGE_NAME := bioinf-fi/ontmet
TAG        := latest
OUTPUT     := ontmet.tar.gz
ENGINE     ?= docker
DATA_DIR   ?= $(PWD)
ARGS       ?=

.PHONY: all build squash import clean run-igv igv run-shell help tools

all: build

help: ## Show this help message
	@echo "Available targets:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "Variables:"
	@echo "  ENGINE=docker|podman  Container engine (default: docker)"
	@echo "  DATA_DIR=<path>       Data directory to mount (default: current dir)"
	@echo "  ARGS=<build-args>     Additional build arguments for version control"
	@echo ""
	@echo "Container Engine Selection:"
	@echo "  make build ENGINE=docker              Use Docker (default)"
	@echo "  make build ENGINE=podman              Use Podman instead"
	@echo ""
	@echo "Build Examples:"
	@echo "  make build                            Build with default versions"
	@echo "  make build ARGS='--build-arg MINIMAP2_VERSION=2.31'"
	@echo "  make build ARGS='--build-arg PYTHON_VERSION=3.9'"
	@echo ""
	@echo "Usage Examples:"
	@echo "  make run-shell                        Interactive shell"
	@echo "  make run-shell DATA_DIR=/path/to/data Mount custom data directory"
	@echo "  make igv                              Launch IGV GUI (Linux X11 only)"
	@echo "  make tools                            List all available tools"

build: ## Build the genomics container image with parallel optimization
	$(ENGINE) build --load -t $(IMAGE_NAME):$(TAG) -f ontmet.Dockerfile $(ARGS) .

squash: ## Export container to compressed tar archive for distribution
	echo "Squashing $(ENGINE) image to $(OUTPUT)"
	@cid=$$($(ENGINE) create $(IMAGE_NAME):$(TAG)); \
	$(ENGINE) export $$cid | gzip > $(OUTPUT); \
	$(ENGINE) rm $$cid; \
	echo "Squashed image written to $(OUTPUT)"

import: ## Import container from compressed tar archive
	gunzip -c $(OUTPUT) | $(ENGINE) import - $(IMAGE_NAME):$(TAG)

# Run container with X11 forwarding for GUI applications like IGV (Linux)
run-igv: ## Run interactive container with GUI support (Linux X11)
	@echo "Starting container with IGV support (Linux X11)..."
	@echo "Make sure to run 'xhost +local:root' first if you get permission errors"
	$(ENGINE) run -it --rm \
		-e DISPLAY=$$DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix:rw \
		-v "$(DATA_DIR)":/data:Z,U \
		$(IMAGE_NAME):$(TAG)

# Launch IGV directly (convenience target)
igv: ## Launch IGV GUI directly (Linux X11)
	@echo "Launching IGV in container (Linux X11)..."
	@echo "Make sure X11 forwarding is enabled: xhost +local:root"
	$(ENGINE) run -it --rm \
		-e DISPLAY=$$DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix:rw \
		-v "$(DATA_DIR)":/data:Z,U \
		$(IMAGE_NAME):$(TAG) \
		bash -l -c "igv"

# Run interactive shell (no X11)
run-shell: ## Run interactive shell with data mounting
	$(ENGINE) run -it --rm \
		-v "$(DATA_DIR)":/data:Z,U \
		$(IMAGE_NAME):$(TAG)

clean: ## Remove generated tar archive
	rm -f $(OUTPUT)

tools: ## Show all available tools in the container
	$(ENGINE) run --rm $(IMAGE_NAME):$(TAG) bash -l -c "tools"
