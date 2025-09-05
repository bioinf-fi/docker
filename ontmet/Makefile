# Makefile for building and exporting genomics container

IMAGE_NAME := bioinf-fi/ontmet
TAG        := latest
OUTPUT     := ontmet.tar.gz
ENGINE     ?= docker
DATA_DIR   ?= $(PWD)

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
	@echo ""
	@echo "Examples:"
	@echo "  make build                    Build the container"
	@echo "  make run-shell                Interactive shell"
	@echo "  make igv                      Launch IGV GUI"
	@echo "  make run-igv DATA_DIR=/data   Run with custom data dir"

build: ## Build the genomics container image
	$(ENGINE) build --load -t $(IMAGE_NAME):$(TAG) -f ontmet.Dockerfile .

squash: ## Export container to compressed tar archive
	echo "Squashing $(ENGINE) image to $(OUTPUT)"
	@cid=$$($(ENGINE) create $(IMAGE_NAME):$(TAG)); \
	$(ENGINE) export $$cid | gzip > $(OUTPUT); \
	$(ENGINE) rm $$cid; \
	echo "Squashed image written to $(OUTPUT)"

import: ## Import container from compressed tar archive
	gunzip -c $(OUTPUT) | $(ENGINE) import - $(IMAGE_NAME):$(TAG)

# Run container with X11 forwarding for GUI applications like IGV
run-igv: ## Run interactive container with GUI support (for IGV)
	@echo "Starting container with IGV support..."
	@echo "Make sure to run 'xhost +local:root' first if you get permission errors"
	$(ENGINE) run -it --rm \
		-e DISPLAY=$$DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix:rw \
		-v "$(DATA_DIR)":/data \
		$(IMAGE_NAME):$(TAG)

# Launch IGV directly (convenience target)
igv: ## Launch IGV GUI directly
	@echo "Launching IGV in container..."
	@echo "Make sure X11 forwarding is enabled: xhost +local:root"
	$(ENGINE) run -it --rm \
		-e DISPLAY=$$DISPLAY \
		-v /tmp/.X11-unix:/tmp/.X11-unix:rw \
		-v "$(DATA_DIR)":/data \
		$(IMAGE_NAME):$(TAG) \
		bash -l -c "igv"

# Run interactive shell (no X11)
run-shell: ## Run interactive shell with data mounting
	$(ENGINE) run -it --rm \
		-v "$(DATA_DIR)":/data \
		$(IMAGE_NAME):$(TAG)

clean: ## Remove generated tar archive
	rm -f $(OUTPUT)

tools: ## Show all available tools in the container
	$(ENGINE) run --rm $(IMAGE_NAME):$(TAG) bash -l -c "tools"

# macOS-specific targets
run-igv-macos: ## Run container with XQuartz support for macOS
	@echo "Starting container with macOS XQuartz support..."
	@echo "Make sure XQuartz is running and network connections are enabled"
	@IP=$$(ifconfig en0 | grep inet | awk '$$1=="inet" {print $$2}') && \
	echo "Using IP: $$IP" && \
	$(ENGINE) run -it --rm \
		-e DISPLAY=$$IP:0 \
		-v "$(DATA_DIR)":/data \
		$(IMAGE_NAME):$(TAG)

igv-vnc: ## Launch IGV via VNC server (cross-platform)
	@echo "Starting IGV with VNC server..."
	@echo "Connect to localhost:5901 with VNC client after startup"
	$(ENGINE) run -it --rm \
		-p 5901:5901 \
		-v "$(DATA_DIR)":/data \
		$(IMAGE_NAME):$(TAG) \
		bash -c "apt update && apt install -y tightvncserver xfce4 && vncserver :1 -geometry 1024x768 -depth 24 && echo 'VNC server started. Connect to localhost:5901' && DISPLAY=:1 igv"