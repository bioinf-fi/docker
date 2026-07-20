.PHONY: help build run shell gui gui-simple test test-x11 versions paper clean

CONTAINER_RUNTIME ?= $(shell command -v docker 2>/dev/null || command -v podman 2>/dev/null)
IMAGE_NAME ?= scrna-workshop
IMAGE_TAG ?= latest
FULL_IMAGE := $(IMAGE_NAME):$(IMAGE_TAG)
DATA_DIR ?= $(CURDIR)
HOST_UID := $(shell id -u)
HOST_GID := $(shell id -g)

ifeq ($(strip $(CONTAINER_RUNTIME)),)
$(error Neither Docker nor Podman was found. Install one, then try again.)
endif

define RUN_BASE
$(CONTAINER_RUNTIME) run --rm -it --init \
	--user $(HOST_UID):$(HOST_GID) \
	-e HOME=/tmp/workshop-home \
	-e TERM="$${TERM:-xterm-256color}" \
	-v "$(DATA_DIR):/data:Z,U"
endef

help:
	@echo "scRNA-seq workshop container"
	@echo ""
	@echo "Container runtime: $(CONTAINER_RUNTIME)"
	@echo "Workshop folder:  $(DATA_DIR) -> /data"
	@echo ""
	@echo "make build       Build the image"
	@echo "make run         Start an interactive terminal"
	@echo "make gui         Start a terminal with authenticated X11 forwarding"
	@echo "make gui-simple  Start X11 using relaxed local access"
	@echo "make test-x11    Check the host's X11 setup"
	@echo "make test        Verify all required commands in the image"
	@echo "make versions    Print installed tool versions"
	@echo "make paper       Download the workshop paper into this folder"
	@echo "make clean       Remove the local image"
	@echo ""
	@echo "Override the runtime if needed, for example: make run CONTAINER_RUNTIME=podman"

build:
	$(CONTAINER_RUNTIME) build -t $(FULL_IMAGE) -f Dockerfile .

run:
	$(RUN_BASE) $(FULL_IMAGE)

shell: run

gui:
	@if [ -z "$${DISPLAY:-}" ]; then \
		echo "DISPLAY is not set. Run 'make test-x11' for diagnostics."; \
		exit 1; \
	fi
	@if ! command -v xauth >/dev/null 2>&1; then \
		echo "xauth is required for authenticated mode; use 'make gui-simple' or install xauth."; \
		exit 1; \
	fi
	@xauth_file=$$(mktemp); \
	trap 'rm -f "$$xauth_file"' EXIT INT TERM; \
	xauth_entries=$$(xauth nlist "$$DISPLAY" 2>/dev/null || true); \
	if [ -z "$$xauth_entries" ]; then \
		echo "No X11 authorization entry was found for $$DISPLAY."; \
		echo "Try 'make gui-simple'."; \
		exit 1; \
	fi; \
	printf '%s\n' "$$xauth_entries" | sed -e 's/^..../ffff/' | xauth -f "$$xauth_file" nmerge -; \
	$(RUN_BASE) \
		-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
		-v "$$xauth_file:/tmp/.Xauthority:ro" \
		-e DISPLAY="$$DISPLAY" \
		-e XAUTHORITY=/tmp/.Xauthority \
		--network=host \
		--security-opt label=disable \
		$(FULL_IMAGE)

gui-simple:
	@if [ -z "$${DISPLAY:-}" ]; then \
		echo "DISPLAY is not set. Run 'make test-x11' for diagnostics."; \
		exit 1; \
	fi
	@set -e; \
	xhost +local: >/dev/null; \
	trap 'xhost -local: >/dev/null 2>&1 || true' EXIT INT TERM; \
	$(RUN_BASE) \
		-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
		-e DISPLAY="$$DISPLAY" \
		--network=host \
		--security-opt label=disable \
		$(FULL_IMAGE)

test:
	$(CONTAINER_RUNTIME) run --rm $(FULL_IMAGE) /bin/bash -lc '\
		for tool in STAR hisat2 hisat2-build samtools featureCounts prefetch fasterq-dump xpdf; do \
			command -v "$$tool" >/dev/null || exit 1; \
		done; \
		show-tools'

test-x11:
	@bash ./scripts/test-x11

versions:
	$(CONTAINER_RUNTIME) run --rm $(FULL_IMAGE) show-tools

paper:
	$(CONTAINER_RUNTIME) run --rm \
		--user $(HOST_UID):$(HOST_GID) \
		-e HOME=/tmp/workshop-home \
		-v "$(DATA_DIR):/data:Z" \
		$(FULL_IMAGE) download-paper /data/paper.pdf

clean:
	$(CONTAINER_RUNTIME) image rm $(FULL_IMAGE)
