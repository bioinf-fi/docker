.PHONY: help build run shell notebook gui gui-simple test test-x11 versions paper clean

CONTAINER_RUNTIME ?= $(shell command -v docker 2>/dev/null || command -v podman 2>/dev/null)
IMAGE_NAME ?= scrna-workshop
IMAGE_TAG ?= latest
FULL_IMAGE := $(IMAGE_NAME):$(IMAGE_TAG)
DATA_DIR ?= $(CURDIR)
HOST_UID := $(shell id -u)
HOST_GID := $(shell id -g)
NOTEBOOK_PORT ?= 8888

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
	@echo "make notebook    Start Jupyter Lab with the embedded workshop notebook"
	@echo "make gui         Start a terminal with authenticated X11 forwarding"
	@echo "make gui-simple  Start X11 using relaxed local access"
	@echo "make test-x11    Check the host's X11 setup"
	@echo "make test        Verify all required commands in the image"
	@echo "make versions    Print installed tool versions"
	@echo "make paper       Download the workshop paper into this folder"
	@echo "make clean       Remove the local image"
	@echo ""
	@echo "Override the runtime if needed, for example: make run CONTAINER_RUNTIME=podman"
	@echo "Override the notebook port if needed: make notebook NOTEBOOK_PORT=8889"

build:
	$(CONTAINER_RUNTIME) build -t $(FULL_IMAGE) -f Dockerfile .

run:
	$(RUN_BASE) $(FULL_IMAGE)

shell: run

notebook:
	@echo "Jupyter Lab will open http://127.0.0.1:$(NOTEBOOK_PORT)/lab/tree/pydeseq2_notebook.ipynb"
	@echo "Embedded notebooks are seeded into $(DATA_DIR)/notebooks when missing."
	$(RUN_BASE) \
		-p $(NOTEBOOK_PORT):8888 \
		-e JUPYTER_PORT=8888 \
		$(FULL_IMAGE) start-notebook

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
		for tool in STAR hisat2 hisat2-build samtools featureCounts prefetch fasterq-dump xpdf jupyter; do \
			command -v "$$tool" >/dev/null || exit 1; \
		done; \
		jupyter kernelspec list | grep -F pydeseq2 >/dev/null; \
		python -c 'import pydeseq2, pandas, numpy, matplotlib, seaborn, scipy, scanpy, decoupler, pertpy'; \
		test -f /opt/workshop/notebooks/pydeseq2_notebook.ipynb; \
		test -f /opt/workshop/notebooks/metadata.tsv; \
		test -f /opt/workshop/notebooks/read_counts_all_s1.txt; \
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
