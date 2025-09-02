# Makefile for building and exporting genomics container

IMAGE_NAME := bioinf-fi/ontmet
TAG        := latest
OUTPUT     := ontmet.tar.gz

# Tools can be switched (docker/podman)
ENGINE     ?= podman

.PHONY: all podman docker clean

all: podman

podman: ENGINE := podman
podman: build squash

docker: ENGINE := docker
docker: build squash

# 1. Build image
build:
	$(ENGINE) build -t $(IMAGE_NAME):$(TAG) -f ontmet.Dockerfile

# 2. Squash to single layer (export+import)
squash:
	@cid=$$($(ENGINE) create $(IMAGE_NAME):$(TAG)); \
	$(ENGINE) export $$cid | gzip > $(OUTPUT); \
	$(ENGINE) rm $$cid; \
	echo "Squashed image written to $(OUTPUT)"

# Optional: Load back into local engine (single-layer import)
import:
	gunzip -c $(OUTPUT) | $(ENGINE) import - $(IMAGE_NAME):$(TAG)

clean:
	rm -f $(OUTPUT)