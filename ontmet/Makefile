# Makefile for building and exporting genomics container

IMAGE_NAME := bioinf-fi/ontmet
TAG        := latest
OUTPUT     := ontmet.tar.gz
ENGINE     ?= docker

.PHONY: all build squash import clean

all: build

build:
	$(ENGINE) build --load -t $(IMAGE_NAME):$(TAG) -f ontmet.Dockerfile .

squash:
	echo "Squashing $(ENGINE) image to $(OUTPUT)"
	@cid=$$($(ENGINE) create $(IMAGE_NAME):$(TAG)); \
	$(ENGINE) export $$cid | gzip > $(OUTPUT); \
	$(ENGINE) rm $$cid; \
	echo "Squashed image written to $(OUTPUT)"

import:
	gunzip -c $(OUTPUT) | $(ENGINE) import - $(IMAGE_NAME):$(TAG)

clean:
	rm -f $(OUTPUT)