#
# Download data via the curl downloader.
#

# The URL to be downloaded
URL ?= https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz

# Destination file.
FILE ?= refs/$(notdir ${URL})

# Additions flags for curl
FLAGS ?= 

# Makefile customizations.
SHELL := bash
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# General usage information.
usage:
	@echo "#"
	@echo "# curl.mk: download data"
	@echo "#"
	@echo "# URL=${URL}"
	@echo "# FILE=${FILE}"
	@echo "#"
	@echo "# make run|test|clean"
	@echo "#"

# Download the file with aria2c.
${FILE}:
	mkdir -p $(dir ${FILE})
	curl ${FLAGS} ${URL} > ${FILE}

# Run the download.
run: ${FILE}
	@ls -lh ${FILE}

# Remove the file.
run!:
	rm -f ${FILE}

clean: run!

# Installation instructions
install:
	@echo "# no installation required"

.PHONY:: usage run run! clean install
