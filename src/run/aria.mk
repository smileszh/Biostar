#
# Download data via the aria2c downloader.
#

# The URL to be downloaded
URL ?= https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz

# The directory to store the data.
DIR ?= data

# Destination file.
FILE ?= ${DIR}/$(notdir ${URL})

# Aria2c flags.
FLAGS = -x 5 -c --summary-interval=10

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# General usage information.
usage:
	@echo "#"
	@echo "# aria.mk: downloads data"
	@echo "#"
	@echo "# URL=${URL}"
	@echo "# FILE=${FILE}"
	@echo "#"
	@echo "#  make run"
	@echo "#"

# Download the file with aria2c.
${FILE}:
	# Create the directory.
	mkdir -p $(dir ${FILE})

	# Download the file.
	aria2c ${FLAGS} -o ${FILE} ${URL}

run: ${FILE}
	@ls -lh $<

# Remove the file.
clean:
	rm -f ${FILE}

run!: clean

# A shortcut to default run
test: clean run

# Installation instructions
install:
	@echo "micromamba install aria2"

.PHONY:: usage run run! install test
