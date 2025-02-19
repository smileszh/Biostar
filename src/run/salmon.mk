#
# Generates alignments with salmon
#

# FASTQ read pair.
R1 ?= reads/read_1.fastq

# Second in pair.
R2 ?=

# Set MODE to SE if R2 is empty otherwise PE
ifeq ($(R2),)
	MODE ?= SE
else
	MODE ?= PE
endif

# Normally it should be a transcriptome file!
REF ?= genome.fa

# The name of the sample that is processed
SAMPLE = $(basename $(notdir ${R1}))

# The name of the program to check
PROG = salmon

# Name of the output directory.
PREFIX ?= salmon/${SAMPLE}

# The name of the quantification file.
QUANT ?= ${PREFIX}/quant.sf

# The directory that holds the index.
IDX_DIR = $(dir ${REF})/idx

# The name of the index
IDX ?= ${IDX_DIR}/$(notdir ${REF}).salmon

# The index file name
IDX_FILE ?= ${IDX}/info.json

# Number of CPUS
NCPU ?= 2

# Additional salmon options.
SALMON_FLAGS = -l A

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage:
	@echo "#"
	@echo "# salmon.mk: classify reads using salmon "
	@echo "#"
	@echo "# MODE=${MODE}"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# SAMPLE=${SAMPLE}"
	@echo "#"
	@echo "# make index align"
	@echo "#"

${REF}:
	@if [ ! -f ${REF} ]; then
		echo "# reference not found: REF=${REF}";
		exit -1
	fi

${R1}:
	@if [ ! -f ${R1} ]; then
		echo "# read not found: REF=${R1}";
		exit -1
	fi

# Index the reference.
${IDX}: ${REF}
	@command -v ${PROG} > /dev/null 2>&1 || (echo "# ${PROG} not installed" && exit 1)
	mkdir -p $(dir ${IDX})
	salmon index --threads ${NCPU} -t ${REF} -i ${IDX}

# Generate the index.
index: ${IDX}
	@ls -lh ${IDX_FILE}

# Remove the index.
index!: ${IDX}
	rm -rf ${IDX_FILE}

# Paired end mode classification.
ifeq ($(MODE), PE)
CMD = salmon quant -q ${SALMON_FLAGS} -i ${IDX}  -1 ${R1} -2 ${R2} --threads ${NCPU} -o ${PREFIX}
else
CMD = salmon quant -q ${SALMON_FLAGS} -i ${IDX}  -r ${R1} --threads ${NCPU} -o ${PREFIX}
endif

# Generate the quantification matrix.
${QUANT}: ${R1} ${R2}
	@command -v ${PROG} > /dev/null 2>&1 || (echo "# ${PROG} not installed" && exit 1)
	@if [ ! -f ${IDX_FILE} ]; then
		echo "# salmon index missing: IDX=${IDX}";
		exit -1
	fi
	mkdir -p ${PREFIX}
	${CMD}

# Target to trigger the alignment.
run: ${QUANT}
	@ls -lh ${QUANT}

# Remove the BAM file.
run!:
	rm -rf ${QUANT}

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_aligner MOD=src/run/salmon.mk

# Targets that are not files.
.PHONY: run run! install usage index test

# Install required software.
install:
	@echo micromamba install salmon

