#
# Main testing makefile
#

TEST_DIR = testrun

# The follow
N = 10000

SRR = SRR1553425

ACC = AF086833
REF = refs/${ACC}.fa
R1 = reads/${SRR}_1.fastq
R2 = reads/${SRR}_2.fastq
BAM = bam/${SRR}
VCF = vcf/${SRR}

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

VERSION = 2024.10.25

usage:
	@echo "#"
	@echo "# Testing the Biostar Toolbox"
	@echo "#"
	@echo "# TEST_DIR=${TEST_DIR}"
	@echo "#"
	@echo "# VERSION=${VERSION}"
	@echo "#"
	@echo "# Use the Source, Luke!"
	@echo "#"
	@echo "# make test         # tests main modules"
	@echo "# make recipes      # tests the recipes"
	@echo "# make workflows    # tests the workflows"
	@echo "# make deepvariant  # tests the deepvariant module (unix only)"
	@echo "#"
	@echo "# make all          # runs all tests"


# Read names.
R1 = reads/${SRR}_1.fastq
R2 = reads/${SRR}_2.fastq

# BAM file name
BAM = bam/${SRR}.bam

setup:
	mkdir -p ${TEST_DIR} && cd ${TEST_DIR} && ln -sf ../src .

# Module tests.
modules: setup
	cd ${TEST_DIR}
	make -f src/run/datasets.mk run
	make -f src/run/sra.mk run SRR=${SRR} N=${N}
	make -f src/run/fastp.mk run! run P1=reads/${SRR}_1.fastq P2=reads/${SRR}_2.fastq
	make -f src/run/genbank.mk fasta ACC=${ACC}

	# Test paired end run.
	make -f src/run/bwa.mk test
	make -f src/run/bowtie2.mk test
	make -f src/run/minimap2.mk test
	make -f src/run/hisat2.mk test
	make -f src/run/wiggle.mk test
	make -f src/run/bcftools.mk test
	make -f src/run/freebayes.mk test
	micromamba run -n salmon_env make -f src/run/salmon.mk test
	make -f src/run/aria.mk test

# This can only be run on UNIX
deepvariant: setup
	cd ${TEST_DIR}
	make -f src/run/deepvariant.mk test

# Test the workflows.
recipes: setup
	cd ${TEST_DIR}
	make -f src/recipes/short-read-alignments.mk run
	make -f src/recipes/variant-calling.mk run
	make -f src/recipes/rnaseq-with-hisat.mk run
	micromamba run -n salmon_env make -f src/recipes/rnaseq-with-salmon.mk run

# Test the workflows.
workflows: setup
	cd ${TEST_DIR}
	make -f src/workflows/airway.mk design run
	make -f src/workflows/presenilin.mk design run
	make -f src/workflows/snpcall.mk run

# Test all functionality.
all: modules salmon recipes workflows

test: modules

test!:
	rm -rf ${TEST_DIR}/bam ${TEST_DIR}/reads ${TEST_DIR}/vcf ${TEST_DIR}/refs ${TEST_DIR}/trim

install:
	@echo "bash src/setup/init-stats.sh"

.PHONY: usage test test! all modules salmon deepvariant recipes workflows
