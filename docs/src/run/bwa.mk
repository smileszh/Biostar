#
# Generate alignments with bwa
#

# The reference genome.
REF ?= refs/genome.fa

# The directory that holds the index.
IDX_DIR = $(dir ${REF})bwaIndex

# The name of the index
IDX ?= ${IDX_DIR}/$(notdir ${REF})

# A file in the index directory.
IDX_FILE ?= ${IDX}.ann

# Number of CPUS
NCPU ?= 2

# Additional flags to pass to BWA.
FLAGS ?= -t ${NCPU}

# First in pair.
R1 ?= reads/reads1.fq

# Second in pair. Runs in single end mode if not defined.
R2 ?=

# The alignment file.
BAM ?= bam/alignment.bam

# Set the values for the read groups.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA

# Build the read groups tag.
RG ?= '@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:${PL}'

# Makefile customizations.
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage:
	@echo "#"
	@echo "# bwa.mk: align reads using BWA"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# IDX=${IDX} (optional)"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make index|run|test|clean"
	@echo "#"

# Read 1 must exist.
${R1}:
	@echo "#"
	@echo "# Read 1 file not found: R1=${R1}"
	@echo "#"
	@exit -1

# If R2 is set, it must exist.
ifneq (${R2},)
${R2}:
	@echo "#"
	@echo "# Read 2 file not found: R2=${R2}"
	@echo "#"
	@exit -1
endif

# Bail out if reference is not found
${REF}:
	@echo "#"
	@echo "# Reference genome not found: REF=${REF}";
	@echo "#"
	@exit -1

# Index the reference genome.
${IDX_FILE}: ${REF}
	@mkdir -p $(dir $@)
	bwa index -p ${IDX} ${REF}

# Create the index.
index: ${IDX_FILE}
	@echo "# bwa index: ${IDX}"

# Remove the index.
index!:
	rm -rf ${IDX_FILE}

# Paired end alignment.
${BAM}: ${R1} ${R2}
	@if [ ! -f ${IDX_FILE} ]; then
		echo "# bwa index not found: IDX=${IDX}";
		exit -1
	fi
	mkdir -p $(dir $@)
	bwa mem ${FLAGS} -R ${RG} ${IDX} ${R1} ${R2} | samtools sort -@ ${NCPU} > ${BAM}

# Create the BAM index file.
${BAM}.bai: ${BAM}
	samtools index ${BAM}

# Generate the alignment.
align: ${BAM}.bai
	@ls -lh ${BAM}

# Generate the alignment.
align!:
	rm -f ${BAM} ${BAM}.bai

# Create a wiggle coverage file.
coverage: ${BAM}.bai ${REF}
	@make -f src/run/coverage.mk BAM=${BAM} REF=${REF} run

# A synonym for coverage.
wiggle: coverage

# Running creates the wiggle file as well.
run: align coverage

# Remove the BAM file.
run!:
	rm -rf ${BAM} ${BAM}.bai

clean: run!

# Run the test.
test:
	make -f src/run/tester.mk test_aligner MOD=src/run/bwa.mk

# The name of the stats file.
STATS = $(basename ${BAM}).stats

# Generate stats on the alignme
${STATS}: ${BAM}.bai
	samtools flagstat ${BAM} > ${STATS}

# Trigger the statistics generation.
stats: ${STATS}
	@echo "# ${STATS}"
	@cat ${STATS}

# Install required software.
install::
	@echo micromamba install bwa samtools



# Targets that are not files.
.PHONY: usage index align run run! test wiggle install
