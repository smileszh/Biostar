#
# Generates alignments with minimap2
#

# First in pair.
R1 ?= reads/read1.fq

# Second in pair. Leave empty for single end library.
R2 ?=

# Number of CPUS to use during alignment.
NCPU ?= 2

# Number of CPUS to use during sorting.
NSORT ?= 2

# Additional minimap2 options. Example: -x asm20 --secondary=no --sam-hit-only
FLAGS =

# Reference genome.
REF ?= genome.fa

# The alignment file.
BAM ?= bam/minimap2.bam

# Set the values for the read groups.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA
RG ?= "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:${PL}"

# Makefile customizations.
SHELL := bash
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage::
	@echo "#"
	@echo "# minimap2.mk: align read using minimap2 "
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make run"
	@echo "#"

# Read 1 must exist.
${R1}:
	@echo "# Read 1 file not found: R1=${R1}"
	@exit -1

# If R2 is set, it must exist.
ifneq (${R2},)
${R2}:
	@echo "# Read 2 file not found: R2=${R2}"
	@exit -1
endif

# Bail out if reference is not found
${REF}:
	@echo "#"
	@echo "# Error! Missing reference genome: REF=${REF}";
	@echo "#"
	@exit -1

# Generate the index.
index: ${REF}
	@echo "# The minimap2 index is generated during alignment!"

# Remove the index.
index!:
	@echo "# The minimap2 index is generated during alignment!"

# Paired end alignment.
${BAM}: ${REF} ${R1} ${R2}
	mkdir -p $(dir $@)
	minimap2 -a --MD -t ${NCPU} ${FLAGS} -R ${RG} ${REF} ${R1} ${R2} |\
		 samtools sort -@ ${NSORT} > ${BAM}

# Create the BAM index file.
${BAM}.bai: ${BAM}
	samtools index ${BAM}

# Generate the alignment.
align: ${BAM}.bai
	@ls -lh ${BAM}

# Create a wiggle coverage file.
coverage: align
	make -f src/run/coverage.mk BAM=${BAM} REF=${REF} run

# Running creates the wiggle file as well.
run: coverage

# Remove the BAM file.
run!:
	rm -rf ${BAM} ${BAM}.bai
	make -f src/run/coverage.mk BAM=${BAM} run!

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_aligner MOD=src/run/minimap2.mk

# Install required software.
install::
	@echo micromamba install minimap2 samtools

# Targets that are not files.
.PHONY: align align! install usage index test
