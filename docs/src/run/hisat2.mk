#
# Generate alignments with hisat2
#

# Number of CPUS
NCPU ?= 2

# Additional flags to pass to HISAT.
HISAT2_FLAGS ?= --threads ${NCPU} --sensitive

# Read group information.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA

# Set the read groups.
RG ?= --rg-id ${ID} --rg SM:${SM} --rg LB:${LB} --rg PL:${PL}

# First in pair.
R1 ?= reads/read1.fq

# Second in pair.
R2 ?=

# The reference genome.
REF ?= refs/genome.fa

# The directory that holds the index.
IDX_DIR = $(dir ${REF})hisat2Index

# The name of the index
IDX ?= ${IDX_DIR}/$(notdir ${REF})

# A file in the index directory.
IDX_FILE ?= ${IDX}.1.ht2

# The alignment file.
BAM ?= bam/hisat2.bam

# Set the values for the read groups.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage:
	@echo "#"
	@echo "# hisat2.mk: align reads using hisat2"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# IDX=${IDX}"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make index run"
	@echo "#"

# Reads must exist.
${R1}:
	@echo "# file not found: R1=${R1}"
	@exit -1

${REF}:
	@echo "# file not found: REF=${REF}"
	@exit -1

# Build the index for the reference genome.
${IDX_FILE}: 
	@if [ ! -f ${REF} ]; then
		echo "# file not found: REF=${REF}";
		exit -1
	fi
	mkdir -p $(dir $@)

	# Unzip the reference if it is gzipped.
	if [ "$(suffix ${REF})" = ".gz" ]; then
		# Unzip if the file is only found in compressed format.
		[ ! -f "$${REF%.*}" ] && gunzip -k "$$REF"
		hisat2-build --threads ${NCPU} $(basename ${REF}) ${IDX};
	else
		hisat2-build --threads ${NCPU} ${REF} ${IDX};
	fi

# Create the index.
index: ${IDX_FILE}
	@echo "# hisat2 index: ${IDX}"

# Remove the index.
index!:
	rm -rf ${IDX_FILE}

# If R2 is empty, use the single-end command.
ifeq (${R2},)
CMD = hisat2 ${HISAT2_FLAGS} --threads ${NCPU} ${RG} -x ${IDX} -U ${R1}
else
CMD = hisat2 ${HISAT2_FLAGS} --threads ${NCPU} ${RG} -x ${IDX} -1 ${R1} -2 ${R2}
endif

# Perform the alignment.
${BAM}: ${R1} ${R2}
	@if [ ! -f ${IDX_FILE} ]; then
		echo "# hisat2 index not found: IDX=${IDX}";
		exit -1
	fi
	@mkdir -p $(dir $@)
	${CMD} | samtools sort -@ ${NCPU} > $@

# Create the BAM index file.
${BAM}.bai: ${BAM}
	samtools index ${BAM}

# Generate the alignment.
align: ${BAM}.bai
	@ls -lh ${BAM}

# Create a wiggle coverage file.
wiggle: align ${REF}
	@make -f src/run/wiggle.mk BAM=${BAM} REF=${REF} run

# Display the BAM file path.
run: wiggle

# Remove the BAM file.
run!:
	rm -rf ${BAM} ${BAM}.bai 

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_aligner MOD=src/run/hisat2.mk

# Install required software.
install::
	@echo micromamba install hisat2 samtools



# Targets that are not files.
.PHONY: align install usage index test


