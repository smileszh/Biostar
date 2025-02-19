#
# Generate alignments with star
#
# When using lots of threads you may need to increase ulimit with:
#
# ulimit -n 10000
#

# Number of CPUS
NCPU ?= 2

# Additional flags to pass to HISAT.
STAR_FLAGS ?= --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within

# Read groups.
RG ?= --rg-id ${ID} --rg SM:${SM} --rg LB:${LB} --rg PL:${PL}

# Sam filter flags to filter the BAM file before sorting.
SAM_FLAGS ?=

# First in pair.
R1 ?= reads/read1.fq

# Second in pair.
R2 ?=

# Set MODE to SE if R2 is empty otherwise PE
ifeq ($(R2),)
	MODE ?= SE
else
	MODE ?= PE
endif

# The reference genome.
REF ?= refs/genome.fa

# The directory that holds the index.
IDX_DIR = $(dir ${REF})/idx

# The name of the index
IDX ?= ${IDX_DIR}/$(notdir ${REF}).star

# A file in the index directory.
IDX_FILE ?= ${IDX}/SAindex

# The alignment file.
BAM ?= bam/star.bam

# Star needs a prefix to run.
PREFIX = $(basename ${BAM})

# Unsorted BAM file
BAM_TMP ?=${PREFIX}Aligned.sortedByCoord.out.bam

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
usage::
	@echo "#"
	@echo "# star.mk: align reads using STAR"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# IDX=${IDX}"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make index run"
	@echo "#"

# Build the index for the reference genome.
${IDX_FILE}:
	@if [ ! -f ${REF} ]; then
		echo "# file not found: REF=${REF}";
		exit -1
	fi
	mkdir -p $(dir $@)

	GL=$$(awk '$$0 ~ ">" {next} {len+=length($$0)} END {print len}' $(REF))

	# Calculate genomeSAindexNbases
	SA=$$(awk -v len=$$GL 'BEGIN{print int(log(len)/log(2)/2 - 1)}')

	# Cap the value at 14
	@if [ $$SA -gt 14 ]; then
		SA=14
	fi

	if [ "$(suffix ${REF})" = ".gz" ]; then
		# Unzip if it does not exist unzipped
		[ ! -f "$${REF%.*}" ] && gunzip -k "$$REF"
		STAR --runThreadN ${NCPU} --runMode genomeGenerate --genomeFastaFiles $(basename ${REF}) --genomeDir ${IDX} --genomeSAindexNbases $$SA
	else
		STAR --runThreadN ${NCPU} --runMode genomeGenerate --genomeFastaFiles ${REF} --genomeDir ${IDX} --genomeSAindexNbases $$SA
	fi



# Create the index.
index: ${IDX_FILE}
	@echo "# STAR index: ${IDX}"

# Remove the index.
index!:
	rm -rf ${IDX_FILE}

# Select the command by alignment mode.
ifeq ($(MODE), PE)
CMD = STAR --genomeDir ${IDX} --runThreadN ${NCPU} --readFilesIn ${R1} ${R2} --outFileNamePrefix ${PREFIX} ${STAR_FLAGS}
else
CMD = STAR --genomeDir ${IDX} --runThreadN ${NCPU} --readFilesIn ${R1} --outFileNamePrefix ${PREFIX} ${STAR_FLAGS}
endif

# Perform the alignment.
${BAM}: ${R1} ${R2}
	@if [ ! -f ${IDX_FILE} ]; then
		echo "# STAR index not found: IDX=${IDX}";
		exit -1
	fi
	@mkdir -p $(dir $@)
	${CMD}
	mv -f ${BAM_TMP} ${BAM}

# Create the BAM index file.
${BAM}.bai: ${BAM}
	samtools index ${BAM}

# Generate the alignment.
align: ${BAM}.bai
	@ls -lh ${BAM}

# Create a wiggle coverage file.
wiggle: align
	@make -f src/run/wiggle.mk BAM=${BAM} REF=${REF} run

# Display the BAM file path.
run: wiggle

# Remove the BAM file.
run!:
	rm -rf ${BAM} ${BAM}.bai ${BAM_TMP}

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_aligner MOD=src/run/star.mk

# Install required software.
install::
	@echo micromamba install star

# Reads must exist.
${R1}:
	@echo "# file not found: R1=${R1}"
	@exit -1

# Targets that are not files.
.PHONY: align install usage index test


