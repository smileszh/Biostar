#
# Generates a bigwig coverage from a bam file
#

# The reference genome.
REF ?= refs/genome.fa

# The name of the samtools faidx for the reference genome.
FAI ?= $(REF).fai

# Name of the BAM file.
BAM ?= bam/results.bam

# The bedgraph file.
BGR ?= ${BAM:.bam=.bedgraph}

# The wiggle file. Swap extension.
WIG ?= $(BAM:.bam=.bw)

# Number of CPUS
NCPU ?= 2

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage::
	@echo "#"
	@echo "# wiggle.mk: create BigWig from BAM"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# WIG=${WIG}"
	@echo "#"
	@echo "# make run"
	@echo "#"

# Checks that reference exists.
${REF}:
	echo "# Reference not found: ${REF}";
	exit 1

# Checks that BAM file exists.
${BAM}:
	echo "# BAM file not found: ${BAM}";
	exit 1

# Create the fasta index file.
${FAI}: ${REF}
	samtools faidx ${REF}

# Creates the bigwig file.
${WIG}: ${FAI} ${BAM}
	mkdir -p $(dir ${WIG})
	LC_ALL=C; bedtools genomecov -ibam  ${BAM} -split -bg  | sort -k1,1 -k2,2n > ${BGR}
	bedGraphToBigWig ${BGR} ${FAI} ${WIG}
	rm -f ${BGR}

# Generate the wiggle file.
run: ${WIG}
	@ls -lh ${WIG}

# Remove the wiggle file.
run!:
	rm -rf ${WIG}

test:
	# Get the reference genome.
	make -f src/run/bwa.mk test

# Install required software.
install::
	@echo micromamba install bedtools ucsc-bedgraphtobigwig

# Targets that are not files.
.PHONY: usage run run! test

