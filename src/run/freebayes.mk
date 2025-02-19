	#
# Generates SNP calls with freebaye
#


# Number of CPUS
NCPU ?= 2

# The alignment file.
BAM ?= bam/results.bam

# The reference genome.
REF ?= refs/genome.fa

# Name of the VCF with all variants.
VCF ?= vcf/$(notdir $(basename ${BAM})).freebayes.vcf.gz

# Additional flags passed to freebayes
BAYES_FLAGS = --pvar 0.5

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage::
	@echo "#"
	@echo "# freebayes.mk: calls variants using freebayes"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make run"
	@echo "#"

# Read 1 must exist.
${REF}:
	@echo "# Reference file not found: REF=${REF}"
	@exit -1

# Read 1 must exist.
${BAM}:
	@echo "# BAM alignment file not found: BAM=${BAM}"
	@exit -1

# Call SNPs with freebayes.
${VCF}: ${BAM} ${REF}
	mkdir -p $(dir $@)
	freebayes ${BAYES_FLAGS} -f ${REF} ${BAM} | bcftools norm -f ${REF} -d all -O z  > ${VCF}

# The VCF index file.
${VCF}.tbi: ${VCF}
	bcftools index -t -f $<

# The main action.
run:: ${VCF}.tbi
	@ls -lh ${VCF}

# Undo the main action.
run!::
	rm -rf ${VCF}

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_caller MOD=src/run/bwa.mk CALL=src/run/freebayes.mk

# Print installation instructions.
install::
	@echo micromamba install bcftools freebayes

# Targets that are not files.
.PHONY: run run! install usage index test

