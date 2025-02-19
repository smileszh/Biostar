#
# Generates SNP calls with bcftools.
#

# A root to derive output default names from.
SRR = SRR1553425

# Number of CPUS
NCPU ?= 2

# Genbank accession number.
ACC ?= AF086833

# The reference genome.
REF ?= refs/${ACC}.fa

# The alignment file.
BAM ?= bam/${SRR}.bam

# The variant file.
VCF ?= vcf/$(notdir $(basename ${BAM})).vcf.gz

# Additional bcf flags for pileup annotation.
PILE_FLAGS =  -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# Additional bcf flags for calling process.
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage::
	@echo "#"
	@echo "# bcftools.mk: calls variants using bcftools"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make run"
	@echo "#"

${VCF}: ${BAM} ${REF}
	mkdir -p $(dir $@)
	bcftools mpileup ${PILE_FLAGS} -O u -f ${REF} ${BAM} | \
		bcftools call ${CALL_FLAGS} -mv -O u | \
		bcftools norm -f ${REF} -d all -O u | \
		bcftools sort -O z > ${VCF}

${VCF}.tbi: ${VCF}
	bcftools index -t -f $<

run:: ${VCF}.tbi
	@ls -lh ${VCF}

run!::
	rm -rf ${VCF} ${VCF}.tbi

install::
	@echo micromamba install bcftools

# Test the entire pipeline.
test:
	make -f src/run/tester.mk test_caller MOD=src/run/bwa.mk CALL=src/run/bcftools.mk

# Targets that are not files.
.PHONY: run run! install usage index test
