# Number of CPUS
NCPU ?= 2

# The reference genome.
REF ?= refs/genome.fa

# The alignment file.
BAM ?= bam/results.bam

# The variant file.
VCF ?= vcf/$(notdir $(basename ${BAM})).deepvariant.vcf.gz

# The temporary intermediate results.
TMP ?= tmp/$(notdir ${VCF})

# The model type in deep variant.
MODEL ?= WGS

# The directory to mount
MNT ?= /deep

# Directory to mount in the singularity container.
MNT_PARAM ?= $(shell pwd):${MNT}

# Deepvariant singularity image
SIF ?= deepvariant_1.5.0.sif

# Deepvariant singularity image URL
SIF_URL ?= docker://google/deepvariant:1.5.0

# The deepvariant command line run
CMD ?= singularity run -B ${MNT_PARAM} ${SIF}

# Additional bcf flags for calling process.
CALL_FLAGS ?=

# Example:
# CALL_FLAGS = --regions chr1:1-1000000

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage:
	@echo "#"
	@echo "# deepvariant.mk: call variants using Google Deepvariant"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make run"
	@echo "#"

${SIF}:
	singularity pull ${SIF_URL}

${REF}.fai: ${REF}
	samtools faidx ${REF}

# Generate the variant calls.
${VCF}: ${SIF} ${BAM} ${REF} ${REF}.fai
	mkdir -p $(dir $@)
	${CMD} \
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=${MODEL} \
	--ref ${MNT}/${REF} \
	--reads ${MNT}/${BAM} \
	--output_vcf ${MNT}/${VCF}  \
	--num_shards ${NCPU} \
	--intermediate_results_dir ${MNT}/${TMP} \
	${CALL_FLAGS}

# Create the VCF index
${VCF}.tbi: ${VCF}
	bcftools index -t -f $<

TOOLS = singularity samtools

check:
	@$(foreach tool,$(TOOLS),\
		if ! command -v $(tool) > /dev/null; then \
			echo "# Error: $(tool) is not installed"; \
			exit 1; \
		fi;)

# Generating a VCF file.
run: check ${VCF}.tbi
	@ls -lh ${VCF}

run!:
	rm -rf ${VCF} ${VCF}.tbi

# Test the entire pipeline.
test: check
	make -f src/run/tester.mk test_caller MOD=src/run/bwa.mk CALL=src/run/deepvariant.mk

# Installation instructions.
install:
	@echo micromamba install singularity
	@echo singularity pull ${SIF_URL}

.PHONY: usage vcf test install
