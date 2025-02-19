#
# Tests the pipeline
#

# Default values for testing.
SRR ?= SRR1553425
ACC ?= AF086833
REF ?= refs/${ACC}.fa

# May plug in other modules for same testing.
MOD ?= src/run/bwa.mk
CALL ?= src/run/bcftools.mk

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print the usage.
usage:
	@echo "#"
	@echo "# Tests various pipelines"
	@echo "#"
	@echo "# make test_aligner MOD=src/run/bwa.mk"
	@echo "#"

# Test an alignment pipeline.
test_aligner:
	# Get the reference genome.
	make -f src/run/genbank.mk ACC=${ACC} REF=${REF} fasta

	# Get the FASTQ reads.
	make -f src/run/sra.mk SRR=${SRR} run

	# Single end run.
	make -f ${MOD} \
			REF=refs/${ACC}.fa \
			R1=reads/${SRR}_1.fastq \
			BAM=bam/${SRR}.bam \
			index run! run

	# Paired end run.
	make -f ${MOD} \
			REF=refs/${ACC}.fa \
			R1=reads/${SRR}_1.fastq \
			R2=reads/${SRR}_2.fastq \
			BAM=bam/${SRR}.bam \
			index run! run

# Test an alignment pipeline.
test_caller: test_aligner
	# Run variant caller
	make -f ${CALL} \
			REF=refs/${ACC}.fa \
			BAM=bam/${SRR}.bam \
			VCF=bam/${SRR}.vcf.gz \
			run! run

# Test an alignment pipeline.
test_wiggle:
	# Get the reference genome.
	make -f src/run/genbank.mk ACC=${ACC} REF=${REF} fasta

	# Get the FASTQ reads.
	make -f src/run/sra.mk SRR=${SRR} run

