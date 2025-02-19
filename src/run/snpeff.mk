#
# Annotating SNP effects with SNPeff.
#

# Where snpEff will store the database
IDX = idx/snpEff

# Genbank accession number.
REF = refs/genome.fa

# The GenBank annotation file.
GFF ?= refs/annotations.gff

# Existing Variant calls.
VCF ?= vcf/variants.vcf.gz

# The prefix for the annotated files
RES = results

# Annotated variants.
EFF ?= ${RES}/snpeff.vcf.gz

# The SNPeff html report.
HTML ?= ${RES}/snpeff.html

# The SNPeff CXV report.
CSV ?= ${RES}/snpeff.csv

# The label for the SNPeff database.
LABEL = genome

# The SNPeff database.
SNPEFF_DB ?= ${IDX}/${LABEL}/snpEffectPredictor.bin

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# General usage information.
usage::
	@echo "#"
	@echo "# snpeff.mk: annotate variants calls with SNPeff."
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# GFF=${GFF}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# LABEL=${LABEL} # Used by snpEff"
	@echo "#"
	@echo "# EFF=${EFF}"
	@echo "# HTML=${HTML}"
	@echo "# CSV=${CSV}"
	@echo "#"
	@echo "# make build|run|clean"
	@echo "#"

# Check the reference file.
${REF}:
	@echo "#"
	@echo "# REF not found: ${REF}"
	@echo "#"
	exit -1

# Check the GFF file.
${GFF}:
	@echo "#"
	@echo "# GFF not found: ${GFF}"
	@echo "#"
	exit -1

# Building a custom snpEff database snpeff needs the files in specific folders.
# 1. copy the genbank to its location.
# 2. append entry to current genome to the config.
# 3. Build the snpEff database.

${SNPEFF_DB}: ${GFF} ${REF}	
	mkdir -p ${IDX}/${LABEL}

	# Copy the files to the snpEff folder.
	cp -f ${REF} ${IDX}/${LABEL}/sequences.fa
	cp -f ${GFF} ${IDX}/${LABEL}/genes.gff

	# Make the configuration file.
	echo "${LABEL}.genome : ${LABEL}" >	snpeff.config

	# Build the database.
	snpEff build -dataDir ${IDX} -v ${LABEL}

# Create the annotated VCF file.
${EFF} ${HTML} ${CSV}: ${SNPEFF_DB}
	mkdir -p $(dir ${EFF})
	snpEff ann -csvStats ${CSV} -s ${HTML} -dataDir ${IDX} -v ${LABEL} ${VCF} | bcftools view -O z -o ${EFF}
	bcftools index ${EFF}

# Creates the SNPeff database
index: ${SNPEFF_DB}
	@ls -lh ${SNPEFF_DB}

# Another name for index
build: index

# Removes the SNPeff database
index!:
	rm -f ${SNPEFF_DB}

# Runs the SNPeff tool.
run: ${EFF} ${HTML} ${CSV}
	@ls -lh ${EFF}
	@ls -lh ${HTML}
	@ls -lh ${CSV}

# Removes the SNPeff output
clean:
	rm -f ${ANN_VCF} ${ANN_HTML}

# Another name cleanup
run!: clean

# Shows the usage information.
install::
	@echo micromamba install snpeff

.PHONY: run
