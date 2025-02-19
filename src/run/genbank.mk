#
# Downloads NCBI data vie Entrez API
#

# Accession number at NCBI
ACC ?= AF086833

# The accession as a fasta file.
REF ?= refs/${ACC}.fa

# The accession as a genbank file.
GBK ?= refs/${ACC}.gb

# The accession as a gff file.
GFF ?= refs/${ACC}.gff

# Check the make version.
ifeq ($(origin .RECIPEPREFIX), undefined)
  $(error "### Error! Please use GNU Make 4.0 or later ###")
endif

# Makefile customizations.
.DELETE_ON_ERROR:
.ONESHELL:
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage::
	@echo "#"
	@echo "# genbank.mk: download sequences from GenBank"
	@echo "#"
	@echo "# ACC=${ACC}"
	@echo "# REF=${REF}"
	@echo "# GBK=${GBK}"
	@echo "# GFF=${GFF}"
	@echo "#"
	@echo "# make fasta|gff|genbank|all"
	@echo "#"

# Obtain a sequence from NCBI.
${REF}:
	mkdir -p $(dir $@)
	bio fetch ${ACC} --format fasta > $@

${GBK}:
	mkdir -p $(dir $@)
	bio fetch ${ACC} > $@

${GFF}:
	mkdir -p $(dir $@)
	bio fetch ${ACC} --format gff > $@

# Download FASTA file.
fasta::  ${REF}
	@ls -lh ${REF}

# Download GenBank file.
genbank::  ${GBK}
	@ls -lh ${GBK}

# Download GFF file.
gff::  ${GFF}
	@ls -lh ${GFF}

# Remove FASTA file.
fasta!::
	rm -rf ${REF}

# Remove GenBank file.
genbank!::
	rm -rf ${GBK}

# Remove GFF file.
gff!::
	rm -rf ${GFF}

# The run target is FASTA
run:: fasta

# Undo run target
run!:: fasta!

# Create all three outputs
all:: fasta genbank gff

test:
	make -f src/run/genbank.mk fasta! fasta ACC=${ACC} REF=${REF}
	make -f src/run/genbank.mk gff! gff ACC=${ACC} GFF=${GFF}
	make -f src/run/genbank.mk genbank! genbank ACC=${ACC} GBK=${GBK}

# Installation instructions
install::
	@echo pip install bio --upgrade

# Targets that are not files.
.PHONY: usage install run run! test fasta fasta! genbank genbank! gff gff!
