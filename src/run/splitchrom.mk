#
# Splits a BAM file by a chromosome
#

# The chromosome to split on.
CHROM ?= chr1

# The alignment file.
BAM ?= bam/aligment.bam

# The output file. name
OUT ?= $(basename ${BAM}).${CHROM}.bam

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# The first target is always the help.
usage::
	@echo "#"
	@echo "# splitchrom.mk: extracts a chromosome from a BAM file"
	@echo "#"
	@echo "# CHROM=${CHROM}"
	@echo "# BAM=${BAM}"
	@echo "# OUT=${OUT}"
	@echo "#"
	@echo "# make run"
	@echo "#"


# Show all the chromosomes in a BAM file
show:
	@samtools view -H  ${BAM} | perl -lne '/SN:(\S+)/ and print "$$1"'

# Generate the names of the chromosomes from the BAM file.
${OUT}: ${BAM}
	mkdir -p $(dir ${OUT})
	samtools view -b ${BAM} ${CHROM} > ${OUT}

# The index file for the BAM file.
${OUT}.bai: ${OUT}
	samtools index ${OUT}

# Split the BAM file into multiple BAM files.
run: ${OUT}.bai
	@ls -lh $<

run!:
	rm -f ${OUT} ${OUT}.bai

install::
	@echo micromamba install samtools

# Targets that are not files.
.PHONY: usage run
