
#
# Runs the ivar package.
#

# Command line flags for ivar.
IVAR_FLAGS = -t 0.5 -m 1

# The accession number
ACC=AF086833

# The SRR number
SRR=SRR1553425

# The alignment file.
BAM ?= bam/${SRR}.bam

# The reference genome.
REF ?= refs/${ACC}.fa

# The reference as GFF.
GFF ?= refs/${ACC}.gff

# The prefix for the output files.
PREFIX = ivar/${ACC}.consensus

# The ivar consensus file.
CONS = ${PREFIX}.fa

# The ivar variants file.
VARS = ${PREFIX}.tsv

# Makefile customizations.
SHELL := bash
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

usage::
> @echo "#"
> @echo "# ivar.mk: runs the ivar suite"
> @echo "#"
> @echo "# REF=${REF}"
> @echo "# BAM=${BAM}"
> @echo "#"
> @echo "# make cons vars"
> @echo "#"

# Generates the ivar consensus sequences.
${CONS}: ${BAM}
> mkdir -p $(dir $@)
> samtools mpileup ${BAM} -aa -A -d 0 -Q 20 | ivar consensus -p ${PREFIX} ${IVAR_FLAGS}

# Generates the ivar variant table.
${VARS}: ${BAM} ${REF} ${GFF}
> mkdir -p $(dir $@)
> samtools mpileup ${BAM} -aa -A -B -d 0 -Q 20 | ivar variants -p ${PREFIX} ${IVAR_FLAGS} -r ${REF} -g ${GFF}

# Prints the consensus file location.
cons: ${CONS}
> @ls -lh ${CONS}

# Remove the consensus.
cons!:
> @rm -rf ${CONS}

# Prints the variant file location.
vars: ${VARS}
> @ls -lh ${VARS}

# Remove the variants.
vars!:
>@rm -rf ${VARS}

test:
> make -f src/run/genbank.mk gff! gff ACC=${ACC} GFF=${GFF}
> make -f src/run/bwa.mk BAM=${BAM} REF=${REF} index align
> make -f src/run/ivar.mk BAM=${BAM} REF=${REF} cons! cons vars! vars

# Install required software.
install::
>@echo micromamba install ivar samtools

# Targets that are not files.
.PHONY: cons cons! vars vars! usage
