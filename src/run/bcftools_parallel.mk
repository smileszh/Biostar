#
# Runs SNP calling benchmarks
#

# The number of compute cores to use.
NCPU ?= 7

# Derived filenames
REF ?= refs/genome.fa
BAM ?= bam/alignment.bam
VCF ?= vcf/variants.vcf.gz

# How to label each split BAM file.
LABEL ?= parts

# The file that contains the parts to split by.
PARTS ?= tmp/${LABEL}.txt

# The directory that stores the split BAM files.
DIR ?= $(dir ${PARTS})

# Apply Makefile customizations.
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print some usage information
usage:
	@echo "#"
	@echo "# Parallel VCF calling"
	@echo "#"
	@echo "# LABEL=${LABEL}"
	@echo "# PARTS=${PARTS} "
	@echo "#"
	@echo "# NCPU=${NCPU} "
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make run"
	@echo "#"
	@echo "# Use the source, Luke!"
	@echo "#"


# Generate the parts file.
${PARTS}: ${BAM}

	# Create the parts directory
	mkdir -p $(dir ${PARTS})

	# Find the chromosomal names in the BAM.
	make -f src/run/splitchrom.mk show BAM=${BAM} > ${PARTS}

	# Keep only the first N lines of the parts file.
	cat ${PARTS} | head -25 > ${PARTS}.tmp && mv ${PARTS}.tmp ${PARTS}

	# Split the file by chromosome.
	cat ${PARTS} | parallel -v --eta -j ${NCPU} \
		make -f src/run/splitchrom.mk run \
			BAM=${BAM} \
			OUT=tmp/${LABEL}_{}.bam \
			CHROM={}

# Generate the split file
split: ${PARTS}
	@ls -lh ${PARTS}

# Call snps on the parts file
${VCF}: ${PARTS}
	mkdir -p tmp $(dir ${VCF})

	# Call variants on each chromosome separately.
	cat ${PARTS}  | parallel -v --eta -j ${NCPU} --lb \
		make -f src/run/bcftools.mk run BAM=tmp/${LABEL}_{}.bam  REF=${REF} VCF=tmp/${LABEL}_{}.vcf.gz

	# List the output VCF files.
	cat ${PARTS} | parallel -j 1 echo tmp/${LABEL}_{}.vcf.gz > tmp/${LABEL}_files.txt

	# Merge the files listed in the input
	bcftools concat -f tmp/${LABEL}_files.txt -o ${VCF}

	# Index the resulting VCF file.
	bcftools index -t -f ${VCF}

	# List the variant file
	@ls -l ${VCF}

# Trigger the variant calling.
run: ${VCF}
	@ls -lh ${VCF}


# Remove the files that have been created.
clean:
	rm -rf tmp
	rm -f ${PARTS} ${VCF}


# Not file targets.
.PHONY: design genome index fastq align vcf parts all clean get
