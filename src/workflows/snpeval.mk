#
# A makefile to evaluate SNP calls.
#

# Use a subset of the chromosome for testing.
CHR = chr1

# The chromosome FASTA file
REF_URL = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${CHR}.fa.gz

# The reference genome.
REF = refs/${CHR}.fa.gz

# The rtg index for the reference genome.
SDF = refs/sdf/${CHR}

# The high confidence SNP calls for the 2021 GIAB dataset.
GIAB_2021_VCF = vcf/GIAB-2021.${CHR}.vcf.gz

# The high confidence regions for the 2021 GIAB dataset.
GIAB_2021_BED = vcf/GIAB-2021.bed

# The high confidence regions for the 2017 GIAB dataset.
GIAB_2017_BED = vcf/GIAB-2017.bed

# The high confidence SNP calls for the 2017 GIAB dataset.
GIAB_2017_VCF = vcf/GIAB-2017.${CHR}.vcf.gz

# The Coriell index (the published paper)
CORIEL_VCF = vcf/DEEP-2017.${CHR}.vcf.gz

# The deep variant calls made in 2023.
DEEP_VCF = vcf/DEEP-2023.${CHR}.vcf.gz

# GATK calls performed with the Biostar Workflows.
GATK_VCF = vcf/GATK-2023.${CHR}.vcf.gz

# BCFTOOLS calls performed with the Biostar Workflows.
BCF_VCF = vcf/BCFT-2023.${CHR}.vcf.gz

# What target to compare it against
TARGET = ${GIAB_2021_VCF}

# Additional flags for rtg vcfeval.
FLAGS =

# Limit the evaluation to the high confidence regions.
FLAGS = -e ${GIAB_2021_BED}

#
# Apply Makefile customizations.
#
SHELL := bash
.ONESHELL:
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-print-directory
.SHELLFLAGS := -eu -o pipefail -c

# Print some usage information
usage:
	@echo "#"
	@echo "# Evaluate SNP calls"
	@echo "#"
	@echo "# TARGET = ${TARGET}"
	@echo "# FLAGS  = ${FLAGS}"
	@echo "#"
	@echo "# Usage: make data index eval "
	@echo "#"

# Get the reference genome.
${REF}:
	mkdir -p $(dir ${REF})
	curl ${URL} > ${REF}

# Trigger the reference download.
refs:${REF}
	@ls -lh ${REF}

# Index genome for rtg
${SDF}: ${REF}
	rtg format -o ${SDF} ${REF}

# Create the SDF index for rtg vcfeval.
index: ${SDF}
	@ls -lh ${SDF}

# Obtain the data
data:
	curl http://data.biostarhandbook.com/vcf/snpeval.tar.gz | tar zxvf -

# Run the evaluation.
eval: ${SDF}

	@printf "#\n# LIMIT=${FLAGS}\n#\n"

	@printf "#\n# ${TARGET} vs GIAB-2017\n#\n"
	@rtg vcfeval -b ${TARGET} -c ${GIAB_2017_VCF} -t ${SDF} -e ${GIAB_2017_BED} -o results/GIAB-2021-vs-GIAB-2017

	@printf "#\n# ${TARGET} vs DEEP-2017\n#\n"
	@rtg vcfeval -b ${TARGET} -c ${CORIEL_VCF} -t ${SDF} ${FLAGS} -o results/GIAB-2021-vs-DEEP-2017

	@printf "#\n# ${TARGET} vs DEEP-2023\n#\n"
	@rtg vcfeval -b ${TARGET} -c ${DEEP_VCF} -t ${SDF} ${FLAGS} -o results/GIAB-2021-vs-DEEP-2023

	@printf "#\n# ${TARGET} vs GATK-2023\n#\n"
	@rtg vcfeval -b ${TARGET} -c ${GATK_VCF} -t ${SDF} ${FLAGS} -o results/GIAB-2021-vs-GATK-2023

	@printf "#\n# ${TARGET} vs BCFTOOLS-2023\n#\n"
	@rtg vcfeval -b ${TARGET} -c ${BCF_VCF} -t ${SDF} ${FLAGS} -o results/GIAB-2021-vs-BCFT-2023

clean:
	rm -rf results

.PHONY: usage refs index eval
