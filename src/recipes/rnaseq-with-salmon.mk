#
# Biostar Workflows: https://www.biostarhandbook.com/
# 
# Chapter: RNA-Seq with salmon
#
# https://www.biostarhandbook.com/books/workflows/rnaseq/rnaseq-using-salmon/
#
# This code performs an RNA-Seq analysis using the salmon classifier.
# 
# Uses the bioinformatics toolbox: bio code
#


# The URL for the data.
DATA_URL = http://data.biostarhandbook.com/data/uhr-hbr.tar.gz

# The downloaded data file.
DATA_FILE = $(notdir ${DATA_URL})

# Genome reference.
REF = refs/chr22.transcripts.fa

# The design file
DESIGN = design.csv

# Final combined counts in CSV format.
COUNTS = res/counts-salmon.csv

# Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Flags passed to parallel.
FLAGS = --eta --lb --header : --colsep ,

# Number of CPUS to use
NCPU = 4

# Print usage
usage:
	@echo "#"
	@echo "# Biostar Workflows: RNA-Seq with salmon" 
	@echo "#"
	@echo "# https://www.biostarhandbook.com/books/workflows/rnaseq/rnaseq-using-salmon/"
	@echo "#"
	
	@echo "# REF=${REF}"
	@echo "# DESIGN=${DESIGN}"
	@echo "# COUNTS=${COUNTS}"
	@echo "#"

	@echo "# make data   (downloads the FASTQ and FASTA files)"
	@echo "# make index  (generates the salmon index)"
	@echo "# make design (generates the ${DESIGN} file)"
	@echo "# make align  (creates salmon count files)"
	@echo "# make count  (combines the counts into a single file)"
	@echo "#"

	@echo "# make all    (runs all steps)"
	@echo "#"

# Generate the design file.
${DESIGN}:
	@cat << EOF > ${DESIGN}
	sample,group
	HBR_1,HBR
	HBR_2,HBR
	HBR_3,HBR
	UHR_1,UHR
	UHR_2,UHR
	UHR_3,UHR
	EOF

# Show the design file.
design: ${DESIGN}
	@ls -lh ${DESIGN}

# Download the sequencing data and references.
${DATA_FILE}:
	# Download the sequencing data
	wget -nc  ${DATA_URL}
	tar xzvf uhr-hbr.tar.gz

# Trigger the data download
data: ${DATA_FILE}
	@ls -lh ${DATA_FILE}

# Generate the salmon index
# Run it in an enviroment that has salmon installed.
index: ${DESIGN}
	# Uses the salmon_env enviroment that has salmon installed.
	make -f src/run/salmon.mk index  REF=${REF} NCPU=${NCPU}

# Run salmon on the data
align:
	# Uses the salmon_env enviroment that has salmon installed. 
	cat ${DESIGN} | parallel ${FLAGS} \
		make -f src/run/salmon.mk \
		REF=${REF} \
		BAM=bam/{sample}.bam  \
		R1=reads/{sample}_R1.fq \
		NCPU=${NCPU} \
		SAMPLE={sample} \
		run

# Combine the counts into a single file
${COUNTS}:
	# Make the directory name for the counts.
	mkdir -p $(dir $@)

	# Combine the outputs into a single file.
	micromamba run -n stats Rscript src/r/combine_salmon.r -d ${DESIGN} -o ${COUNTS}

# Print the counts
counts: ${COUNTS}
	@ls -lh ${COUNTS}

all: data index align counts

.PHONY: usage data align counts run
