#
# Biostar Workflows: https://www.biostarhandbook.com/
#

# Chromosome Y of the human genome.
ENSEMBL_URL = http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz

# The name of the local file for the data.
ENSEMBL_FILE = refs/chrY.fa.gz

# GenBank accession.
GENBANK_ACC=AF086833

# The name of the fasta file
GENBANK_FASTA = refs/${GENBANK_ACC}.fa

# For NCBI assembly ids
NCBI_ASSEMBLY_ID = GCA_000005845

# We found out the file name structure only after we downloaded it.
NCBI_ASSEMBLY_FASTA = ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna

usage:
	@echo "#"
	@echo "# make ensembl genbank assembly"
	@echo "#"

# Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Downlaoding Ensembl files.
# Ensembl provides gzip files instead of bg-zipped FASTA. We need to recompress with bgzip.
${ENSEMBL_FILE}:

	curl ${ENSEMBL_URL} > ${ENSEMBL_FILE}

	make -f src/run/bgzip.mk run FILE=${ENSEMBL_FILE}

ensembl: ${ENSEMBL_FILE}
	@ls -lh $<

# Accessing GENBANk files.
${GENBANK_FASTA}:
	mkdir -p $(dir $@)
	bio fetch ${GENBANK_ACC} -format fasta > $@

# Download the GENBANK fasta.
genbank: ${GENBANK_FASTA}
	@ls -lh $<

# The assembly fasta file generation.
${NCBI_ASSEMBLY_FASTA}:

	datasets download genome accession ${NCBI_ASSEMBLY_ID}  --filename ncbi_dataset.zip

	# Unzip while forcing overwrite
	unzip -o ncbi_dataset.zip

# Download the NCBI genome assembly
assembly: ${NCBI_ASSEMBLY_FASTA}
	@ls -lh $<

# Run all the tests.
run: ensembl genbank assembly

run!:
	rm -rf ${ENSEMBL_FILE} ${GENBANK_FASTA} ncbi_dataset.zip ncbi_dataset

.PHONY: usage ensembl genbank assembly run
