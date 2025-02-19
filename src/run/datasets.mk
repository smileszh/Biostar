#
# Downloads Genome data from NCBI.
#

# Accession number
ACC ?= GCF_000006945.2

#
# Datasets will add an assembly name for the genome.
#
# You can find the assembly name with: 
#
#  make info ACC=GCF_000840245.1 | grep assembly_name
#
NAME ?= ASM694v2

# The genome file name as downloaded.
GENOME_FA ?= ncbi_dataset/data/${ACC}*/${ACC}*_genomic.fna

# The GFF file name as downloaded.
GENOME_GFF ?= ncbi_dataset/data/${ACC}*/genomic.gff

# The GTF file name as downloaded.
GENOME_GTF ?= ncbi_dataset/data/${ACC}*/genomic.gtf

# The final names for the reference and GFF files.
REF = refs/${ACC}.fa
GFF = refs/${ACC}.gff
GTF = refs/${ACC}.gtf

# Which files to download.
INCLUDE ?= genome,gff3,gtf

# Makefile customizations.
SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Print usage information.
usage:
	@echo "#"
	@echo "# datasets.mk: download NCBI genomes"
	@echo "#"
	@echo "# ACC=${ACC}"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# GFF=${GFF}"
	@echo "# GTF=${GTF}"
	@echo "#"
	@echo "# make info|name|run|clean "
	@echo "#"

# Print the summary information on the genome.
info:
	@datasets summary genome accession ${ACC} | jq

# Finds the assembly name for an accession
name:
	@datasets summary genome accession ${ACC} | jq '.reports[] | {accession: .paired_accession, assembly_name: .assembly_info.assembly_name}'

# This is how we make the genome.
${GENOME_FA}:
	
	# Download the genome.
	datasets download genome accession ${ACC} --include ${INCLUDE} --filename ncbi_dataset.zip
	
	# Never overwrite files when unzipping!
	unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	
	# Clean up the zip file.
	rm -f ncbi_dataset.zip

# We don't know for sure which files exist hence we check for each one.
# It is not as robust is it could be but we blame NCBI for the wonky process.
${REF}: ${GENOME_FA} 

	# Make the reference directories
	mkdir -p $(dir ${REF}) $(dir ${GFF})
	cp ${GENOME_FA} ${REF}

	# If the GFF file exists copy it to destination.
	if [ -f ${GENOME_GFF} ]; then
		cp ${GENOME_GFF} ${GFF}
	fi

	# If the GTF file exists copy it to destination.
	if [ -f ${GENOME_GTF} ]; then
		cp ${GENOME_GTF} ${GTF}
	fi

# Downloads the genome file.
run: ${REF}
	ls -lh ${REF} ${GFF} ${GTF}

# Triggers the download if the files does not exist.
stats: ${REF}
	seqkit stats ${REF}

# Cleanup the downloaded files.
clean:
	rm -rf ncbi_dataset/data/${ACC}*
	rm -f md5sum.txt ncbi_dataset.zip ${REF} ${GFF} ${GTF}

# A shortcut to clean.
run!: clean

# Test the tool
test: clean run

# Mark the targets that do not create files.
.PHONY: usage info name run stats clean run! test