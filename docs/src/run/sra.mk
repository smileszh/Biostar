#
# Downloads sequencing reads from SRA.
#

# The directory that stores FASTQ reads that we operate on.
DIR ?= reads

# SRR number (sequencing run from the Ebola outbreak data of 2014)
SRR ?= SRR1553425

# The name of the unpacked reads
P1 ?= ${DIR}/${SRR}_1.fastq
P2 ?= ${DIR}/${SRR}_2.fastq

# The name of the reads
# You can rename reads to be more descriptive.
R1 ?= ${P1}
R2 ?= ${P2}

# How many reads to download (N=ALL downloads everything).
N ?= 10000

# Makefile customizations.
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information.
usage::
	@echo "#"
	@echo "# sra.mk: downloads FASTQ reads from SRA"
	@echo "#"
	@echo "# SRR=${SRR}"
	@echo "# N=${N} (use N=ALL to download all reads)"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "#"
	@echo "# make run|test|aria|clean"
	@echo "#"
	

# Determine if we download all reads or just a subset.
ifeq ($(N), ALL)
FLAGS ?= -F --split-files
else
FLAGS ?= -F --split-files -X ${N}
endif

# Obtain the reads from SRA.
${R1}:
	# Create the directory.
	mkdir -p ${DIR}

	# Download the reads.
	fastq-dump ${FLAGS} -O ${DIR} ${SRR}

	# Rename the first in pair if final name is different.
	if [ "${P1}" != "${R1}" ]; then mv -f ${P1} ${R1}; fi

	# Rename the second pair if exists and is different.
	if [ -f ${P2} ] && [ "${P2}" != "${R2}" ]; then mv -f ${P2} ${R2}; fi

# List the data. We don't know if the reads are paired or not.
run: ${R1}
	@if [ -f ${R2} ]; then 
		@ls -lh ${R1} ${R2}
	else
		@ls -lh ${R1}
	fi

# Removes the SRA files.
clean:
	rm -f ${P1} ${P2} ${R1} ${R2}

# A synonym for clean.
run!: clean

# Download via aria2c command line too.
# The process may be more reliable than fastq-dump.
# Will not rename files to R1 and R2!
aria:
	# Extact URLs from the search output
	# Then use aria2c to download the reads for the SRR number.
	bio search ${SRR} | jq -r '.[].fastq_url[]' | \
		parallel -j 1 --lb make -f src/run/aria.mk URL={} DIR=${DIR} run

	# Shows the resulting files.
	ls -l ${DIR}

# Run the test suite.
test: clean run

install::
	@echo micromamba install sra-tools

# Targets that are not files.
.PHONY: usage run run! test install
