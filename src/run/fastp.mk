#
# Trims FASTQ files and runs FASTQC on them.
#

# The input read pairs.
P1 ?= reads/read1.fq
P2 ?=

# Extract the extension from P1
EXT1 ?= $(suffix ${P1})
EXT2 ?= $(suffix ${P2})

# Built the output name for read pairs
R1 = $(basename ${P1}).trimmed${EXT1}

# Single end mode if P2 not set
ifeq ($(P2),)
R2 =
else
R2 = $(basename ${P2}).trimmed${EXT2}
endif

# Fastp reporting files.
FASTP_JSON = $(basename ${R1}).fastp.json
FASTP_HTML = $(basename ${R1}).fastp.html

# Number of CPUS to use.
NCPU ?= 2

# The adapter sequence for trimming.
ADAPTER ?= AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage
usage:
	@echo "#"
	@echo "# fastp.mk: run fastp trimming on FASTQ files"
	@echo "#"
	@echo "# P1=${P1}"
	@echo "# P2=${P2}"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "#"
	@echo "# make run|clean "
	@echo "#"

# List the FLAGS used as trimming options
FASTP_FLAGS ?= --adapter_sequence ${ADAPTER} \
               --cut_right --cut_right_window_size 4 \
               --cut_right_mean_quality 30 \
               --length_required 50 \
                -j ${FASTP_JSON} -h ${FASTP_HTML}

# IF.
ifeq ($(P2),)
CMD = fastp -i ${P1} -o ${R1} -w ${NCPU} ${FASTP_FLAGS}
else
CMD = fastp -i ${P1} -I ${P2} -o ${R1} -O ${R2} -w ${NCPU} ${FASTP_FLAGS}
endif

# Read 1 must exist.
${P1}:
	@echo "# Input read 1 file not found: P1=${P1}"
	@exit -1

# If R2 is set, it must exist.
ifneq (${P2},)
${P2}:
	@echo "# Input read 2 file not found: P2=${P2}"
	@exit -1
endif

# Perform the trimming
${R1} ${R2}: ${P1} ${P2}
	 mkdir -p $(dir ${R1})
	 ${CMD}

# Run the trimming
run: ${R1} ${R2}
	 @ls -lh ${R1} ${R2}

# Removes the trimmed files.
run!:
	 rm -f ${R1} ${R2}

# Synonym for run!
clean: run!

# Test the module
test:
	make -f src/run/sra.mk run SRR=SRR1553425 N=1000
	make -f src/run/fastp.mk run! run P1=reads/SRR1553425_1.fastq run
	make -f src/run/fastp.mk run! run P1=reads/SRR1553425_1.fastq P2=reads/SRR1553425_2.fastq run

# Installation instuctions
install:
	@echo micromamba install fastp

# Targets that are not valid files.
.PHONY: usage run run! clean test install
