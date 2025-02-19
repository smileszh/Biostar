#
# Biostar Workflows: https://www.biostarhandbook.com/
#

# The SRR number.
SRR=SRR10971381

# How many reads to download at once (must be first parameter of the script)
N = 1000000

# Original fastq reads.
P1 = reads/${SRR}_1.fastq
P2 = reads/${SRR}_2.fastq

# Trimmed fastq reads.
R1 = trim/${SRR}_1.fastq
R2 = trim/${SRR}_2.fastq

# How much memory in GigaBytes to use.
MEM = 2

# Memory in bytes.
MEM_BYTES = $(shell echo $$((${MEM} * 1024 * 1024 * 1024)))

# The number of CPUS
NCPU = 4

# The name of the assembly output directory.
ASM = asm

# The name of the contigs after the assembly.
CONTIGS = ${ASM}/final.contigs.fa

# The name of the blast database
BLAST_DB = ~/db/ref_viruses_rep_genomes

# The name of a blast database file.
BLAST_DB_FILE = ${BLAST_DB}.ndb

# The assembly blast results file.
BLAST_RESULT = blast-${N}.txt

# Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

usage:
	@echo "#"
	@echo "# make fastq trim assemble blast"
	@echo "#"

# Download the FASTQ files.
fastq:
	make -f src/run/sra.mk run SRR=${SRR} N=${N} MODE=PE

trim:
	make -f src/run/fastp.mk run P1=${P1} P2=${P2} R1=${R1} R2=${R2} 

# Run the assembler.
${CONTIGS}: ${R1} ${R2}
	rm -rf asm
	@echo "#"
	@echo "# Memory: ${MEM_GB} GBytes "
	@echo "# Threads: ${NCPU}"
	@echo "#"
	megahit -t ${NCPU} -m ${MEM_BYTES} -1 ${R1} -2 ${R2} -o ${ASM}

# Display the assembly file.
assemble: ${CONTIGS}
	@ls -lh $<

# Display the assembly file.
assemble!: ${CONTIGS}
	rm -f ${CONTIGS}

# Create the blast database.
${BLAST_DB_FILE}:
	mkdir -p $(dir $@)
	cd $(dir $@) && update_blastdb.pl --source aws --decompress ref_viruses_rep_genomes

# Align contigs against all viral genomes. Sort by longest alignment length.
${BLAST_RESULT}: ${CONTIGS} ${BLAST_DB_FILE}
	blastn -db ${BLAST_DB} \
		   -query ${CONTIGS} \
		   -outfmt "6 qacc qlen length qcovs sacc staxid ssciname" | \
		    sort -k 4 -rn > ${BLAST_RESULT}

# Show the blast result.
blast: ${BLAST_RESULT}
	@ls -lh $<

# Remove the blast result.
blast!:
	rm -rf ${BLAST_RESULT}

install:
	@echo mamba install megahit=1.2.8 blast

.PHONY: usage fastq trim assemble blast install
