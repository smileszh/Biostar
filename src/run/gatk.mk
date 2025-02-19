#
# Generates SNP calls with gatk4
#

# A root to derive output default names from.
SRR=SRR1553425

# Number of CPUS
NCPU ?= 2

# Accession number
ACC = AF086833

# Reducing the location
LOC = AF086833.2

# The reference genome.
REF ?= refs/${ACC}.fa

# Reference dictionary.
DICT = $(basename ${REF}).dict

# The known sites in VCF format.
SITES ?= vcf/${SRR}.knownsites.vcf.gz

# The alignment file.
BAM ?= bam/${SRR}.bam

# Mark duplicates BAM file.
DUP = $(basename ${BAM}).markdup.bam

# Recalibration table.
TAB = $(basename ${BAM}).recal.txt

# Recalibrated BAM file.
RCB = $(basename ${BAM}).recal.bam

# The directory for the variant files.
VCF_DIR = vcf

# The variant file.
VCF ?= ${VCF_DIR}/$(notdir $(basename ${BAM})).gatk.vcf.gz

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

JAVA_OPTS = '-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=4'

# The first target is always the help.
usage::
	@echo "#"
	@echo "# gatk.mk: call variants using GATK4"
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# BAM=${BAM}"
	@echo "# LOC=${LOC}"
	@echo "# SITES=${SITES}"
	@echo "#"
	@echo "# DUP=${DUP}"
	@echo "# TAB=${TAB}"
	@echo "# RCB=${RCB}"
	@echo "#"
	@echo "# VCF=${VCF} "
	@echo "# "
	@echo "# make mark calibrate apply vcf"
	@echo "#"

# Generate the FASTA sequence dictionary.
${DICT}: ${REF}
	gatk --java-options ${JAVA_OPTS} CreateSequenceDictionary --REFERENCE ${REF}

# Shortcut to FASTA sequence dictionary.
dict: ${DICT}

# Mark duplicates.
${DUP}: ${BAM} ${DICT}
	gatk MarkDuplicates -I ${BAM} -O ${DUP} -M ${DUP}.metrics.txt

# Shortcut to mark duplicates.
mark: ${DUP}
	ls -l $<

# Generate the recalibration table.
${TAB}: ${DUP} ${SITES}
	gatk BaseRecalibrator -R ${REF} -I ${BAM} --known-sites ${SITES} -O ${TAB}

# Compute the calibration table.
calibrate: ${TAB}
	ls -l $<

# Remove the recalibration table.
calibrate!:
	rm -rf ${TAB}

# Apply the recalibration table.
${RCB}: ${TAB}
	gatk ApplyBQSR -R ${REF} -I ${BAM} --bqsr-recal-file ${TAB} -O ${RCB}

# Apply a recalibration onto the BAM file.
apply: ${RCB}
	ls -l $<

# Remove recalibrated BAM file.
apply!:
	rm -rf ${RCB}

# Call variants.
${VCF}: ${BAM} ${DICT}
	mkdir -p $(dir $@)
	gatk HaplotypeCaller --java-options ${JAVA_OPTS} -I ${BAM} -R ${REF} -L ${LOC} --output ${VCF}

# Shortcut to call variants.
vcf: ${VCF}
	@ls -lh ${VCF}

# Delete variant call.
vcf!:
	rm -rf ${VCF}

# Test the entire pipeline.
test:
	# Get the reference genome.
	make -f src/run/genbank.mk ACC=${ACC} fasta
	# Get the FASTQ reads.

	make -f src/run/sra.mk SRR=${SRR} get
	# Align the FASTQ reads.

	make -f src/run/bwa.mk BAM=${BAM} REF=${REF} index align
	# Create known sites with bcftools

	make -f src/run/bcftools.mk BAM=${BAM} REF=${REF} VCF=${SITES} vcf
	# Call the variants.

	make -f src/run/gatk.mk TARGET=${TARGET} VCF=${VCF} BAM=${BAM} REF=${REF} mark calibrate apply vcf! vcf

install::
	@echo micromamba install gatk4 test

# Targets that are not files.
.PHONY: vcf vcf! usage test mark calibrate apply
