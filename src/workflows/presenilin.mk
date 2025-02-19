#
# Presenilin RNA-Seq in the Biostar Workflows
#
# http://www.biostarhandbook.com
#

# Single end mode.
MODE = SE

# How many reads to download N=ALL for all.
N = 4000000

# The number of CPUs to use.
NCPU = 4

# The URL for the CDNA
URL = http://ftp.ensembl.org/pub/release-107/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz

# The name of the CDNA file locally.
REF = ~/refs/drerio/Danio_rerio.cdna.fa

# The name of the ensembl database to use for transcript to gene mapping.
ENSEMBL_DB ?= drerio_gene_ensembl

# The transcript to gene mapping file.
TX2GENE ?= ~/refs/tx2gene/drerio_tx2gene.csv

# Change the design matrix name.
DESIGN ?= presen.csv

# The name of the counts file.
COUNTS ?= res/psen-counts.csv

# Include the airway RNA-Seq pipeline.
include src/workflows/airway.mk

# Apply Makefile customizations.
.DELETE_ON_ERROR:
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

design:
	@cat << EOF > ${DESIGN}
	run,group,sample,celltype
	SRR8530750,WT,WT_1,1
	SRR8530751,WT,WT_2,2
	SRR8530752,WT,WT_3,3
	SRR8530753,WT,WT_4,4
	SRR8530754,Q96,Q96_1,5
	SRR8530755,Q96,Q96_2,6
	SRR8530756,Q96,Q96_3,7
	SRR8530757,Q96,Q96_4,8
	EOF
	@echo "# Created presenilin design: ${DESIGN}"

# Run alignments for all samples in single end mode
align: ${DESIGN}
	# Run the alignment on each sample.
	cat ${DESIGN} | parallel ${FLAGS} \
	                micromamba run -n salmon_env make -f src/run/salmon.mk run \
	                REF=${REF} \
	                SAMPLE={sample} \
	                R1=reads/{run}_1.fastq



