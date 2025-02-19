#
# Streamline RefGenie initialization
#

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Default configuration
CONF = ~/refs/refgenie.yaml

usage::
	@echo "#"
	@echo "# refgenie.mk: init pull"
	@echo "#"

# Initialize refgenie
init:
	mkdir -p $(dir ${CONF})
	refgenie init -c ${CONF}
	echo "export REFGENIE=${CONF}" >> ~/.bashrc

# Pull human genome
hg38:
	refgenie pull hg38/fasta -c ${CONF} --no-overwrite
	refgenie pull hg38/bwa_index -c ${CONF} --no-overwrite

install:
	@echo pip install refgenie

