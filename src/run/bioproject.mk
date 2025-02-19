#
# Downloads NCBI run information based on a bioproject number
#

# Project number
ID ?= PRJNA257197

# Project runinfo file.
OUT ?= ${ID}.csv

# Makefile customizations.
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# General usage information.
usage::
	@echo "#"
	@echo "# bioproject.mk: downloads runinfo for an SRA bioproject"
	@echo "#"
	@echo "# ID=PRJNA257197"
	@echo "# OUT=PRJNA257197.csv"
	@echo "#"
	@echo "# make run|clean "
	@echo "#"

# Project run information.
${OUT}:
	mkdir -p $(dir $@)
	bio search ${ID} --header --csv > $@

# Target to download all the data.
run::  ${OUT}
	@ls -lh ${OUT}

# Remove bioprject
run!::
	rm -rf ${OUT}

# For backward compatibility.
get: run
clean: run!

# Installation instructions
install::
	@echo "pip install bio --upgrade"
