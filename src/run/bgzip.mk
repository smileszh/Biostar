#
# Recompress a gzip file with bgzip.
#

# The file to compress
FILE ?= source.gz

# Temporary file.
TMP ?= $(basename ${FILE}).tmpfile.000.gz

# Number of CPUs to use
NCPU ?= 4

# Makefile customizations.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

usage:
	@echo "#"
	@echo "# make run FILE=source.gz"
	@echo "#"

# Indicate that the file is missing.
${FILE}:
	@echo "#"
	@echo "# File not found: FILE=${FILE}"
	@echo "#"
	@echo "# Please specify an existig file to re-compress."
	@echo "#"
	@echo "# make run FILE=source.gz"
	@echo "#"
	@exit 1

# Recompress the file.
run: ${FILE}
	@gunzip -c ${FILE} | bgzip -@ ${NCPU} -c > ${TMP}
	@mv -f ${TMP} ${FILE}
	@ls -lh ${FILE}

# Remove the target.
clean:
	rm -f ${TMP}

.PHONY: usage run run! clean
