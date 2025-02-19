#
# Download data via rsync
#
# You might ask yourself Why have a makefile for what seems to be be trivial command?
#
# In this case to help us remember the correct commands and to make it fit with
# the rest of the workflow building process.
#
# The rsync rules everyone always forgets:
#
#   - trailing slash on URL means copy the contents of the directory
#   - no trailing slash on URL means copy the directory itself
#   - a slash on the destination directory has no effect
#

# The remote location we wish to download
URL ?= rsync://hgdownload.soe.ucsc.edu/goldenPath/eboVir3/bigZips

# The destination that the download will be placed in.
DEST ?= refs

# Makefile customizations.
SHELL := bash
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# General usage information.
usage::
> @echo ""
> @echo "# rsync.mk: downloads data via rsync protocol"
> @echo ""
> @echo "# make run URL=? FILE=?"
> @echo ""

# rsync will automatically detect updated files.
run::
> mkdir -p ${DEST}
> rsync --times --copy-links --recursive -vz -P ${URL} ${DEST}
> find ${DEST}/*

run!::
> @echo "# cannot undo an rsync download"

# Installation instructions
install::
> @echo "# no installation required"

