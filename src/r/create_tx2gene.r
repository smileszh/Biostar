#!/usr/bin/env Rscript

#
# Obtain gene names for transcripts with the biomaRt package.
#
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html
#
# Command line install: mamba install bioconductor-biomart
#

OUTPUT_FILE = "tx2gene.csv"

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Create the maoping of the command line arguments.
option_list = list(
  make_option(c("-d", "--dataset"), action="store", default='mmusculus_gene_ensembl', dest="dataset", type='character',
              help="the dataset under study [%default]"),

  make_option(c("-o", "--out"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the output file name. [%default]"),

  make_option(c("-s", "--show"), action="store_true", default=F, dest="show_all",
          help="list the valid databases")
)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# Load the biomart manager
library(biomaRt)

# The biomaRt dataset name.
dataset <- opt$dataset

# The ouput file name.
output_file = opt$output_file

# Show the available datasets.
if (opt$show_all) {
    # Connect to the database.
    mart <- useEnsembl(biomart = "ensembl")

    ## list the available datasets in this Mart
    listDatasets(mart=mart)

} else {

	# Print output
	cat("# Create tx2gene mapping\n")
	cat("# Dataset:", dataset, "\n")

	# Make a connection to the validated dataset.
	cat("# Connecting to ensembl.\n")
	mart <- useEnsembl(dataset=dataset, biomart = 'ensembl')

	# The first column must match the feature id used during quantification.
	attributes <- c(
	  "ensembl_transcript_id",
	  "ensembl_gene_id",
	  "transcript_length",
	  "external_gene_name"
	)

	cat("# Submitting the query.\n")

	# Perform the query.
	data <- biomaRt::getBM(attributes = attributes, mart = mart)

	# Save the data into a file.
	write.csv(data, output_file, row.names=FALSE, quote=FALSE)

	# Inform the user.
	cat("# Output:", output_file, "\n")
}
