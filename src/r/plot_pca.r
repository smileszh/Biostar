#!/usr/bin/env Rscript

#
# Differential expression analysis with the DESeq2 package.
#
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#
# Command line install: mamba install bioconductor-deseq2
#

# How much to move the label
NUDGE_Y = 1

DESIGN_FILE = "design.csv"

# The name of the file that contains the counts.
COUNTS_FILE = "counts.csv"

# The output file.
OUTPUT_FILE = "pca.pdf"

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Command line options.
option_list = list(

  make_option(c("-d", "--design_file"), action="store", default=DESIGN_FILE, dest="design_file", type='character',
              help="the design file for the samples [%default]"),

  make_option(c("-c", "--counts"), action="store", default=COUNTS_FILE, dest="counts_file", type='character',
              help="input counts file  [%default]"),

  make_option(c("-n", "--nudge"), action="store", default=NUDGE_Y, dest="nudge_y", type='double',
              help="how much to move the label [%default]"),

  make_option(c("-o", "--out"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the  results file  [%default]"),

  make_option(c("-f", "--factor_name"), action="store", default="group", dest="factor_name", type='character',
              help="the factor column name in the design file [%default]"),
  
  make_option(c("-s", "--sample_name"), action="store", default="sample", dest="sample_name", type='character',
              help="the sample column name in the design file [%default]")
  
  )

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))


# The name of the file that contains the counts.
counts_file = opt$counts_file

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = opt$design_file

# The output file.
output_file = opt$output_file

# The output file.
nudge_y = opt$nudge_y

# The name of the factor column
factor_name = opt$factor_name

# The name of the sample column
sample_name = opt$sample_name

# Inform the user.
cat("# Generating PCA plot", "\n")
cat("# Design:", design_file, "\n")
cat("# Counts:", counts_file, "\n")

# Read the design file.
inp <- read.csv(design_file, strip.white = T, stringsAsFactors = F)

# Sample column name.
cat("# Sample column:", sample_name, "\n")

if (!sample_name %in% colnames(inp)) {
  stop("# Sample column: ", sample_name, " not found in design file (use -s).")
}

# Factor column name.
cat("# Factor column:", factor_name, "\n")
if (!factor_name %in% colnames(inp)) {
  stop("# Factor column: ", factor_name, " not found in design file (use -f).")
}

# Create a new data frame to work with.
df = data.frame(
  sample = inp[, sample_name],
  groups = inp[, factor_name]
)

# The default groups come from the design
groups = unique(df$groups)

# Bioconductor packages used in the script.
bio.p <- c("DESeq2", "tibble", "dplyr")

# Print information on the groups
for (grp in groups) {
  n_grp <- sum(df$groups == grp)
  cat("# Group", grp, "has", n_grp, "samples.\n")
}

# Load Bioconductor packages
cat("# Initializing ", bio.p, "...")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
)
cat(" done\n")

# Turn conditions into factors.
df$groups = factor(df$groups)

# The first level should correspond to the first entry in the file
df$groups = relevel(df$groups, toString(df$groups[1]))

# Isolate the sample names.
sample_names = df$sample

# Read the data from input file.
count_df = read.csv(counts_file, strip.white = T, stringsAsFactors = F, check.names=F)

# Created rounded integers for the count data
count_data = as.matrix(round(count_df[, sample_names]))

# Sanity check.
if (nrow(count_data) == 0) {
  cat("# Warning: The count data has no rows.\n")
  quit()
}

# Convert counts to integer.
mode(count_data) <- "integer"

# Initialize the DESEq2 dataset.
dds = DESeqDataSetFromMatrix(
  countData = count_data,
  colData = df,
  design = ~groups,
)

# Find the normalization factors.
dds <- estimateSizeFactors(dds, quiet = TRUE)

# Count the number of valid genes.
nsub = sum( rowMeans( counts(dds, normalized=T)) > 5 )

# Generate the variance stabilization.
vsd <- vst(dds, nsub = nsub)

# Load the plotting library.
suppressPackageStartupMessages(library(ggplot2))

# Create the PCA plot.
pdf(output_file)

# The title of the plot.
title = paste("PCA plot")

# Here we may need to move the labels around a bit.
nudge <- position_nudge(y = nudge_y)

# Generate the PCA plot.
plotPCA(vsd, intgroup=c("groups")) +
  geom_text(aes(label = colnames(vsd)), position=nudge, size = 2.5) +
  ggtitle(title)

# Close the device.
tmp = dev.off()

# Message to the user.
cat("# PCA plot:", output_file, "\n")
