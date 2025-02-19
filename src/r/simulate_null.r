#!/usr/bin/env Rscript

# The simulated design file.
DESIGN_FILE = "design.csv"

# The simulated counts file.
COUNTS_FILE = "counts.csv"

# The number of replicates.
N_REPS = 5

# Minimum count
MIN_COUNT = 10

# Default input counts to select from.
INFILE = "src/data/barton_counts.csv"

# Detect R version and stop if under 4
if (as.numeric(R.version$major) < 4) {
  stop("R version 4 or higher is required")
}

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Command line options.
option_list = list(
  
  make_option(c("-d", "--design_file"), action="store", default=DESIGN_FILE, dest="design_file", type='character',
              help="simulated design file [%default]"),
  
  make_option(c("-o", "--output"), action="store", default=COUNTS_FILE, dest="counts_file", type='character',
              help="simulated counts  [%default]"),
  
  make_option(c("-r", "--nreps"), action="store", default=N_REPS, dest="n_reps", type='integer',
              help="the number of replicates [%default]"),
  
  make_option(c("-i", "--infile"), action="store", default=INFILE, dest="infile", type='character',
              help="input file with counts [%default]")
)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# Set the values from the options.
design_file = opt$design_file
counts_file = opt$counts_file
n_reps = opt$n_reps
infile = opt$infile

# Read the input data.
counts <- read.csv(infile, stringsAsFactors = F)

# Gene names come from the first column.
names <- counts[,1]

# Isolate the WT count data alone.
counts <- counts[ ,2:49]

# Build a new dataframe.
df = data.frame(name=names, FDR=1)

# The subset that we select.
set = counts[ ,sample(1:48, 2*n_reps)]

# Combine the counts into a data frame
df <- cbind(df, set)

# Generating A1, A2 ... and B1, B2 ... as sample names.
samples = colnames(set)

# Generate the groups: A and B
groups = c(rep("A", n_reps), rep("B", n_reps))

# The data column names
header = c("name", "FDR", samples)

# Set the column header on the subset.
colnames(df) <- header

# Write the simulated counts
write.csv(df, counts_file, quote = F, row.names = F)

# Create a design matrix.
ds = data.frame(sample=samples, group=groups)

# Save the design file.
write.csv(ds, design_file, quote = F, row.names = F)

# The number of genes
n_genes = nrow(counts)

# Number high expression genes.
n_high = sum(rowMeans(counts) >= MIN_COUNT)

# Inform the user on the actions.
cat("# Generating null data", "\n")
cat("# Input:", infile, "\n")
cat("# Total rows:", n_genes, "\n")
cat("# Above minimum:", n_high,"\n")
cat("# Replicates:", n_reps, "\n")
cat("# Design:", design_file,"\n")
cat("# Counts:", counts_file,"\n")

