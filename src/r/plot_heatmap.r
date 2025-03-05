#!/usr/bin/env Rscript

#
# Create heat map from a differential expression count table.
#

# The name of the file that contains the DE results.
COUNT_FILE = "edger.csv"

# The name of the design file.
DESIGN_FILE = "design.csv"

# The name of the output file.
OUTPUT_FILE = "heatmap.pdf"

# FDR cutoff.
MIN_FDR = 0.05

# Plot width
WIDTH = 12

# Plot height.
HEIGHT = 13

# Set the margins
MARGINS = c(9, 12)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Command line options.
option_list = list(

  make_option(c("-d", "--design_file"), action="store", default=DESIGN_FILE, dest="design_file", type='character',
              help="the design file for the samples [%default]"),

  make_option(c("-c", "--counts"), action="store", default=COUNT_FILE, dest="counts_file", type='character',
              help="input counts file  [%default]"),

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

# The name of the factor column
factor_name = opt$factor_name

# The name of the sample column
sample_name = opt$sample_name

# Packages used in the script.
bio.p <- c("gplots", "tibble", "dplyr", "tools")

cat("# Initializing ", bio.p, "...")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
)
cat(" done\n")

# Inform the user.
cat("# Tool: Create heatmap", "\n")
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
groups = unique(df$group)

# Print information on the groups
for (grp in groups) {
  n_grp <- sum(df$groups == grp)
  cat("# Group", grp, "has", n_grp, "samples.\n")
}

# Read normalized counts from the standard input.
count_data = read.csv(counts_file, header=T, as.is=TRUE)

# Subset data for values under a treshold.
count_data = subset(count_data, count_data$FDR <= MIN_FDR)

# Sanity check.
if (nrow(count_data) == 0) {
  cat("# Warning: The count data has no rows that pass the FDR cutoff.\n")
  quit()
}

# Isolate the sample names.
sample_names = df$sample

# Use the gene column if available.
if( "gene" %in% colnames(count_data)) {
  cat("# Found gene name column.\n")
  row_names = count_data$gene
} else {
  # The heatmap row names will be taken from the first column.
  row_names = count_data[, 1, drop=T]
}

# Check that columns are present.
miss = sample_names[!(sample_names %in% names(count_data))]
if (length(miss)>0) {
  stop("# Sample not found in counts: ", miss)
}

# Extract the samples
count_matrix = as.matrix(count_data[, sample_names])

# Adds a little noise to each element to avoid the
# clustering function failure on zero variance rows.
count_matrix = jitter(count_matrix, factor = 1, amount = 0.00001)

# Normalize each row to a z-score
# Zscore scaling with a trick I read somewhere else.
count_zscores = t(scale(t(count_matrix)))

# Normalize each row to a z-score
# count_zscores = NULL
# for (i in 1 : nrow(count_matrix)) {
#   row = count_matrix[i,]
#   zrow = (row - mean(row)) / sd(row)
#   count_zscores = rbind(count_zscores, zrow)
# }

# Set the row names on the zscores.
row.names(count_zscores) = row_names

# Turn the data into a matrix for heatmap2.
count_zscores = as.matrix(count_zscores)

# Set the color palette.
col = greenred

# Create a PDF device
pdf(output_file, width = WIDTH, height = HEIGHT)

heatmap.2(count_zscores, col=col, density.info="none", Colv=NULL,
    dendrogram="row", trace="none", margins=MARGINS, lhei=LHEI)

# Close the device.
dev.off()

# Inform the user.
cat("# Output:", output_file, "\n")


