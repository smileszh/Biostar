#!/usr/bin/env Rscript

#
# Create volcano from a differential expression count table.
#

# The name of the file that contains the DE results.
COUNT_FILE = "edger.csv"

# The name of the design file.
DESIGN_FILE = "design.csv"

# The name of the output file.
OUTPUT_FILE = "volcano.pdf"

# Significant variable
p = "FDR"

# Significant cutoff.
p_t = 0.05

# log2FoldChange cutoff
min_log2FoldChange = 1

# Plot width
WIDTH = 8

# Plot height.
HEIGHT = 8

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
bio.p <- c("gplots", "tibble", "dplyr", "tools","tidyverse")

cat("# Initializing ", bio.p, "...")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
)
cat(" done\n")

# Inform the user.
cat("# Tool: Create volcano", "\n")
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
  # The volcano row names will be taken from the first column.
  row_names = count_data[, 1, drop=T]
}

# Check that columns are present.
miss = sample_names[!(sample_names %in% names(count_data))]
if (length(miss)>0) {
  stop("# Sample not found in counts: ", miss)
}

# Extract the samples
coln = c("gene", "log2FoldChange", "PValue", "PAdj", "FDR")
deg = count_data[, coln]

# Screening for differentially expressed genes

k1 <- (deg[, p] < p_t)&(deg$log2FoldChange < -min_log2FoldChange)
k2 <- (deg[, p] < p_t)&(deg$log2FoldChange > min_log2FoldChange)
deg = mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))

# Differentially expressed gene statistics
up <- sum(deg$change == "up")
down <- sum(deg$change == "down")
stable <- sum(deg$change == "stable")

cat("# There are", up, "up-regulated genes.\n")
cat("# There are", down, "down-regulated genes.\n")
cat("# There are", stable, "unchanged genes.\n")

# Plot volcano

ggplot(data = deg, aes(x = log2FoldChange, y = -log10(.data[[p]]))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-min_log2FoldChange,min_log2FoldChange),lty=4,col="black",linewidth=0.8) +
  geom_hline(yintercept = -log10(p_t),lty=4,col="black",linewidth=0.8) +
  theme_bw()

ggsave(output_file, width = WIDTH, height = HEIGHT)

# Inform the user.
cat("# Output:", output_file, "\n")

