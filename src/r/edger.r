#!/usr/bin/env Rscript

#
# Differential expression analysis with the edgeR package.
# 
# https://bioconductor.org/packages/release/bioc/html/edgeR.html
#
# Command line install: mamba install bioconductor-edger r-tibble r-dplyr r-optparse
#

# The input counts file.
COUNTS_FILE = "counts.csv"

# The input design file.
DESIGN_FILE="design.csv"

# The default output file.
OUTPUT_FILE = "edger.csv"

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
              help="input count file  [%default]"),

  make_option(c("-o", "--out"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the output file [%default]"),

  make_option(c("-m", "--method"), action="store", default="glm", dest="method", type='character',
              help="the method glm/classic [%default]"),

  make_option(c("-f", "--factor_name"), action="store", default="group", dest="factor_name", type='character',
              help="the factor column name in the design file [%default]"),
  
  make_option(c("-s", "--sample_name"), action="store", default="sample", dest="sample_name", type='character',
              help="the sample column name in the design file [%default]")

)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# The name of the file that contains the counts.
counts_file = opt$counts_file

# The sample file is in CSV format and must have the headers "sample" and "group".
design_file = opt$design_file

# The final result file.
output_file = opt$output_file

# The final result file.
method = opt$method

# The name of the factor column
factor_name = opt$factor_name

# The name of the sample column
sample_name = opt$sample_name

# Bioconductor packages used in the script.
pkgs <- c("edgeR", "tibble", "dplyr", "tools")

# Load required packages
cat("# Initializing", pkgs, "...")
status = suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)
cat(" done\n")

# Inform the user.
cat("# Tool: edgeR", "\n")
cat("# Design:", design_file, "\n")
cat("# Counts:", counts_file, "\n")

# Read the sample file.
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
groups = unique(df$groups )

# Print the factors
cat("# Factors:", groups, "\n")

# Check that the lenght of the groups is 2.
if (length(groups) != 2) {
  stop("# The number of factors must be 2.")
}

# Print information on the groups
for (grp in groups) {
  n_grp <- sum(df$groups == grp)
  cat("# Group", grp, "has", n_grp, "samples.\n")
}

# Turn conditions into factors.
df$groups = factor(df$groups)

# The first level should correspond to the first entry in the file
# Required later when building a model.
df$groups = relevel(df$groups, toString(df$groups[1]))

# Isolate the sample names.
sample_names <- df$sample

# Sample names for the two groups.
col_names_A <- sample_names [df$group == groups[[1]]]
col_names_B <- sample_names [df$group == groups[[2]]]

# Read the data from the counts file.
counts_df = read.csv(counts_file, strip.white = T, stringsAsFactors = F)

# How many rows at the start.
n_start = nrow(counts_df)

# Created rounded integers for the count data
counts_mat = round(counts_df[, sample_names])

# Set the rownames on the matrix.
row.names(counts_mat) = counts_df[,1]

# Column names that should be carried over. Skip sample names and FDR (if exists)
other_cols = names(counts_df)[!names(counts_df) %in% c(sample_names, "FDR")]

# The starting columns of the dataframe
start_df = counts_df[, other_cols, drop=F]

# Set the rownames as well.
row.names(start_df) = start_df[,1]

# Group selector
group <- df$groups

#
# The steps below will perform the edgeR analysis.
#
# Select the GLM method
if (method=='glm') {

  cat("# Method: glm", "\n")
  
  # Following the recommended edger vignette here verbatim
  dge <- DGEList(counts=counts_mat, group=group)
  
  # Filter by expression.
  keep <- filterByExpr(dge)
  
  # Filter the expression
  dge <- dge[keep, ,keep.lib.sizes=FALSE]
  
  # The design matrix.
  design <- model.matrix(~group)
  
  # Calculate the normalization factors
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit,coef=2)
  etp = topTags(qlf, n=Inf)
  
} else {
 
  # 
  # This is the classic mode of running EdgeR
  #
  cat("# Method: classic", "\n")
  # Creates a DGEList object from a table of counts and group.
  dge <- DGEList(counts=counts_mat, group=group)
  
  # Filter by expression.
  keep <- filterByExpr(dge)
  
  # Filter the expression
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  # Maximizes the negative binomial conditional common likelihood to estimate a common dispersion value across all genes.
  dge <- estimateCommonDisp(dge)
  
  # Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
  dge <- estimateTagwiseDisp(dge)
  
  # Compute genewise exact tests for differences in the means between the groups.
  etx <- exactTest(dge)
  
  # Extracts the most differentially expressed genes.
  etp <- topTags(etx, n=Inf)
} 
  
#
# From here on out is about formatting the output.
#

# Get the scale of the data
scale = dge$samples$lib.size * dge$samples$norm.factors

# Get the normalized counts
normed = round(t(t(counts_mat)/scale) * mean(scale))

# Turn the results into a data frame.
data = merge(start_df, etp$table, by="row.names")

# Add the rownames
row.names(data) = data[,1]

# Get rid of extra column it gained
data = data[, 2:ncol(data)]

# Create column placeholders.
data$baseMean = 1
data$baseMeanA = 1
data$baseMeanB = 1
data$foldChange = 2 ^ data$logFC
data$log2FoldChange=data$logFC
data$falsePos = 1

# Compute the adjusted p-value
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Create a merged output that contains the normalized counts.
total <- merge(data, normed, by='row.names')

# Get rid of extra column it gained
total = total[, 2:ncol(total)]

# Sort again by P-value.
total = arrange(total, PValue)

# Compute the false discovery counts on the sorted table.
total$falsePos = 1:nrow(total) * total$FDR

# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])
total$baseMean = total$baseMeanA + total$baseMeanB

# Round the numbers to make them look better
#total$foldChange = round(total$foldChange, 3)
#total$FDR = round(total$FDR, 4)
#total$PAdj = round(total$PAdj, 4)
#total$logCPM = round(total$logCPM, 1)
#total$log2FoldChange = round(total$log2FoldChange, 3)
total$baseMean = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB = round(total$baseMeanB, 1)
total$falsePos = round(total$falsePos, 0)

# The total numbers
n_total = nrow(total)

# Count the number of significant results.
fdr_count = sum(total$FDR < 0.05)
fdr_frac = round(fdr_count/n_total, 3) * 100

pval_count = sum(total$PValue < 0.05)
pval_frac = round(pval_count/n_total, 3) * 100

# Reorganize columns names to make more sense.
new_cols = c( other_cols, 
             "baseMean","baseMeanA","baseMeanB",
             "foldChange", "log2FoldChange",
             "PValue","PAdj", 
             "FDR","falsePos", col_names_A, col_names_B)

# Rearrange the data frame with new column order.
total = total[, new_cols]

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 2)
total$PValue = formatC(total$PValue, format = "e", digits = 2)

# Write the results to the standard output.
write.csv(total, file=output_file, row.names = F, quote = F)

# Compute more stats.
n_removed = n_start - n_total

# Nicer formatting
pval_count = sprintf("%4d", pval_count)
pval_frac = sprintf("%.2f", pval_frac)
fdr_count = sprintf("%4d", fdr_count)
fdr_frac = sprintf("%.2f", fdr_frac)

# Inform the user.
cat("# Input:", n_start, "rows\n")
cat("# Removed:", n_removed, "rows\n")
cat("# Fitted:", n_total, "rows\n")
cat("# Significant PVal:", pval_count, "(", pval_frac, "%)\n")
cat("# Significant FDRs:", fdr_count, "(", fdr_frac, "%)\n")
cat("# Results:", output_file, "\n")
