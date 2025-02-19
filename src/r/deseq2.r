#!/usr/bin/env Rscript

#
# Differential expression analysis with the DESeq2 package.
#
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#
# Command line install: mamba install bioconductor-deseq2
#

DESIGN_FILE = "design.csv"

# The name of the file that contains the counts.
COUNTS_FILE = "counts.csv"

# The output file.
OUTPUT_FILE = "deseq2.csv"

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

# Inform the user.
cat("# Running DESeq2", "\n")
cat("# Design:", design_file, "\n")
cat("# Counts:", counts_file, "\n")

# Read the design file.
inp <- read.csv(design_file, strip.white = T, stringsAsFactors = F)

# Check the sample column name.
cat("# Sample column:", sample_name, "\n")

if (!sample_name %in% colnames(inp)) {
  stop("# Sample column: ", sample_name, " not found in design file (use -s).")
}

# Check the factor column name.
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

# Create a new data frame to work with.
df = data.frame(
  sample = inp[, sample_name],
  groups = inp[, factor_name]
)

# The default groups come from the design
groups = unique(df$groups )

# Print information on the groups
for (grp in groups) {
  n_grp <- sum(df$groups == grp)
  cat("# Group", grp, "has", n_grp, "samples.\n")
}

# Turn conditions into factors.
df$group = factor(df$group)

# The first level should correspond to the first entry in the file
# Required later when building a model.
df$group = relevel(df$group, toString(df$group[1]))

# Isolate the sample names.
sample_names = df$sample

# Read the data from input file.
count_df = read.csv(counts_file, strip.white = T, stringsAsFactors = F)

# Created rounded integers for the count data
count_data = round(count_df[, sample_names, drop=F])

# Column names should be carried over. Skip sample names and FDR
other_cols = names(count_df)[!names(count_df) %in% c(sample_names, "FDR")]

# Carry over all additional columns not in sample names.
row_data = count_df[, other_cols, drop=F]

#
# Running DESeq2
#

# Bioconductor packages used in the script.
bio.p <- c("DESeq2", "tibble", "dplyr", "tools")

# Load Bioconductor packages
cat("# Initializing ", bio.p, "...")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
)
cat(" done\n")

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(
  countData = count_data,
  rowData = row_data,
  colData = df,
  design = ~group
)

# Number of rows before
n_start = nrow(dds)

# Remove rows with no counts
keep <- rowSums(counts(dds)) > 3

# At least 3 samples with a count of 10 or higher
# keep <- rowSums(counts(dds) >= 10) >= 3

# Apply the filtering
dds <- dds[keep, ]

# Run deseq
dds = DESeq(dds)

# Get the normalized counts.
normed = counts(dds, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

# Format the results.
res = results(dds)

#
# The rest of the code is about formatting the output dataframe.
#

# The first columns that carry over from the input.
start =  dplyr::as_tibble(rowData(dds))[, other_cols]

# The main tibble
data = bind_cols(start, as_tibble(res))

# Create the foldChange column.
data$foldChange = 2 ^ data$log2FoldChange

# Rename the columns.
data =  dplyr::rename(data, PValue=pvalue, FDR=padj)

# Set all NA values to 1.
data$FDR[is.na(data$FDR)] <- 1

# Create a real adjusted pvalue
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Create the additional columns that we wish to present.
data$baseMeanA = 1
data$baseMeanB = 1

# Combine the normed matrix with the results
total <- bind_cols(data, normed)

# Sort again by P-value.
total = arrange(total, FDR)

# Compute the false discovery counts on the sorted table.
total$falsePos = 1:nrow(total) * total$FDR

# Sample names for condition A and B
col_names_A <- sample_names [df$group == groups[[1]]]
col_names_B <- sample_names [df$group == groups[[2]]]

#col_names_A = unlist(df %>% filter(group==groups[1]) %>% select(sample))
#col_names_B = unlist(df %>% filter(group==groups[2]) %>% select(sample))

# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])

# Bringing some sanity to numbers. Round columns to fewer digits.
total$foldChange = round(total$foldChange, 3)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
total$FDR = round(total$FDR, 4)
total$falsePos = round(total$falsePos, 0)

# Reorganize columns names to make more sense.
new_cols = c(other_cols, 
             "baseMean","baseMeanA","baseMeanB",
             "foldChange", "log2FoldChange",
             "lfcSE","stat","PValue","PAdj", 
             "FDR","falsePos", col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]

# Count the number of significant results.
n_total = nrow(total)
fdr_count = sum(total$FDR < 0.05, na.rm=TRUE)
fdr_frac = round(fdr_count/n_total, 3) * 100

pval_count = sum(total$PValue < 0.05, na.rm=TRUE)
pval_frac = round(pval_count/n_total, 3) * 100

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)

# Write the results to the standard output.
write.csv(total, file=output_file, row.names = F, quote = F)

# Inform the user.
cat("# Input:", n_start, "rows\n")
cat("# Removed:", n_start-n_total, "rows\n")
cat("# Fitted:", n_total, "rows\n")
cat("# Significant PVal:", pval_count, "(", pval_frac, "%)\n")
cat("# Significant FDRs:", fdr_count, "(", fdr_frac, "%)\n")
cat("# Results:", output_file, "\n")
