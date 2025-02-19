#!/usr/bin/env Rscript

#
# Compares the overlap between two DE runs.
#
# Reads data under an FRD cutoff and prints the true positive, false positive etc
#

# The name of the first file.
FILE1 = "counts.csv"

# The name of the second file.
FILE2 = "edger.csv"

# Check and install packages
all.p <- c("optparse", "tibble", "dplyr")
mis.p <- all.p[!(all.p %in% installed.packages()[, "Package"])]
if (length(mis.p)) {
  install.packages(mis.p, repos = "https://repo.miserver.it.umich.edu/cran/")
}

# Load core packages
status = suppressPackageStartupMessages(
  lapply(all.p, require, character.only = TRUE)
)

# Command line options.
option_list = list(

  make_option(c("-a", "--file1"), action = "store", default = FILE1, dest = "file1", type = 'character',
              help = "the query file [%default]"),

  make_option(c("-b", "--file2"), action = "store", default = FILE2, dest = "file2", type = 'character',
              help = "the target file [%default]"),

  make_option(c("-c", "--colname"), action = "store", default = "name", dest = "col_name", type = 'character',
              help = "the column name to match on [%default]"),

  make_option(c("-p", "--pval_name"), action = "store", default = "FDR", dest = "pval_name", type = 'character',
              help = "the column name for the pvalue [%default]"),

  make_option(c("-t", "--pval_level"), action = "store", default = 0.05, dest = "pval_level", type = 'double',
              help = "the pvalue cutoff [%default]"),

  my_option <- make_option( c("-s", "--show_common"), action="store_true", dest = "show_common", type = "logical",  default = F,
                            help = "prints common genes"),

  make_option(c("-o", "--out"), action = "store", default = "summary.csv", dest = "out_file", type = 'character',
              help = "the name of the output file [%default]")

)


# Parse the options.
opt = parse_args(OptionParser(option_list = option_list))

# The name of the first file.
file1 = opt$file1

# The name of the second file.
file2 = opt$file2

# The name of the column to find identities over
col_name = opt$col_name

# The pvalue threshold to operate on.
pval_name = opt$pval_name

# The pvalue threshold to operate on.
pval_level = opt$pval_level

# The name of the output file.
output_file = opt$out_file

show_common = opt$show_common

# Read the sample file.
data1 <- read.csv(file1, stringsAsFactors = F, strip.white = T, row.names = NULL)

# Read the sample file.
data2 <- read.csv(file2, stringsAsFactors = F, strip.white = T, row.names = NULL)

# Finds the differentially expressed names
de1 <- subset(data1, data1[, pval_name] <= pval_level)
de2 <- subset(data2, data2[, pval_name] <= pval_level)

# Extract the gene names only
name1 = de1[, col_name]
name2 = de2[, col_name]

# Subset only in de1 and de2
only_de1 = de1[!(name1 %in% name2),]
only_de2 = de2[!(name2 %in% name1),]

# Totals
tot1 = length(name1)
tot2 = length(name2)

# Intersect: common elements.
common <- intersect(name1, name2)


if (show_common) {
  # Print the common elements
  cat("gene\n")
  cat(common, sep = "\n")
} else {
  # Report the differences
  cat("# Tool: evaluate_results.r", "\n")
  cat("#", tot1, "in", file1, "\n")
  cat("#", tot2, "in", file2, "\n")
  cat("#", length(common), "found in both\n")
  cat("#", nrow(only_de1), "found only in", file1, "\n")
  cat("#", nrow(only_de2), "found only in", file2, "\n")
  cat("# Summary:", output_file, "\n")
}


# Generate the summary matrix
con <- file(output_file, open = "wt")

# Write elements only in file 1
writeLines(
  paste0("#\n", "# Found only in ", file1, "\n"), con
)

write.csv(only_de1, con, quote = F, row.names = F)

# Write elements only in file 2
writeLines(
  paste0("#\n", "# Found only in ", file2, "\n"), con
)
write.csv(only_de2, con, quote = F, row.names = F)

close(con)

