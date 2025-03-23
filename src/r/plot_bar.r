#!/usr/bin/env Rscript
#
# Create barplot from variant group.
#

# The name of the data file,
DATA_FILE = "data.csv"

# The name of the output file.
OUTPUT_FILE = "bar.pdf"

# Plot width
WIDTH = 7

# Plot height.
HEIGHT = 7

# Set the margins
MARGINS = c(9, 12)

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Command line options.
option_list = list(
  
  make_option(c("-d", "--data_file"), action="store", default=DATA_FILE, dest="data_file", type='character',
              help="the data file for the group [%default]"),
  
  make_option(c("-o", "--out"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the  results file  [%default]")
)


# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# The name of the file that contains the data.
data_file = opt$data_file

# The output file.
output_file = opt$output_file

# Packages used in the script.
bio.p <- c("ggpubr","tidyverse")

cat("# Initializing ", bio.p, "...")
status = suppressPackageStartupMessages(
  lapply(bio.p, library, character.only = TRUE)
)
cat(" done\n")

# Read data from the standard input.
data = read.csv(data_file, header=T, as.is=TRUE)

# Ensure that group columns are factor.
data$group <- factor(data$group)

# Automatically generates all two-by-two comparison combinations
group_levels <- levels(data$group)
comparisons_list <- combn(group_levels, 2, simplify = FALSE)

# Plot bar plot.

ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_sd", 
               width = 0.2, size = 0.5) +
  geom_jitter(aes(color = sample), width = 0.1, size = 1, alpha = 0.5, show.legend = FALSE) +
  stat_compare_means(
    method = "t.test",  
    comparisons = comparisons_list, 
    label = "p.signif", 
    tip.length = 0.03,
    size = 6,
    vjust = 0.5
  ) +
  scale_fill_manual(values = c("#1F77B4", "#E15759")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(size = 8),
        axis.line = element_line(color = "black", size = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(output_file, width = WIDTH, height = HEIGHT, units = "cm")

# Inform the user.
cat("# Output:", output_file, "\n")

