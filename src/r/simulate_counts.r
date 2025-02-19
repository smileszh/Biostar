#!/usr/bin/env Rscript

#
# Read count simulator with PROPER
#
# Wu H, Wang C, Wu Z (2014). "PROPER: Comprehensive Power Evaluation for Differential Expression using RNA-seq." Bioinformatics.
#
# https://bioconductor.org/packages/release/bioc/html/PROPER.html
#
# BiocManager::install("PROPER")
#

# The simulated design file.
DESIGN_FILE = "design.csv"

# The simulated counts file.
COUNTS_FILE = "counts.csv"

# The number of simulated genes.
N_GENES = 20000

# The number of replicates.
N_REPS = 3

# The minimum number of reads to consider a gene DE.
MIN_COUNT = 10

# Set the model name (see the documentation for PROPER)
# The name of the data that generated the simulation mode (different error models)
# Valid values from high dispersion to low 1-4:  bottomly, maqc, gilad, cheung,
E_LEVEL = 1

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

  make_option(c("-n", "--ngenes"), action="store", default=N_GENES, dest="n_genes", type='integer',
              help="the  results file  [%default]"),

  make_option(c("-r", "--nreps"), action="store", default=N_REPS, dest="n_reps", type='integer',
              help="the number of replicates [%default]"),

  make_option(c("-e", "--error"), action="store", default=E_LEVEL, dest="e_level", type='integer',
              help="the dispersion level (0-4) [%default]"),

  make_option(c("-s", "--seed"), action="store", default=0, dest="n_seed", type='integer',
              help="set the seed for the simulation [%default]")
)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# Set the values from the options.
design_file = opt$design_file
counts_file = opt$counts_file
n_genes = opt$n_genes
n_reps = opt$n_reps
e_level = opt$e_level
n_seed= opt$n_seed

if (e_level == 0) {
  e_name = "maqc"
} else if (e_level == 1) {
  e_name = "bottomly"
} else if (e_level == 2) {
    e_name = "gilad"
} else {
    e_name = "cheung"
}

# Load required Bioconductor packages.
pkgs <- c("PROPER")
cat("# Initializing ", pkgs, "...")
status = suppressPackageStartupMessages(
  lapply(pkgs, library, character.only = TRUE)
)
cat(" done\n")

# Set the seed. Random seed if not specified.
if (n_seed == 0) {
    seed <- sample.int(100,size=1)
} else {
    seed <- n_seed
}

# Generate the simulation options.
opts = RNAseq.SimOptions.2grp(ngenes=n_genes, sim.seed=seed,
       lBaselineExpr=e_name,
       lOD=e_name,
)

# Run the simulation.
sim = simRNAseq(opts, n1=n_reps, n2=n_reps)

# The simulation matrix.
mat <- sim$counts

# Generating A1, A2 ... and B1, B2 ... as sample names.
s1 = paste(rep("A", n_reps), 1:n_reps, sep="")
s2 = paste(rep("B", n_reps), 1:n_reps, sep="")

# Combine the sample names.
samples = c(s1, s2)

# Generate the groups: A and B
groups = c(rep("A", n_reps), rep("B", n_reps))

# Set the amples names on the simulation.
colnames(mat) <- samples

# Create gene names for the simulated counts.
names = paste(rep("GENE-"), 1:n_genes, sep="")

# Create the dataframe. All genes have NO differential expression.
df = data.frame(name=names, state="NO", FDR=1)

# Differentially expressed genes get a YES value in the state column.
df$state[sim$DEid] = "YES"

# Set the FDR on the DE genes to 1.
df$FDR[sim$DEid] = 0

# How many genes changes.
n_changed = length(sim$DEid)

# Combine the metadata with the counts.
df = cbind(df, mat)

# Set the FDR back to 1 for rows with insufficient readcounts.
df[rowMeans(mat) < MIN_COUNT,]$FDR = 1

# Number high expression genes.
n_high = nrow(df[rowMeans(mat) >= MIN_COUNT,])

# The number of DE genes with data
n_detect = nrow(df[df$FDR<1, ])

# Sort to put low FDR first.
df = df[order(df$FDR),]

# Save the simulation matrix.
write.csv(df, counts_file, quote = F, row.names = F)

# Create a designe matrix.
ds = data.frame(sample=samples, group=groups)

# Save the design file.
write.csv(ds, design_file, quote = F, row.names = F)

# Inform the user on the actions.
cat("# PROspective Power Evaluation for RNAseq", "\n")
cat("# Error level: ", e_level, " (", e_name, ")\n", sep="")
cat("# All genes:", n_genes, "\n")
cat("# Genes with data:", n_high,"\n")
cat("# Genes that changed:", n_changed,"\n")
cat("# Changes we can detect:", n_detect,"\n")
cat("# Replicates:", n_reps, "\n")
cat("# Design:", design_file,"\n")
cat("# Counts:", counts_file,"\n")
