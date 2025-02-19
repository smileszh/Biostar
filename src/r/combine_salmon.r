#!/usr/bin/env Rscript

#
# Combines salmon classification output with the tximport package.
#
# Also adds gene names to the columns.
#
# https://bioconductor.org/packages/release/bioc/html/tximport.html
#
# The program can also be run at the command line:
#
# Rscript combine_salmon.r -h
#

# The design file.
DESIGN_FILE = "design.csv"

# The file mapping transcripts to genes (NULL means don't use)
TX2GENE = NULL

# The output counts file.
OUTPUT_FILE = "counts.csv"

# The directory that salmon produced the results in.
DATA_DIR = 'salmon'

# Summarize over genes. Needs a TX2GENE file as well.
GENE_OUTPUT = F

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}


# Command line options.
option_list = list(

  make_option(c("-d", "--design_file"), action="store", default=DESIGN_FILE, dest="design_file", type='character',
              help="the design file for the samples [%default]"),

  make_option(c("-t", "--tx2gene"), action="store", default=TX2GENE, dest="txt2gene_file", type='character',
              help="the file with the gene2transcript mappings [%default]"),

  make_option(c("-o", "--out"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the output file name. [%default]"),

  make_option(c("-G", "--gene_count"), action="store_true", default=GENE_OUTPUT, dest="gene_output",
          help="summarize over genes"),

  make_option(c("-D", "--data_dir"), action="store", default=DATA_DIR, dest="data_dir", type='character',
              help="the data directory for the sample [%default]")

)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# Load the packages
library(tximport)

# The directory where the counts for each sample are located.
data_dir = opt$data_dir

# The design file must have a column named "sample"
design_file = opt$design_file

# The name of the file that contains transcript to gene mapping.
tx2gene_file = opt$txt2gene_file

# The output file name.
output_file = opt$output_file

# Keep the output at transcript level.
gene_output = opt$gene_output

# Create transcripts
txOut = !gene_output

# What software created the mappings.
method = "salmon"

# Should the method attempt to add gene names.
has_gene_names = !is.null(tx2gene_file)

# Inform the user.
cat("# Combine salmon results\n")
cat("# Sample: ", design_file, "\n")
cat("# Data dir: ", data_dir, "\n")


# Read the sample file
sample_data = read.csv(design_file, stringsAsFactors=F)

# Isolate the sample names.
sample_names = sample_data$sample

# Generate the file names that contain the quantification data.
files = file.path(data_dir,  sample_names, "quant.sf")

# We have a transcipt to gene name mapping.
if (has_gene_names) {
  
  cat("# Gene names:", tx2gene_file, "\n")
  
  # Read the files
  tx2gene = read.csv(tx2gene_file)
  
  # Back fill empty externals names with gene ids. Better then empty!
  empty = tx2gene$external_gene_name == ''
  tx2gene$external_gene_name[empty]=tx2gene$ensembl_gene_id[empty]
}

# Summarize the quantification data.
if (gene_output && has_gene_names) {
  
  # Print information
  cat("# Summarizing over genes\n:")
  
  # We have mapping data.
  tx = tximport(files, type=method, txOut=F, tx2gene=tx2gene, ignoreTxVersion=T)
  
} else {
  
  # Transcript level summary.
  tx = tximport(files, type=method, txOut=T)

}

# The lenght of each feature
len = tx$length[,1]

# Transform counts into a dataframe.
df = data.frame(tx$counts)

# Set the column names.
colnames(df) = sample_names

# Add lenghts to the counts
df$length = round(len,1)

# Unique names used during classification.
df$name = rownames(df)

# Default column order.
cols = c("name", "length", sample_names)

# Add the gene name to each transcript
if (has_gene_names){
  
  # Remove version number from the transcript ids.
  targets = sub("\\..*", "", df$name)
  
  # The default gene names are the targets
  df$gene = targets
  
  # Different columns need to be matched depending on quantification type.
  if (gene_output) {
    # Indices for the matching gene ids.
    indices = match(targets, tx2gene$ensembl_gene_id)

  } else {
    
    # Indices for the matching transcript ids.
    indices = match(targets, tx2gene$ensembl_transcript_id)
    
  }

  # Remove missing matches.
  #indices<-indices[!is.na(indices)]
  
  # Add the gene name column by slicing via indices.
  df$gene = tx2gene$external_gene_name[indices]
  
  df$gene [is.na(df$gene)] = df$name[is.na(df$gene)]
  
  # Missing values
  #miss.idx = is.null(df$gene == '')
  
  # Backfill the missing indices with the names.
  #df$gene[miss.idx] = df$name[miss.idx]

  # Build the new column order.
  cols = c("name", "gene", "length", sample_names)
  
} 

# Reorganize the columns.
out = df[, cols]

# Need to add gene lenght information

# Save the  resulting summarized counts.
write.csv(out, file = output_file,  row.names = FALSE, quote = FALSE)

# Inform the user.
cat("# Results: ", output_file, "\n")
