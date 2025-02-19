#!/usr/bin/env Rscript

#
# Adds gene name to a featurecounts output file.
#

# The design file must have a column named "sample"
COUNTS_FILE = "counts.txt"

# The name of the file that contains transcript to gene mapping.
TX2GENE_FILE = ""

# The output file name.
OUTPUT_FILE= "counts.csv"

# Install the optparse package automatically.
if(!require(optparse, quietly=T)){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library(optparse)
}

# Command line options.
option_list = list(

  make_option(c("-c", "--counts"), action="store", default=COUNTS_FILE, dest="counts_file", type='character',
              help="input count file  [%default]"),

  make_option(c("-o", "--output_file"), action="store", default=OUTPUT_FILE, dest="output_file", type='character',
              help="the design file for the samples [%default]"),

  make_option(c("-t", "--tx2gene"), action="store", default=TX2GENE_FILE, dest="tx2gene_file", type='character',
              help="the transcript to gene mapping [%default]")

)

# Parse the options.
opt = parse_args(OptionParser(option_list=option_list))

# The name of the file that contains the counts.
input_file = opt$counts_file

# The final result file.
output_file = opt$output_file

# The transcript to gene mapping file.
tx2gene_file = opt$tx2gene_file

# Inform the user.
cat("# Reformating featurecounts.\n")
cat("# Input:", input_file, "\n")

# Read the input data
data = read.table(input_file, stringsAsFactors = FALSE, header=T, sep="\t")

# Remove certain columns.
data = subset(data, select = -c(2, 3, 4, 5) )

# Remove file related information from the column
headers = sub("\\.bam", "", colnames(data))
headers = sub("bam\\.", "", headers)

# Attach the new names
colnames(data) <- headers

# Add columns with specific names
data$name = data$Geneid
data$gene = data$Geneid

targets =  sub("\\..*", "", data$name)

if (tx2gene_file != '') {
  cat("# Tx2gene:", tx2gene_file, "\n")

  # Read the transcript to gene mapping file.
  tx2gene = read.csv(tx2gene_file)

  # Find the indices where the targets match to the mapping.
  indices = match(targets, tx2gene$ensembl_gene_id)

  # Add the gene name column.
  data$gene = tx2gene$external_gene_name[indices]
  
  data$gene[is.na(data$gene)] = targets[is.na(data$gene)]
  
  data$gene[data$gene ==""] = targets[data$gene ==""]
  
} else {
  data$gene = targets
}

# List the desired column order.
cols <- c("name", "gene", headers[3:length(headers)])

# Slice the data.
out = data[,cols]

# Save the  resulting summarized counts.
write.csv(out, file=output_file,  row.names = FALSE, quote = FALSE)

cat("# Output:", output_file, "\n")


