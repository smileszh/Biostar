#
# Biostar Workflows: https://www.biostarhandbook.com/
#
# Installs R based required packages into RStudio
#


# The list of packages that are already installed.
found = installed.packages()[,"Package"]

# Required core packages.
all.pkg <- c("optparse", "tibble", "dplyr")

# Missing core packages.
mis.pkg <- all.pkg[!(all.pkg %in% found )]

cat("# Installing core packages\n")

# Install missing core packages.
if(length(mis.pkg)) {
  install.packages(mis.pkg, repos = "https://repo.miserver.it.umich.edu/cran/")
}

# Initialize the BioConductor installer
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://repo.miserver.it.umich.edu/cran/")
}

# The list of packages that are already installed.
found = installed.packages()[,"Package"]

# Required BioConductor packages
bio.pkg <- c("edgeR", "DESeq2", "PROPER",
             "tximport", "biomaRt", "GenomeInfoDbData",
             "goseq")

# Missing BioConductor packages
bio.miss <- bio.pkg[!(bio.pkg %in% found)]

cat("# Installing BioConductor packages\n")

# Install missing BioConductor packages
if(length(bio.miss)) {
  BiocManager::install(bio.miss)
}


