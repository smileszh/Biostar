#
# Biostar Workflows: https://www.biostarhandbook.com/
#
# Install R based packages into a conda environment.
#
# Bash strict mode
set -ue

# Set the python version
PY='python=3.12'

echo "#"
echo "# Creating the salmon_env environment"
echo "#"

# Create a salmon environment
micromamba create -yn salmon_env salmon make parallel ${PY}

# Set the statistics package.
ENV=stats

echo "#"
echo "# Creating the stats environment"
echo "#"

# Create the stats environment.
micromamba create -yn ${ENV} ${PY}

# Install general R packages.
micromamba install -n ${ENV} -y r-optparse r-tibble r-dplyr r-gplots r-pacman python

# Install Bioconductor based packages.
micromamba install -n ${ENV} -y bioconductor-proper bioconductor-biomart bioconductor-edger bioconductor-deseq2 bioconductor-tximport bioconductor-goseq

# Install bio into the current environment
micromamba run -n ${ENV} pip install bio --upgrade

# Run the doctor to verify the packages
echo "#"
echo "# Verifying the installation"
echo "#"
micromamba run -n ${ENV} Rscript src/setup/doctor.r

