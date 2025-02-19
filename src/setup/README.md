## Installing in R or RStudio

You have different options to install the required R packages.

If you are using RStudio set the working directory to the location of the code then execute:
   
```r
source("src/setup/init-rstudio.r")
```

Or you can run the following from the command line:

```bash
Rscript src/setup/init-rstudio.r"
```

## Installing into a conda environment

We recommend that you create a separate conda environment for the statistical analysis.

```bash
# Create new environment
bash src/setup/init-stats.sh
```




