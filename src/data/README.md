# Datasets 

Count datasets from various publications.

## golden

The counts for the Golden Snitch dataset in the RNA-Seq by Example book.

## barton

Data from the experiments described in:

* [Statistical models for RNA-seq data derived from a two-condition 48-replicate experiment][1]

[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4754627/

* [How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?][2]

[2]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/

* [Optimization of an RNA-Seq Differential Gene Expression Analysis Depending on Biological Replicate Number and Library Size]

[3]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5817962/

Data locations:

* [barton48_WT_raw][wt]
* [barton48_KO_raw][wt]


[wt]: https://figshare.com/articles/dataset/Wild_type_yeast_gene_read_counts_from_48_replicate_experiment/1425503
[ko]: https://figshare.com/articles/dataset/SNF2_knock_out_yeast_gene_read_counts_from_a_48_replicate_experiment/1425502

## Parameters for simulate_counts.r

* cheung: parameters from Cheung data. They are measurements from unrelated individuals,
so the dispersions are large.
* gilad: parameters from Gilad data. They are for Human liver sample comparisons between
male and female. This dataset has moderate dispersions.
* bottomly: parameters from Bottomly data. They compare two strains of inbred mice and the
within group dispersions are small.
* maqc: MAQC data which are technical replicates. There are no biological variation from the
replicates, so the dispersions are close to 0.
* GE.human: a vector of aggregated counts from all samples in Cheung data. This is used for
generating marginal expression distribution when provide sequnencing depth.
* pbmc: data from gene expression microarray data. This is for an example of using historical
data to estimate effect size
