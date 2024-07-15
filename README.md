# Co-occurrence analysis of AMR genes

## Introduction 

Short read sequencing makes de novo genome assembly difficult at the site of repeated sequences, which cannot be bridged if the reads are shorter than the repeat length. The inability to traverse repetitive regions makes de novo genome assembly difficult at the site of repeated sequences, which cannot be bridged if the reads are shorter than the length of the repeat. 

![Example Plot](figures/the-short-read-sequencing-workflow.png)

We used co-occurrence analysis to overcome these limitations of short-read sequencing. Co-occurrence analysis allows for the comparison of genetic interactions across genomes. It can show which genes are carried together hinting at carriage on the same or different replicons. We ran a probabilistic analysis of gene co-occurrence to understand the genetic location of AMR genes relative to each other in short-read genomic data. 

## Methodology

We used SPAdes v3.11.1 to assemble 772 genomes of Escherichia coli and Klebsiella pneumoniae previously sequenced and isolated from blood and stool samples from Blantyre, Malawi. [Abricate](https://github.com/tseemann/abricate) v0.0.9 was used to screen for AMR genes against the Resfinder database (60% minimum length and 90% percentage identity).

We then used the [cooccur package](https://griffithdan.github.io/pages/code_and_data/cooccur.html) v1.3 in R v4.3.1 to employ a probabilistic model to determine whether each pair of genes had an observed co-occurrence which was significantly different from the expected co-occurrence, revealing genetic interactions among AMR genes. If a pair of genes were observed to co-occur significantly less than expected by chance their relationship was termed negative, if they were observed to co-occur significantly more it was termed positive and any pair of genes which were observed to co-occur with no significant difference from their expected value were termed random. A significant difference was defined as having a p-value ≤ 0.05. 

The final figure will look like this:

![Example Plot](figures/co-occurrence-probablisitic-heatmap-displaying-relationships-between-select-AMR-genes.png)

## How to run a co-occurrence analysis

You will find the complete code including instructions about how to run it [here](https://another-goodman.github.io/co-occurrence-analysis/Co-occurrence-analysis-of-AMR-genes.html).
