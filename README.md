
# transGI: Nonparametric transformation method of transcriptome based on prior gene interactions

The interaction between genes can be represented as a scale-free
directed network. A common practice is to infer the strength of signals
for genes and their downstream functions based on gene expression
levels. Although this approach is simple, it does not consider prior
information on gene-gene interactions.

To incorporate prior information on gene interactions into the model, a
non-parametric transformation of transcript expression is a worthwhile
consideration. transGI has compiled several methods and pairs of gene
interactions from previously published literature, allowing users to
quantitatively construct gene interaction networks de novo.

Nonparametric transformation method of transcriptome based on prior gene
interactions has the following potential advantages:

1.  It allows users to construct a gene interaction network for single
    samples, by comparing it with a reference network or selecting
    high-variance features from it, quantifying the perturbation effects
    of gene expression in specific tissues or central networks that
    promote biological events;
2.  Non-parametric transformation methods have good cross-dataset
    performance and can robustly remove batch effects between data sets.
    For microarrays and RNA-seq, due to differences in platforms,
    technologies, quantification methods, and sensitivities between the
    two, the effectiveness of nonparametric transformation methods is
    comparable to traditional methods such as quantile regression;
3.  Nonparametric transformation methods rely on fewer assumptions.

Currently, the functions supported by transGI include:

1.  Using transformation methods from previously publication like
    deltarank or pairwise calculation methods to construct single-sample
    networks on prior gene interaction pairs.
2.  Providing a gene enrichment method that uses prior gene interaction
    pairs as background gene sets. Simply put, it is a weighted gene
    enrichment method that uses the node degree of each prior gene
    interaction pair as weight.
3.  Providing a series of methods for selecting features with high
    variance or correlation with certain indicators, as well as
    visualizing gene interaction networks.

## Installation

You can install the development version of transGI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Vinnish-A/transGI")
```

## Example

## Future Plan
