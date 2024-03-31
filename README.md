
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

Use the `transGI` function to transform transcriptome data.

``` r
library(transGI)

testMat = system.file('extdata', 'inputMatTest.csv', package = 'transGI') |>
  read.csv(row.names = 'symbol') |>
  as.matrix() |> 
  _[1:2000, ]
res = transGI(testMat, 'deltarank', 'reactome', nThreads_ = 4)
```

Use `enrichGraph` to calculate the enrichment results of genes in the
background gene interaction network.

``` r
symbols = readLines('material/symbols.csv')

db = read_gmt('material/h.all.v2023.2.Hs.symbols.gmt')
res_enrich = enrichGraph(symbols, db_ = db, pathways_ = 'all', bgNet_ = 'reactome')

autoBar(res_enrich)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" style="display: block; margin: auto;" />

## Future Plan

Currently, transGI has implemented a series of functions including
feature engineering, feature selection, gene enrichment, and simple
visualization. In the future, we plan to add the following features to
transGI:

1.  Infer new potential gene pairs for interaction based on perturbation
    effects in experimental and control groups
2.  For networks of interactions between two or more genes, add
    functions such as network comparison and key node selection
3.  Add more methods for constructing single-sample interaction networks
    and calculating perturbation effects.

## Citation

The calculation method for deltarank is obtained from (Chen et al.
2020), and the pairwise calculation method is obtained from (Zhu et al.
2022)

If you use transGI in your research, please cite:

    @Manual{,
      title = {transGI: Nonparametric transformation method of transcriptome based on prior gene interactions},
      author = {Zhihao Xu},
      year = {2024},
      note = {R package version 0.1.0},
      url = {https://github.com/Vinnish-A/transGI},
    }

## reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-10.1093/bib/bbaa268" class="csl-entry">

Chen, Yuanyuan, Yu Gu, Zixi Hu, and Xiao Sun. 2020.
“<span class="nocase">Sample-specific perturbation of gene interactions
identifies breast cancer subtypes</span>.” *Briefings in Bioinformatics*
22 (4): bbaa268. <https://doi.org/10.1093/bib/bbaa268>.

</div>

<div id="ref-10.1093/bib/bbac344" class="csl-entry">

Zhu, Sujie, Weikaixin Kong, Jie Zhu, Liting Huang, Shixin Wang, Suzhen
Bi, and Zhengwei Xie. 2022. “<span class="nocase">The genetic
algorithm-aided three-stage ensemble learning method identified a robust
survival risk score in patients with glioma</span>.” *Briefings in
Bioinformatics* 23 (5): bbac344. <https://doi.org/10.1093/bib/bbac344>.

</div>

</div>
