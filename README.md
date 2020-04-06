
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ass

<!-- badges: start -->

<!-- badges: end -->

The goal of ass is to
…

## Installation

<!-- You can install the released version of ass from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("ass") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
Sys.setenv(GITHUB_PAT = "16946cf3c89e74cd07363d0f660559e4339d730b")
remotes::install_github("biomedbigdata/as_simulator")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ass, quietly = T)
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which, which.max, which.min
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:data.table':
#> 
#>     first, second
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:data.table':
#> 
#>     shift
## basic example code
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- #summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

``` r
library(ggbio, quietly = T)
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> Need specific help about ggbio? try mailing 
#>  the maintainer or visit http://tengfei.github.com/ggbio/
#> 
#> Attaching package: 'ggbio'
#> The following objects are masked from 'package:ggplot2':
#> 
#>     geom_bar, geom_rect, geom_segment, ggsave, stat_bin, stat_identity,
#>     xlim
gff = rtracklayer::import('data/splicing_variants_ENSG00000122257.gff3')
exons = gff[gff$type == 'exon']
autoplot(split(exons, exons$transcript_id))
#> Constructing graphics...
#> Warning: `quo_expr()` is deprecated as of rlang 0.2.0.
#> Please use `quo_squash()` instead.
#> This warning is displayed once per session.
```

<img src="man/figures/README-visualization-1.png" width="100%" />

``` r
read.csv('data/event_annotation_ENSG00000122257.tsv')
#>    event_annotation                                     variant
#> 1                es                          ENSG00000122257_es
#> 2                a3                          ENSG00000122257_a3
#> 3               mee                         ENSG00000122257_mee
#> 4               mee                    ENSG00000122257_template
#> 5               afe ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 6               afe                    ENSG00000122257_template
#> 7               ale ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 8               ale                    ENSG00000122257_template
#> 9               mee ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 10              mee                    ENSG00000122257_template
#> 11              mes ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 12               es ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 13               ir ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 14               a3 ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#> 15               a5 ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5
#>                                       template     genomic_start
#> 1                     ENSG00000122257_template          24559505
#> 2                     ENSG00000122257_template          24567792
#> 3                     ENSG00000122257_template          24555821
#> 4                          ENSG00000122257_mee          24556308
#> 5                     ENSG00000122257_template          24537693
#> 6  ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5          24539556
#> 7                     ENSG00000122257_template          24570876
#> 8  ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5          24568745
#> 9                     ENSG00000122257_template          24555821
#> 10 ENSG00000122257_ale,afe,mee,ir,mes,a3,es,a5          24556308
#> 11                    ENSG00000122257_template 24561824,24563199
#> 12                    ENSG00000122257_template          24559505
#> 13                    ENSG00000122257_template          24549402
#> 14                    ENSG00000122257_template          24567792
#> 15                    ENSG00000122257_template          24563645
#>          genomic_end transcriptomic_start transcriptomic_end
#> 1           24559805                 4760               5060
#> 2           24567810                 6293               6311
#> 3           24555917                 4506               4602
#> 4           24556447                 4506               4645
#> 5           24537933                    1                241
#> 6           24540792                    1               1237
#> 7           24572863                 7966               9953
#> 8           24570499                 6395               8149
#> 9           24555917                 7101               7197
#> 10          24556447                 4506               4645
#> 11 24562161,24563501            5165,5503          5502,5805
#> 12          24559805                 4760               5060
#> 13          24552992                  799               4389
#> 14          24567810                 6293               6311
#> 15          24563664                 5806               5825
```

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
