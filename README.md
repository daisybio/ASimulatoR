
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ass

<!-- badges: start -->

<!-- badges: end -->

The goal of ass is to simulate RNA-seq reads with alternative splicing
events.

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

Please note that we use a custon Version of Polyester, that is available
at <https://github.com/quirinmanz/polyester>

## Example

### Creating exon supersets

``` r
suppressMessages(library(ass))

# create exon superset for genes on chromosome 21 of ensembl release 99
gtf_file = system.file('data', 'Homo_sapiens.GRCh38.99.21.gtf', package = 'ass')
# by default the produed superset will be saved as .rda file into the same directory
exon_superset = get_exon_supersets(gtf_file)
#> importing gtf...
#> finished importing gtf
#> 
#> creating superset...
#> finished creating superset
#> 
#> saving superset...
#> finished saving superset
#> 
exon_superset[[1]]
#> GRanges object with 32 ranges and 10 metadata columns:
#>        seqnames            ranges strand |       source        type       score
#>           <Rle>         <IRanges>  <Rle> |  <character> <character> <character>
#>    [1]       21 41879270-41879482      - | as_simulator        exon           .
#>    [2]       21 41878996-41879140      - | as_simulator        exon           .
#>    [3]       21 41878695-41878860      - | as_simulator        exon           .
#>    [4]       21 41871494-41871621      - | as_simulator        exon           .
#>    [5]       21 41867300-41867377      - | as_simulator        exon           .
#>    ...      ...               ...    ... .          ...         ...         ...
#>   [28]       21 41815705-41815836      - | as_simulator        exon           .
#>   [29]       21 41810154-41814962      - | as_simulator        exon           .
#>   [30]       21 41804534-41804723      - | as_simulator        exon           .
#>   [31]       21 41802712-41803118      - | as_simulator        exon           .
#>   [32]       21 41798225-41801722      - | as_simulator        exon           .
#>              phase         gene_id            transcript_id  template
#>        <character>     <character>              <character> <logical>
#>    [1]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>    [2]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>    [3]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>    [4]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>    [5]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>    ...         ...             ...                      ...       ...
#>   [28]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>   [29]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>   [30]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>   [31]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>   [32]           . ENSG00000141956 ENSG00000141956_template      TRUE
#>        gene_exon_number  tr_start    tr_end
#>               <integer> <integer> <integer>
#>    [1]                1         1       213
#>    [2]                2       214       358
#>    [3]                3       359       524
#>    [4]                4       525       652
#>    [5]                5       653       730
#>    ...              ...       ...       ...
#>   [28]               28      3969      4100
#>   [29]               29      4101      8909
#>   [30]               30      8910      9099
#>   [31]               31      9100      9506
#>   [32]               32      9507     13004
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### Simulating Alternative Splicing

``` r
# now we create splice variants from 9 exon supersets
max_genes = 9

# to create one variant per gene we set prob_as_freq to TRUE and then produce one of each AS events as well as one variant containing every event
probs_as_freq = T
event_freq = 
    setNames(rep(1/9, 9),
             c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee', 'es,ir,mes,a3,a5,afe,ale,mee'))

# we will create a small experiment with just one samplesper group
num_reps = c(1,1)

# we can use the previously created superset to simulate splice variants from, since it is saved in the same directory as the gtf
# if no superset is found, a new one will be created
simulate_alternative_splicing(input_dir = system.file('data', package = 'ass'),
                              event_probs = event_freq, 
                              outdir = 'simulation', 
                              probs_as_freq = probs_as_freq, 
                              max_genes = max_genes,
                              num_reps = num_reps)
#> found the following fasta files: 21.fa
#> note that splice variants will only be constructed from chromosomes that have a corresponding fasta file
#> 
#> loading superset...
#> finished loading superset
#> 
#> create splicing variants and annotation...
#> finished creating splicing variants and annotation
#> 
#> exporting gtf for read simulation...
#> finished exporting gtf
#> 
#> exporting event_annotation...
#> finished exporting event_annotation...
#> 
#> start simulation with polyester:
#> parsing gtf and sequences...
#> done parsing
#> start sequencing... (1m reads per iteration)
#> sample_01: overall 38753 reads
#> sample_01: iteration 01
#> sample_01: fragments generated
#> sample_01: write read pairs
#> sample_02: overall 41504 reads
#> sample_02: iteration 01
#> sample_02: fragments generated
#> sample_02: write read pairs
#> finished sequencing
```

### Check AS Event Annotation

``` r
event_anno = read.csv('simulation/event_annotation.tsv', sep = '\t')
event_anno[event_anno$event_annotation == 'es',]
#>    event_annotation                                     variant
#> 1                es                          ENSG00000183844_es
#> 19               es ENSG00000273840_es,ir,mes,a3,a5,afe,ale,mee
#>                    template genomic_start genomic_end transcriptomic_start
#> 1  ENSG00000183844_template      41347013    41347100                 2988
#> 19 ENSG00000273840_template      10491540    10491644                  225
#>    transcriptomic_end
#> 1                3075
#> 19                329
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- #summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

### Visualize Splice Variants

``` r
# to visualitze the splice variants we will use ggbio
suppressMessages(library(ggbio))
# firstly, we load the newly created gtf file 
gtf = rtracklayer::import('simulation/splicing_variants.gtf')
# the gene id of the variant with all events
gene_id = gtf$gene_id[grep('es,ir,mes,a3,a5,afe,ale,mee', gtf$transcript_id, fixed = T)[1]]
exons = gtf[gtf$type == 'exon' & gtf$gene_id == gene_id]
suppressWarnings(ggbio::autoplot(split(exons, exons$transcript_id)))
#> Constructing graphics...
```

<img src="man/figures/README-visualization-1.png" width="100%" />

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
