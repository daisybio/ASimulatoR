
# ASimulatoR

The goal of ASimulatoR is to simulate RNA-seq reads with alternative
splicing events. The alternative splicing events are well documented and
the true origin of each read is used for exon and juntion coverage via a
modified version of the bioconductor polyester package.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("biomedbigdata/ASimulatoR")
```

Please note that we use a custom version of Polyester, that is available
at <https://github.com/biomedbigdata/polyester>

## Example

### Creating exon supersets

Firstly, we create exon supersets by joining all exons of a gene from a
gtf file. These supersets are then used to create splice variants. Since
all exons from one gene are used to create the exon superset, you may
find that the term exon superset is used analogously to gene.

``` r
suppressMessages(library(ASimulatoR))

# create exon superset for genes on chromosome 21 of ensembl release 99
gtf_file = system.file('extdata', 'Homo_sapiens.GRCh38.99.21.gtf', package = 'ASimulatoR')

# by default the produced superset will be saved as .rda file into the same directory
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
exon_superset[[1]][1:5, ]
#> GRanges object with 5 ranges and 10 metadata columns:
#>       seqnames            ranges strand |      source        type       score       phase         gene_id            transcript_id  template gene_exon_number  tr_start    tr_end
#>          <Rle>         <IRanges>  <Rle> | <character> <character> <character> <character>     <character>              <character> <logical>        <integer> <integer> <integer>
#>   [1]       21 41879270-41879482      - |  ASimulatoR        exon           .           . ENSG00000141956 ENSG00000141956_template      TRUE                1         1       213
#>   [2]       21 41878996-41879140      - |  ASimulatoR        exon           .           . ENSG00000141956 ENSG00000141956_template      TRUE                2       214       358
#>   [3]       21 41878695-41878860      - |  ASimulatoR        exon           .           . ENSG00000141956 ENSG00000141956_template      TRUE                3       359       524
#>   [4]       21 41871494-41871621      - |  ASimulatoR        exon           .           . ENSG00000141956 ENSG00000141956_template      TRUE                4       525       652
#>   [5]       21 41867300-41867377      - |  ASimulatoR        exon           .           . ENSG00000141956 ENSG00000141956_template      TRUE                5       653       730
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### Simulating Alternative Splicing

You can find more information about the main function of this package at
the end of the page.

This simulator supports eight different AS events:

| es           | mes                    | ir               | a3                                  | a5                               | mee                      | afe                    | ale                   |
| ------------ | ---------------------- | ---------------- | ----------------------------------- | -------------------------------- | ------------------------ | ---------------------- | --------------------- |
| exon skiping | multiple exon skipping | intron retention | alternative 3’/acceptor splice site | alternative 5’/donor splice site | mutually exclusive exons | alternative first exon | alternative last exon |

``` r
# define your input_dir, where the annotation gtf (or the exon supersets if you have already created them) and the genome fasta files are located
# here we will use the example data
input_dir = system.file('extdata', package = 'ASimulatoR')

# define, how many groups and samples per group you analyze. Here we create a small experiment with two groups with one sample per group:
num_reps = c(1,1)

# define your outdir with NO slash
outdir = 'simulation'

# define the number of genes you want to work with. If you want all exons, do not specify this parameter or set it to NULL
# here we create splice variants from 9 exon supersets:
max_genes = 9
```

You could define the distribution of the events by probability or
relative frequency.

  - Probability: For each superset we create an event with the
    probability mentioned in `event_prob`.
  - Frequency: Set `probs_as_freq = T`. The exon supersets are
    partitioned corresponding to the `event_prob` parameter.

<!-- end list -->

``` r
# in this example we use relative frequencies
# here we produce eight variants with one of each AS events as well as one variant containing every event
# if probs_as_freq was FALSE, a random number would be drawn for each event-superset combination and only if it was smaller than 1/9 the AS event would be created
probs_as_freq = T
event_freq = 
    setNames(rep(1/9, 9),
             c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee', 'es,ir,mes,a3,a5,afe,ale,mee'))



# we use the previously created superset to simulate splice variants from, since it is saved in the same directory as the gtf
# if no superset is found, a new one will be created
simulate_alternative_splicing(input_dir = input_dir,
                              event_probs = event_freq,
                              outdir = outdir, 
                              probs_as_freq = probs_as_freq, 
                              max_genes = max_genes,
                              num_reps = num_reps)
#> found the following fasta files: 21.fa
#> note that splice variants will only be constructed from chromosomes that have a corresponding fasta file
#> 
#> loading superset...
#> finished loading superset
#> 
#> assign variants to supersets...
#> create splicing variants and annotation. This may take a while...
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
#> sample_01: overall 91608 reads
#> sample_01: iteration 01
#> sample_01: fragments generated
#> sample_01: write read pairs
#> sample_02: overall 78660 reads
#> sample_02: iteration 01
#> sample_02: fragments generated
#> sample_02: write read pairs
#> finished sequencing
```

### Visualize Splice Variants

``` r
# to visualize the splice variants we will use ggbio
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

``` r

# have a look at the event annotation
event_anno = read.csv('simulation/event_annotation.tsv', sep = '\t')
event_anno[grepl(gene_id, event_anno$template) | grepl(gene_id, event_anno$variant), ]
#>    event_annotation                                     variant                                    template                                                                    genomic_start                                                                      genomic_end                         transcriptomic_start                           transcriptomic_end
#> 1               afe ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18460572                                                                         18460731                                            1                                          160
#> 2               afe                    ENSG00000154646_template ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                                                                         18485799                                                                         18485879                                            1                                           81
#> 3               ale ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18269116                                                                         18270124                                         6950                                         7958
#> 4               ale                    ENSG00000154646_template ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                                                                         18275197                                                                         18275336                                         3018                                         3157
#> 5               mee ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18379283                                                                         18379318                                          792                                          827
#> 6               mee                    ENSG00000154646_template ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                                                                         18372193                                                                         18372324                                          786                                          917
#> 7               mes ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template 18332084,18329169,18326432,18315146,18312945,18297734,18294603,18294270,18281040 18332173,18329294,18326572,18315256,18313077,18297829,18294652,18294444,18281221 1818,1908,2034,2175,2286,2419,2515,2565,2740 1907,2033,2174,2285,2418,2514,2564,2739,2921
#> 8                es ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18352903                                                                         18353052                                         1275                                         1424
#> 9                ir ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18359864                                                                         18365139                                         1223                                         6498
#> 10               a3 ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18383728                                                                         18383778                                          589                                          639
#> 11               a5 ENSG00000154646_es,ir,mes,a3,a5,afe,ale,mee                    ENSG00000154646_template                                                                         18398199                                                                         18398220                                          499                                          520
```

## `simulate_alternative_splicing`: Simulate RNA-seq experiment with splicing variants

### Description

Firstly, exon supersets are created by joining all exons of a gene from
a gtf file. Next, splicing variants are created with documentation and
event annotation based on the users input. Finally, fastq files
containing RNA-seq reads from the splice variants and the real exon and
junction coverage are created using a modified version of the polyester
R package available on <https://github.com/quirinmanz/polyester>.

### Usage

``` r
simulate_alternative_splicing(input_dir, event_probs, outdir, ncores = 1L, ...)
```

### Arguments

| Argument      | Description                                                                                                                                                                                                                                                                                     |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input_dir`   | Character path to directory containing the gtf file from which splice variants are created and genome fasta files with one file per chromosome i.e. <chr_name>.fa passed to polyester                                                                                                           |
| `event_probs` | Named list/vector containing numerics corresponding to the probabilites to create the event (combination). If `probs_as_freq` is `TRUE` `event_probs` correspond to the relative frequency of occurences for the event(combination) and in this case the sum of all frequencies has to be \<=1. |
| `outdir`      | character, path to folder where simulated reads and all annotations should be written, with *no* slash at the end. By default, reads are written to current working directory.                                                                                                                  |
| `ncores`      | the number of cores to be utilized for parallel generation of splice variant creation and read simulation.                                                                                                                                                                                      |
| `...`         | any of several other arguments that can be used to add nuance to the simulation and splice variant creation. See details.                                                                                                                                                                       |

### Details

Reads are simulated from a GTF file which is produced by
`create_splicing_variants_and_annotation` plus DNA sequences.

Several optional parameters can be passed to this function to adjust the
simulation. For polyester parameters refer to `simulate_experiment` from
the polyester R package:

  - `write_gff` : Additionally to the gtf file containing the splice
    variants, a gff3 file with the same content will be printed to the
    outdir. Default `TRUE`

  - `max_genes` : The maximum number of genes/exon supersets to be
    included in the process of splice variant creation. Default `NULL`
    which means that all available exon supersets will be used.

  - `exon_junction_coverage` : Should the real coverage of exons,
    junctions and retained introns be written into a additional file.
    Default `TRUE`

  - `multi_events_per_exon` : Should it be possible to have more than
    one AS event at the same exon if multiple variants are created for
    the same exon superset? \!If this option is set to `TRUE` , there
    may occur unforeseen AS events that are not documented in the
    event\_annotation file\!. Default `FALSE`

  - `probs_as_freq` : Should `event_probs` be treated as relative
    frequencies instead of probabilities? Default `FALSE`

  - `save_exon_superset` : Should the exon supersets be saved to .rda
    file? Default `TRUE`

Parameters passed to polyester that we assigned different defaults to
than in `simulate_experiment`:

  - `fold_changes` : Currently, ASimulatoR introduces random isoform
    switches. Those can be retraced in the sim\_tx\_info.txt file
    written by polyester. We plan on improving this in the future.

  - `strand_specific` : Strand-specific simulation (1st read forward
    strand, 2nd read reverse strand with respect to transcript
    sequence). Default `TRUE` .

  - `meanmodel` : `reads_per_transcripts` as a function of transcript
    length. Default `TRUE` .

  - `verbose` : Should progress messages be printed during the
    sequencing process? Default `TRUE` .

  - `exon_junction_coverage` : Should the coverage of exons, junctions
    and retained introns be determined? Default `TRUE` .

  - `exon_junction_table` : If `exon_junction_coverage=TRUE` a
    `data.table` produced by `create_splicing_variants_and_annotation`
    to determine exon and intron coverage.

### Value

No return, but simulated reads, a simulation info file, an alternative
splicing event annotation and exon and junction coverages are written to
`outdir` .

### References

Alyssa C. Frazee, Andrew E. Jaffe, Ben Langmead, Jeffrey T. Leek,
Polyester: simulating RNA-seq datasets with differential transcript
expression, Bioinformatics, Volume 31, Issue 17, 1 September 2015, Pages
2778–2784, <https://doi.org/10.1093/bioinformatics/btv272>

### License

ASimulatoR Copyright (C) 2020 Manz, Quirin

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <https://www.gnu.org/licenses/>.
