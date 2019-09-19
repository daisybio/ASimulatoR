library(rtracklayer)
library(polyester)
library(data.table)
source('functions.R')

### comments
#
# TODO: report nr_events (parameter), 
# coordinates in global/transcript, 
# psi values, 
# look at abundance,
# exon_usage
# 
###


### parameters ----
# general parameters
seed <- 19L
nr_cores <- 12L
set.seed(seed)

# as parameters
gtf_path <- 'ensembl/Homo_sapiens.GRCh38.97.gtf'
nr_events <- list('ES' = 20L, 'MES' = 10L, 'IR' = 10L, 'A3' = 10L, 'A5' = 10L, 'ATSS' = 5L, 'ATTS' = 5L, 'MEE' = 10L, 'NO_AS' = 20L) 

# polyester parameters
outdir <- 'outdir'
seq_path <- 'ensembl/Homo_sapiens.GRCh38.97.fa'
err_rate = 0L

### importing gtf ----
print('importing gtf...')
gtf <- import(gtf_path)
print('finished importing gtf')


### creating premRNA and splice variants----
# TODO: make one data.table containing the eventannotation
# event_annotation <- data.table(
#   variant_1 = character(),
#   affected_genomic_region_1 = numeric(),
#   variant_2 = character(),
#   affected_genomic_region_2 = numeric(),
#   event_annotation = character(),
#   inclusion_reads = integer(),
#   exclusion_reads = integer()
# )
#
# TODO: make one data.table containing the exon resolution
# exon_annotation <- data.table(
#   genomic_coordinates = numeric(),
#   transcript_id = character(),
#   readcounts = integer()
# )

# determine how many genes to choose to generate splice variants from
nr_genes <- sum(unlist(nr_events))

# take only exons which are protein coding to generate hypothetical premRNAs
gtf_exons <- gtf[gtf$transcript_biotype == 'protein_coding' & gtf$type == 'exon']

# create the splice variants for every event
print('creating splice variants...')
#TODO: check if there is a method for parallel and side-effects
splicing_variants <- #mc
  lapply(names(nr_events), function(event_names){
  result <- GRanges()
  event_names <- strsplit(event_names, ',', fixed = T)[[1]]
  count <- 0L
  
  # certain as-events need a certain length, otherwise it is not possible to construct this event
  cutoff <- ifelse(('MEE' %in% event_names| 'MES' %in% event_names), 3L, ifelse(('ES' %in% event_names), 2L, 1L))
  
  # if the hypothetical premRNA is not long enough, another gene is drawn randomly
  while(count < nr_events[[event_names]]){
    drawn_gene_id <- sample(unique(gtf_exons$gene_id), 1L)
    drawn_gene <- gtf_exons[gtf_exons$gene_id == drawn_gene_id]
    gtf_exons <<- gtf_exons[gtf_exons$gene_id != drawn_gene_id]
    drawn_hypPreMRNA <- create_hyp_premRNA(drawn_gene)
    if(length(drawn_hypPreMRNA) > cutoff){
      count <- count + 1L
      result <- c(result, construct_splice_variants(drawn_hypPreMRNA, event_names))
    }
  }
  result
})#, mc.cores = nr_cores)
names(splicing_variants) <- names(nr_events)
print('finished creating splice variants')


### simulate ----
# combine all splice variants to create custom gtf for the read-simulation
variants_gtf <- Reduce(c, splicing_variants)
variants_gtf_path <- file.path(outdir, 'variants.gtf')
nr_transcripts <- length(unique(variants_gtf$transcript_id))
export(variants_gtf, variants_gtf_path)

# create a fold_change matrix
# TODO: use as parameters for polyester, especially for the abundance etc.
fold_change_mat <- matrix(c(rep(2, 5), rep(1, nr_transcripts - 5), 
                            rep(1, 5), rep(4, 5), rep(1, nr_transcripts - 10)), 
                          nrow = nr_transcripts)

# use polyester to simulate reads
simulate_experiment(gtf = variants_gtf_path, seqpath = seq_path,
                    fold_changes = fold_change_mat, outdir = outdir, 
                    error_rate = err_rate, reads_per_transcript = rep(300, nr_transcripts))


out_files <- list.files(outdir)
fasta_files <- file.path(outdir, out_files[endsWith(out_files, '.fasta')])

# create fastq files from fasta files
# TODO: do more while going through
lapply(fasta_files, fastqFromFasta)


### Notes ----
# TODO: 1. combine all exons to the hypothetical pre-mRNA using the genomic ranges package!
# TODO: 2. introduce as-events from this hypothetical pre-mRNA
# TODO: 2.1. keep and express hypothetical pre-mRNA
# TODO: 3. write gtf to file
# TODO: 4. give gtf and fasta files to polyester
# TODO: 5. use rccp to modify fasta to fastq if not fast enough