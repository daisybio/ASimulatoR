require(rtracklayer)
require(polyester)
require(data.table)
source('functions.R')

### comments
#
# TODO: report nr_events (parameter), 
# coordinates in global/transcript, 
# psi values, 
# look at abundance,
# exon_usage
# 
# rccp to modify fasta to fastq if not fast enough
###


### parameters ----
# general parameters
seed <- 19L
nr_cores <- 20L
set.seed(seed)

# as parameters
gtf_path <- '/nfs/scratch/ensembl/Homo_sapiens.GRCh38.97.gtf'
nr_events <- list('ES' = 20L, 'MES' = 10L, 'IR' = 10L, 'A3' = 10L, 'A5' = 10L, 'ATSS' = 5L, 'ATTS' = 5L, 'MEE' = 10L, 'NO_AS' = 20L) 

# polyester parameters
outdir <- 'outdir'
seq_path <- '/nfs/scratch/ensembl/Homo_sapiens.GRCh38.97.fa'
err_rate = 0L

### importing gtf ----
message('importing gtf...')
gtf <- import(gtf_path)
message('finished importing gtf')
message('')


### creating premRNA and splice variants----
# TODO: make one data.table containing the eventannotation
# event_annotation <- data.table(
#   variant_1 = character(),
#   affected_genomic_region_1 = integer(),
#   affected_transcriptomic_region_1 = integer(),
#   variant_2 = character(),
#   affected_genomic_region_2 = integer(),
#   affected_transcriptomic_region_2 = integer(),
#   strand = character(),
#   event_annotation = character(),
#   inclusion_reads = integer(),
#   exclusion_reads = integer()
# )

# TODO: make one data.table containing the exon resolution, maybe add to gtf?
# exon_annotation <- data.table(
#   exon_id = character(),
#   genomic_coordinates = numeric(),
#   transcript_id = character(),
#   readcounts = integer()
# )

# determine how many genes to choose to generate splice variants from
nr_genes <- sum(unlist(nr_events))

# take only exons which are protein coding to generate hypothetical premRNAs
gtf_exons <- gtf[gtf$transcript_biotype == 'protein_coding' & gtf$type == 'exon']

# create the splice variants for every event
message('creating splice variants...')
#TODO: check if there is a method for parallel and side-effects
#TODO: foreach and no global operator, create blacklist with used and small genes
variants_event_annotation <- #mc
  sapply(names(nr_events), function(event_names){
  result <- list()
  event_names <- strsplit(event_names, ',', fixed = T)[[1]]
  count <- 0L
  
  # as-events need a certain length, otherwise it is not possible to construct this event
  cutoff <- ifelse(('MEE' %in% event_names | 'MES' %in% event_names), 3L, ifelse(('ES' %in% event_names), 2L, 1L))
  
  # if the hypothetical premRNA is not long enough, another gene is drawn randomly
  while(count < nr_events[[event_names]]){
    drawn_gene_id <- sample(unique(gtf_exons$gene_id), 1L)
    drawn_gene <- gtf_exons[gtf_exons$gene_id == drawn_gene_id]
    gtf_exons <<- gtf_exons[gtf_exons$gene_id != drawn_gene_id]
    negative_strand <- runValue(strand(drawn_gene) == '-')
    drawn_hypPreMRNA <- create_hyp_premRNA(drawn_gene, negative_strand)
    if(length(drawn_hypPreMRNA) > cutoff){
      count <- count + 1L
      result[[count]] <- construct_splice_variants(drawn_hypPreMRNA, event_names, negative_strand)
    }
  }
  result
})#, mc.cores = nr_cores)
variants_event_annotation <- unlist(unlist(variants_event_annotation, F), F)
variants_ind <- endsWith(names(variants_event_annotation), "variant")
message('finished creating splice variants')
message('')


### simulate ----
# combine all splice variants to create custom gtf for the read-simulation
variants_gtf <- Reduce(c, variants_event_annotation[variants_ind])
event_annotation <- rbindlist(variants_event_annotation[!variants_ind])
exon_resolution <- as.data.table(variants_gtf)
exon_resolution <- exon_resolution[,which(unlist(lapply(exon_resolution, function(x)!all(is.na(x))))),with=F]
# variants_gtf <- sortSeqlevels(variants_gtf)
# variants_gtf <- sort(variants_gtf)
# mcols(variants_gtf) <- data.frame(lapply(mcols(variants_gtf), as.character))
variants_gtf_path <- file.path(outdir, 'variants.gtf')
nr_transcripts <- length(unique(variants_gtf$transcript_id))
export(variants_gtf, variants_gtf_path)

# create a fold_change matrix
# TODO: use as parameters for polyester, especially for the abundance etc.
fold_change_mat <- matrix(c(rep(2, 5), rep(1, nr_transcripts - 5), 
                            rep(1, 5), rep(4, 5), rep(1, nr_transcripts - 10)), 
                          nrow = nr_transcripts)

message('starting simulation with polyester...')
# use polyester to simulate reads
simulate_experiment(gtf = variants_gtf_path, seqpath = seq_path,
                    fold_changes = fold_change_mat, outdir = outdir, 
                    error_rate = err_rate, reads_per_transcript = rep(300, nr_transcripts))
                    
message('done simulating reads with polyester')
message('')

### make statistics ----
# create fastq files from fasta files
# TODO: do more while going through fasta files
out_files <- list.files(outdir)
fasta_files <- file.path(outdir, out_files[endsWith(out_files, '.fasta')])
message('creating fastq-files and statistics...')
event_annotation_and_exon_resolution <- unlist(mclapply(fasta_files, check_reads_and_convert_format, event_annotation, exon_resolution, mc.cores = nr_cores), recursive = F)
event_annotation_and_exon_resolution <- split(event_annotation_and_exon_resolution, names(event_annotation_and_exon_resolution))
event_annotation <- Reduce(merge, event_annotation_and_exon_resolution$event_annotation)
exon_resolution <- Reduce(merge, event_annotation_and_exon_resolution$exons)
exon_resolution <- exon_resolution[order(transcript_id)]
#tmp <- check_reads_and_convert_format(fasta_files[1], event_annotation, exon_dt)
#TODO: merge multiple data tables
fwrite(x = event_annotation, file = file.path(outdir, "event_annotation.csv"), quote = F)
fwrite(x = exon_resolution, file = file.path(outdir, "exon_resolution.csv"), quote = F)
message('done creating fastq-files and statistics')
message('')
