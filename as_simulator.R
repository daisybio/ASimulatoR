require(rtracklayer)
require(Biostrings)
require(polyester)
require(data.table)
source('functions.R')

### comments
# rccp to modify fasta to fastq if not fast enough, but first see if it works with Biostrings
###


### parameters ----
# general parameters
seed <- 19L
nr_cores <- 30L
set.seed(seed)

# simulator parameters
gtf_path <- 'ensembl/Homo_sapiens.GRCh38.97.gtf'
nr_events <- setNames(c(1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9), 
                      c('ES', 'MES', 'IR', 'A3', 'A5', 'AFE', 'ALE', 'MEE', 'NO_AS'))
min_exon_number <- 4L
max_exon_number <- 30L
nr_genes <- sum(nr_events)
events_as_fractions <- nr_genes == 1

# polyester parameters
seq_path <- 'ensembl/Homo_sapiens.GRCh38.97.fa'
valid_chromosomes <- sub('.fa', '', list.files('ensembl/Homo_sapiens.GRCh38.97.fa'))
outdir <- paste('outdir', min_exon_number, max_exon_number, sep = '_')
dir.create(outdir, F)
num_reps <- c(10,10)
err_rate = 0L
read_length <- 150L
strand_specific <- T
paired <- T
reportCoverage <- F


### importing/creating exon superset gtf ----
if (file.exists(paste0(gtf_path, '.exon_superset.rda'))) {
  message('loading superset...')
  load(paste0(gtf_path, '.exon_superset.rda'))
  message('finished loading superset')
} else {
  exon_supersets <- create_exon_supersets(gtf_path, nr_cores)
}

### create the splice variants for every event ----
message('creating splice variants...')

# kick out small genes, big genes and genes with small exons
message(paste('filtering out genes that are not n the chromosomes', paste(valid_chromosomes, collapse = ','), 'or with less exons than', min_exon_number, 'or more than', max_exon_number))
#gene_lengths <- sapply(exon_supersets, length)
filter_out <- sapply(exon_supersets, function(gene){
  return(length(gene) < min_exon_number | length(gene) > max_exon_number | any(width(gene) <= 3 | !runValue(seqnames(gene)) %in% valid_chromosomes))
})

exon_supersets <- sample(exon_supersets[!filter_out])

#gene_lengths <- gene_lengths[gene_lengths >= min_exon_number & gene_lengths <= max_exon_number]
if (events_as_fractions) {
  nr_events <- nr_events * length(exon_supersets)
} else if (length(exon_supersets) < nr_genes) {
  stop('not enough genes left to construct the nr of indicated events')
} 
nr_events <- floor(nr_events)
nr_genes <- sum(nr_events)
if (nr_genes < length(exon_supersets)) {
  exon_supersets <- exon_supersets[1:nr_genes]
}

# create splicing variants
# event_names <- strsplit(names(nr_events), ',', fixed = T)

#message('creating splice variants...')
exon_supersets <- split(exon_supersets, rep(names(nr_events), nr_events))
#TODO: something went wrong here!

variants <- mclapply(names(exon_supersets), 
                                      function(event_name) 
                                        mclapply(exon_supersets[[event_name]], 
                                                 construct_splice_variants, 
                                                 strsplit(event_name, ',', fixed = T)[[1]]
                                                 , mc.cores = max(floor((nr_cores - length(exon_supersets))/length(exon_supersets)), 1), mc.allow.recursive = T)
                                      , mc.cores = min(c(nr_cores, length(exon_supersets))), mc.allow.recursive = T)
rm(exon_supersets)
variants <- unlist(unlist(variants, F), F)
variants_ind <- endsWith(names(variants), 'variant')
message('finished creating splice variants')
message('')

# min_lengths <- sort(setNames(sapply(event_names, function(event_names) 
#   ifelse(('MEE' %in% event_names | 'MES' %in% event_names), 4L, 
#          ifelse(('ES' %in% event_names | 'AFE' %in% event_names | 'ALE' %in% event_names), 3L, 
#                ifelse('IR' %in% event_names | 'A3' %in% event_names | 'A5' %in% event_names, 2L, 1L)))), names(nr_events)), T)



# mee and mes need minimum length 4
# es and afe and ale need minimum length 3
# a3 and a5 and ir need length 2
# no as need length 1
# variants_event_annotation <- #mc
#   sapply(names(nr_events), function(event_names){
#     result <- list()
#     event_names <- strsplit(event_names, ',', fixed = T)[[1]]
#     count <- 0L
#     
#     # as-events need a certain length, otherwise it is not possible to construct this event
#     cutoff <- ifelse(('MEE' %in% event_names | 'MES' %in% event_names), 3L, ifelse(('ES' %in% event_names | 'AFE' %in% event_names | 'ALE' %in% event_names), 2L, 1L))
#     
#     # if the hypothetical premRNA is not long enough, another gene is drawn randomly
#     while (count < nr_events[[event_names]]) {
#       drawn_gene_id <- draw_one_sample_safely(unique(gtf_exons$gene_id))
#       drawn_gene <- gtf_exons[gtf_exons$gene_id == drawn_gene_id]
#       gtf_exons <<- gtf_exons[gtf_exons$gene_id != drawn_gene_id]
#       negative_strand <- runValue(strand(drawn_gene) == '-')
#       drawn_hypPreMRNA <- create_hyp_premRNA(drawn_gene, negative_strand)
#       if (length(drawn_hypPreMRNA) > cutoff) {
#         count <- count + 1L
#         result[[count]] <- construct_splice_variants(drawn_hypPreMRNA, event_names, negative_strand)
#       }
#     }
#     result
#   })#, mc.cores = nr_cores)



### simulate ----
# combine all splice variants to create custom gtf for the read-simulation
event_annotation <- rbindlist(variants[!variants_ind])
variants <- rbindlist(variants[variants_ind])
variants_gtf_path <- file.path(outdir, 'variants.gtf')
message('exporting gtf for read simulation...')
export(variants[variants$type != 'junction'], variants_gtf_path)
message('finished exporting gtf')

#variants_gtf <- variants[variants$type != 'junction']

#exon_junction_resolution <- as.data.table(variants)
#exon_junction_resolution <- exon_junction_resolution[, which(unlist(lapply(exon_junction_resolution, function(x)!all(is.na(x))))), with = F]
# variants_gtf <- sortSeqlevels(variants_gtf)
# variants_gtf <- sort(variants_gtf)
# mcols(variants_gtf) <- data.frame(lapply(mcols(variants_gtf), as.character))
nr_transcripts <- length(unique(variants$transcript_id))

# create a fold_change matrix
# TODO: use as parameters for polyester, especially for the abundance etc.
fold_change_mat <- matrix(c(rep(2, 5), rep(1, nr_transcripts - 5), 
                            rep(1, 5), rep(4, 5), rep(1, nr_transcripts - 10)), 
                          nrow = nr_transcripts)

message('starting simulation with polyester...')
# use polyester to simulate reads
simulate_experiment(gtf = variants_gtf_path, 
                    seqpath = seq_path,
                    outdir = outdir, 
                    num_reps = num_reps,
                    fold_changes = fold_change_mat, 
                    reads_per_transcript = rep(300L, nr_transcripts), 
                    paired = paired,
                    reportCoverage = reportCoverage,
                    readlen = read_length,
                    error_rate = err_rate, 
                    strand_specific = strand_specific,
                    seed = seed,
                    gzip = T)

message('done simulating reads with polyester')
message('')

### make statistics ----
# create fastq files from fasta files
out_files <- list.files(outdir)
fasta_files <- file.path(outdir, out_files[grepl('.*?\\.fasta(.gz)?$', out_files)])
message('creating fastq-files and statistics...')
#setkey(exon_junction_resolution, exon_id)
#tmp <- check_reads_and_convert_format(fasta_files[1], exon_junction_resolution, read_length)

exon_junction_counts <- 
  mclapply(fasta_files, check_reads_and_convert_format, exon_junction_resolution, read_length, mc.cores = nr_cores)
exon_junction_counts <- exon_junction_counts[!sapply(exon_junction_counts, is.null)]
exon_junction_counts <- Reduce(function(x, y) merge(x, y, by = 'exon_id', all = T), exon_junction_counts)
sample_names <- names(exon_junction_counts)['exon_id' != names(exon_junction_counts)]
exon_junction_resolution[exon_junction_counts, on = 'exon_id', get('sample_names') := mget(paste0('i.', sample_names))]
#event_annotation_and_exon_junction_resolution <- 
#  split(event_annotation_and_exon_junction_resolution, names(event_annotation_and_exon_junction_resolution))
#event_annotation <- Reduce(merge, event_annotation_and_exon_junction_resolution$event_annotation)
#exon_junction_resolution_merged <- Reduce(merge, exon_junction_resolution)
##Reduce(merge, event_annotation_and_exon_junction_resolution$exons)
#exon_junction_resolution <- exon_junction_resolution[order(seqnames, transcript_id)]

# write annotations to files
fwrite(x = event_annotation, file = file.path(outdir, 'event_annotation.csv'), quote = F)
fwrite(x = exon_junction_resolution, file = file.path(outdir, 'exon_junction_resolution.csv'), quote = F)
message('done creating fastq-files and statistics')
message('')
