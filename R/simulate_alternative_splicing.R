check_parameters <- function(extras) {
  if (is.null(extras$exon_junction_coverage)) {
    extras$exon_junction_coverage <- T
  }
  if (is.null(extras$error_rate)) {
    extras$error_rate <- 0
  }
  if (is.null(extras$strand_specific)) {
    extras$strand_specific <- T
  }
  if (is.null(extras$gzip)) {
    extras$gzip = F
  }
  if (is.null(extras$shuffle)) {
    extras$shuffle = T
  }
  if (is.null(extras$fastq)) {
    extras$fastq = T
  }
  if (is.null(extras$readlen)) {
    extras$readlen <- 150
  }
  if (is.null(extras$ncores)) {
    extras$ncores <- 1
  }

  return(extras)
}


simulate_alternative_splicing <-
  function(gtf_path,
           event_probs,
           fold_changes = NULL,
           seqpath,
           outdir,
           ...)
  {
    # require(rtracklayer)
    # require(data.table)
    # require(polyester)
    # require(foreach)
    extras <- check_parameters(list(...))

    # Store the current random number generator to restore at the end
    # Changing to L'Ecuyer-CMRG allows for reproducibility for parallel runs
    # See ?parallel::mc.parallel for more information
    old_rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    if (is.null(extras$seed)) {
      extras$seed <- 142 # allows any run to be reproducible
    }
    set.seed(extras$seed)

    # prep output directory
    outdir = gsub(' ', '\\\\ ', outdir)
    if (.Platform$OS.type == 'windows') {
      shell(paste('mkdir', outdir))
    } else {
      system(paste('mkdir -p', outdir))
    }

    #TODO: check event_probs

    ### create the splice variants for every event ----
    # extras$max_genes <- 100
    # extras$ncores <- 40
    variants <- create_splicing_variants_and_annotation(gtf_path,
                                                        extras$max_genes,
                                                        seqpath,
                                                        event_probs,
                                                        outdir,
                                                        extras$ncores)

    #TODO: make the transcript expression
    nr_transcripts <- length(unique(variants$transcript_id))
    if (is.null(fold_changes)) {
      extras$fold_changes <-
        matrix(c(
          rep(2, 2),
          rep(1, nr_transcripts - 2),
          rep(1, 2),
          rep(4, 2),
          rep(1, nr_transcripts - 4)
        ),
        nrow = nr_transcripts)
    } else {
      extras$fold_changes <- fold_changes
    }
    if (is.null(extras$reads_per_transcript)) {
      extras$reads_per_transcript <- rep(300, nr_transcripts)
    } else if (length(extras$reads_per_transcript) == 1) {
      extras$reads_per_transcript <-
        rep(extras$reads_per_transcript, nr_transcripts)
    }
    if (extras$exon_junction_coverage) {
      extras$exon_junction_table <- variants[, ID := .I]
    } else {
      extras$exon_junction_table <- NULL
    }
    extras$gtf <- file.path(outdir, 'splicing_variants.gtf')
    extras$seqpath <- seqpath
    extras$outdir <- outdir

    ### simulate with polyester----
    #TODO: Problem with the exon_id which is not present
    do.call(simulate_experiment, extras)

    ### make statistics ----
    # create_exon_junction_resolution(outdir,
    #                                 variants_annotation$variants,
    #                                 variants_annotation$event_annotation,
    #                                 extras$readlen,
    #                                 extras$nr_cores)
    #

    # Restore whatever RNG the user had set before running this function
    RNGkind(old_rng[1], old_rng[2])

  }


### debugging ----
library(rtracklayer)
library(data.table)
library(devtools)
install(pkg = 'polyester', quick = T)
library(polyester)

ncores <- 40
data.table::setDTthreads(ncores)

gtf_path <-
  '/nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99.gtf'
event_probs <-
  setNames(
    list(0.1, 0.1, 0.1, 0.4, 0.1, 0.1, 0.1, 0.1, 0.1),
    c(
      'es',
      'mes',
      'ir',
      'a3',
      'a5',
      'afe',
      'ale',
      'mee',
      'ale,afe,mee,ir,mes,a3,es,a5'
    )
  )
seqpath <-
  '/nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99.fa'
outdir <- '/nfs/home/students/ga89koc/hiwi/as_simulator/outdir_small'

source('R/create_splicing_variants_and_annotation.R')
# exon_junction_resolution <- create_splicing_variants_and_annotation(
#   gtf_path = gtf_path,
#   max_genes = 100,
#   seqpath = seqpath,
#   event_probs = event_probs,
#   outdir = outdir,
#   ncores = ncores
# )
