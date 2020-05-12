.check_parameters <- function(args, input_dir) {
  # check input dir
  if (!dir.exists(input_dir)) {
    stop(sprintf("could not find input directory %s", input_dir))
  } else {
    gtfs <- list.files(input_dir, pattern = '\\.gtf$|\\.gff$')
    if (length(gtfs) == 0) {
      stop(sprintf("could not find gtf/gff in input directory %s", input_dir))
    } else {
      args$gtf_path <- file.path(input_dir, gtfs[1])
      if (length(gtfs) > 1) {
        warning(sprintf("found more than one gtf/gff file in input directory. using %s...", args$gtf_path))
      }
    }
    fastas <- list.files(input_dir, pattern = '\\.fa$')
    if (length(fastas) == 0)
      stop(sprintf("could not find fasta files in input directory %s", input_dir))
    args$valid_chromosomes <- sub('.fa', '', fastas)
    message(sprintf("found the following fasta files: %s", paste(fastas, collapse = ", ")))
    message("note that splice variants will only be constructed from chromosomes that have a corresponding fasta file")
    message('')
  }
  # check other params
  
  if (is.null(args$exon_junction_coverage))
    args$exon_junction_coverage <- T
  else
    stopifnot(is.logical(args$exon_junction_coverage))
  if (is.null(args$error_rate))
    args$error_rate <- 0
  else
    stopifnot(is.numeric(args$error_rate))
  if (is.null(args$strand_specific))
    args$strand_specific <- T
  else
    stopifnot(is.logical(args$strand_specific))
  if (is.null(args$shuffle))
    args$shuffle <- T
  stopifnot(is.logical(args$shuffle))
  if (is.null(args$fastq))
    args$fastq <- T
  else
    stopifnot(is.logical(args$fastq))
  if (is.null(args$write_gff))
    args$write_gff <- T
  else
    stopifnot(is.logical(args$write_gff))
  if (is.null(args$multi_events_per_exon))
    args$multi_events_per_exon <- F
  else
    stopifnot(is.logical(args$multi_events_per_exon))
  if (is.null(args$probs_as_freq))
    args$probs_as_freq <- F
  else
    stopifnot(is.logical(args$probs_as_freq))
  if (is.null(args$verbose))
    args$verbose <- T
  else
    stopifnot(is.logical(args$verbose))
  if (!is.null(args$seq_depth)) {
    stopifnot(is.numeric(args$seq_depth))
  }
  if (is.null(args$save_exon_superset)) {
    args$save_exon_superset <- T
  } else {
    stopifnot(is.logical(args$save_exon_superset))
  }

  return(args)
}

.check_event_probs <- function(event_probs, probs_as_freq){
  as_names <- c('a3', 'a5', 'ir', 'es', 'ale', 'afe', 'mee', 'mes')
  if (is.list(event_probs))
    event_probs <- unlist(event_probs)
  if (is.numeric(event_probs)) {
    if (any(event_probs < 0 | event_probs > 1) || is.null(names(event_probs)))
      stop('Event probabilites have to be provided as named list/vector and each entry must be a probability.')
    if (probs_as_freq && (sum(event_probs) > 1))
      stop('If probabilites should be used as relative frequencies the sum of all entries cannot be greater than one.')
    names(event_probs) <- tolower(names(event_probs))
    split_names <- unique(unlist(strsplit(names(event_probs), ',', fixed = T)))
    match.arg(split_names, as_names, T)
    return(event_probs)
  } else stop('Event probabilites have to be provided as named list/vector and each entry must be a probability.')
}

#' Simulate RNA-seq experiment with splicing variants 
#'
#' Firstly, exon supersets are created by joining all exons of a gene from a gtf file.
#' Next, splicing variants are created with documentation and event annotation based on the users input.
#' Finally, fastq files containing RNA-seq reads from the splice variants and the real exon and junction coverage 
#' are created using a modified version of the polyester R package available on https://github.com/quirinmanz/polyester.
#'
#' @param input_dir Character path to directory containing the gtf file 
#' from which splice variants are created and genome fasta files passed to polyester.
#' @param event_probs Named list/vector containing numerics corresponding
#'  to the probabilites to create the event(combination). 
#'  If \code{probs_as_freq} is \code{TRUE} \code{event_probs} correspond 
#'  to the relative frequency of occurences for the event (combination) and 
#'  in this case the sum of all frequencies has to be <=1.
#' @param outdir character, path to folder where simulated reads and all
#'   annotations should be written, with *no* slash at the end. By default,
#'   reads are written to current working directory.
#' @param ncores the number of cores to be utilized for parallel generation
#'   of splice variant creation and read simulation.
#' @param ... any of several other arguments that can be used to add nuance
#'   to the simulation and splice variant creation. See details.
#'
#' @details Reads are simulated from a GTF file which is produced by 
#'   \code{\link{create_splicing_variants_and_annotation}} plus DNA
#'   sequences.
#'
#'   Several optional parameters can be passed to this function to adjust the
#'   simulation. For polyester parameters refer to \code{\link{simulate_experiment}}:
#'   
#'   \itemize{
#'   \item \code{write_gff}: Additionally to the gtf file containing the splice variants,
#'   a gff3 file with the same content will be printed to the outdir. 
#'   Default \code{TRUE}
#'   \item \code{max_genes}: The maximum number of genes/exon supersets to be included 
#'   in the process of splice variant creation. 
#'   Default \code{NULL} which means that all available exon supersets will be used.
#'   \item \code{exon_junction_coverage}: Should the real coverage of exons, junctions 
#'   and retained introns be written into a additional file.
#'   Default \code{TRUE}
#'   \item \code{multi_events_per_exon}: Should it be possible to have more than one AS event 
#'   at the same exon if multiple variants are created for the same exon superset?
#'   !If this option is set to \code{TRUE}, there may occur unforeseen AS events 
#'   that are not documented in the event_annotation file!.
#'   Default \code{FALSE}
#'   \item \code{probs_as_freq}: Should \code{event_probs} be treated as relative frequencies instead of probabilities?
#'   Default \code{FALSE}
#'   \item \code{save_exon_superset}: Should the exon supersets be saved to .rda file?
#'   Default \code{TRUE}
#'   }
#'
#'   Parameters passed to polyester that we assigned different defaults to than in \code{\link{simulate_experiment}}:
#'   \itemize{
#'   \item \code{fold_changes}: Currently, ass introduces random isoform switches. 
#'   Those can be retraced in the sim_tx_info.txt file written by polyester.
#'   We plan on improving this in the future.
#'   \item \code{strand_specific}: Strand-specific simulation (1st read forward strand,
#'   2nd read reverse strand with respect to transcript sequence). Default \code{TRUE}.
#'   \item \code{meanmodel}: \code{reads_per_transcripts} as a function of transcript length. Default \code{TRUE}.
#'   \item \code{verbose}: Should progress messages be printed during the sequencing process? Default \code{TRUE}.
#'   \item \code{exon_junction_coverage}: Should the coverage of exons, junctions and retained introns be determined? Default \code{TRUE}.
#'   \item \code{exon_junction_table}: If \code{exon_junction_coverage=TRUE} a \code{data.table} produced by \code{\link{create_splicing_variants_and_annotation}} 
#'   to determine exon and intron coverage.
#'   }
#'
#' @references
#'   Alyssa C. Frazee, Andrew E. Jaffe, Ben Langmead, Jeffrey T. Leek, 
#'   Polyester: simulating RNA-seq datasets with differential transcript expression, 
#'   Bioinformatics, Volume 31, Issue 17, 1 September 2015, Pages 2778–2784, 
#'   https://doi.org/10.1093/bioinformatics/btv272
#'
#'
#' @return No return, but simulated reads, a simulation info file,
#'   an alternative splicing event annotation and exon and junction coverages are written
#'   to \code{outdir}.
#' 
#' @export
#' @import data.table
#' @importFrom stats runif
#' @importFrom polyester simulate_experiment
#' @importFrom parallel mclapply

simulate_alternative_splicing <-
  function(input_dir,
           event_probs,
           outdir,
           ncores = 1L,
           ...)
  {
    args <- .check_parameters(list(...), input_dir)
    
    event_probs <- .check_event_probs(event_probs, args$probs_as_freq)

    # Store the current random number generator to restore at the end
    # Changing to L'Ecuyer-CMRG allows for reproducibility for parallel runs
    # See ?parallel::mc.parallel for more information
    old_rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    if (is.null(args$seed)) {
      args$seed <- 142 # allows any run to be reproducible
    }
    set.seed(args$seed)
    data.table::setDTthreads(ncores)
    args$ncores <- ncores

    # prep output directory
    outdir <- gsub(' ', '\\\\ ', outdir)
    if (.Platform$OS.type == 'windows') {
      shell(paste('mkdir', outdir))
    } else {
      system(paste('mkdir -p', outdir))
    }

    ### create the splice variants for every event ----
    args$exon_junction_table <- create_splicing_variants_and_annotation(
      args$gtf_path,
      args$valid_chromosomes,
      event_probs,
      outdir,
      args$ncores,
      args$write_gff,
      args$max_genes,
      args$exon_junction_coverage,
      args$multi_events_per_exon,
      args$probs_as_freq,
      args$save_exon_superset
    )

    
    tr_per_gene <- args$exon_junction_table[, .(nr_transcripts = length(unique(transcript_id))), by = gene_id]$nr_transcripts
    args$fold_changes <- do.call(rbind, lapply(tr_per_gene, function(nr_tr){
      if (nr_tr == 1)
        return(matrix(rep(1, 2), ncol = ifelse(is.null(args$num_reps), 2, length(args$num_reps))))
      else {
        nr_groups = ifelse(is.null(args$num_reps), 2, length(args$num_reps))
        if (runif(1) < 0.5) 
          return(matrix(rep(1, nr_groups * nr_tr), nrow = nr_tr))
        m = matrix(c(2, rep(1, nr_tr - 1)))
        for (i in 2:nr_groups) {
          c2 = sample(m[,1])
          while (identical(m[,1], c2)) {
            c2 <- sample(m[,1])
          }
          m = cbind(m, c2)
        }
        return(m)
      }
    }))
    
    args$meanmodel = T
    args$gtf <- file.path(outdir, 'splicing_variants.gtf')
    args$seqpath <- input_dir
    args$outdir <- outdir

    ### simulate with polyester----
    message('start simulation with polyester:')
    do.call(simulate_experiment, args)

    # Restore whatever RNG the user had set before running this function
    RNGkind(old_rng[1], old_rng[2])
    invisible(NULL)
  }


### debugging ----
# max_genes = 16
# multi_events_per_exon = T
# prob_as_freq = T
# seq_depth = 2e06
# params = list(
#   ncores = 4,
#   input_dir = '../ensembl_data/Homo_sapiens.GRCh38.99/',
#   event_probs =
#     setNames(c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
#              c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')),
#   max_genes = max_genes,
#   outdir = sprintf(
#     '../ensembl_data/Homo_sapiens.GRCh38.99_out/maxGenes%d_SeqDepth%d_multiEventsPerExon%s_eventsAsFreq%s',
#     max_genes,
#     seq_depth,
#     multi_events_per_exon,
#     prob_as_freq
#   ),
#   seq_depth = seq_depth,
#   multi_events_per_exon = multi_events_per_exon,
#   prob_as_freq = prob_as_freq
# )
#
# do.call(simulate_alternative_splicing, params)
