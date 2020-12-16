.check_input_dir <- function(input_dir) {
  params <- list()
  # check input dir
  if (!dir.exists(input_dir)) {
    stop(sprintf("could not find input directory %s", input_dir))
  } else {
    gtfs <- list.files(input_dir, pattern = '\\.gtf$|\\.gff3?$')
    if (length(gtfs) == 0) {
      stop(sprintf("could not find gtf/gff in input directory %s", input_dir))
    } else {
      params$gtf_path <- file.path(input_dir, gtfs[1])
      if (length(gtfs) > 1) {
        warning(sprintf("found more than one gtf/gff file in input directory. using %s...", params$gtf_path))
      }
    }
    fastas <- list.files(input_dir, pattern = '\\.fa$')
    if (length(fastas) == 0)
      stop(sprintf("could not find fasta files in input directory %s", input_dir))
    params$valid_chromosomes <- sub('.fa', '', fastas)
    message(sprintf("found the following fasta files: %s", paste(fastas, collapse = ", ")))
    message("note that splice variants will only be constructed from chromosomes that have a corresponding fasta file")
    message('')
    return(params)
  }
}

.check_parameters <- function(params) {
  # check other arguments
  if (is.null(params$exon_junction_coverage))
    params$exon_junction_coverage <- T
  else
    stopifnot(is.logical(params$exon_junction_coverage))
  if (is.null(params$strand_specific))
    params$strand_specific <- T
  else
    stopifnot(is.logical(params$strand_specific))
  if (is.null(params$shuffle))
    params$shuffle <- T
  stopifnot(is.logical(params$shuffle))
  if (is.null(params$fastq))
    params$fastq <- T
  else
    stopifnot(is.logical(params$fastq))
  if (is.null(params$write_gff))
    params$write_gff <- T
  else
    stopifnot(is.logical(params$write_gff))
  if (is.null(params$multi_events_per_exon))
    params$multi_events_per_exon <- F
  else
    stopifnot(is.logical(params$multi_events_per_exon))
  if (is.null(params$probs_as_freq))
    params$probs_as_freq <- F
  else
    stopifnot(is.logical(params$probs_as_freq))
  if (is.null(params$verbose))
    params$verbose <- T
  else
    stopifnot(is.logical(params$verbose))
  if (!is.null(params$seq_depth)) {
    stopifnot(is.numeric(params$seq_depth))
  }
  if (is.null(params$save_exon_superset)) {
    params$save_exon_superset <- T
  } else {
    stopifnot(is.logical(params$save_exon_superset))
  }
  if (!is.null(params$novel_variants)){
    stopifnot(params$novel_variants >= 0 && params$novel_variants <= 1)
  }

  return(params)
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
#' Firstly, exon supersets are created by joining all exons of a gene from a gtf/gff file.
#' Next, splicing variants are created with documentation and event annotation based on the users input.
#' Finally, fastq files containing RNA-seq reads from the splice variants and the real exon and junction coverage 
#' are created using a modified version of the polyester R package available on https://github.com/biomedbigdata/polyester.
#'
#' @param input_dir Character path to directory containing the gtf/gff file 
#' from which splice variants are created and genome fasta files with 
#' one file per chromosome i.e. <chr_name>.fa passed to polyester.
#' @param outdir character, path to folder where simulated reads and all
#'   annotations should be written, with *no* slash at the end. By default,
#'   reads are written to current working directory.
#' @param event_probs Named list/vector containing numerics corresponding
#'  to the probabilites to create the event(-combination). 
#'  If \code{probs_as_freq} is \code{TRUE} \code{event_probs} correspond 
#'  to the relative frequency of occurences for the event (combination) and 
#'  in this case the sum of all frequencies has to be <=1.
#'  No default, must not be \code{NULL}, except if \code{preset} is given.
#' @param preset if you want to use preset parameters one of 
#' 'event_partition', 'experiment_bias', 'event_combination_2'.
#' Check \code{?\link{presets}} for more information
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
#'   simulation. For polyester parameters refer to \code{\link{polyester::simulate_experiment}}:
#'   
#'   \itemize{
#'   \item \code{novel_variants}: Numeric value between 0 and 1 indicating the percentage 
#'   of splicing variants that will not be written to an additional gtf file splicing_variants_novel.gtf.
#'   \item \code{write_gff}: Additionally to the gtf file containing the splice variants,
#'   a gff3 file with the same content will be printed to the outdir. 
#'   Default \code{TRUE}
#'   \item \code{max_genes}: The maximum number of genes/exon supersets to be included 
#'   in the process of splice variant creation. 
#'   Default \code{NULL} which means that all available exon supersets will be used.
#'   **This is a computation heavy default and you might want to adjust it!**
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
#'   \item \code{fold_changes}: Currently, ASimulatoR introduces random isoform switches. 
#'   Those can be retraced in the sim_tx_info.txt file written by polyester.
#'   We plan on improving this in the future.
#'   \item \code{strand_specific}: Strand-specific simulation (1st read forward strand,
#'   2nd read reverse strand with respect to transcript sequence). Default \code{TRUE}.
#'   \item \code{meanmodel}: \code{reads_per_transcripts} as a function of transcript length. Always \code{TRUE} in ASimulatoR.
#'   \item \code{frag_GC_bias}: A sample-specific GC content bias on the fragment level. Currently not supported in ASimulatoR: always 'none'.
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
#' @importFrom utils data
#' @importFrom polyester simulate_experiment
#' @importFrom pbmcapply pbmclapply

simulate_alternative_splicing <-
  function(input_dir,
           outdir,
           event_probs = NULL,
           preset = NULL,
           ncores = 1L,
           ...)
  {
    # check parameters and compatibility
    params <- .check_input_dir(input_dir)
    
    preset_res <- assign_preset(preset, event_probs)
    
    extras <- list(...)
    
    extras_presets <- .check_parameters(c(extras, 
                                          preset_res[!(names(preset_res) %in% names(extras))]))
    params <- c(extras_presets, 
                params)
    
    # means that we have preset event_probs
    if (is.null(event_probs)) event_probs <- params$event_probs
    event_probs <- .check_event_probs(event_probs, params$probs_as_freq)

    # Store the current random number generator to restore at the end
    # Changing to L'Ecuyer-CMRG allows for reproducibility for parallel runs
    # See ?parallel::mc.parallel for more information
    old_rng <- RNGkind()
    RNGkind("L'Ecuyer-CMRG")
    if (is.null(params$seed)) {
      params$seed <- 142 # allows any run to be reproducible
    }
    set.seed(params$seed)
    data.table::setDTthreads(ncores)
    message(sprintf('set data.table threads to %s', data.table::getDTthreads()))
    params$ncores <- ncores

    # prep output directory
    outdir <- gsub(' ', '\\\\ ', outdir)
    if (.Platform$OS.type == 'windows') {
      shell(paste('mkdir', outdir))
    } else {
      system(paste('mkdir -p', outdir))
    }

    ### create the splice variants for every event ----
    params$exon_junction_table <- create_splicing_variants_and_annotation(
      params$gtf_path,
      params$valid_chromosomes,
      event_probs,
      outdir,
      params$ncores,
      params$write_gff,
      params$max_genes,
      params$exon_junction_coverage,
      params$multi_events_per_exon,
      params$probs_as_freq,
      params$save_exon_superset,
      params$novel_variants
    )

    
    tr_per_gene <- params$exon_junction_table[, .(nr_transcripts = length(unique(transcript_id))), by = gene_id]$nr_transcripts
    params$fold_changes <- do.call(rbind, lapply(tr_per_gene, function(nr_tr){
      if (nr_tr == 1)
        return(matrix(rep(1, 2), ncol = ifelse(is.null(params$num_reps), 2, length(params$num_reps))))
      else {
        nr_groups = ifelse(is.null(params$num_reps), 2, length(params$num_reps))
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
    
    params$meanmodel <- T
    params$frag_GC_bias <- NULL
    params$gtf <- file.path(outdir, 'splicing_variants.gtf')
    params$seqpath <- input_dir
    params$outdir <- outdir

    ### simulate with polyester----
    message('start simulation with polyester:')
    do.call(simulate_experiment, params)

    # Restore whatever RNG the user had set before running this function
    RNGkind(old_rng[1], old_rng[2])
    invisible(NULL)
  }


# ## debugging ----
# ## set parameters
# input = '../ensembl_data/Homo_sapiens.GRCh38.99'
# output = '../ensembl_data/Homo_sapiens.GRCh38.99_out/experiment_bias'
# 
# ## run simulator
# # library(ASimulatoR)
# simulate_alternative_splicing(input_dir = input, outdir = output, preset = 'experiment_bias', max_genes = 100, seq_depth = 2e6)
