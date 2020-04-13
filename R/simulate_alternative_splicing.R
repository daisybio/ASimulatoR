.check_parameters <- function(args, input_dir) {
  #TODO: check if they are the right class?
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
  if (is.null(args$exon_junction_coverage)) args$exon_junction_coverage <- T
  else stopifnot(is.logical(args$exon_junction_coverage))
  if (is.null(args$error_rate)) args$error_rate <- 0
  else stopifnot(is.numeric(args$error_rate))
  if (is.null(args$strand_specific))
    args$strand_specific <- T
  else stopifnot(is.logical(args$strand_specific))
  if (is.null(args$gzip))
    args$gzip <- T
else stopifnot(is.logical(args$gzip))
  if (is.null(args$shuffle))
    args$shuffle <- T
  stopifnot(is.logical(args$shuffle))
  if (is.null(args$fastq))
    args$fastq <- T
  else stopifnot(is.logical(args$fastq))
  if (is.null(args$readlen)) 
    args$readlen <- 150
   else stopifnot(is.numeric(args$readlen))
  if (is.null(args$write_gff)) 
    args$write_gff <- F
    else stopifnot(is.logical(args$write_gff))
  if (is.null(args$multi_events_per_exon)) 
    args$multi_events_per_exon <- F
  else stopifnot(is.logical(args$multi_events_per_exon))
  if (is.null(args$probs_as_freq)) 
    args$probs_as_freq <- F
  else stopifnot(is.logical(args$probs_as_freq))
  if (is.null(args$verbose)) 
    args$verbose <- T
  else stopifnot(is.logical(args$verbose))
  if (!is.null(args$seq_depth)){
    stopifnot(is.numeric(args$seq_depth))
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

#' simulate RNA-seq experiment with splicing variants 
#'
#' @param input_dir directory containing the gtf file from which splice variants are created
#' and fasta files passed to polyester
#' @param event_probs 
#' @param outdir character, path to folder where simulated reads and all
#'   annotations should be written, with *no* slash at the end. By default,
#'   reads are written to current working directory.
#' @param ncores the number of cores to be utilized for parallel generation
#'   of splice variant creation and read simulation.
#' @param ... any of several other arguments that can be used to add nuance
#'   to the simulation. See details.
#'
#' @details Reads can either be simulated from a FASTA file of transcripts
#'   (provided with the \code{fasta} argument) or from a GTF file plus DNA
#'   sequences (provided with the \code{gtf} and \code{seqpath} arguments).
#'   Simulating from a GTF file and DNA sequences may be a bit slower: it took
#'   about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, X,
#'   and Y in hg19.
#'
#'   Several optional parameters can be passed to this function to adjust the
#'   simulation. The options are:
#'
#'   \itemize{
#'   \item \code{readlen}: read length. Default 100.
#'   \item \code{lib_sizes}: Library size factors for the biological replicates.
#'   \code{lib_sizes} should have length equal to the total number of
#'   replicates in the experiment, i.e., \code{sum(num_reps)}. For each
#'   replicate, once the number of reads to simulate from each transcript for
#'   that replicate is known, all read numbers across all transcripts from that
#'   replicate are multiplied by the corresponding entry in \code{lib_sizes}.
#'   \item \code{distr} One of 'normal', 'empirical', or 'custom', which
#'   specifies the distribution from which to draw RNA fragment lengths. If
#'   'normal', draw fragment lengths from a normal distribution. You can provide
#'   the mean of that normal distribution with \code{fraglen} (defaults to 250)
#'   and the standard deviation of that normal distribution with \code{fragsd}
#'   (defaults to 25). You can provide a single number for each, or a
#'   vector with length equal to the total number of samples.
#'   If 'empirical', draw fragment lengths
#'   from a fragment length distribution estimated from a real data set. If
#'   'custom', draw fragment lengths from a custom distribution, which you can
#'   provide as the \code{custdens} argument. \code{custdens} should be a
#'   density fitted using \code{\link{logspline}}.
#'   \item \code{error_model}: The error model can be one of:
#'     \itemize{
#'     \item \code{'uniform'}: errors are distributed uniformly across reads.
#'     You can also provide an \code{'error_rate'} parameter, giving the overall
#'     probability of making a sequencing error at any given nucleotide. This
#'     error rate defaults to 0.005.
#'     \item \code{'illumina4'} or \code{'illumina5'}: Empirical error models.
#'     See \code{?add_platform_error} for more information.
#'     \item \code{'custom'}: A custom error model you've estimated from an
#'     RNA-seq data set using \code{GemErr}. See \code{?add_platform_error}
#'     for more info. You will need to provide both \code{model_path} and
#'     \code{model_prefix} if using a custom error model. \code{model_path} is
#'     the output folder you provided to \code{build_error_model.py}. This path
#'     should contain either two files suffixed _mate1 and _mate2, or a file
#'     suffixed _single. \code{model_prefix} is the 'prefix' argument you
#'     provided to \code{build_error_model.py} and is whatever comes before the
#'     _mate1/_mate2 or _single files in \code{model_path}.
#'     }
#'   \item \code{bias} One of 'none', 'rnaf', or 'cdnaf'. 'none'
#'   represents uniform fragment selection (every possible fragment in a
#'   transcript has equal probability of being in the experiment); 'rnaf'
#'   represents positional bias that arises in protocols using RNA
#'   fragmentation, and 'cdnaf' represents positional bias arising in protocols
#'   that use cDNA fragmentation (Li and Jiang 2012). Using the 'rnaf' model,
#'   coverage is higher in the middle of the transcript and lower at both ends,
#'   and in the 'cdnaf' model, coverage increases toward the 3' end of the
#'   transcript. The probability models used come from Supplementary Figure S3
#'   of Li and Jiang (2012). Defaults to 'none' if you don't provide this.
#'   \item \code{gcbias} list indicating which samples to add GC bias to, and
#'   from which models. Should be the same length as \code{sum(num_reps)};
#'   entries can be either numeric or of class \code{loess}. A numeric entry of
#'   0 indicates no GC bias. Numeric entries 1 through 7 correspond to the
#'   7 empirical GC models that ship with Polyester, estimated from GEUVADIS
#'   HapMap samples NA06985, NA12144, NA12776, NA18858, NA20542, NA20772,
#'   and NA20815, respectively. The code used to derive the empirical GC models
#'   is available at
#'   \url{https://github.com/alyssafrazee/polyester/blob/master/make_gc_bias.R}.
#'   A loess entry should be a loess prediction model
#'   that takes a GC content percent value (between 0 and 1) a transcript's
#'   deviation from overall mean read count based on that GC value. Counts for
#'   each replicate will be adjusted based on the GC bias model specified for
#'   it. Numeric and loess entries can be mixed. By default, no bias is
#'   included.
#'   \item \code{frag_GC_bias} Either a matrix of dimensions 101 x \code{sum(num_reps)}
#'   or 'none'. The default is 'none'. If specified, the matrix contains the probabilities
#'   (a number in the range [0,1]) that a fragment will appear in the output given its GC content.
#'   The first row corresponds to a fragment with GC content of 0 percent, the second row
#'   1 percent, the third row 2 percent, etc., and the last row 100 percent.
#'   The columns correspond to different probabilites for each sample. Internally,
#'   a coin flip (a Bernoulli trial) determines if each fragment is kept, depending on its GC content.
#'   Note that the final library size will depend on the elements of the matrix, and
#'   it might make sense to scale up the \code{lib_size} of the samples with low
#'   probabilites in the matrix in the range of the transcriptome GC content distribution.
#'   Note that the \code{count_matrix} written to \code{outdir} contains the counts before
#'   applying fragment GC bias.
#'   \item \code{strand_specific} defaults to \code{FALSE}, which means fragments are
#'   generated with equal probability from both strands of the transcript sequence.
#'   set to \code{TRUE} for strand-specific simulation (1st read forward strand,
#'   2nd read reverse strand with respect to transcript sequence).
#'   \item \code{meanmodel}: set to TRUE if you'd like to set
#'   \code{reads_per_transcripts} as a function of transcript length. We
#'   fit a linear model regressing transcript abundance on transcript length,
#'   and setting \code{meanmodel=TRUE} means we will use transcript lengths
#'   to draw transcript abundance based on that linear model. You can see our
#'   modeling code at \url{http://htmlpreview.github.io/?https://github.com/alyssafrazee/polyester_code/blob/master/length_simulation.html}
#'   \item \code{write_info}: set to FALSE if you do not want files of
#'   simulation information written to disk. By default, transcript fold
#'   changes and expression status, replicate library sizes and group
#'   identifiers, and an R data object of the counts matrix (before
#'   application of fragment GC bias) are written to \code{outdir}.
#'   \item \code{seed}: specify a seed (e.g. \code{seed=142} or some other
#'   integer) to set before randomly drawing read numbers, for reproducibility.
#'   If one is not provided, a default of 142 will be used.
#'   \item \code{transcriptid}: optional vector of transcript IDs to be written
#'   into \code{sim_info.txt} and used as transcript identifiers in the output
#'   fasta files. Defaults to \code{names(readDNAStringSet(fasta))}. This
#'   option is useful if default names are very long or contain special
#'   characters.
#'   \item \code{gzip}: pass \code{gzip=TRUE} to write gzipped fasta files as
#'   output (by default, fasta output files are not compressed when written to
#'   disk).
#'   \item \code{exononly}: (passed to \code{\link{seq_gtf}}) if \code{TRUE}
#'   (as it is by default), only create transcript sequences from the features
#'   labeled \code{exon} in \code{gtf}.
#'   \item \code{idfield}: (passed to \code{\link{seq_gtf}})in the
#'   \code{attributes} column of \code{gtf}, what is the name of the field
#'   identifying transcripts? Should be character. Default
#'   \code{"transcript_id"}.
#'   \item \code{attrsep}: (passed to \code{\link{seq_gtf}}) in the
#'   \code{attributes} column of \code{gtf}, how are attributes separated?
#'   Default \code{"; "}.
#'   \item \code{shuffle}: should the reads be shuffled before written to file?
#'   Default \code{TRUE}.
#'   }
#'
#' @references
#'   't Hoen PA, et al (2013): Reproducibility of high-throughput mRNA and
#'   small RNA sequencing across laboratories. Nature Biotechnology 31(11):
#'   1015-1022.
#'
#'   Li W and Jiang T (2012): Transcriptome assembly and isoform expression
#'   level estimation from biased RNA-Seq reads. Bioinformatics 28(22):
#'   2914-2921.
#'
#'   McElroy KE, Luciani F and Thomas T (2012): GemSIM: general,
#'   error-model based simulator of next-generation sequencing data. BMC
#'   Genomics 13(1), 74.
#'
#' @return No return, but simulated reads, a simulation info file,
#'   an alternative splicing event annotation and exon and junction coverages are written
#'   to \code{outdir}. Note that reads are written out transcript by transcript
#'   and so need to be shuffled if used as input to quantification algorithms such
#'   as eXpress or Salmon.
#'
#' @export

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
      args$probs_as_freq
    )

    
    #TODO: make the transcript expression
    # browser()
    nr_transcripts <- length(unique(args$exon_junction_table$transcript_id))
    # tr_lengths <- args$exon_junction_table[type == 'exon', .(tr_length = sum(width)), by = transcript_id]
    # b0 = -3.0158
    # b1 = 0.8688
    # sigma = 4.152
    # logmus = b0 + b1 * log2(tr_lengths$tr_length) + rnorm(tr_lengths$tr_length, 
    #                                                     0, sigma)
    # reads_per_transcript = 2^logmus - 1
    # reads_per_transcript = pmax(reads_per_transcript, 1e-06)
    # basemeans = matrix(c(reads_per_transcript, reads_per_transcript), 
    #        nrow = length(reads_per_transcript))
    # size = basemeans/3
    # NB(as.matrix(basemeans), as.matrix(size))
    if (is.null(args$fold_changes)) {
      args$fold_changes <-
        matrix(c(
          rep(2, 2),
          rep(1, nr_transcripts - 2),
          rep(1, 2),
          rep(4, 2),
          rep(1, nr_transcripts - 4)
        ),
        nrow = nr_transcripts)
    }
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
# params = list(
#   ncores = 4,
#   input_dir = '../ensembl_data/Homo_sapiens.GRCh38.99/',
#   event_probs =
#     setNames(c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
#              c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')),
#   max_genes = max_genes,
#   outdir = sprintf(
#     'outdir/maxGenes%d_multiEventsPerExon%s_eventsAsFreq%s',
#     max_genes,
#     multi_events_per_exon,
#     prob_as_freq
#   ),
#   gzip = T,
#   write_gff = T,
#   verbose = T,
#   multi_events_per_exon = multi_events_per_exon,
#   prob_as_freq = prob_as_freq
# )
#
# do.call(simulate_alternative_splicing, params)
