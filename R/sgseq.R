.check_coverage <- function(mates, transcript, readlen){
  fromStart <- mates[1] == 'mate1Start:1'
  if (fromStart) {
    #TODO: check if readlength is right with polyester
    reads <- IRanges(c(1, max(transcript$tr_end) - readlen + 1L), c(readlen, max(transcript$tr_end)))
  } else {
    coords <- as.numeric(unlist(strsplit(substring(mates, 7), '-', T)))
    reads <- IRanges(coords[c(1,3)], coords[c(2,4)])
  }
  exon_ids <- transcript[type == 'exon' & IRanges(tr_start, tr_end) %over% reads]$exon_id
  jct_ids <- transcript[type == 'junction' & IRanges(tr_start, tr_end) %within% reads]$exon_id
  return(c(exon_ids, jct_ids))
}

#' Internal sequencing function
#'
#' This internal function actually does the sequencing in the following steps:
#' generate the fragments, get the reads from the fragments, add any errors,
#' and then write the reads. It does this in 1 million reads chunks to
#' balance speed and memory usage. This is not intended to be called directly;
#' instead it is meant to be called via \code{\link{simulate_experiment}},
#' \code{\link{simulate_experiment_countmat}}, or \code{\link{simulate_experiment_empirical}}.
#'
#' @param readmat a N by M matrix describing the full counts for the
#'   experiment; each row is an experiment and each column is a sample in the
#'   experiment.
#' @param transcripts a \code{DNAStringSet} object containing the transcripts
#'   from which reads are simulated
#' @param paired a boolean determining whether or not to simulated paired-end
#'   or single-end reads
#' @param outdir the output directory where the simulated reads will be written
#' @param extras a list of extra options, generated internally by
#'   \code{simulate_experiment}, \code{simulate_experiment_countmat}, or
#'   \code{simulate_experiment_empirical}
#' @param reportCoverage whether to write out coverage information to
#'   \code{sample_coverages.rda} file in the \code{outdir}.
#'   defaults to \code{FALSE}
#' @param ncores the number of cores to be utilized for parallel generation
#'   of simulated reads. Note that if more than one core is specified,
#'   the code will parallelize by replicate, so if num_reps == 1, this
#'   will not be utilized.
#'
#' @return no return value. This function is called for its side effect of
#'   generating simulated reads and writing them to file.
#' @importFrom parallel mclapply
#' @import data.table
sgseq = function(readmat, transcripts, paired, outdir, extras, reportCoverage=FALSE, ncores=1L){
  #report theoretically perfect coverage if reportCoverage=TRUE, will write a file
  if(reportCoverage==TRUE){
    templates = unique(transcripts)
    coverage_matrices = list()
    for(i in 1:length(templates)){coverage_matrices = c(coverage_matrices, list(matrix(0, ncol=dim(readmat)[2], width(templates)[i])))}
    names(coverage_matrices) = names(templates)
  }

  message('start sequencing... (1m reads per iteration)')

  exon_junction_counts <- parallel::mclapply(seq_len(ncol(readmat)), function(i) {
    sample_name = sprintf("sample_%02d", i)

    message(sprintf('%s: overall %d reads', sample_name, sum(readmat[,i])))
    ##$ begin small chunk regarding fragment GC bias or not
    if (is.matrix(extras$frag_GC_bias)) {
      frag_GC_bias <- extras$frag_GC_bias[,i]
    } else {
      frag_GC_bias <- 'none'
    }
    ### end
    tObj = rep(transcripts, times=readmat[,i])
    iterations = ceiling(length(tObj) / 1e6L)
    offset = 1L
    region_counts <- list()
    for(iteration in seq_len(iterations)) {
      message(sprintf('%s: iteration %02d', sample_name, iteration))
      tSubset = tObj[offset:min(offset+999999L, length(tObj))] ## corrected value of integer added to offset to avoid duplicating reads
      tFrags = generate_fragments(tSubset, extras$fraglen[i], extras$fragsd[i],
                                  extras$readlen, extras$distr, extras$custdens,
                                  extras$bias, frag_GC_bias, extras$exon_junction_table)
      message(sprintf('%s: fragments generated', sample_name))
      if (extras$exon_junction_coverage) {
        covered_IDs = table(tFrags$coveredIDs)
        covered_IDs = data.table(
          ID = names(covered_IDs),
          count_ = as.vector(covered_IDs)
        )
        region_counts[[iteration]] = covered_IDs
      }
      tFrags = tFrags$tFrags
      if (!extras$strand_specific) {
        #reverse_complement some of those fragments
        tFrags = reverse_complement(tFrags)
      }

      #get reads from fragments
      reads = get_reads(tFrags, extras$readlen, paired)
      if(reportCoverage==TRUE){
            read_info = unique(names(reads))
            read_info_split = strsplit(read_info, ";mate1:|;mate2:")
            read_info_matrix = matrix(unlist(read_info_split), ncol=3, byrow=T)
            for(j in 1:dim(read_info_matrix)[1]){
                  read = read_info_matrix[j,]
                  target = which(names(coverage_matrices)==read[1])
                  # ML: changing these to strsplit (str_split requires stringr depends or imports)
                  coords1 = unlist(strsplit(read[2], "-"))
                  coords2 = unlist(strsplit(read[3], "-"))
                  coverage_matrices[[target]][coords1[1]:coords1[2],i]=coverage_matrices[[target]][coords1[1]:coords1[2],i]+1
                  coverage_matrices[[target]][coords2[1]:coords2[2],i]=coverage_matrices[[target]][coords2[1]:coords2[2],i]+1
                  save(coverage_matrices, file=file.path(outdir, 'sample_coverages.rda') )
            }
      }

      #add sequencing error
      if(extras$error_model == 'uniform'){
          errReads = add_error(reads, extras$error_rate)
      }else if(extras$error_model == 'custom'){
          errReads = add_platform_error(reads, 'custom', paired, extras$path)
      }else{
          errReads = add_platform_error(reads, extras$error_model, paired)
      }

      #write read pairs
      message(sprintf('%s: write read pairs', sample_name))
      write_reads(errReads, readlen=extras$readlen,
          fname=file.path(outdir, sample_name), paired=paired,
          gzip=extras$gzip, offset=offset, shuffle = extras$shuffle, fastq = extras$fastq)
      offset = offset + 1e6L
    }
    if(extras$exon_junction_coverage){
      region_counts = rbindlist(region_counts)
      region_counts = region_counts[, (sum = sum(.SD[[1]])), by = ID]
      names(region_counts)[2] <- sample_name
      return(region_counts)
    }
  }, mc.cores = ifelse(ncol(readmat) >= ncores, ncores, ncol(readmat)))
  if (extras$exon_junction_coverage) {
    exon_junction_counts <- Reduce(function(x, y)
      merge(x, y, by = 'ID', all = T), exon_junction_counts)
    sample_names <- names(exon_junction_counts)['ID' != names(exon_junction_counts)]
    extras$exon_junction_table[, ID := as.character(ID)]
    extras$exon_junction_table[exon_junction_counts, on = 'ID', get('sample_names') := mget(paste0('i.', sample_names))]
    extras$exon_junction_table$ID <- NULL
    fwrite(x = extras$exon_junction_table, file = file.path(outdir, 'exon_junction_coverage.tsv'), quote = F, sep = '\t')
  }
  message('finished sequencing')
  invisible(NULL)
}

