.get_transcriptomic_coord <- function(widths){
  tr_end <- cumsum(widths)
  tr_start <- c(1L, (tr_end + 1L)[-length(widths)])
  list(tr_start = tr_start, tr_end = as.integer(tr_end))
}

.add_transcript_and_junction_lines <- function(variant){
  # create transcript line
  tr <- range(variant)
  S4Vectors::mcols(tr) <- S4Vectors::mcols(variant[1])
  tr$type <- 'transcript'
  # create junction lines
  junctions <- GenomicRanges::gaps(variant, start = min(IRanges::start(variant))) + 1
  if (length(junctions) > 0)
    S4Vectors::mcols(junctions) <- S4Vectors::DataFrame(
      source = "as_simulator",
      type = "junction",
      score = ".",
      phase = ".",
      gene_id = variant$gene_id[1],
      transcript_id = variant$transcript_id[1],
      template = variant$template[1],
      gene_exon_number = 0L,
      tr_start = variant$tr_end[-length(variant)],
      tr_end = variant$tr_start[-1]
    )
  # add retained intron
  if ((!variant$template[1]) && grepl('ir', variant$transcript_id[1], fixed = t)) {
    variant <- c(parent.frame()$ri, variant)
  }
  # merge everything
  variant <- sort(c(variant, junctions), decreasing = (S4Vectors::runValue(GenomicRanges::strand(variant)) == '-'))
  variant <- c(tr, variant)
  variant$gene_exon_number[variant$type != 'exon'] <- NA
  variant[1]$tr_start <- 1L
  variant[1]$tr_end <- variant[length(variant)]$tr_end
  variant
}

.add_gene_line <- function(orig_template, template){
  g <- range(orig_template)
  S4Vectors::mcols(g) <- S4Vectors::mcols(template[1])
  g$type <- 'gene'
  template <- c(g, template)
  template[1]$transcript_id <- NA
  template[1]$template <- NA
  template[1]$tr_start <- 0L
  template[1]$tr_end <- 0L
  template
}


.assign_events <- function(v, min_nr_exons_per_event) {
  # done with all assignments
  if (length(min_nr_exons_per_event) == 0) return(list())
  # assign event and split vector, try until it works
  # TODO: assign greedy after maximum number of attemtps exceeded
  while (T) {
    result <- list()
    event <- draw_one_sample_safely(min_nr_exons_per_event)
    # does not work
    if (event > length(v)) {
      browser()
      stop('Error in assign events: Event length is bigger than vector.')
    }
    event_start <- draw_one_sample_safely(1:(length(v) - event + 1))
    vs <-
      list(v1 = v[1:event_start], v2 = v[(event_start + event - 1):length(v)])
    result[[names(event)]] <-
      v[event_start:(event_start + event - 1)]
    min_nr_exons_per_event_tmp <-
      min_nr_exons_per_event[names(min_nr_exons_per_event) != names(event)]
    events_first <-
      sample(c(TRUE, FALSE),
             length(min_nr_exons_per_event_tmp),
             replace = TRUE)
    if (sum(min_nr_exons_per_event_tmp[events_first]) - (sum(events_first) - 1L) > length(vs$v1) ||
        sum(min_nr_exons_per_event_tmp[!events_first]) - (sum(!events_first) - 1L) > length(vs$v2))
      next
    result_v1 <-
      .assign_events(vs$v1, min_nr_exons_per_event_tmp[events_first])
    result_v2 <-
      .assign_events(vs$v2, min_nr_exons_per_event_tmp[!events_first])
    if ((length(result_v1) + length(result_v2)) == length(min_nr_exons_per_event_tmp))
      return(c(result, result_v1, result_v2))
  }
}

#' Internal function to create alternative splicing events
#'
#' This is not intended to be called directly;
#' instead it is meant to be called via \code{\link{simulate_alternative_cplicing}}
#'
#' @param gtf_path
#' @param valid_chromosomes
#' @param event_probs
#' @param outdir
#' @param ncores the number of cores to be utilized for parallel generation
#'   of splice variants.
#' @param write_gff
#' @param max_genes
#' @param exon_junction_coverage
#'
#' @return if \code{exon_junction_coverage = TRUE} the exons, junctions and retained introns as data table in gtf style
#'
create_splicing_variants_and_annotation <-
  function(gtf_path,
           valid_chromosomes,
           event_probs,
           outdir,
           ncores,
           write_gff,
           max_genes = NULL,
           exon_junction_coverage = T,
           multi_events_per_exon = F,
           probs_as_freq = F) {

    ### create exon_superset ----
    exon_supersets <- get_exon_supersets(gtf_path, valid_chromosomes, ncores)

    ### assign as events to supersets ----
    gene_lengths <- sapply(exon_supersets, length)
    nr_genes <- min(sum(gene_lengths > 1), max_genes)
    if (probs_as_freq) {
      construct_all <- sapply(names(event_probs), function(event) {
        rep(names(event_probs) == event, floor(event_probs[[event]] * nr_genes))
      })
      construct_all <- c(construct_all, rep(F, length(event_probs)*nr_genes - length(construct_all)))
    } else {
      construct_all <- runif(nr_genes * length(event_probs)) < event_probs
    }
    construct_all_list <- split(construct_all, ceiling(seq_along(construct_all)/length(event_probs)))
    min_nr_exons_per_event <-
      setNames(c(2, 2, 2, 3, 3, 3, 4, 4),
               c('a3', 'a5', 'ir', 'es', 'ale', 'afe', 'mee', 'mes'))
    min_nr_exons <- sapply(construct_all_list, function(construct){
      construct <- names(event_probs)[construct]
      get_min_nr_exons(min_nr_exons_per_event, construct, multi_events_per_exon)
    })
    drawn_genes <- character(nr_genes)
    for (i in 1:nr_genes) {
      if (min_nr_exons[i] > 0) {
        drawn_genes[i] <-
          draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
        gene_lengths[[drawn_genes[i]]] <- -1
      } else {
        #TODO: what to do if we did not draw an as event?
        drawn_genes[i] <-
          draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
        gene_lengths[[drawn_genes[i]]] <- -1
      }
    }

    ### create splice variants and annotation ----
    message('create splicing variants and annotation...')
    #TODO: is random number generation ok?
    all_variants_and_event_annotation <-
     parallel::mclapply(1:nr_genes, function(i) {
        construct <- names(event_probs)[construct_all_list[[i]]]
        orig_template <- exon_supersets[[drawn_genes[i]]]
        if (length(construct) == 0) {
          return(list(
            variants = .add_gene_line(
              orig_template,
              .add_transcript_and_junction_lines(orig_template)
            )
          ))
        } else {
          neg_strand <-
            S4Vectors::runValue(GenomicRanges::strand(orig_template)) == '-'
          gene_id <- orig_template$gene_id[1]
          event_combs <- strsplit(construct, ',', T)
          available_exons <- length(orig_template)
          exon_vector <-
            setNames(rep(T, available_exons), seq_len(available_exons))

          if (!multi_events_per_exon) {
            res <-
              construct_variant(
                unique(unlist(event_combs)),
                exon_vector,
                min_nr_exons_per_event,
                available_exons,
                orig_template,
                multi_events_per_exon,
                neg_strand,
                min_nr_exons[i]
              )
            exon_vector <- res$exon_vector
            event_exons <- res$event_exons
          }
          template <- orig_template[exon_vector]
          S4Vectors::mcols(template)[c('tr_start', 'tr_end')] <-
            .get_transcriptomic_coord(IRanges::width(template))


          variants_and_event_annotation <-
            lapply(event_combs, function(comb) {
              ### create splice variants ----
              if (multi_events_per_exon) {
                res <-
                  construct_variant(
                    comb,
                    exon_vector,
                    min_nr_exons_per_event,
                    available_exons,
                    orig_template,
                    multi_events_per_exon,
                    neg_strand
                  )
                exon_vector <- res$exon_vector
                event_exons <- res$event_exons
                variant <- res$variant
              } else {
                variant <-
                  orig_template[apply_exon_events(exon_vector, comb, event_exons)]
              }
              
              ### implement a3, a5 and ir
              if ('a3' %in% comb) {
                new_a3 <- res$new_a3
                a3_index <- res$a3_index
                if (neg_strand)
                  end(variant[variant$gene_exon_number == a3_index]) <- new_a3
                else
                  start(variant[variant$gene_exon_number == a3_index]) <- new_a3
              }
              if ('a5' %in% comb) {
                new_a5 <- res$new_a5
                a5_index <- res$a5_index
                if (neg_strand)
                  start(variant[variant$gene_exon_number == a5_index]) <- new_a5
                else
                  end(variant[variant$gene_exon_number == a5_index]) <- new_a5
              }
              if ('ir' %in% comb) {
                ir_indices <- as.integer(names(event_exons[['ir']]))
                IRanges::ranges(variant[variant$gene_exon_number == ir_indices[1]]) <-
                  range(IRanges::ranges(variant[variant$gene_exon_number %in% ir_indices]))
                variant[variant$gene_exon_number == ir_indices[2]] <-
                  NULL
                retaining_exon <- variant[variant$gene_exon_number == ir_indices[1]]
                flanking_exons <-
                  template[template$gene_exon_number %in% as.integer(names(event_exons$ir))]
                ri <-
                  GenomicRanges::gaps(flanking_exons, start = min(start(flanking_exons)))
                S4Vectors::mcols(ri) <-
                  S4Vectors::mcols(retaining_exon)
                ri$tr_start <- retaining_exon$tr_start +
                  ifelse(neg_strand,
                         (end(retaining_exon) - end(ri)),
                         start(ri) - start(retaining_exon))
                ri$tr_end <- ri$tr_start + IRanges::width(ri) - 1L
                ri$type <- 'ri'
                ri$transcript_id <-
                  sprintf('%s_%s', gene_id, paste(comb, collapse = ","))
              }

              variant$transcript_id <-
                sprintf('%s_%s', gene_id, paste(comb, collapse = ","))
              variant$template <- F
              S4Vectors::mcols(variant)[c('tr_start', 'tr_end')] <-
                .get_transcriptomic_coord(IRanges::width(variant))

              ### create event annotation ----

              if (multi_events_per_exon && res$two_variants) {
                event_annotation <-
                  get_event_annotation(
                    comb,
                    event_exons,
                    exon_vector,
                    neg_strand,
                    res$template_new,
                    variant
                  )
                variant <- c(.add_transcript_and_junction_lines(res$template_new), .add_transcript_and_junction_lines(variant))
              } else {
                event_annotation <-
                  get_event_annotation(
                    comb,
                    event_exons,
                    exon_vector,
                    neg_strand,
                    template,
                    variant
                  )
                variant <- .add_transcript_and_junction_lines(variant)
              }

              ### add transcript line ----
              list(event_annotation = event_annotation,
                   variant = variant)
            })

          variants_and_event_annotation <-
            unlist(variants_and_event_annotation, recursive = F)
          event_annos <-
            'event_annotation' == names(variants_and_event_annotation)
          event_annotation <-
            data.table::rbindlist(variants_and_event_annotation[event_annos])
          variants <-
            unlist(as(variants_and_event_annotation[!event_annos], "GRangesList"),
                   use.names = F)

          return(list(
            event_annotation = event_annotation,
            variants = c(
              .add_gene_line(
                orig_template,
                .add_transcript_and_junction_lines(template)
              ),
              variants
            )
          ))
        }
      }, mc.cores = ncores)
    all_variants_and_event_annotation <- unlist(all_variants_and_event_annotation, recursive = F)
    event_annos <- endsWith(names(all_variants_and_event_annotation), 'event_annotation')
    event_annotation <- data.table::rbindlist(all_variants_and_event_annotation[event_annos])
    variants <- unlist(as(all_variants_and_event_annotation[!event_annos], "GRangesList"), use.names = F)
    message('finished creating splicing variants and annotation')
    message('')

    ### exporting ----
    message('exporting gtf for read simulation...')
    rtracklayer::export(variants[variants$type != 'junction' & variants$type != 'ri'],
                        file.path(outdir, 'splicing_variants.gtf'))
    if (write_gff)
      rtracklayer::export(variants[variants$type != 'junction' & variants$type != 'ri'],
                          file.path(outdir, 'splicing_variants.gff3'))
    message('finished exporting gtf')
    message('')

    ### write annotations to file ----
    message('exporting event_annotation...')
    data.table::fwrite(
      x = event_annotation,
      file = file.path(outdir, 'event_annotation.tsv'),
      quote = F,
      sep = '\t'
    )
    message('finished exporting event_annotation...')
    message('')

    if (exon_junction_coverage) return(data.table::as.data.table(variants)[type != 'gene' & type != 'transcript'])
    else invisible(NULL)
  }
