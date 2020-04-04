.draw_one_sample_safely <- function(v) {
  if (length(v) == 1)
    return(v[1])
  else
    return(sample(v, 1))
}

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
    event <- .draw_one_sample_safely(min_nr_exons_per_event)
    # does not work
    if (event > length(v)) {
      browser()
      stop('Error in assign events: Event length is bigger than vector.')
    }
    event_start <- .draw_one_sample_safely(1:(length(v) - event + 1))
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

.apply_exon_events <- function(v, comb, event_exons) {
  if ('afe' %in% comb)
    v[c(1, 2)] <- v[c(2, 1)]
  if ('ale' %in% comb)
    v[c(length(v) - 1, length(v))] <- v[c(length(v), length(v) - 1)]
  if ('mee' %in% comb)
    v[names(event_exons[['mee']][c(2, 3)])] <-
      v[names(event_exons[['mee']][c(3, 2)])]
  if ('es' %in% comb)
    v[names(event_exons[['es']][2])] <- F
  if ('mes' %in% comb)
    v[names(event_exons[['mes']][-c(1, length(event_exons[['mes']]))])] <-
      F
  return(v)
}

.get_min_nr_exons <-
  function(min_nr_exons_per_event,
           construct,
           multi_events_per_exon = F) {
    event_combs <- strsplit(construct, ',', T)
    if (length(event_combs) == 0)
      return(0)
    if (multi_events_per_exon) {
      return(max(sapply(event_combs, function(all_events) {
        sum(min_nr_exons_per_event[all_events]) - (length(all_events) - 1L)
      })))
    } else {
      all_events <- unique(unlist(event_combs))
      return(sum(min_nr_exons_per_event[all_events]) - (length(all_events) - 1L))
    }
  }

.get_event_annotation <-
  function(comb,
           event_exons,
           exon_vector,
           neg_strand,
           template,
           variant) {

    event_annotation <- data.table::data.table(
      event_annotation = character(),
      variant = character(),
      template = character(),
      genomic_start = integer(),
      genomic_end = integer(),
      transcriptomic_start = integer(),
      transcriptomic_end = integer()
    )

    #TODO: wrap in function
    if ('afe' %in% comb) {
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'afe',
          variant$transcript_id[1],
          template$transcript_id[1],
          start(variant[1]),
          end(variant[1]),
          variant[1]$tr_start,
          variant[1]$tr_end
        ),
        list(
          'afe',
          template$transcript_id[1],
          variant$transcript_id[1],
          start(template[1]),
          end(template[1]),
          template[1]$tr_start,
          template[1]$tr_end
        )
      ))
    }
    if ('ale' %in% comb) {
      variant_length <- length(variant)
      template_length <- length(template)
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'ale',
          variant$transcript_id[1],
          template$transcript_id[1],
          start(variant[variant_length]),
          end(variant[variant_length]),
          variant[variant_length]$tr_start,
          variant[variant_length]$tr_end
        ),
        list(
          'ale',
          template$transcript_id[1],
          variant$transcript_id[1],
          start(template[template_length]),
          end(template[template_length]),
          template[template_length]$tr_start,
          template[template_length]$tr_end
        )
      ))
    }
    if ('mee' %in% comb) {
      incl_exons <-
        list(template = as.integer(names(event_exons$mee)[2]),
             variant = as.integer(names(event_exons$mee)[3]))
      if (exon_vector[names(event_exons$mee)[3]]) incl_exons[c(1,2)] <- incl_exons[c(2,1)]
      variant_exon <- variant[variant$gene_exon_number == incl_exons$variant]
      template_exon <- template[template$gene_exon_number == incl_exons$template]
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'mee',
          variant$transcript_id[1],
          template$transcript_id[1],
          start(variant_exon),
          end(variant_exon),
          variant_exon$tr_start,
          variant_exon$tr_end
        ),
        list(
          'mee',
          template$transcript_id[1],
          variant$transcript_id[1],
          start(template_exon),
          end(template_exon),
          template_exon$tr_start,
          template_exon$tr_end
        )
      ))
    }
    if ('mes' %in% comb) {
      skipped_exons <- template[template$gene_exon_number %in% as.integer(names((event_exons$mes[-c(1, length(event_exons$mes))])))]
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'mes',
          variant$transcript_id[1],
          template$transcript_id[1],
          paste(start(skipped_exons), collapse = ','),
          paste(end(skipped_exons), collapse = ','),
          paste(skipped_exons$tr_start, collapse = ','),
          paste(skipped_exons$tr_end, collapse = ',')
        )
      ))
    }
    if ('es' %in% comb) {
      skipped_exon <- template[template$gene_exon_number == as.integer(names(event_exons$es[2]))]
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'es',
          variant$transcript_id[1],
          template$transcript_id[1],
          start(skipped_exon),
          end(skipped_exon),
          skipped_exon$tr_start,
          skipped_exon$tr_end
        )
      ))
    }
    if ('ir' %in% comb) {
      ri <- parent.frame()$ri
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'ir',
          variant$transcript_id[1],
          template$transcript_id[1],
          start(ri),
          end(ri),
          ri$tr_start,
          ri$tr_end
        )
      ))
    }
    if ('a3' %in% comb) {
      new_a3 <- parent.frame()$new_a3
      alt_exon <- template[template$gene_exon_number == parent.frame()$a3_index]
      old_a3 <- ifelse(neg_strand, end(alt_exon), start(alt_exon))
      skipped_bases <- abs(old_a3 - new_a3)
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'a3',
          variant$transcript_id[1],
          template$transcript_id[1],
          ifelse(neg_strand, (new_a3 + 1L), old_a3),
          ifelse(neg_strand, old_a3, (new_a3 - 1L)),
          alt_exon$tr_start,
          alt_exon$tr_start + skipped_bases - 1L
        )
      ))
    }
    if ('a5' %in% comb) {
      new_a5 <- parent.frame()$new_a5
      alt_exon <- template[template$gene_exon_number == parent.frame()$a5_index]
      old_a5 <- ifelse(neg_strand, start(alt_exon), end(alt_exon))
      skipped_bases <- abs(old_a5 - new_a5)
      event_annotation <- data.table::rbindlist(list(
        event_annotation,
        list(
          'a5',
          variant$transcript_id[1],
          template$transcript_id[1],
          ifelse(neg_strand, old_a5, (new_a5 + 1)),
          ifelse(neg_strand, (new_a5 - 1L), old_a5),
          alt_exon$tr_start,
          alt_exon$tr_start + skipped_bases - 1L
        )
      ))
    }
    event_annotation
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
           max_genes,
           exon_junction_coverage,
           multi_events_per_exon) {

    ### create exon_superset ----
    exon_supersets <- get_exon_supersets(gtf_path, valid_chromosomes, ncores)

    ### assign as events to supersets ----
    gene_lengths <- sapply(exon_supersets, length)
    nr_genes <- min(sum(gene_lengths > 1), max_genes)
    construct_all <- runif(nr_genes * length(event_probs)) < event_probs
    construct_all_list <- split(construct_all, ceiling(seq_along(construct_all)/length(event_probs)))
    min_nr_exons_per_event <-
      setNames(c(2, 2, 2, 3, 3, 3, 4, 4),
               c('a3', 'a5', 'ir', 'es', 'ale', 'afe', 'mee', 'mes'))
    min_nr_exons <- sapply(construct_all_list, function(construct){
      construct <- names(event_probs)[construct]
      .get_min_nr_exons(min_nr_exons_per_event, construct, multi_events_per_exon)
    })
    drawn_genes <- character(nr_genes)
    for (i in 1:nr_genes) {
      if (min_nr_exons[i] > 0) {
        drawn_genes[i] <-
          .draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
        gene_lengths[[drawn_genes[i]]] <- -1
      } else {
        #TODO: what to do if we did not draw an as evemt?
        drawn_genes[i] <-
          .draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
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
          all_events <- unique(unlist(event_combs))
          exon_vector <-
            setNames(rep(T, available_exons), seq_len(available_exons))

          #TODO: wrap in function
          if (!multi_events_per_exon) {
            exon_vector_tmp <- exon_vector
            all_events_tmp <- all_events
            min_nr_exons_per_event_tmp <- min_nr_exons_per_event
            if ('afe' %in% all_events) {
              exon_vector[.draw_one_sample_safely(c(1, 2))] <- F
              exon_vector_tmp <- exon_vector_tmp[-c(1, 2)]
              all_events_tmp <-
                all_events_tmp[all_events_tmp != 'afe']
            }
            if ('ale' %in% all_events) {
              exon_vector[.draw_one_sample_safely(c((length(exon_vector) - 1), length(exon_vector)))] <-
                F
              exon_vector_tmp <-
                exon_vector_tmp[-c(length(exon_vector_tmp) - 1, length(exon_vector_tmp))]
              all_events_tmp <-
                all_events_tmp[all_events_tmp != 'ale']
            }
            if ('mes' %in% all_events) {
              mee_max_skipped_exons <- available_exons - min_nr_exons[i] + 2
              mee_nr_skipped_exons <-
                .draw_one_sample_safely(2:mee_max_skipped_exons)
              min_nr_exons_per_event_tmp['mes'] <-
                mee_nr_skipped_exons + 2
            }

            # assign the chunks with as events to the exon_vector
            event_exons <-
              .assign_events(exon_vector_tmp, min_nr_exons_per_event_tmp[all_events_tmp])

            if ('mee' %in% all_events) {
              exon_vector[.draw_one_sample_safely(names(event_exons$mee)[2:3])] <-
                F
            }

            if ('a3' %in% all_events) {
              a3_index <- as.integer(names(event_exons[['a3']][2]))
              new_a3 <-
                .draw_one_sample_safely((start(orig_template)[a3_index] + 1):(end(orig_template)[a3_index] - 1))
            }

            if ('a5' %in% all_events) {
              a5_index <- as.integer(names(event_exons[['a5']][1]))
              if (!is.null(event_exons[['a3']]) &&
                  a5_index == a3_index) {
                new_a5 <- ifelse(neg_strand,
                                 .draw_one_sample_safely((start(
                                   orig_template
                                 )[a5_index] + 1):(new_a3 - 1)),
                                 .draw_one_sample_safely((new_a3 + 1):(end(
                                   orig_template
                                 )[a5_index] - 1)))
              }
              else
                new_a5 <-
                  .draw_one_sample_safely((start(orig_template)[a5_index] + 1):(end(orig_template)[a5_index] - 1))
            }
          }
          template <- orig_template[exon_vector]
          S4Vectors::mcols(template)[c('tr_start', 'tr_end')] <-
            .get_transcriptomic_coord(IRanges::width(template))


          variants_and_event_annotation <-
            lapply(event_combs, function(comb) {
              ### create splice variants ----
              #TODO: wrap in function
              if (multi_events_per_exon) {
                exon_vector_tmp <- exon_vector
                exon_vector_new <- exon_vector
                comb_tmp <- comb
                min_nr_exons_per_event_tmp <- min_nr_exons_per_event
                two_variants <- F
                if ('afe' %in% comb) {
                  two_variants <- T
                  exon_vector_new[.draw_one_sample_safely(c(1, 2))] <- F
                  exon_vector_tmp <- exon_vector_tmp[-c(1, 2)]
                  comb_tmp <-
                    comb_tmp[comb_tmp != 'afe']
                }
                if ('ale' %in% comb) {
                  two_variants <- T
                  exon_vector_new[.draw_one_sample_safely(c((length(exon_vector_new) - 1), length(exon_vector_new)))] <-
                    F
                  exon_vector_tmp <-
                    exon_vector_tmp[-c(length(exon_vector_tmp) - 1, length(exon_vector_tmp))]
                  comb_tmp <-
                    comb_tmp[comb_tmp != 'ale']
                }
                if ('mes' %in% comb) {
                  mee_max_skipped_exons <- available_exons - .get_min_nr_exons(min_nr_exons_per_event, paste(comb, collapse = ','), multi_events_per_exon) + 2L
                  mee_nr_skipped_exons <-
                    .draw_one_sample_safely(2:mee_max_skipped_exons)
                  min_nr_exons_per_event_tmp['mes'] <-
                    mee_nr_skipped_exons + 2
                }

                # assign the chunks with as events to the exon_vector
                event_exons <-
                  .assign_events(exon_vector_tmp, min_nr_exons_per_event_tmp[comb_tmp])

                if ('mee' %in% comb) {
                  two_variants <- T
                  exon_vector_new[.draw_one_sample_safely(names(event_exons$mee)[2:3])] <-
                    F
                }

                if ('a3' %in% comb) {
                  a3_index <- as.integer(names(event_exons[['a3']][2]))
                  new_a3 <-
                    .draw_one_sample_safely((start(orig_template)[a3_index] + 1):(end(orig_template)[a3_index] - 1))
                }

                if ('a5' %in% comb) {
                  a5_index <- as.integer(names(event_exons[['a5']][1]))
                  if (!is.null(event_exons[['a3']]) &&
                      a5_index == a3_index) {
                    new_a5 <- ifelse(neg_strand,
                                     .draw_one_sample_safely((start(
                                       orig_template
                                     )[a5_index] + 1):(new_a3 - 1)),
                                     .draw_one_sample_safely((new_a3 + 1):(end(
                                       orig_template
                                     )[a5_index] - 1)))
                  }
                  else
                    new_a5 <-
                      .draw_one_sample_safely((start(orig_template)[a5_index] + 1):(end(orig_template)[a5_index] - 1))
                }
                variant <-
                  orig_template[.apply_exon_events(exon_vector_new, comb, event_exons)]
                if (two_variants) {
                  template_new <- orig_template[exon_vector_new]
                  template_new$transcript_id <-
                    sprintf('%s_%s_template', gene_id, paste(comb, collapse = ","))
                  S4Vectors::mcols(template_new)[c('tr_start', 'tr_end')] <-
                    .get_transcriptomic_coord(IRanges::width(template_new))
                }
              } else {
                variant <-
                  orig_template[.apply_exon_events(exon_vector, comb, event_exons)]
              }
              if ('a3' %in% comb) {
                if (neg_strand)
                  end(variant[variant$gene_exon_number == a3_index]) <- new_a3
                else
                  start(variant[variant$gene_exon_number == a3_index]) <- new_a3
              }
              if ('a5' %in% comb) {
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

              if (multi_events_per_exon && two_variants) {
                event_annotation <-
                  .get_event_annotation(
                    comb,
                    event_exons,
                    exon_vector_new,
                    neg_strand,
                    template_new,
                    variant
                  )
                variant <- c(.add_transcript_and_junction_lines(template_new), .add_transcript_and_junction_lines(variant))
              } else {
                event_annotation <-
                  .get_event_annotation(
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
