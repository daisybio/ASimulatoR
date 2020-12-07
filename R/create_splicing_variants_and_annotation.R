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
  junctions <- IRanges::gaps(variant, start = min(BiocGenerics::start(variant))) + 1
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
  variant <- sort(c(variant, junctions), decreasing = (S4Vectors::runValue(BiocGenerics::strand(variant)) == '-'))
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
  template[1]$tr_start <- NA
  template[1]$tr_end <- NA
  template
}

#' Internal function to create alternative splicing events and annotation
#'
#' This is not intended to be called directly;
#' instead it is meant to be called via \code{\link{simulate_alternative_splicing}}
#'
#' @param gtf_path path to the gtf/gff file from which splice variants are created
#' @param valid_chromosomes character vector. Only from these chromosomes splice variants are created. 
#' When used from \code{\link{create_splicing_variants_and_annotation}} chromosomes for which fasta files exist are used.
#' @param event_probs Named list/vector containing numerics corresponding
#'  to the probabilites to create the event(combination). 
#'  If \code{probs_as_freq} is \code{TRUE} \code{event_probs} correspond 
#'  to the relative frequency of occurences for the event(combination) and 
#'  in this case the sum of all frequencies has to be <=1.
#' @param outdir character, path to folder where simulated reads and all
#'   annotations should be written, with *no* slash at the end. By default,
#'   reads are written to current working directory.
#' @param ncores the number of cores to be utilized for parallel generation
#'   of splice variants.
#' @param novel_variants Numeric value between 0 and 1 indicating the percentage 
#'   of splicing variants that will not be written to an additional gtf file splicing_variants_novel.gtf.
#' @param write_gff Additionally to the gtf file containing the splice variants,
#'   a gff3 file with the same content will be printed to the outdir. 
#'   Default \code{TRUE}
#' @param max_genes The maximum number of genes/exon supersets to be included 
#'   in the process of splice variant creation. 
#' @param exon_junction_coverage Should the real coverage of exons, junctions 
#'   and retained introns should be written into a additional file?
#'   Default \code{TRUE}, which means, that this function returns an exon junction table
#' @param multi_events_per_exon Should it be possible to have more than one AS event 
#'   at the same exon if multiple variants are created for the same exon superset?
#'   !If this option is set to \code{TRUE}, there may occur unforeseen AS events 
#'   that are not documented in the event_annotation file!.
#'   Default \code{FALSE}
#' @param probs_as_freq Should \code{event_probs} be treated as relative frequencies instead of probabilities?
#'   Default \code{FALSE}
#' @param save_exon_superset should the exon_supersets be saved to .rda file?
#' Default \code{TRUE}
#'
#' @return if \code{exon_junction_coverage = TRUE} the exons, junctions and retained introns as data table in gtf style
#' @import rtracklayer
#' @importFrom stats runif setNames
#' @importFrom methods as
#' @importFrom pbmcapply pbmclapply
create_splicing_variants_and_annotation <-
  function(gtf_path,
           valid_chromosomes,
           event_probs,
           outdir,
           ncores,
           write_gff = F,
           max_genes = NULL,
           exon_junction_coverage = T,
           multi_events_per_exon = F,
           probs_as_freq = F,
           save_exon_superset = T,
           novel_variants = NULL) {

    ### create exon_superset ----
    exon_supersets <- get_exon_supersets(gtf_path, ncores, save_exon_superset)
    
    ### filter for supersets on valid chromosomes ---
    superset_chr <- sapply(exon_supersets, function(e) S4Vectors::runValue(GenomeInfoDb::seqnames(e)))
    exon_supersets <- exon_supersets[superset_chr %in% valid_chromosomes]
    if (length(exon_supersets) == 0){
      stop('no exon superset found on valid chromosomes. Do the fasta file names match the gtf/gff chromosome names?')
    }
    
    
    ### assign splicing variants with as events to supersets ----
    message('assign variants to supersets...')
    
    gene_lengths <- sapply(exon_supersets, length)
    total_gene_lengths <- cumsum(rev(table(gene_lengths)))
    nr_genes <- min(sum(gene_lengths > 1), max_genes)
    too_few_supersets <- T
    while (too_few_supersets){
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
      total_min_nr_exons <- cumsum(rev(table(min_nr_exons)))
      for (x in names(total_min_nr_exons)[names(total_min_nr_exons) != '0']) {
        if(total_min_nr_exons[[x]] > total_gene_lengths[[x]]){
            missing_genes <- total_min_nr_exons[[x]] - total_gene_lengths[[x]]
            nr_genes <- nr_genes - ceiling((nr_genes / total_min_nr_exons[[x]]) * missing_genes)
            message(sprintf('%d supersets with at least %s exons found but %d needed. Reducing total number of genes to %d.',
                            total_gene_lengths[[x]],
                            x,
                            total_min_nr_exons[[x]],
                            nr_genes))
            too_few_supersets <- T
            break
        }
        too_few_supersets <- F
      }
    }
    decr <- order(min_nr_exons, decreasing = T)
    min_nr_exons <- min_nr_exons[decr]
    construct_all_list <- construct_all_list[decr]
    drawn_genes <- character(nr_genes)
    for (i in 1:nr_genes) {
      # if (min_nr_exons[i] > 0) {
        drawn_genes[i] <-
          draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
        gene_lengths[[drawn_genes[i]]] <- -1
      # } else {
      #   #TODO: what to do if we did not draw an as event?
      #   drawn_genes[i] <-
      #     draw_one_sample_safely(names(gene_lengths)[gene_lengths >= min_nr_exons[i]])
      #   gene_lengths[[drawn_genes[i]]] <- -1
      # }
    }
    
    ### create splice variants and annotation ----
    message('create splicing variants and annotation. This may take a while...')
    all_variants_and_event_annotation <-
      pbmcapply::pbmclapply(1:nr_genes, function(i) {
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
            S4Vectors::runValue(BiocGenerics::strand(orig_template)) == '-'
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
          } else {
            mee_events <- c('afe', 'ale', 'mee')
            mee_events <- mee_events[mee_events %in% unique(unlist(event_combs))]
            # if no mee event should be constructed every event_comb is made from the same template vector
            # if there are mee_events we make sure that they are at the beginning and end of the exon vector 
            # and afterwards apply the normal function for the rest
            if (length(mee_events) != 0) {
              res_tmp <- construct_variant(mee_events, exon_vector, min_nr_exons_per_event, available_exons, orig_template, F, neg_strand, min_nr_exons[i], assign_mees_only = T)
              mee_exons <- res_tmp$event_exons$mee
              res_list <- lapply(event_combs, function(comb){
                construct_variant(comb, res_tmp$exon_vector, min_nr_exons_per_event, available_exons, orig_template, multi_events_per_exon, neg_strand,
                                  get_min_nr_exons(min_nr_exons_per_event, c(comb, mee_events)), mee_events, mee_exons)
              })
              exon_vector <- res_tmp$exon_vector
            } else {
              res_list <- lapply(event_combs, function(comb){
                construct_variant(comb, exon_vector, min_nr_exons_per_event, available_exons, orig_template, multi_events_per_exon, neg_strand, 
                                  get_min_nr_exons(min_nr_exons_per_event, comb, F))
              })
            }
          }
          template <- orig_template[exon_vector]
          S4Vectors::mcols(template)[c('tr_start', 'tr_end')] <-
            .get_transcriptomic_coord(IRanges::width(template))


          variants_and_event_annotation <-
            lapply(1:length(event_combs), function(comb_ind) {
              comb <- event_combs[[comb_ind]]
              
              ### create splice variants ----
              if (multi_events_per_exon) {
                res <- res_list[[comb_ind]]
              #   res <-
              #     construct_variant(
              #       comb,
              #       exon_vector,
              #       min_nr_exons_per_event,
              #       available_exons,
              #       orig_template,
              #       multi_events_per_exon,
              #       neg_strand,
              #       get_min_nr_exons(min_nr_exons_per_event, c(comb, mee_events), multi_events_per_exon)
              #     )
              #   # exon_vector <- res$exon_vector
              #   event_exons <- res$event_exons
              #   # variant <- res$variant
              }# else {
              event_exons <- res$event_exons
                variant <-
                  orig_template[apply_exon_events(exon_vector, comb, event_exons)]
              #}
              
              ### apply a3, a5 and ir
              if ('a3' %in% comb) {
                new_a3 <- res$new_a3
                a3_index <- res$a3_index
                if (neg_strand)
                  BiocGenerics::end(variant[variant$gene_exon_number == a3_index]) <- new_a3
                else
                  BiocGenerics::start(variant[variant$gene_exon_number == a3_index]) <- new_a3
              }
              if ('a5' %in% comb) {
                new_a5 <- res$new_a5
                a5_index <- res$a5_index
                if (neg_strand)
                  BiocGenerics::start(variant[variant$gene_exon_number == a5_index]) <- new_a5
                else
                  BiocGenerics::end(variant[variant$gene_exon_number == a5_index]) <- new_a5
              }
              if ('ir' %in% comb) {
                ir_indices <- event_exons[['ir']]
                IRanges::ranges(variant[variant$gene_exon_number == ir_indices[1]]) <-
                  range(IRanges::ranges(variant[variant$gene_exon_number %in% ir_indices]))
                variant[variant$gene_exon_number == ir_indices[2]] <-
                  NULL
                retaining_exon <- variant[variant$gene_exon_number == ir_indices[1]]
                flanking_exons <-
                  template[template$gene_exon_number %in% ir_indices]
                ri <-
                  GenomicRanges::gaps(flanking_exons, start = min(BiocGenerics::start(flanking_exons)))
                S4Vectors::mcols(ri) <-
                  S4Vectors::mcols(retaining_exon)
                ri$tr_start <- retaining_exon$tr_start +
                  ifelse(neg_strand,
                         (BiocGenerics::end(retaining_exon) - BiocGenerics::end(ri)),
                         BiocGenerics::start(ri) - BiocGenerics::start(retaining_exon))
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

              # if (multi_events_per_exon && res$two_variants) {
              #   event_annotation <-
              #     get_event_annotation(
              #       comb,
              #       event_exons,
              #       exon_vector,
              #       neg_strand,
              #       res$template_new,
              #       variant
              #     )
              #   variant <- c(.add_transcript_and_junction_lines(res$template_new), .add_transcript_and_junction_lines(variant))
              # } else {
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
              # }

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
            unlist(methods::as(variants_and_event_annotation[!event_annos], "GRangesList"),
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
      }, mc.cores = ncores, mc.set.seed = T)
    all_variants_and_event_annotation <- unlist(all_variants_and_event_annotation, recursive = F)
    event_annos <- endsWith(names(all_variants_and_event_annotation), 'event_annotation')
    event_annotation <- data.table::rbindlist(all_variants_and_event_annotation[event_annos])
    variants <- unlist(methods::as(all_variants_and_event_annotation[!event_annos], "GRangesList"), use.names = F)
    message('finished creating splicing variants and annotation')
    message('')

    ### exporting ----
    message('exporting gtf for read simulation...')

    gene_idx <- variants$type == "gene"
    transcript_idx <- variants$type == "transcript"
    exon_idx <- variants$type == "exon"
    rtracklayer::export(variants[gene_idx | transcript_idx | exon_idx],
                        file.path(outdir, 'splicing_variants.gtf'))
    if(!is.null(novel_variants)) {
      variant_ids <-
        variants[transcript_idx]$transcript_id[!endsWith(variants[transcript_idx]$transcript_id, 'template')]
      remove_from_novel <- 
        sample(variant_ids , floor(length(variant_ids)*novel_variants))
      rtracklayer::export(variants[(gene_idx | transcript_idx | exon_idx) 
                                   & (!variants$transcript_id %in% remove_from_novel)],
                          file.path(outdir, 'splicing_variants_novel.gtf'))
    }
    if (write_gff) {
      # create Id and Parent field before exporting
      variants$Parent <- character(length(variants))
      variants$ID <- character(length(variants))
      variants$ID[gene_idx] <- variants$gene_id[gene_idx]
      variants$Parent[gene_idx] <- NA
      variants$ID[transcript_idx] <- variants$transcript_id[transcript_idx]
      variants$Parent[transcript_idx] <- variants$gene_id[transcript_idx]
      variants$ID[exon_idx] <- sprintf('%s_x%s', variants$transcript_id[exon_idx], variants$gene_exon_number[exon_idx])
      variants$Parent[exon_idx] <- variants$transcript_id[exon_idx]
      # export gff3
      rtracklayer::export(variants[gene_idx | transcript_idx | exon_idx],
                          file.path(outdir, 'splicing_variants.gff3'))
      if(!is.null(novel_variants)) {
        rtracklayer::export(variants[(gene_idx | transcript_idx | exon_idx) 
                                     & (!variants$transcript_id %in% remove_from_novel)],
                            file.path(outdir, 'splicing_variants_novel.gff3'))
      }
      variants$ID <- NULL
      variants$Parent <- NULL
    }
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

    if (exon_junction_coverage) return(data.table::as.data.table(variants[!(gene_idx | transcript_idx)]))
    else invisible(NULL)
  }
