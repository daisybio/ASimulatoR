#internal function to determine details about the splice variant(s)
construct_variant <-
  function(comb,
           exon_vector,
           min_nr_exons_per_event,
           available_exons,
           orig_template,
           multi_events_per_exon,
           neg_strand,
           min_nr_exons = NULL) {
    
  result = list()
  exon_vector_tmp <- exon_vector
  comb_tmp <- comb
  min_nr_exons_per_event_tmp <- min_nr_exons_per_event
  two_variants <- F
  if ('afe' %in% comb) {
    two_variants <- T
    exon_vector[draw_one_sample_safely(c(1, 2))] <- F
    exon_vector_tmp <- exon_vector_tmp[-c(1, 2)]
    comb_tmp <-
      comb_tmp[comb_tmp != 'afe']
  }
  if ('ale' %in% comb) {
    two_variants <- T
    exon_vector[draw_one_sample_safely(c((length(exon_vector) - 1), length(exon_vector)))] <-
      F
    exon_vector_tmp <-
      exon_vector_tmp[-c(length(exon_vector_tmp) - 1, length(exon_vector_tmp))]
    comb_tmp <-
      comb_tmp[comb_tmp != 'ale']
  }
  if ('mes' %in% comb) {
    mee_max_skipped_exons <-
      available_exons - ifelse(
        multi_events_per_exon,
        get_min_nr_exons(
          min_nr_exons_per_event,
          paste(comb, collapse = ','),
          multi_events_per_exon
        ) + 2L,
        min_nr_exons)
    mee_nr_skipped_exons <-
      draw_one_sample_safely(2:mee_max_skipped_exons)
    min_nr_exons_per_event_tmp['mes'] <-
      mee_nr_skipped_exons + 2
  }
  
  # assign the chunks with as events to the exon_vector
  event_exons <-
    .assign_events(exon_vector_tmp, min_nr_exons_per_event_tmp[comb_tmp])
  
  if ('mee' %in% comb) {
    two_variants <- T
    exon_vector[draw_one_sample_safely(names(event_exons$mee)[2:3])] <-
      F
  }
  
  if ('a3' %in% comb) {
    result$a3_index <- as.integer(names(event_exons[['a3']][2]))
    result$new_a3 <-
      draw_one_sample_safely((start(orig_template)[result$a3_index] + 1):(end(orig_template)[result$a3_index] - 1))
  }
  
  if ('a5' %in% comb) {
    result$a5_index <- as.integer(names(event_exons[['a5']][1]))
    if (!is.null(event_exons[['a3']]) &&
        result$a5_index == result$a3_index) {
      result$new_a5 <- ifelse(neg_strand,
                       draw_one_sample_safely((start(
                         orig_template
                       )[result$a5_index] + 1):(result$new_a3 - 1)),
                       draw_one_sample_safely((result$new_a3 + 1):(end(
                         orig_template
                       )[result$a5_index] - 1)))
    }
    else
      result$new_a5 <-
        draw_one_sample_safely((start(orig_template)[result$a5_index] + 1):(end(orig_template)[result$a5_index] - 1))
  }
  if (multi_events_per_exon) {
    result$variant <-
      orig_template[apply_exon_events(exon_vector, comb, event_exons)]
    if (two_variants) {
      result$template_new <- orig_template[exon_vector]
      result$template_new$transcript_id <-
        sprintf('%s_%s_template', orig_template$gene_id[1], paste(comb, collapse = ","))
      S4Vectors::mcols(result$template_new)[c('tr_start', 'tr_end')] <-
        .get_transcriptomic_coord(IRanges::width(result$template_new))
    }
  }
  result$two_variants <- two_variants
  result$exon_vector <- exon_vector
  result$event_exons <- event_exons
  return(result)
}