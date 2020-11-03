.assign_events <- function(v, min_nr_exons_per_event, max_attempts = 10) {
  # if (sum(min_nr_exons_per_event) - (length(min_nr_exons_per_event) - 1L) > length(v)) 
  #   stop('Error in assign events: Event lengths are bigger than vector. This is never supposed to happen.')
  # done with all assignments
  if (length(min_nr_exons_per_event) == 0) return(list())
  # assign event and split vector, try until it works
  # assign greedy after maximum number of attempts exceeded
  counter <- 0
  while (T) {
    counter <- counter + 1
    result <- list()
    event <- draw_one_sample_safely(min_nr_exons_per_event)
    # does not work
    if (event > length(v)) {
      stop('Error in assign events: Event length is bigger than vector. This is never supposed to happen.')
    }
    event_start <- ifelse(counter > max_attempts, 
                          1,
                          draw_one_sample_safely(1:(length(v) - event + 1)))
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
    if (sum(min_nr_exons_per_event_tmp[events_first]) - sum(events_first) + 1L > length(vs$v1) ||
        sum(min_nr_exons_per_event_tmp[!events_first]) - sum(!events_first) + 1L > length(vs$v2))
      next
    result_v1 <-
      .assign_events(vs$v1, min_nr_exons_per_event_tmp[events_first])
    result_v2 <-
      .assign_events(vs$v2, min_nr_exons_per_event_tmp[!events_first])
    if ((length(result_v1) + length(result_v2)) == length(min_nr_exons_per_event_tmp))
      return(c(result, result_v1, result_v2))
  }
}

#internal function to determine details about the splice variant(s)
construct_variant <-
  function(comb,
           exon_vector,
           min_nr_exons_per_event,
           available_exons,
           orig_template,
           multi_events_per_exon,
           neg_strand,
           min_nr_exons,
           mee_events = NULL,
           mee_exons = NULL,
           assign_mees_only = F) {
    
  result = list()
  exon_vector_tmp <- exon_vector
  comb_tmp <- comb
  min_nr_exons_per_event_tmp <- min_nr_exons_per_event
  # if multi_events_per_exon == T means that the mee events have already been asigned, we have to remove the affected exons
  if (multi_events_per_exon) {
    if ('afe' %in% mee_events) {
      exon_vector_tmp <- exon_vector_tmp[-c(1, 2)]
      comb_tmp <-
        comb_tmp[comb_tmp != 'afe']
    }
    if ('ale' %in% mee_events){
      exon_vector_tmp <-
        exon_vector_tmp[-c(length(exon_vector_tmp) - 1, length(exon_vector_tmp))]
      comb_tmp <-
        comb_tmp[comb_tmp != 'ale']
    }
  } else {
    if ('afe' %in% comb) {
      exon_vector[draw_one_sample_safely(c(1, 2))] <- F
      exon_vector_tmp <- exon_vector_tmp[-c(1, 2)]
      comb_tmp <-
        comb_tmp[comb_tmp != 'afe']
    }
    if ('ale' %in% comb) {
      exon_vector[draw_one_sample_safely(c((length(exon_vector) - 1), length(exon_vector)))] <-
        F
      exon_vector_tmp <-
        exon_vector_tmp[-c(length(exon_vector_tmp) - 1, length(exon_vector_tmp))]
      comb_tmp <-
        comb_tmp[comb_tmp != 'ale']
    }
  }
  if ('mes' %in% comb) {
    mes_max_skipped_exons <-
      available_exons + 2L - min_nr_exons
    mes_nr_skipped_exons <-
      draw_one_sample_safely(2:mes_max_skipped_exons)
    min_nr_exons_per_event_tmp['mes'] <-
      mes_nr_skipped_exons + 2
  }
  
  # assign the chunks with as events to the exon_vector
  # if multi_events_per_exon, make sure it is compatible with mee events
  if (multi_events_per_exon) {
    if ('mee' %in% mee_events) {
      comb_tmp <-
        comb_tmp[comb_tmp != 'mee']
      if (mee_exons[1] == names(exon_vector_tmp)[1]) 
        exon_vector_tmp <- exon_vector_tmp[-(1:(min_nr_exons_per_event_tmp[['mee']] - 1))]
      else 
        exon_vector_tmp <- exon_vector_tmp[-((length(exon_vector_tmp) - min_nr_exons_per_event_tmp[['mee']] + 2):length(exon_vector_tmp))]
    }
  }
  if (assign_mees_only & 'mee' %in% comb) {
    possible_ind <- list(names(exon_vector_tmp)[1:min_nr_exons_per_event_tmp[['mee']]], names(exon_vector_tmp)[(length(exon_vector_tmp) + 1 - min_nr_exons_per_event_tmp[['mee']]):length(exon_vector_tmp)])
    event_exons <- list(mee = as.integer(unlist(draw_one_sample_safely(possible_ind))))
  } else {
    event_exons <- .assign_events(as.integer(names(exon_vector_tmp)), min_nr_exons_per_event_tmp[comb_tmp])
  }
  
  if ('mee' %in% comb) {
    if (multi_events_per_exon) 
      event_exons$mee <- mee_exons
    else 
      exon_vector[draw_one_sample_safely(event_exons$mee[2:3])] <- F
  }
  
  
  if ('a3' %in% comb) {
    result$a3_index <- event_exons[['a3']][2]
    result$new_a3 <-
      draw_one_sample_safely((BiocGenerics::start(orig_template)[result$a3_index] + 1):(end(orig_template)[result$a3_index] - 1))
  }
  
  if ('a5' %in% comb) {
    result$a5_index <- event_exons[['a5']][1]
    if (!is.null(event_exons[['a3']]) &&
        result$a5_index == result$a3_index) {
      result$new_a5 <- ifelse(neg_strand,
                       draw_one_sample_safely((BiocGenerics::start(
                         orig_template
                       )[result$a5_index] + 1):(result$new_a3 - 1)),
                       draw_one_sample_safely((result$new_a3 + 1):(end(
                         orig_template
                       )[result$a5_index] - 1)))
    }
    else
      result$new_a5 <-
        draw_one_sample_safely((BiocGenerics::start(orig_template)[result$a5_index] + 1):(end(orig_template)[result$a5_index] - 1))
  }
  result$exon_vector <- exon_vector
  result$event_exons <- event_exons
  return(result)
}