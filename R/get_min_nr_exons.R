# internal function to determine the minimal number of exons to create certain event combinations
get_min_nr_exons <-
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