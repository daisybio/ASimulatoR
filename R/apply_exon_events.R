# internal function to determine the exons, that are not used, based on the template
apply_exon_events <- function(v, comb, event_exons) {
  if ('afe' %in% comb)
    v[c(1, 2)] <- v[c(2, 1)]
  if ('ale' %in% comb)
    v[c(length(v) - 1, length(v))] <- v[c(length(v), length(v) - 1)]
  if ('mee' %in% comb)
    v[event_exons[['mee']][c(2, 3)]] <-
      v[event_exons[['mee']][c(3, 2)]]
  if ('es' %in% comb)
    v[event_exons[['es']][2]] <- F
  if ('mes' %in% comb)
    v[event_exons[['mes']][-c(1, length(event_exons[['mes']]))]] <-
      F
  return(v)
}