#internal function to assign preset if present
assign_preset <- function(preset, event_probs){
  change_event_probs <- is.null(event_probs)
  if (is.null(preset) & change_event_probs) stop('you have to supply at least one of "preset", "event_probs"')
  if (is.null(preset)) return(list())
  stopifnot(is.character(preset) & (length(preset) == 1))
  utils::data(presets, package = 'ASimulatoR', envir = environment())
  preset <- match.arg(preset, names(presets))
  return(presets[[preset]])
}

# # this code was used to assign presets
# 
# presets <- list()
# error_rate <- 0.001
# readlen <- 76
# seq_depth <- 5e7
# num_reps <- c(1,1)
# as_events <- c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')
# event_probs <- rep(1/(length(as_events) + 1), length(as_events))
# names(event_probs) <- as_events
# presets$event_partition <- list(
#   event_probs = event_probs,
#   seq_depth = seq_depth,
#   max_genes = 20000,
#   error_rate = error_rate,
#   readlen = readlen,
#   probs_as_freq = T,
#   num_reps = num_reps
# )
# presets$rna_seq_experiment <- list(
#   event_probs = event_probs,
#   seq_depth = seq_depth,
#   max_genes = 20000,
#   error_model = 'illumina5',
#   pcr_rate = 0.001,
#   bias = 'cdnaf',
#   distr = 'empirical',
#   adapter_contamination = T,
#   readlen = readlen,
#   probs_as_freq = T,
#   num_reps = num_reps
# )
# as_combs <- combn(as_events, 2, FUN = function(...) paste(..., collapse = ','))
# event_probs <- rep(1/(length(as_combs) + 1), length(as_combs))
# presets$event_combination_2 <- list(
#   event_probs = event_probs,
#   seq_depth = seq_depth,
#   max_genes = 10000,
#   error_rate = error_rate,
#   readlen = readlen,
#   multi_events_per_exon = T,
#   num_reps = num_reps
# )
# save(presets, file = 'data/presets.rda')

#' Parameter presets for the \code{\link{simulate_alternative_splicing}} function
#'
#' @name presets
#' @docType data
#' @author Quirin Manz \email{quirin.manz@@tum.de}
#' @usage data(presets)
#' @keywords data
"presets"