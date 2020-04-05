draw_one_sample_safely <- function(v) {
  if (length(v) == 1)
    return(v[1])
  else
    return(sample(v, 1))
}