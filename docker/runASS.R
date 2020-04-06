library(ass)
args = commandArgs(trailingOnly = TRUE)
input = args[1] # /nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99.fa
output = args[2]
max_genes = 125
multi_events_per_exon = T
prob_as_freq = T
params = list(
  ncores = 1,
  input_dir = input,
  event_probs = 
    setNames(c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125),
             c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')),
  max_genes = max_genes,
  outdir = sprintf(
    '%s/maxGenes%d_multiEventsPerExon%s_eventsAsFreq%s',
    output,
    max_genes,
    multi_events_per_exon,
    prob_as_freq
  ),
  multi_events_per_exon = multi_events_per_exon,
  prob_as_freq = prob_as_freq
)

do.call(simulate_alternative_splicing, params)