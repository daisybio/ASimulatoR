#Rscript runASimulatoR /path/to/input_folder/ /path/to/output_folder/


### set parameters ----
args = commandArgs(trailingOnly = TRUE)

input = args[1] #path to the input folder. Note: contains the chromosomes fasta file and the genome annotation (in gtf or gff3 format)
output = args[2] #path to the output folder
seed = 142 #seed for a reproducible simulation
ncores = 8 #number of cores for the parallel simulation

multi_events_per_exon = T #T - allow multiple events per exon; F - allow only one event per exon
probs_as_freq = F #T - use fixed frequencies (should be summed up to 1); F - use probabilities

error_rate = 0.001 #sequencing error rate. 0.001 = 0.1%
readlen = 76 #read length
max_genes = NULL #number of genes. NULL - use all compatible exon supersets
seq_depth = 2e06 #sequencing depth i.e. the number of reads per sample
num_reps = c(1,0) #number of samples and groups. E.g., c(1,1) - two groups with 1 sample each 

#distribution of events. One event per gene; equal distribution
as_events = c('es','ir','a3','a5','mes','mee','ale','afe')
event_probs = rep(1/(length(as_events) + 1), length(as_events))
names(event_probs) = as_events


##other examples---
#one event per gene; custom distribution
#event_probs = setNames(c(0.060, 0.078, 0.098, 0.049, 0.022, 0.004), c('a5', 'a3', 'es', 'ir', 'mes', 'mee')) 

#two events per gene; equal distribution
#as_events = c('es','ir','a3','a5','mes','mee','ale','afe')
#as_combs = combn(as_events, 2, FUN = function(...) paste(..., collapse = ','))
#event_probs = rep(1/(length(as_combs) + 1), length(as_combs))
#names(event_probs) = as_combs

#custom combinations of events; custom distributio
#event_probs = setNames(c(0.060, 0.078, 0.098), c('a5,a3', 'a3,es,ir', 'es,mes'))


#name of the output folder
outdir = sprintf(
    '%s/%s_maxGenes%d_SeqDepth%g_errRate%g_readlen%d_multiEventsPerExon%s_probsAsFreq%s',
    output,
        gsub('[ ]+?', '-', Sys.time()),
    ifelse(is.null(max_genes), 0, max_genes),
    seq_depth,
    error_rate,
    readlen,
    multi_events_per_exon,
    probs_as_freq
  )

params = list(
  seed = seed,
  ncores = ncores,
  input_dir = input,
  event_probs = event_probs,
  outdir = outdir,
  seq_depth = seq_depth,
  max_genes = max_genes,
  error_rate = error_rate,
  readlen = readlen,
  multi_events_per_exon = multi_events_per_exon,
  probs_as_freq = probs_as_freq,
  num_reps = num_reps
)


### run simulator ----
library(ASimulatoR)
do.call(simulate_alternative_splicing, params) #simulate RNA-Seq datasets


### copy this script to output folder for reproducibility ----
rscript = sub('--file=', '', commandArgs()[4], fixed = T)
bn_rscript = basename(rscript)
file.copy(rscript, file.path(outdir, bn_rscript))
