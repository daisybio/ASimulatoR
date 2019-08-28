library(rtracklayer)
library(polyester)
library(data.table)
source("functions.R")

### importing gtf ----
gtf_path <- 'ensembl/Homo_sapiens.GRCh38.97.gtf'

message("importing gtf...")
gtf <- import(gtf_path)

# TODO: set seed
# TODO: check frameshift
# TODO: filter for genes with more than 3 exons
nr_genes <- 100
gtf <- gtf[gtf$gene_id %in% sample(unique(gtf$gene_id), nr_genes)]

gtf_exons <- gtf[gtf$type == "exon"]
gtf_genes <- split(x = gtf_exons, f = gtf_exons$gene_id)

### creating premRNA ----
nr_cores <- 16
gtf_premRNA <- GRangesList(mclapply(gtf_genes, create_premRNA, mc.cores = nr_cores))

# TODO: report nr_events (paramter), 
# coordinates in global/transcript, 
# psi values, 
# look at abundance,
# exon_usage


# # make data.table from gtf ----
# gtf_dt <- setDT(as.data.frame(gtf))

### simulate ----
# make gtf smaller for an example run
allowed_genes <- c("ENSG00000063177", "ENSG00000198755")
small_gtf <- gtf[gtf$gene_id %in% allowed_genes]
nr_transcripts <- length(small_gtf[small_gtf$type == "transcript"]$transcript_id)
small_gtf_path <- 'exampleRun/example.gtf'
export(small_gtf, small_gtf_path)
fold_change_mat <- matrix(c(rep(2, 5), rep(1, nr_transcripts - 5), 
                            rep(1, 5), rep(4, 5), rep(1, nr_transcripts - 10)), 
                          nrow = nr_transcripts)
out_dir <- 'exampleRun'
seq_path <- 'ensembl/fasta'
err_rate = 0
simulate_experiment(gtf = small_gtf_path, seqpath = seq_path,
                    fold_changes = fold_change_mat, outdir = out_dir, 
                    error_rate = err_rate, reads_per_transcript = rep(300, nr_transcripts))
out_files <- list.files(out_dir)
fasta_files <- file.path(out_dir, out_files[endsWith(out_files, ".fasta")])

fastqFromFasta <- function(in_file){
  if(!class(in_file) == "character"){
    stop("in_file argument cannot be converted to a file")
  }else{
    in_con <- file(in_file, open = "rt")
  }
  message(paste("creating fastq file for", in_file))
  out_file <- sub(".fasta$", ".fastq", in_file)
  out_con <- file(out_file, open = "wt")
  while( length( one_line <- readLines( in_con , 1 ) ) > 0 ){
    if(startsWith(one_line, ">")){
      one_line <- sub(pattern = "^>", replacement = "@", x = one_line)
    }else{
      one_line <- paste(one_line, "+", paste(rep("~", nchar(one_line)), collapse = ""), sep = "\n")
    }
    writeLines( one_line, out_con )
  }
  close( in_con ) ; close( out_con )
}

lapply(fasta_files, fastqFromFasta)


### Notes ----
# TODO: 1. combine all exons to the hypothetical pre-mRNA using the genomic ranges package!
# TODO: 2. introduce as-events from this hypothetical pre-mRNA
# TODO: 2.1. keep and express hypothetical pre-mRNA
# TODO: 3. write gtf to file
# TODO: 4. give gtf and fasta files to polyester
# TODO: 5. use rccp to modify fasta to fastq if not fast enough