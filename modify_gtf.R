library(rtracklayer)
library(polyester)
library(data.table)

### gtf 

gtf <- import('ensembl/Homo_sapiens.GRCh38.97.gtf')

# make gtf smaller for an example run
y_gtf <- gtf[seqnames(gtf) == "Y"]
small_y_gtf <- y_gtf[y_gtf$transcript_id %in% sample(unique(y_gtf$transcript_id), 100)]
export(small_y_gtf, 'ensembl/Homo_sapiens.GRCh38.97.chrY.small.gtf')

# make data.frame from gtf
# gtf_dt <- setDT(as.data.frame(gtf))

# simulate with DE 
fold_change_mat <- matrix(c(4,4,1,1,1,1,1,1,4,4,1,1), nrow=6)
out_dir <- 'readsY'
err_rate = 0
simulate_experiment(gtf = 'ensembl/Homo_sapiens.GRCh38.97.chrY.small.gtf', seqpath = 'ensembl/fasta',
                    fold_changes = fold_change_mat, outdir = out_dir, error_rate = err_rate, strand_specific = T)
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



# highest value fastq
# TODO: 1. combine all exons to the hypothetical pre-mRNA genomic ranges!
# TODO: 2. delete some exons from this hypothetical pre-mRNA
# TODO: 2.1. keep and express hypothetical pre-mRNA
# TODO: 3. write gtf to file
# TODO: 4. give gtf and fasta files to polyester
# TODO: 5. rccp to modify fasta to fastq