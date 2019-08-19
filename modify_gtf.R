library(rtracklayer)
library(polyester)
library(data.table)

### gtf manipulation ----

gtf <- import('ensembl/Homo_sapiens.GRCh38.97.gtf')
gtf_exons <- gtf[gtf$type == "exon"]
gtf_genes <- split(x = gtf_exons, f = gtf_exons$gene_id)
gtf_premRNA <- lapply(gtf_genes, function(gene){
  trs <- split(gene, gene$transcript_id)
  hyp_premRNA <- Reduce(union, trs)
  mdata <- mcols(gene)[1,]
  gene_id <- mdata$gene_id
  mdata[!as.vector(is.na(mdata))] <- NA
  mdata$gene_id <- gene_id
  mdata$type <- "exon"
  mdata$transcript_id <- paste0(gene_id, "_preMRNA")
  # TODO add metadata to preMRNA
  # TODO check hyp. premRNA once again
})

# TODO: introduce splicing


# # make data.table from gtf ----
# gtf_dt <- setDT(as.data.frame(gtf))

### simulate ----
# make gtf smaller for an example run take a few random transcripts from the Y-chromosome
y_gtf <- gtf[seqnames(gtf) == "Y"]
nr_transcripts <- 100
small_y_gtf <- y_gtf[y_gtf$transcript_id %in% sample(unique(y_gtf$transcript_id), nr_transcripts)]
small_y_gtf_path <- 'ensembl/Homo_sapiens.GRCh38.97.chrY.small.gtf'
export(small_y_gtf, small_y_gtf_path)
fold_change_mat <- matrix(c(rep(2, 5), rep(1, nr_transcripts - 5), 
                            rep(1, 5), rep(4, 5), rep(1, nr_transcripts - 10)), 
                          nrow = nr_transcripts)
out_dir <- 'readsY'
err_rate = 0
simulate_experiment(gtf = small_y_gtf_path, seqpath = 'ensembl/fasta',
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