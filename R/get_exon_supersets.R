#' Create or load exon_supersets from gtf/gff file
#'
#' @param gtf_path path to the gtf/gff file from which exon supersets are created
#' @param ncores number of cores used to generate exon supersets in parallel.
#' Default \code{1}: no parallel computing.
#' @param save should the exon_supersets be saved to .rda file?
#' Default \code{TRUE}
#'
#' @export
#' @return exon supersets generated from gtf/gff file
#' @import rtracklayer
#' @importFrom pbmcapply pbmclapply
get_exon_supersets <-
  function(gtf_path, ncores=1L, save=TRUE) {
    exon_supersets_path <- sprintf('%s.exon_superset.rda', gtf_path)
    if (file.exists(exon_supersets_path)) {
      message('loading superset...')
      load(exon_supersets_path)
      message('finished loading superset')
      message('')
    } else {
      message('importing gtf/gff...')
      exon_supersets <- rtracklayer::import(gtf_path)
      message('finished importing gtf/gff')
      
      if (grepl('\\.gff3?$', gtf_path)){
        message('Found GFF: inferring gene_id for each exon...')
        exon_supersets$gene_id <- sub(pattern = "^gene:", replacement = "", x = exon_supersets$gene_id)
        tr_ids <- unique(unlist(exon_supersets[exon_supersets$type == 'exon']$Parent))
        tr2gene <- setNames(unlist(exon_supersets[exon_supersets$ID %in% tr_ids]$Parent), 
                            unlist(exon_supersets[exon_supersets$ID %in% tr_ids]$ID))
        exon_supersets[exon_supersets$type == 'exon']$gene_id <- 
          tr2gene[unlist(exon_supersets[exon_supersets$type == 'exon']$Parent)]
      }
      message('')
      
      exon_supersets <-
        exon_supersets[exon_supersets$type == 'exon']
      S4Vectors::mcols(exon_supersets) <- S4Vectors::mcols(exon_supersets)["gene_id"]
      exon_supersets <-
        split(exon_supersets, exon_supersets$gene_id)
      
      message('creating superset...')
      gc()
      exon_supersets <-
        pbmcapply::pbmclapply(exon_supersets, function(gene) {
          template <- GenomicRanges::reduce(gene)
          neg_strand <-
            S4Vectors::runValue(BiocGenerics::strand(gene)) == '-'
          if (neg_strand)
            template <- rev(template)
          gene_id <- gene$gene_id[1]
          transcript_id <- sprintf('%s_template', gene_id)
          S4Vectors::mcols(template) <-
            cbind(
              S4Vectors::DataFrame(
                source = "ASimulatoR",
                type = "exon",
                score = ".",
                phase = ".",
                gene_id = gene_id,
                transcript_id = transcript_id,
                template = T,
                gene_exon_number = 1L:length(template)
              ),
              .get_transcriptomic_coord(IRanges::width(template))
            )
          template
        }, mc.cores = ncores)
      message('finished creating superset')
      message('')
      if (save) {
        message('saving superset...')
        save(exon_supersets, file = exon_supersets_path)
        message('finished saving superset')
        message('')
      }
    }
    return(exon_supersets)
  }
