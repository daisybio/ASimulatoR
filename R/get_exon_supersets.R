#' Title
#'
#' @param gtf_path
#' @param valid_chromosomes
#' @param ncores
#'
#' @return
#'
#' @examples
get_exon_supersets <-
  function(gtf_path, valid_chromosomes, ncores) {
    exon_supersets_path <- sprintf('%s..exon_superset.rda')
    if (file.exists(exon_supersets_path)) {
      message('loading superset...')
      load(exon_supersets_path)
      exon_supersets <-
        exon_supersets[as.character(sapply(exon_supersets, function(g)
          S4Vectors::runValue(GenomeInfoDb::seqnames(g))))
          %in% valid_chromosomes]
      message('finished loading superset')
      message('')
    } else {
      message('importing gtf...')
      exon_supersets <- rtracklayer::import(gtf_path)
      message('finished importing gtf')
      message('')

      exon_supersets <-
        exon_supersets[GenomeInfoDb::seqnames(exon_supersets) %in% valid_chromosomes &
                         exon_supersets$type == 'exon']
      exon_supersets <-
        split(exon_supersets, exon_supersets$gene_id)

      message('creating superset...')

      exon_supersets <-
        parallel::mclapply(exon_supersets, function(gene) {
          template <- GenomicRanges::reduce(gene)
          neg_strand <-
            S4Vectors::runValue(GenomicRanges::strand(gene)) == '-'
          if (neg_strand)
            template <- rev(template)
          gene_id <- gene$gene_id[1]
          transcript_id <- sprintf('%s_template', gene_id)
          S4Vectors::mcols(template) <-
            cbind(
              S4Vectors::DataFrame(
                source = "as_simulator",
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

      message('saving superset...')
      save(exon_supersets, file = exon_supersets_path)
      message('finished saving superset')
      message('')
    }
    return(exon_supersets)
  }
