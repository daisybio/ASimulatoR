draw_one_sample_safely <- function(v){
  ifelse(length(v) == 1, v[1], sample(v, 1))
}

get_not_first <- function(n, first){
  ifelse(first, draw_one_sample_safely((2:n)), draw_one_sample_safely((1:(n - 1))))
}

get_transcriptomic_coord <- function(transcript, negative_strand){
  if(is.null(transcript) | !(class(transcript) == 'GRanges')) 
    stop('input must be of class GRanges')
  if(negative_strand) transcript <- transcript[length(transcript):1]
  mcols(transcript)$"TR_end" <- cumsum(width(ranges(transcript)))
  mcols(transcript)$"TR_start" <- c(1, (transcript$'TR_end' + 1)[-length(transcript)])
  sort(transcript)
}

set_variant <- function(splice_variant, ase, start, end, tr_start, tr_end, negative_strand){
  if(is.null(splice_variant) | !(class(splice_variant) == 'GRanges')) 
    stop('input must be of class GRanges')
  mcols(splice_variant)[[paste(ase, 'genomic_start', sep = '_')]] <- start
  mcols(splice_variant)[[paste(ase, 'genomic_end', sep = '_')]] <- end
  mcols(splice_variant)[[paste(ase, 'transcriptomic_start', sep = '_')]] <- tr_start
  mcols(splice_variant)[[paste(ase, 'transcriptomic_end', sep = '_')]] <- tr_end
  mcols(splice_variant)$CS <- F
  return(get_transcriptomic_coord(splice_variant, negative_strand))
}

create_hyp_premRNA <- function(gene, negative_strand){
  if(is.null(gene) | !class(gene) == 'GRanges') stop('input has to be GRanges object')
  hyp_premRNA <- reduce(gene)
  mdata <- mcols(gene)[1,]
  if(!'gene_id' %in% names(mdata)) stop('no gene_id found for input')
  gene_id <- mdata$gene_id
  mdata[!as.vector(is.na(mdata))] <- NA
  mdata$gene_id <- gene_id
  mdata$type <- 'exon'
  mdata$transcript_id <- paste0(gene_id, '_preMRNA')
  mdata$CS <- T
  mcols(hyp_premRNA) <- mdata
  hyp_premRNA <- get_transcriptomic_coord(hyp_premRNA, negative_strand)
  return(hyp_premRNA)
}

construct_ES_variant <- function(premRNA, negative_strand){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 2){
    exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 1)))
    splice_variant <- premRNA[-exon_ind]
    return(set_variant(splice_variant, 'ES', start(ranges(premRNA[exon_ind])), end(ranges(premRNA[exon_ind])),
                       premRNA[exon_ind]$TR_start, premRNA[exon_ind]$TR_end, negative_strand))
  } else {
    stop('not long enough to introduce ES')
  }
}

construct_MES_variant <- function(premRNA, negative_strand){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 3){
    first_exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 2)))
    nr_exons <- draw_one_sample_safely((2:(length(premRNA) - first_exon_ind)))
    last_exon_ind <- (first_exon_ind + nr_exons - 1)
    splice_variant <- premRNA[-(first_exon_ind:last_exon_ind)]
    return(set_variant(splice_variant, 'MES', start(ranges(premRNA[first_exon_ind])), end(ranges(premRNA[(first_exon_ind + nr_exons - 1)])),
                        ifelse(negative_strand, premRNA[last_exon_ind]$TR_start, premRNA[first_exon_ind]$TR_start), 
                        ifelse(negative_strand, premRNA[first_exon_ind]$TR_end, premRNA[last_exon_ind]$TR_end),
                         negative_strand))
  } else {
    stop('not long enough to introduce MES')
  }
}

construct_IR_variant <- function(premRNA, negative_strand){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 1){
    exon_ind <- draw_one_sample_safely((1:(length(premRNA) - 1)))
    splice_variant <- premRNA[-(exon_ind + 1)]
    end(ranges(splice_variant[exon_ind])) <- end(ranges(premRNA[exon_ind + 1]))
    tr_start <- (ifelse(negative_strand, premRNA[exon_ind + 1L]$TR_end, premRNA[exon_ind]$TR_end) + 1L)
    #browser() # TODO: stuck here
    return(set_variant(splice_variant, 'IR', (end(ranges(premRNA[exon_ind])) + 1L), (start(ranges(premRNA[exon_ind + 1])) - 1L), 
                       tr_start, (tr_start + start(ranges(premRNA[exon_ind + 1L])) - end(ranges(premRNA[exon_ind])) - 2), negative_strand))
  } else {
    stop('not long enough to introduce IR')
  }
}

# TODO: check if exon is longer than 3
construct_A5_A3_variant <- function(premRNA, negative_strand, a5){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 1){
    exon_ind <- get_not_first(length(premRNA), (negative_strand == a5))
    splice_variant <- premRNA
    new_site <- draw_one_sample_safely(seq((start(ranges(splice_variant[exon_ind])) + 3), (end(ranges(splice_variant[exon_ind]))), 3))
    if((negative_strand == a5)){
      start(ranges(splice_variant[exon_ind])) <- new_site
      cutout <- IRanges(start(ranges(premRNA[exon_ind])), (new_site - 1))
    } else {
      end(ranges(splice_variant[exon_ind])) <- new_site - 1
      cutout <- IRanges(new_site, end(ranges(premRNA[exon_ind])))
    }
    tr_first <- ifelse(a5, premRNA[exon_ind]$TR_end, premRNA[exon_ind]$TR_start)
    tr_second <- ifelse(a5, (tr_first - width(cutout) + 1L), (tr_first + width(cutout) - 1L))
    return(set_variant(splice_variant, ifelse(a5, 'A5', 'A3'), start(cutout), end(cutout), 
                       ifelse(a5, tr_second, tr_first), ifelse(a5, tr_first, tr_second), negative_strand))
  } else {
    stop('not long enough to introduce A5/A3')
  }
}

construct_ATSS_ATTS_variant <- function(premRNA, negative_strand, atss){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 1){
    exon_ind <- get_not_first(length(premRNA), (negative_strand == atss))
    if((negative_strand == atss)){
      splice_variant <- premRNA[-(exon_ind:length(premRNA))]
      cutout <- IRanges(start(ranges(premRNA[exon_ind])), end(ranges(premRNA[length(premRNA)])))
    } else {
      splice_variant <- premRNA[-(1:exon_ind)]
      cutout <- IRanges(start(ranges(premRNA[1])), end(ranges(premRNA[exon_ind])))
    }
    return(set_variant(splice_variant, ifelse(atss, 'ATSS', 'ATTS'), start(cutout), end(cutout),
                      ifelse(atss, 1L, premRNA[exon_ind]$TR_start), ifelse(atss, premRNA[exon_ind]$TR_end, max(premRNA$TR_end)), negative_strand))
  } else {
    stop('not long enough to introduce ATSS/ATTS')
  }
}

construct_MEE_variant <- function(premRNA, negative_strand){
  if(is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
  if(length(premRNA) > 3){
    exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 2)))
    splice_variant_2 <- premRNA[-exon_ind]
    mcols(splice_variant_2)$transcript_id <- gsub('variant1', 'variant2', mcols(splice_variant_2)$transcript_id, fixed = T)
    splice_variant_1 <- premRNA[-(exon_ind + 1)]
    tr_start <- ifelse(negative_strand, premRNA[exon_ind + 2L]$TR_end, premRNA[exon_ind - 1L]$TR_end) + 1L
    return(c(set_variant(splice_variant_1, 'MEE', start(ranges(premRNA[exon_ind])), end(ranges(premRNA[exon_ind])),
                         tr_start, (tr_start + width(ranges(premRNA[exon_ind])) - 1L), negative_strand),
             set_variant(splice_variant_2, 'MEE', start(ranges(premRNA[exon_ind + 1])), end(premRNA[exon_ind + 1]),
                         tr_start, (tr_start + width(ranges(premRNA[exon_ind + 1L])) - 1L), negative_strand)))
  } else {
    stop('not long enough to introduce MEE')
  }
}

construct_splice_variants <- function(hypPreMRNA, as_events, negative_strand){
  if(is.null(hypPreMRNA) | !(class(hypPreMRNA) == 'GRanges')) 
    stop('input must be of class GRanges')
  # event_coordinates <- c('ES_start', 'ES_end', 'MES_start', 'MES_end', 'IR_start', 'IR_end', 'A5_start', 'A5_end', 'A3_start', 'A3_end',
  #                        'ATSS_start', 'ATSS_end', 'ATTS_start', 'ATTS_end', 'MEE_start', 'MEE_end')
  # mcols(hypPreMRNA)[event_coordinates] <- NA
  splice_variant_1 <- sort(hypPreMRNA)
  mcols(splice_variant_1)$transcript_id <- gsub('preMRNA', 'variant1', mcols(splice_variant_1)$transcript_id, fixed = T)
  
  splice_variant_2 <- splice_variant_1
  mcols(splice_variant_2)$transcript_id <- gsub('variant1', 'variant2', mcols(splice_variant_2)$transcript_id, fixed = T)
  
  # event_annotation <- data.table(
  #   cs = character(),
  #   variant = character(),
  #   affected_genomic_region = integer(),
  #   affected_transcriptomic_region = integer(),
  #   event_annotation = character(),
  #   inclusion_reads = integer(),
  #   exclusion_reads = integer()
  # )
  # TODO: check if the as_events are compatible
  if('NO_AS' %in% as_events){
    return(hypPreMRNA)} 
  if('MEE' %in% as_events){
    splice_variant_1 <- construct_MEE_variant(splice_variant_1, negative_strand)
    return(splice_variant_1)
  }
  if('ATSS' %in% as_events){
    splice_variant_2 <- construct_ATSS_ATTS_variant(splice_variant_2, negative_strand, T)
  }else if('ATTS' %in% as_events){
    splice_variant_2 <- construct_ATSS_ATTS_variant(splice_variant_2, negative_strand, F)
  }else if('MES' %in% as_events){
    splice_variant_2 <- construct_MES_variant(splice_variant_2, negative_strand)
  }else if('ES' %in% as_events){
    splice_variant_2 <- construct_ES_variant(splice_variant_2, negative_strand)
  }else if('IR' %in% as_events){
    splice_variant_2 <- construct_IR_variant(splice_variant_2, negative_strand)
  }
  if('A5' %in% as_events){
    splice_variant_2 <- construct_A5_A3_variant(splice_variant_2, negative_strand, T)
  }
  if('A3' %in% as_events){
    splice_variant_2 <- construct_A5_A3_variant(splice_variant_2, negative_strand, F)
  }
  return(list("variants" = c(splice_variant_1, splice_variant_2), 'annotation' = event_annotation))
}

fastqFromFasta <- function(in_file){
  if(!class(in_file) == 'character'){
    stop('in_file argument cannot be converted to a file')
  }else{
    in_con <- file(in_file, open = 'rt')
  }
  message(paste('creating fastq file for', in_file))
  out_file <- sub('.fasta$', '.fastq', in_file)
  out_con <- file(out_file, open = 'wt')
  while( length( one_line <- readLines( in_con , 1 ) ) > 0 ){
    if(startsWith(one_line, '>')){
      one_line <- sub(pattern = '^>', replacement = '@', x = one_line)
    }else{
      one_line <- paste(one_line, '+', paste(rep('~', nchar(one_line)), collapse = ''), sep = '\n')
    }
    writeLines( one_line, out_con )
  }
  close( in_con ) ; close( out_con )
}