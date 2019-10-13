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
  mcols(transcript)$'tr_end' <- cumsum(width(ranges(transcript)))
  mcols(transcript)$'tr_start' <- c(1L, (transcript$'tr_end' + 1L)[-length(transcript)])
  transcript$exon_id <- 1L:length(transcript)
  sort(transcript)
}

set_variant <- function(splice_variant, ase, gen_start, gen_end, tr_start, tr_end, negative_strand){
  if(is.null(splice_variant) | !(class(splice_variant) == 'GRanges')) 
    stop('input must be of class GRanges')
  cs_id <- splice_variant$transcript_id[1]
  substr(cs_id, nchar(cs_id), nchar(cs_id)) <- ifelse(endsWith(cs_id, '2'), '1', '2')
  event_annotation <- data.table(
    cs = cs_id,
    variant = splice_variant$transcript_id[1],
    genomic_start = gen_start,
    genomic_end = gen_end,
    transcriptomic_start = tr_start,
    transcriptomic_end = tr_end,
    event_annotation = ase
  )
  
  # # TODO: do not also keep info in gtf?
  # mcols(splice_variant)[[paste(ase, 'genomic_start', sep = '_')]] <- gen_start
  # mcols(splice_variant)[[paste(ase, 'genomic_end', sep = '_')]] <- gen_end
  # mcols(splice_variant)[[paste(ase, 'transcriptomic_start', sep = '_')]] <- tr_start
  # mcols(splice_variant)[[paste(ase, 'transcriptomic_end', sep = '_')]] <- tr_end
  
  mcols(splice_variant)$CS <- F
  return(list(variant = get_transcriptomic_coord(splice_variant, negative_strand), 
              event_annotation = event_annotation))
}

add_as_event <- function(variant_event_annotation, funct, negative_strand, ...){
  new_variant_event_annotation <- funct(variant_event_annotation$variant, negative_strand, ...)
  new_variant_event_annotation$event_annotation <- rbind(new_variant_event_annotation$event_annotation, variant_event_annotation$event_annotation)
  return(new_variant_event_annotation)
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
                       premRNA[exon_ind]$tr_start, premRNA[exon_ind]$tr_end, negative_strand))
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
                       ifelse(negative_strand, premRNA[last_exon_ind]$tr_start, premRNA[first_exon_ind]$tr_start), 
                       ifelse(negative_strand, premRNA[first_exon_ind]$tr_end, premRNA[last_exon_ind]$tr_end),
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
    tr_start <- (ifelse(negative_strand, premRNA[exon_ind + 1L]$tr_end, premRNA[exon_ind]$tr_end) + 1L)
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
    tr_first <- ifelse(a5, premRNA[exon_ind]$tr_end, premRNA[exon_ind]$tr_start)
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
                       ifelse(atss, 1L, premRNA[exon_ind]$tr_start), ifelse(atss, premRNA[exon_ind]$tr_end, max(premRNA$tr_end)), negative_strand))
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
    tr_start <- ifelse(negative_strand, premRNA[exon_ind + 2L]$tr_end, premRNA[exon_ind - 1L]$tr_end) + 1L
    splice_variant_1 <- set_variant(splice_variant_1, 'MEE', start(ranges(premRNA[exon_ind])), end(ranges(premRNA[exon_ind])),
                                    tr_start, (tr_start + width(ranges(premRNA[exon_ind])) - 1L), negative_strand)
    splice_variant_2 <- set_variant(splice_variant_2, 'MEE', start(ranges(premRNA[exon_ind + 1])), end(premRNA[exon_ind + 1]),
                                    tr_start, (tr_start + width(ranges(premRNA[exon_ind + 1L])) - 1L), negative_strand)
    return(list(variant = c(splice_variant_1$variant, splice_variant_2$variant), 
                event_annotation = rbind(splice_variant_1$event_annotation, splice_variant_2$event_annotation)))
  } else {
    stop('not long enough to introduce MEE')
  }
}

construct_splice_variants <- function(hypPreMRNA, as_events, negative_strand){
  if(is.null(hypPreMRNA) | !(class(hypPreMRNA) == 'GRanges')) 
    stop('input must be of class GRanges')
  
  splice_variant_1 <- sort(hypPreMRNA)
  mcols(splice_variant_1)$transcript_id <- gsub('preMRNA', 'variant1', mcols(splice_variant_1)$transcript_id, fixed = T)
  
  splice_variant_2 <- splice_variant_1
  mcols(splice_variant_2)$transcript_id <- gsub('variant1', 'variant2', mcols(splice_variant_2)$transcript_id, fixed = T)
  
  event_annotation <- data.table(
    cs = character(),
    variant = character(),
    genomic_start = integer(),
    genomic_end = integer(),
    transcriptomic_start = integer(),
    transcriptomic_end = integer(),
    event_annotation = character()
  )
  
  variant_event_annotation <- list(variant = splice_variant_2, event_annotation = event_annotation)
  
  # TODO: check if the as_events are compatible
  if('NO_AS' %in% as_events){
    return(list(variant = hypPreMRNA, event_annotation = event_annotation))} 
  if('MEE' %in% as_events){
    splice_variant_1 <- construct_MEE_variant(splice_variant_1, negative_strand)
    return(splice_variant_1)
  }
  if('ATSS' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_ATSS_ATTS_variant, negative_strand, T)
  }else if('ATTS' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_ATSS_ATTS_variant, negative_strand, F)
  }else if('MES' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_MES_variant, negative_strand)
  }else if('ES' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_ES_variant, negative_strand)
  }else if('IR' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_IR_variant, negative_strand)
  }
  if('A5' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_A5_A3_variant, negative_strand, T)
  }
  if('A3' %in% as_events){
    variant_event_annotation <- add_as_event(variant_event_annotation, construct_A5_A3_variant, negative_strand, F)
  }
  
  variant_event_annotation$variant <- c(splice_variant_1, variant_event_annotation$variant)
  return(variant_event_annotation)
}

check_region_coverage <- function(region_start, region_end, check_start, check_end){
  return(!((region_start > check_end) | (check_start > region_end)))
}

check_junction_coverage <- function(junction_start, junction_end, check_start, check_end){
  return((junction_start >= check_start) & (junction_end <= check_end))
}

check_incl_and_excl_reads <- function(event_annotation, tr_id, read_coordinates, incl, excl){
  if(event_annotation[, sum(variant == tr_id)] != 0){
    event <- event_annotation[variant == tr_id]
    if(event$event_annotation %in% c('MEE', 'IR')){
      cover <- sapply(read_coordinates, function(mate) 
        check_region_coverage(event$transcriptomic_start, event$transcriptomic_end, mate[1], mate[2]))
    } else {
      cover <- sapply(read_coordinates, function(mate)
        check_junction_coverage((event$transcriptomic_start - 1L), event$transcriptomic_start, mate[1], mate[2]))
    }
    if(sum(cover) > 0){
      event_annotation[variant == tr_id][[incl]] <- event_annotation[variant == tr_id][[incl]] + 1L
      if(event$event_annotation == 'MEE')
        event_annotation[cs == tr_id][[excl]] <- event_annotation[cs == tr_id][[excl]] + 1L
    }
  } else if(event_annotation[, sum(cs == tr_id)] != 0){
    event <- event_annotation[cs == tr_id]
    if (event$event_annotation == 'IR'){
      cover <- sapply(read_coordinates, function(mate)
        check_junction_coverage((event$transcriptomic_start - 1L), event$transcriptomic_start, mate[1], mate[2]))
    } else {
      cover <- sapply(read_coordinates, function(mate) 
        check_region_coverage(event$transcriptomic_start, event$transcriptomic_end, mate[1], mate[2]))
    }
    if(sum(cover) > 0){
      event_annotation[cs == tr_id][[excl]] <- event_annotation[cs == tr_id][[excl]] + 1L
    }
  }
  return(event_annotation)
}

check_exon_coverage <- function(exons, tr_id, read_coordinates, sample){
  #tr_exons <- exons[transcript_id == tr_id]
  #browser()
  #ids_2 <- unique(unlist(sapply(read_coordinates, function(mate) subjectHits(findOverlaps(IRanges(mate[1], mate[2]), IRanges(tr_exons$TR_start, tr_exons$TR_end))))))
  ids <- exons[transcript_id == tr_id][unique(unlist(lapply(read_coordinates, function(mate) subjectHits(findOverlaps(IRanges(mate[1], mate[2]), IRanges(tr_start, tr_end))))))]$exon_id
  exons[transcript_id == tr_id & exon_id %in% ids][[sample]] <- exons[transcript_id == tr_id & exon_id %in% ids][[sample]] + 1L
  exons
}

check_reads_and_convert_format <- function(in_file, event_annotation, exons){
  if(!class(in_file) == 'character'){
    stop('in_file argument cannot be converted to a file')
  }else{
    in_con <- file(in_file, open = 'rt')
  }
  splitted_filename <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', basename(in_file)), ' ')[[1]]
  fwreads <- startsWith(splitted_filename[2], '1')
  message(paste('creating fastq file for', in_file))
  out_file <- sub('.fasta$', '.fastq', in_file)
  out_con <- file(out_file, open = 'wt')
  if(fwreads){
    incl <- paste(splitted_filename[1], 'inclusion_reads', sep = '_')
    excl <- paste(splitted_filename[1], 'exclusion_reads', sep = '_')
    event_annotation[[incl]] <- 0L
    event_annotation[[excl]] <- 0L
    exons[[splitted_filename[1]]] <- 0L
  }
  while( length( one_line <- readLines( in_con , 1 ) ) > 0 ){
    if(startsWith(one_line, '>')){
      one_line <- sub(pattern = '^>', replacement = '@', x = one_line)
      if(fwreads){
        splitted_line <- strsplit(one_line, ';', T)[[1]]
        tr_id <- strsplit(splitted_line[1], '/', T)[[1]][2]
        read_coordinates <- lapply(strsplit(unlist(strsplit(splitted_line[-1], ':', T))[-c(1,3)], '-', T), as.numeric)
        event_annotation <- check_incl_and_excl_reads(event_annotation, tr_id, read_coordinates, incl, excl)
        # TODO: exon_resolution
        exons <- check_exon_coverage(exons, tr_id, read_coordinates, splitted_filename[1])
      }
    }else{
      one_line <- paste(one_line, '+', paste(rep('~', nchar(one_line)), collapse = ''), sep = '\n')
    }
    writeLines( one_line, out_con )
  }
  close( in_con ) ; close( out_con )
  return(list(event_annotation = event_annotation, exons = exons))
}
