# TODO: finish constructing function
# TODO: call variants variants  with gsub
# TODO: restructure  construction functions
construct_splice_variants <- function(gene_list, as_events, nr_cores){
  if(is.null(gene_list) | !class(gene_list) == "GRangesList") stop("input has to be GRangesList object")
  as_events <- strsplit(as_events, ",", fixed = T)[[1]]
  new_variants <- GRangesList(mclapply(gene_list, function(gene){
    splice_variant_1 <- create_hyp_premRNA(gene)
    splice_variant_2 <- splice_variant_1
    negative_strand <- runValue(strand(premRNA) == "-")
    # TODO: check if the as_events are compatible
    if("ATSS" %in% as_events){
      splice_variant_2 <- construct_ATSS_variant(splice_variant_2)
      
    } 
    if("ATTS" %in% as_events){
      splice_variant_2 <- construct_ATTS_variant(splice_variant_2)
      
    }
    if("MES" %in% as_events){
      
    }else if("ES" %in% as_events){
      
    }else if("IR" %in% as_events){
      
    }else if("MEE" %in% as_events){
      
    }
    if("A3" %in% as_events){
      
    }
    if("A5" %in% as_events){
      
    }
  }, mc.cores = nr_cores))
  
  
    
}

create_hyp_premRNA <- function(gene){
  if(is.null(gene) | !class(gene) == "GRanges") stop("input has to be GRanges object")
  trs <- split(gene, gene$transcript_id)
  hyp_premRNA <- Reduce(union, trs)
  mdata <- mcols(gene)[1,]
  if(!"gene_id" %in% names(mdata)) stop("no gene_id found for input")
  gene_id <- mdata$gene_id
  mdata[!as.vector(is.na(mdata))] <- NA
  mdata$gene_id <- gene_id
  mdata$type <- "exon"
  mdata$transcript_id <- paste0(gene_id, "_preMRNA")
  mcols(hyp_premRNA) <- mdata
  hyp_premRNA
}

construct_ES_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 2){
    exon_ind <- sample((2:(length(premRNA) - 1)), 1)
    splice_variant <- premRNA[-exon_ind]
    splice_variant
  } else {
    stop("not long enough to intodruce ES")
  }
}

construct_MES_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 3){
    first_exon_ind <- sample((2:(length(premRNA) - 2)), 1)
    nr_exons <- sample((2:(length(premRNA) - first_exon_ind)), 1)
    splice_variant <- premRNA[-(first_exon_ind:(first_exon_ind + nr_exons - 1))]
    splice_variant
  } else {
    stop("not long enough to intodruce MES")
  }
}

construct_IR_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 1){
    exon_ind <- sample((1:(length(premRNA) - 1)), 1)
    splice_variant <- premRNA[-(exon_ind + 1)]
    end(ranges(splice_variant[exon_ind])) <- end(ranges(premRNA[exon_ind + 1]))
    splice_variant
  } else {
    stop("not long enough to intodruce IR")
  }
}

# TODO: check if exon is longer than 3
construct_A5_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 1){
    negative_strand <- runValue(strand(premRNA) == "-")
    exon_ind <- ifelse(negative_strand, sample((2:length(premRNA)), 1), sample((1:(length(premRNA) -1)), 1))
    splice_variant <- premRNA
    new_site <- sample(seq((start(ranges(splice_variant[exon_ind])) + 3), (end(ranges(splice_variant[exon_ind]))), 3), 1)
    if(negative_strand){
      start(ranges(splice_variant[exon_ind])) <- new_site
    } else {
      end(ranges(splice_variant[exon_ind])) <- new_site - 1
    }
    splice_variant
  } else {
    stop("not long enough to intodruce A5")
  }
}

construct_A3_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 1){
    negative_strand <- runValue(strand(premRNA) == "-")
    exon_ind <- ifelse(negative_strand, sample((1:(length(premRNA) -1)), 1), sample((2:length(premRNA)), 1))
    splice_variant <- premRNA
    new_site <- sample(seq((start(ranges(splice_variant[exon_ind])) + 3), (end(ranges(splice_variant[exon_ind]))), 3), 1)
    if(negative_strand){
      end(ranges(splice_variant[exon_ind])) <- new_site - 1
    } else {
      start(ranges(splice_variant[exon_ind])) <- new_site
    }
    splice_variant
  } else {
    stop("not long enough to intodruce A3")
  }
}

construct_ATSS_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 1){
    negative_strand <- runValue(strand(premRNA) == "-")
    exon_ind <- ifelse(negative_strand, sample((2:length(premRNA)), 1), sample((1:(length(premRNA) -1)), 1))
    if(negative_strand){
      splice_variant <- premRNA[-(exon_ind:length(premRNA))]
    } else {
      splice_variant <- premRNA[-(1:exon_ind)]
    }
    splice_variant
  } else {
    stop("not long enough to intodruce ATSS")
  }
}

construct_ATTS_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 1){
    negative_strand <- runValue(strand(premRNA) == "-")
    exon_ind <- ifelse(negative_strand, sample((1:(length(premRNA) -1)), 1), sample((2:length(premRNA)), 1))
    if(negative_strand){
      splice_variant <- premRNA[-(1:exon_ind)]
    } else {
      splice_variant <- premRNA[-(exon_ind:length(premRNA))]
    }
    splice_variant
  } else {
    stop("not long enough to intodruce ATTS")
  }
}

construct_MEE_variant <- function(premRNA){
  if(is.null(premRNA) | !class(premRNA) == "GRanges") stop("input has to be GRanges object")
  if(length(premRNA) > 3){
    possible_exons <- (2:(length(premRNA) - 2))
    exon_ind <- ifelse(length(possible_exons) == 1, possible_exons[1], sample(possible_exons, 1))
    splice_variant_2 <- premRNA[-exon_ind]
    splice_variant_1 <- premRNA[-(exon_ind + 1)]
    GRangesList(splice_variant_1 = splice_variant_1, splice_variant_2 = splice_variant_2)
  } else {
    stop("not long enough to intodruce MEE")
  }
}

add_as_information <- function(variant, event, start, end){
  mdatavariant
}