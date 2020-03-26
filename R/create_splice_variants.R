# ### obsolete functions----
#
#
# draw_one_sample_safely <- function(v){
#   ifelse(length(v) == 1, v[1], sample(v, 1))
# }
#
# get_not_first <- function(n, first){
#   ifelse(first, draw_one_sample_safely((2:n)), draw_one_sample_safely((1:(n - 1))))
# }
#
# get_transcriptomic_coord <- function(transcript, negative_strand){
#   if (is.null(transcript) | !(class(transcript) == 'GRanges'))
#     stop('input must be of class GRanges')
#   if (negative_strand) transcript <- transcript[length(transcript):1]
#   transcript$tr_end <- cumsum(width(ranges(transcript)))
#   transcript$tr_start <- c(1L, (transcript$tr_end + 1L)[-length(transcript)])
#   transcript$exon_id <- paste(transcript$transcript_id, 'exon', (1L:length(transcript)), sep = '_')
#   # TODO: add transcript line
#   sort(transcript)
# }
#
# add_junctions <- function(transcript, negative_strand) {
#   if (length(transcript) == 1) return(transcript)
#   junction_ranges <- IRanges(end(ranges(transcript))[-length(transcript)], start(ranges(transcript))[-1])
#   junctions <- GRanges(rep(runValue(seqnames(transcript)), times = length(junction_ranges)),
#                        junction_ranges,
#                        rep(runValue(strand(transcript)), times = length(junction_ranges)))
#   if (negative_strand) {
#     transcript <- transcript[length(transcript):1]
#     junctions <- junctions[length(junctions):1]
#   }
#   mcols(junctions) <- mcols(transcript)[1:length(junctions), ]
#   mcols(junctions)$tr_start <- transcript$tr_end[-length(transcript)]
#   mcols(junctions)$tr_end <- transcript$tr_start[-1]
#   mcols(junctions)$type <- 'junction'
#   mcols(junctions)$exon_id <- #NULL mcols(junctions)$junction_id <-
#     paste(junctions$transcript_id, 'junction', (1L:length(junctions)), sep = '_')
#   merged_dt <- as.data.table(sort(c(junctions, transcript)))
#   return(merged_dt)
# }
#
#
# set_variant <- function(splice_variant, ase, gen_start, gen_end, tr_start, tr_end, negative_strand){
#   if (is.null(splice_variant) | !(class(splice_variant) == 'GRanges'))
#     stop('input must be of class GRanges')
#   cs_id <- splice_variant$transcript_id[1]
#   substr(cs_id, nchar(cs_id), nchar(cs_id)) <- ifelse(endsWith(cs_id, '2'), '1', '2')
#   event_annotation <- data.table(
#     cs = cs_id,
#     variant = splice_variant$transcript_id[1],
#     genomic_start = gen_start,
#     genomic_end = gen_end,
#     transcriptomic_start = tr_start,
#     transcriptomic_end = tr_end,
#     event_annotation = ase
#   )
#
#   # # TODO: do not also keep info in gtf?
#   # mcols(splice_variant)[[paste(ase, 'genomic_start', sep = '_')]] <- gen_start
#   # mcols(splice_variant)[[paste(ase, 'genomic_end', sep = '_')]] <- gen_end
#   # mcols(splice_variant)[[paste(ase, 'transcriptomic_start', sep = '_')]] <- tr_start
#   # mcols(splice_variant)[[paste(ase, 'transcriptomic_end', sep = '_')]] <- tr_end
#
#   mcols(splice_variant)$CS <- F
#   return(list(variant = get_transcriptomic_coord(splice_variant, negative_strand),
#               event_annotation = event_annotation))
# }
#
# add_as_event <- function(variant_event_annotation, funct, negative_strand, ...){
#   new_variant_event_annotation <- funct(variant_event_annotation$variant, negative_strand, ...)
#   new_variant_event_annotation$event_annotation <- rbind(new_variant_event_annotation$event_annotation,
#                                                          variant_event_annotation$event_annotation)
#   return(new_variant_event_annotation)
# }
#
# create_hyp_premRNA <- function(gene){
#   if (is.null(gene) | !class(gene) == 'GRanges') stop('input has to be GRanges object')
#   negative_strand <- runValue(strand(gene) == '-')
#   hyp_premRNA <- reduce(gene)
#   gene_id <- mcols(gene)$gene_id[1]
#   # TODO: add gene line
#   mcols(hyp_premRNA) <- DataFrame(source = "as_simulator",
#                                   type = "exon",
#                                   score = ".",
#                                   phase = ".",
#                                   gene_id = gene_id,
#                                   transcript_id = paste0(gene_id, '_preMRNA'),
#                                   exon_id = 0L,
#                                   CS = T)
#   hyp_premRNA <- get_transcriptomic_coord(hyp_premRNA, negative_strand)
#   return(hyp_premRNA)
# }
#
# create_exon_supersets <- function(gtf_path, nr_cores){
#   if (file.exists(paste0(gtf_path, '.exon_superset.rda'))) {
#     message('loading superset...')
#     load(paste0(gtf_path, '.exon_superset.rda'))
#     message('finished loading superset')
#   } else {
#     message('importing gtf...')
#     exon_supersets <- import(gtf_path)
#     message('finished importing gtf')
#     message('')
#     #TODO: make references
#     exon_supersets <- exon_supersets[exon_supersets$type == 'exon']
#
#     exon_supersets <- split(exon_supersets, exon_supersets$gene_id)
#
#     message('creating superset...')
#     exon_supersets <- mclapply(exon_supersets, function(gene){create_hyp_premRNA}, mc.cores = nr_cores)
#     message('finished creating superset')
#
#     message('saving...')
#     save(exon_supersets, file = paste0(gtf_path, '.exon_superset.rda'))
#     message('finished saving')
#   }
#   exon_supersets
# }
#
# construct_ES_variant <- function(premRNA, negative_strand){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 2) {
#     exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 1)))
#     splice_variant <- premRNA[-exon_ind]
#     return(set_variant(splice_variant, 'ES', start(ranges(premRNA[exon_ind])), end(ranges(premRNA[exon_ind])),
#                        premRNA[exon_ind]$tr_start, premRNA[exon_ind]$tr_end, negative_strand))
#   } else {
#     stop('not long enough to introduce ES')
#   }
# }
#
# construct_MES_variant <- function(premRNA, negative_strand){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 3) {
#     first_exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 2)))
#     nr_exons <- draw_one_sample_safely((2:(length(premRNA) - first_exon_ind)))
#     last_exon_ind <- (first_exon_ind + nr_exons - 1)
#     splice_variant <- premRNA[-(first_exon_ind:last_exon_ind)]
#     return(set_variant(splice_variant, 'MES', start(ranges(premRNA[first_exon_ind])), end(ranges(premRNA[(first_exon_ind + nr_exons - 1)])),
#                        ifelse(negative_strand, premRNA[last_exon_ind]$tr_start, premRNA[first_exon_ind]$tr_start),
#                        ifelse(negative_strand, premRNA[first_exon_ind]$tr_end, premRNA[last_exon_ind]$tr_end),
#                        negative_strand))
#   } else {
#     stop('not long enough to introduce MES')
#   }
# }
#
# construct_IR_variant <- function(premRNA, negative_strand){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 1) {
#     exon_ind <- draw_one_sample_safely((1:(length(premRNA) - 1)))
#     splice_variant <- premRNA[-(exon_ind + 1)]
#     end(ranges(splice_variant[exon_ind])) <- end(ranges(premRNA[exon_ind + 1]))
#     tr_start <- (ifelse(negative_strand, premRNA[exon_ind + 1L]$tr_end, premRNA[exon_ind]$tr_end) + 1L)
#     return(set_variant(splice_variant, 'IR', (end(ranges(premRNA[exon_ind])) + 1L), (start(ranges(premRNA[exon_ind + 1])) - 1L),
#                        tr_start, (tr_start + start(ranges(premRNA[exon_ind + 1L])) - end(ranges(premRNA[exon_ind])) - 2), negative_strand))
#   } else {
#     stop('not long enough to introduce IR')
#   }
# }
#
# # TODO: check if exon is longer than 3
# construct_A5_A3_variant <- function(premRNA, negative_strand, a5){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 1) {
#     exon_ind <- get_not_first(length(premRNA), (negative_strand == a5))
#     splice_variant <- premRNA
#     new_site <- draw_one_sample_safely(seq((start(ranges(splice_variant[exon_ind])) + 3), (end(ranges(splice_variant[exon_ind]))), 3))
#     if ((negative_strand == a5)) {
#       start(ranges(splice_variant[exon_ind])) <- new_site
#       cutout <- IRanges(start(ranges(premRNA[exon_ind])), (new_site - 1))
#     } else {
#       end(ranges(splice_variant[exon_ind])) <- new_site - 1
#       cutout <- IRanges(new_site, end(ranges(premRNA[exon_ind])))
#     }
#     tr_first <- ifelse(a5, premRNA[exon_ind]$tr_end, premRNA[exon_ind]$tr_start)
#     tr_second <- ifelse(a5, (tr_first - width(cutout) + 1L), (tr_first + width(cutout) - 1L))
#     return(set_variant(splice_variant, ifelse(a5, 'A5', 'A3'), start(cutout), end(cutout),
#                        ifelse(a5, tr_second, tr_first), ifelse(a5, tr_first, tr_second), negative_strand))
#   } else {
#     stop('not long enough to introduce A5/A3')
#   }
# }
#
# construct_AFE_ALE_variant <- function(premRNA, negative_strand, afe){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 2) {
#     if (negative_strand == afe) {
#       return(exclude_one_exon_per_variant(premRNA, negative_strand, (length(premRNA) - 1L)))
#     } else {
#       return(exclude_one_exon_per_variant(premRNA, negative_strand, 1L))
#     }
#     #return(set_variant(splice_variant, ifelse(atss, 'ATSS', 'ATTS'), start(cutout), end(cutout),
#     #                   ifelse(atss, 1L, premRNA[exon_ind]$tr_start), ifelse(atss, premRNA[exon_ind]$tr_end, max(premRNA$tr_end)), negative_strand))
#   } else {
#     stop('not long enough to introduce ATSS/ATTS')
#   }
# }
#
# construct_MEE_variant <- function(premRNA, negative_strand){
#   if (is.null(premRNA) | !class(premRNA) == 'GRanges') stop('input has to be GRanges object')
#   if (length(premRNA) > 3) {
#     exon_ind <- draw_one_sample_safely((2:(length(premRNA) - 2)))
#     return(exclude_one_exon_per_variant(premRNA, negative_strand, exon_ind))
#   } else {
#     stop('not long enough to introduce MEE')
#   }
# }
#
# exclude_one_exon_per_variant <- function(premRNA, negative_strand, exon_ind){
#   as_name <- ifelse((exon_ind > 1 & exon_ind < (length(premRNA) - 1L)), 'MEE',
#                     ifelse((exon_ind == (length(premRNA) - 1L)) == negative_strand, 'AFE', 'ALE'))
#   splice_variant_2 <- premRNA[-exon_ind]
#   splice_variant_2$transcript_id <- gsub('variant1', 'variant2', splice_variant_2$transcript_id, fixed = T)
#   splice_variant_1 <- premRNA[-(exon_ind + 1)]
#   tr_start <- ifelse(as_name == 'AFE', 1L,
#                      ifelse(negative_strand, premRNA[exon_ind + 2L]$tr_end, premRNA[exon_ind - 1L]$tr_end) + 1L)
#   splice_variant_1 <- set_variant(splice_variant_1, as_name, start(ranges(premRNA[exon_ind])), end(ranges(premRNA[exon_ind])),
#                                   tr_start, (tr_start + width(ranges(premRNA[exon_ind])) - 1L), negative_strand)
#   splice_variant_2 <- set_variant(splice_variant_2, as_name, start(ranges(premRNA[exon_ind + 1])), end(premRNA[exon_ind + 1]),
#                                   tr_start, (tr_start + width(ranges(premRNA[exon_ind + 1L])) - 1L), negative_strand)
#   return(list(variant = rbind(add_junctions(splice_variant_1$variant, negative_strand), add_junctions(splice_variant_2$variant, negative_strand)),
#               event_annotation = rbind(splice_variant_1$event_annotation, splice_variant_2$event_annotation)))
# }
#
# construct_splice_variants <- function(exon_vector, events){
#   # if (is.null(hypPreMRNA) | !(class(hypPreMRNA) == 'GRanges'))
#   #   stop('input must be of class GRanges')
#
#   if (events['afe']) {
#     exon_vector[1] <- ifelse(exon_vector[1] == 1, 2, 1)
#   }
#   if (events['ale']) {
#     exon_vector[length(exon_vector)] <-
#       ifelse(exon_vector[length(exon_vector)]
#              == exon_vector[length(exon_vector) - 1] + 1,
#              exon_vector[length(exon_vector)] + 1,
#              exon_vector[length(exon_vector)] - 1)
#   }
#   if (events['mee']) {
#
#   }
#
#   browser()
#   negative_strand <- runValue(strand(hypPreMRNA) == '-')
#
#   #splice_variant_1 <- sort(hypPreMRNA)
#   splice_variant <- hypPreMRNA
#   splice_variant$transcript_id <- gsub('preMRNA', paste(as_events, sep = ','), splice_variant$transcript_id, fixed = T)
#   #splice_variant_1$exon_id <- gsub('preMRNA', 'variant1', splice_variant_1$exon_id, fixed = T)
#
#
#   event_annotation <- data.table(
#     cs = character(),
#     variant = character(),
#     genomic_start = integer(),
#     genomic_end = integer(),
#     transcriptomic_start = integer(),
#     transcriptomic_end = integer(),
#     event_annotation = character()
#   )
#
#   variant_event_annotation <- list(variant = splice_variant, event_annotation = event_annotation)
#
#   # construct_two_variants <- any(c('MEE', 'ALE', 'AFE') %in% as_events)
#   # if (construct_two_variants) {
#   #   splice_variant_2 <- splice_variant
#   #   splice_variant_2$transcript_id <- paste(splice_variant_2$transcript_id, '2', sep = '_')
#   # }
#
#   # TODO: check if the as_events are compatible
#   # if ('CS' %in% as_events) {
#   #   return(list(variant = add_junctions(hypPreMRNA, negative_strand), event_annotation = event_annotation))
#   # } else
#   if ('MEE' %in% as_events) {
#     splice_variant <- construct_MEE_variant(splice_variant, negative_strand)
#     return(splice_variant)
#   } else if ('AFE' %in% as_events) {
#     splice_variant <- construct_AFE_ALE_variant(splice_variant, negative_strand, T)
#     return(splice_variant)
#   } else if ('ALE' %in% as_events) {
#     splice_variant <- construct_AFE_ALE_variant(splice_variant, negative_strand, F)
#     return(splice_variant)
#   }
#   if ('MES' %in% as_events) {
#     variant_event_annotation <- add_as_event(variant_event_annotation, construct_MES_variant, negative_strand)
#   }
#   if ('ES' %in% as_events) {
#     variant_event_annotation <- add_as_event(variant_event_annotation, construct_ES_variant, negative_strand)
#   }
#   if ('IR' %in% as_events) {
#     variant_event_annotation <- add_as_event(variant_event_annotation, construct_IR_variant, negative_strand)
#   }
#   if ('A5' %in% as_events) {
#     variant_event_annotation <- add_as_event(variant_event_annotation, construct_A5_A3_variant, negative_strand, T)
#   }
#   if ('A3' %in% as_events) {
#     variant_event_annotation <- add_as_event(variant_event_annotation, construct_A5_A3_variant, negative_strand, F)
#   }
#
#   variant_event_annotation$variant <- add_junctions(splice_variant, negative_strand)
#   if (construct_two_variants) {
#
#   }
#   return(variant_event_annotation)
# }
