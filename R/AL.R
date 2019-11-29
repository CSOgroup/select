#' Init AL object
#'
#' Create an Alteration Landscape (AL) object
#'
#' @param am = NULL The binary alteration matrix
#' @param samples.on.row = TRUE Are samples on rows?
#' @return An Alteration Landscape (AL) object with the gam 
#' 
#' @examples
#' am = matrix(rbinom(50, 1, 0.1), ncol=5, nrow=10)
#' new.AL(am)
#' 
#' @export
new.AL <- function(am, samples.on.row = TRUE) {
	if(is.null(am)) stop('Input am is NULL')
	al = list()
	al = list('am' = as.matrix(am))
	if(!is.matrix(al$am)) stop('Input am is not a matrix')
	if(!samples.on.row) al$am = t(al$am)
	if(is.null(rownames(al$am))) rownames(al$am) = paste('sample_', 1:nrow(al$am), sep='')
	if(is.null(colnames(al$am))) colnames(al$am) = paste('alteration_', 1:ncol(al$am), sep='')
	al$am = al$am > 0 # binarize!
	al$alterations = list()
	al$samples = list()
	if(is.null(al$alterations$alteration.class)) { al$alterations$alteration.class = rep('alteration', ncol(al$am)); names(al$alterations$alteration.class) = colnames(al$am) }
	if(is.null(al$samples$sample.class)) { al$samples$sample.class = rep('sample', nrow(al$am)); names(al$samples$sample.class) = rownames(al$am) }
	class(al) <- "AL";
	return(al)
}

#' Get sample/alteration blocks
#' 
#' @param al The alteration landscape
#' @return Classification of samples and alterations in blocks.
#' 
get.blocks <- function(al) {
	if(is.null(al$alterations$alteration.class)) { al$alterations$alteration.class = rep('alteration', ncol(al$am)); names(al$alterations$alteration.class) = colnames(al$am) }
	if(is.null(al$samples$sample.class)) { al$samples$sample.class = rep('sample', nrow(al$am)); names(al$samples$sample.class) = rownames(al$am) }
	alteration.class = al$alterations$alteration.class[colnames(al$am)]
	sample.class = al$samples$sample.class[rownames(al$am)]
	feature.blocks = lapply(unique(alteration.class), function(x) which(alteration.class==x))
	names(feature.blocks) = unique(alteration.class)
	sample.blocks = lapply(unique(sample.class), function(x) which(sample.class==x))
	names(sample.blocks) = unique(sample.class)
	missing_rows = which(! 1:nrow(al$am) %in% unlist(sample.blocks))
	missing_cols = which(! 1:ncol(al$am) %in% unlist(feature.blocks))
	return(list('sample.blocks' = sample.blocks, 'alteration.blocks' = feature.blocks, 'missing.samples' = missing_rows, 'missing.alterations' = missing_cols))
}

#' Filter AL
#' 
#' @param al The alteration landscape
#' @param params parameter list
#' @return An Alteration Landscape (AL) object with the filtered am 
#' @export
filter.al <- function(al, params=NULL) {
	if(is.null(params)) return(al)
# We need to perform the permutation considering all the features, but some of them will be ignored afterwards due to low freq count.
# However, some samples might end up with 0 alteration occurrences after removal of such features in the real data, but not in the randomized ones.
# This introduces a bias in the null model. To take this into account, (i) we need to pre-filter samples that will eventually be discarded.
# The alterantive (ii) is to keep all the samples, even if they have all 0s. Method (i) requires the list of final features to be kept. while method (ii) does not require any additional knowledge.
# Solution implemented here: get rid of all features and samples that do not reach minimum occurrence count before generating the null model
	samples_to_consider = 1:nrow(al$am)
	if(!is.null(params$remove.unknown.class.samples)) if(params$remove.unknown.class.samples == TRUE) {	
		temp = which(!is.na(al$samples$sample.class[rownames(al$am)]))
		samples_to_consider = samples_to_consider[which(samples_to_consider %in% temp)]	
	}
	al$am = al$am[samples_to_consider, , drop=FALSE]

	samples_to_consider = 1:nrow(al$am)
	if(!is.null(params$remove.0.samples)) if(params$remove.0.samples == TRUE) {	
		temp = which(rowSums(al$am) > 0)
		samples_to_consider = samples_to_consider[which(samples_to_consider %in% temp)]	
	}
	al$am = al$am[samples_to_consider, , drop=FALSE]

	features_to_consider = 1:ncol(al$am)
	if(!is.null(params$min.feature.support)) {
		temp = which(colSums(al$am) >= params$min.feature.support)
		features_to_consider = features_to_consider[which(features_to_consider %in% temp)]
	}
	if(!is.null(params$min.feature.freq)) {
		temp = which(colSums(al$am) >= (params$min.feature.freq * nrow(al$am)))
		features_to_consider = features_to_consider[which(features_to_consider %in% temp)]
	}
	al$am = al$am[, features_to_consider, drop=FALSE]

	samples_to_consider = 1:nrow(al$am)
	if(!is.null(params$remove.0.samples)) if(params$remove.0.samples == TRUE) {	
		temp = which(rowSums(al$am) > 0)
		samples_to_consider = samples_to_consider[which(samples_to_consider %in% temp)]	
	}
	al$am = al$am[samples_to_consider, , drop=FALSE]

	return(al)
}

#' Build naive veto list by vetoing interactions of CNAs within the same chromosome
#' 
#' @param pool list of CNAs in format AMP.chr..... or DEL.chr.....
#' @return list of pairwise vetos
#' 
#' @export
get_vetos <- function(pool) {
	
	pool.chr = sapply(pool, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
	pool.chr.l = lapply(unique(pool.chr), function(x) names(pool.chr[which(pool.chr==x)]))
	pool.chr.l = pool.chr.l[which(sapply(pool.chr.l, length)>1)]
	pool.chr.c = lapply(pool.chr.l, function(x) { a=t(combn(x, m=2)); a = rbind(a, a[,c(2,1)])})
	vetos = do.call(rbind, pool.chr.c)
	return(vetos)
}

