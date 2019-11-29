#' Core FDR estimation
#' 
#' @param temp list of m matrices (n x n)
#' @return n lists of n vectors of m elements
#' 
estimate_FDR <- function(real_scores, random_scores, invert.scores=FALSE, max.score=NA) {
	# invert.scores: by default (FALSE) assumes that LOWER scores are better
	real = real_scores
	random = random_scores
	if(invert.scores) {
		real = -real
		random = -random
	}
	# print(real)
	t = sort(unique(real), decreasing=FALSE)
	if(!is.na(max.score)) t = t[which(t<=max.score)]
	if(length(t)==0) t = sort(unique(real), decreasing=FALSE)[1]
	# print(t)
	# print((random))
	FPs = matrix(ncol=ncol(random), nrow=length(t))
	TPs = rep(0, length(t))
	FDR_profiles = apply(random, 2, function(x) ecdf(x))
	TPR_profiles = ecdf(real)
	for(y in 1:length(t)) {
		FPs[y,] = sapply(FDR_profiles, function(x) x(t[y]))
		TPs[y] =  TPR_profiles(t[y])
	}
	FPn = FPs * nrow(random)
	TPn = TPs * length(real)
	FPn_median = apply(FPn, 1, function(x) median(x))
	res = data.frame(t, TPn, FPn_median)
	res$FDR = res$FPn_median / res$TPn
	res$N = length(real)
	if(invert.scores) res$t = -res$t
	return(res)
}

#' Interaction type aware FDR estimation
#' 
#' @param temp list of m matrices (n x n)
#' @return n lists of n vectors of m elements
#' 
estimate_FDR_by_type <- function(results, real, random, var='int_type', FDR.thresh=0.1, invert.scores = FALSE) {
	# Divide by type
	classes = unique(results[,var])
	real_scores = lapply(classes, function(x) {
		# print(x)
		real = real[which(results[,var] == x)]
	})
	random_scores = lapply(classes, function(x) {
		# print(x)
		random = random[which(results[,var] == x),,drop=FALSE]
	})
	names(real_scores) = classes
	names(random_scores) = classes

	FDR = lapply(classes, function(x) {
		# print(x)
		estimate_FDR(real_scores[[x]], random_scores[[x]], invert.scores=invert.scores, max.score=NA)
	})
	names(FDR) = classes
	FDR_void = FDR
	FDR = lapply(names(FDR_void),function(x) {a = FDR_void[[x]]; a$name = x; return(a)})
	names(FDR) = names(FDR_void)
	FDR = lapply(FDR, function(x) { a = x$FDR; for(i in (length(a)-1):1) {a[i] = min(a[i], a[i+1])}; x$FDR_ist = x$FDR; x$FDR=a; return(x) })
	# FDR1 = do.call(rbind, FDR)
	thresh = do.call(rbind,lapply(FDR, function(a) {c = a[order(abs(a$FDR_fixed - FDR.thresh))[1],]}))
	return(list('FDR'=FDR, 'thresh'=thresh))
}

#' Get random P.values
#' 
#' @param temp list of m matrices (n x n)
#' @return n lists of n vectors of m elements
#' 
estimate_random_Pval <- function(r_MI, invert.scores = TRUE) {
	# Use invert.scores=TRUE if higher scores are better
	temp_r_MI = r_MI
	if(invert.scores) temp_r_MI = -r_MI # Invert sign to evaluate the correct P(x<=k)
	r_MI_ecdf = apply(temp_r_MI,1,ecdf) # build null model for each feature
	r_MI_p.value = t(sapply(1:nrow(r_MI), function(x) r_MI_ecdf[[x]](temp_r_MI[x,])))
	return(r_MI_p.value)
}

# #' Real Pvalue estimation .. ??
# #' 
# #' @param temp list of m matrices (n x n)
# #' @return n lists of n vectors of m elements
# #' 
# estimate_real_Pval <- function(MI, r_MI, invert.scores = TRUE) {
# 	# Use invert.scores=TRUE if higher scores are better
# 	temp_r_MI = r_MI
# 	if(invert.scores) {
# 		temp_r_MI = -r_MI # Invert sign to evaluate the correct P(x<=k)
# 		MI = -MI
# 	}
# 	r_MI_ecdf = apply(temp_r_MI,1,ecdf) # build null model for each feature
# 	MI_p.value = sapply(1:nrow(r_MI), function(x) r_MI_ecdf[[x]](MI[x]))
# 	return(MI_p.value)
# }

#' Determine optimal P value thresholds to have a 0.1 FDR
#' 
#' @param temp list of m matrices (n x n)
#' @return n lists of n vectors of m elements
#' 
FDR_threshold <- function(FDR) {
	thresh = do.call(rbind,lapply(FDR, function(a) {
			# thresh = lapply(x$FDR, function(a) {
			a1 = a
			a = a[which(a$t <= 0.1),,drop=FALSE]
			if(nrow(a)==0) {
				temp = a1[1,]
				temp['t'] = NA
				temp['TPn'] = NA
				temp['FPn_median'] = NA
				temp['FDR'] = NA
				return(temp)
			}
			# if(nrow(a)==0) return(c())
			b = which(a$FDR <= 0.1)
			if(length(b)>0) {
				best_single_pos = max(b)
			} else {
				b = which(a$FDR > 0.1)
				best_single_pos = b[1]
			}
				return(a[best_single_pos,])
			# c = min(a$FDR[b] - 0.1)

			# c = a[order(abs(a$FDR - 0.1))[1],]
			# })
		}))
	thresh = rbind(thresh, 'TOTAL' = c(NA, sum(thresh$TPn, na.rm=TRUE), sum(thresh$FPn_median,na.rm=TRUE), sum(thresh$FPn_median,na.rm=TRUE) / sum(thresh$TPn, na.rm=TRUE), sum(thresh$N, na.rm=TRUE), 'TOTAL'))
	thresh$t = as.numeric(thresh$t)
	FDR_threshold = thresh
	return(FDR_threshold)
}

#' Perform FDR analysis and determine optimal Pvalue threshold
#' 
#' @param alpi
#' @param r.MI
#' @param column='wMI_p.value'
#' @return FDR estimates, optimal thresholds, annotated alpi table
#' 
FDR_analysis <- function(alpi, r.MI, column='wMI_p.value') {
	r.MI_table = do.call(cbind,lapply(r.MI, function(x) x[cbind(alpi$SFE_1,alpi$SFE_2),drop=FALSE]))
	# print(dim(r.MI_table))
	rownames(r.MI_table) = rownames(alpi)
	r.MI_p.value = estimate_random_Pval(r.MI_table)
	void_global = estimate_FDR_by_type(alpi, alpi[,column], r.MI_p.value)
	FDR = void_global$FDR
	FDR_threshold = FDR_threshold(FDR)
	alpi$FDR = FALSE
	alpi[which(alpi[,column] <= FDR_threshold[alpi$int_type, 't']), 'FDR'] = TRUE
	return(list('FDR'=FDR, 'FDR_threshold'=FDR_threshold, 'alpi'=alpi))
}
