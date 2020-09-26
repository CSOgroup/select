#' Derive TP and FP stats for APC based on FDR scores for a specific threshold
TP.FP.stat <- function(A,th) {
	A = A[which(A$APC > 1e-8),,drop=FALSE]
	t =th
	TN = ecdf(A[!A$wMI_p.value_FDR,,drop=FALSE]$APC)(t)
	FN = ecdf(A[A$wMI_p.value_FDR,,drop=FALSE]$APC)(t)
	TP = 1 - FN
	FP = 1 - TN
	F1 = 2 * TP / (2 * TP + FN + FP)
	MCC = ((TP * TN) - (FP * FN)) / sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
	return(c(t=t, 'TN' = TN, FN = FN, TP = TP, FP = FP, F1 = F1, MCC = MCC))
}

#' Derive TP and FP stats for APC based on FDR scores for a range of thresholds
TP.FP.stats <- function(A,res=1000) {
	A = A[which(A$APC > 1e-8),,drop=FALSE] # Do not check APC values proximal to 0 
	if(nrow(A)<1) return(NULL)
	if(nrow(A[A$wMI_p.value_FDR,])<1) return(NULL)
	if(nrow(A[!A$wMI_p.value_FDR,])<1) return(NULL)
	t1 = min(A$APC) + ((0 : res) * ((max(A$APC)-min(A$APC)) / res))
	t2 = unique(A[A$wMI_p.value_FDR,,drop=FALSE]$APC)
	t = sort(unique(c(t1,t2)))
	AA = data.frame('t' = t)
	AA$TN = ecdf(A[!A$wMI_p.value_FDR,,drop=FALSE]$APC)(AA$t)
	AA$FN = ecdf(A[A$wMI_p.value_FDR,,drop=FALSE]$APC)(AA$t)
	AA$TP = 1 - AA$FN
	AA$FP = 1 - AA$TN
	AA$F1 = 2 * AA$TP / (2 * AA$TP + AA$FN + AA$FP)
	AA$MCC = ((AA$TP * AA$TN) - (AA$FP * AA$FN)) / sqrt((AA$TP + AA$FP)*(AA$TP + AA$FN)*(AA$TN + AA$FP)*(AA$TN + AA$FN))
	return(as.vector(AA))
}

#' Establish which threshold is the best, using the 95\% of FP exclusion criterion
get_thresh.2 <- function(x, A) {
	TP.FP.thresh_v1 = x[x[,'FP'] <= 0.05,,drop=FALSE]
	TP.FP.thresh_v2 = x[x[,'FP'] >= 0.05,,drop=FALSE]
	if(is.null(x)) return(NA) # there are no solutions at all. return NA
	if(is.null(dim(x))) return(NA) # there are no solutions at all. return NA
	if(nrow(x)==0) return(NA) # there are no solutions at all. return NA
	if(nrow(TP.FP.thresh_v1)==0) return(x[nrow(x),]) # there are no good solutions. Set highest th possible 
	if(nrow(TP.FP.thresh_v2)==0) return(TP.FP.thresh_v1[1,]) # All solutions are good. Set threshold to minimum score of good solutions
	# print('here')
	# print(TP.FP.thresh_v1[1,'t'])
	# print(TP.FP.thresh_v2[nrow(TP.FP.thresh_v2),'t'])
	th = (TP.FP.thresh_v1[1,'t'] + TP.FP.thresh_v2[nrow(TP.FP.thresh_v2),'t']) / 2
	th = (TP.FP.thresh_v1[1,'t']) # + TP.FP.thresh_v2[nrow(TP.FP.thresh_v2),'t']) / 2
	th = TP.FP.thresh_v2[nrow(TP.FP.thresh_v2),'t'] # Set threshold at maximum score associated to FP > 0.05
	# print(th)
	# print(TP.FP.stat(A, th))
	return(TP.FP.stat(A, th))
	# return(TP.FP.thresh_v1[1,])
}

#' Calculate ASC score
#' 
#' @param alpi table of pairwise interactions
#' @return a dataframe with stats on thresholds and the selected threshold
#' 
establish_APC_threshold <- function(alpi) {
	if(nrow(alpi) < 25) {
		# If we have less that 25 interactions, the threshold estimation does not make sense.
		stats = NULL
		thresh = rep(NA, 7)
	} else {
		stats = TP.FP.stats(alpi, 10000)
		thresh = get_thresh.2(stats, alpi)
	}
	return(list('stats' = stats, thresh = thresh))

}
