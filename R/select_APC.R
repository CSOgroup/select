#' Convert a list of m matrices (n x n) to n lists of n vectors of m elements
#' 
#' @param temp list of m matrices (n x n)
#' @return n lists of n vectors of m elements
#' 
listmat_to_elemvec <- function(temp) {
	temp1 = lapply(1:nrow(temp[[1]]), function(pos) temp1 = lapply(temp, function(x) x[,pos]))
	names(temp1) = rownames(temp[[1]])
	temp2 = lapply(temp1, function(x) do.call(rbind,x))
	temp3 = temp2 # lapply(temp2, function(x) {a=lapply(1:ncol(x), function(y) x[,y]); names(a) = colnames(x); return(a)})
	return(temp3)
}


background_wMI <- function(r.MI, probs = c(0.25,0.5,0.75)) { # Fast and memory efficiend version
	temp1 = sapply(1:nrow(r.MI[[1]]), function(pos) {
		temp1 = lapply(r.MI, function(x) x[,pos])
		temp2 = do.call(rbind,temp1)
		result = apply(matrixStats:::colQuantiles(temp2, probs = probs),1,max)
		result
	})
	rownames(temp1) = rownames(r.MI[[1]])
	colnames(temp1) = colnames(r.MI[[1]])
	return(temp1)
}

#' Remove interactions with observed overlap strictly similar to expected overlap
#' 
#' @param results table of pairwise interactions
#' @return a dataframe without invariant interactions
#' 
prune_invariants <- function(results) {
	ov_balance = abs(results$overlap - results$r_overlap)
	ov_balance_m1 = abs(results$overlap - ceiling(results$r_overlap))
	ov_balance_p1 = abs(results$overlap - floor(results$r_overlap))
	ov_balance = pmin(ov_balance, ov_balance_m1, ov_balance_p1)
	ov_balance[which(is.na(ov_balance))] = 0
	ov_balance_2 = abs(ov_balance)
	# invariant = results[(which(ov_balance_2<=0.0)),]
	not_invariant = results[(which(ov_balance_2>0)),]
	not_invariant_pairs = rownames(not_invariant)
	results = results[not_invariant_pairs,,drop=FALSE]
	return(results)
}

#' Transform a pairwise table column into tabular form, taking care of NAs
#' 
#' @param lis table of pairwise interactions
#' @param features features to keep in the matrix
#' @param score name of the column to consider
#' @return a 2d table
#' 
lis2mat <- function(lis, features = NULL, score='score') {
	if(nrow(lis)==0) return(NULL)
	if(is.null(features)) features = unique(c(lis$SFE_1, lis$SFE_2))
	MI_vs_null = matrix(NA, nrow=length(features), ncol=length(features))
	rownames(MI_vs_null) = features
	colnames(MI_vs_null) = features
	for(i in 1:nrow(lis)) {
		MI_vs_null[lis[i,'SFE_1'], lis[i,'SFE_2']] = lis[i,score]
		MI_vs_null[lis[i,'SFE_2'], lis[i,'SFE_1']] = lis[i,score]
	}
	# diag(MI_vs_null) = NA
	return(MI_vs_null)
}

#' APC correction. Core version, unaware of any subdivision of the dataset.
#' 
#' @param results table of pairwise interactions
#' @param do.ASC whether calculating APC or ASC
#' @param do.fix_positivity Whether setting negative scores to 0
#' @param use.median Use median instead of average
#' @param column Which column to use as score
#' @return APC scores
APC <- function(results, do.ASC = FALSE, do.fix_positivity = FALSE, do.asinh.before=FALSE, do.asinh.after=FALSE, use.median = TRUE, column='score', out.column='APC', verbose=FALSE) {
	# Build matrix of pairwise scores
	if(do.asinh.before)	results[,column] = asinh(results[,column])
	score_mat = lis2mat(results, NULL, score=column)
	if(is.null(score_mat)) return(NULL)
	diag(score_mat) = NA
	# if required, fix positivity
	if(do.fix_positivity) score_mat = score_mat - min(score_mat, na.rm=TRUE)

	tresults = results[which(!is.na(results[,column])),,drop=FALSE]	
	# APC correction cannot be applied if the background signal cannot be estimated. As measure of estimability, we check whether there are enough alterations participating to the interactions with non-null score. If there is one interaction participating to more that 50% of the interactions we do not apply the APC correction for the specific block.
	estimability = 1 - max(table(c(tresults$SFE_1, tresults$SFE_2)) / nrow(tresults))
	# if(verbose) {
	# 	print(estimability)
	# }
	if(estimability <= 0.5) {
		print('Not correcting with APC as there are not enough interactions to estimate the background distribution. Uncorrected scores will be used.')
		results[,out.column] = results[,column]
		return(results)
	}

	# 'Calculating average background information...'
	tempMI_clean = score_mat[which(!is.na(score_mat))]
	if(length(tempMI_clean)==0) print('Warning! 0 pairs with score different from NA!')
	if(use.median) APC_background = median(tempMI_clean,na.rm=TRUE) else APC_background = mean(tempMI_clean,na.rm=TRUE)

	# Get single feature APC backgrounds (numerators) on single int type blocks
	MIa = apply(score_mat, 1, function(x) {
			if(length(which(!is.na(x)))==0) return(NA)
			x = x[which(!is.na(x))]
			if(use.median) return(median(x, na.rm=TRUE))
			return(mean(x, na.rm = TRUE))
		})

	APC_score1 = rep(MIa, nrow(score_mat))
	dim(APC_score1) = c(nrow(score_mat),nrow(score_mat))
	APC_score2 = t(APC_score1)

	# Put all together
	APC = APC_score1 * APC_score2 / APC_background
	if(do.ASC) APC = APC_score1 + APC_score2 - APC_background
	rownames(APC) = rownames(score_mat)
	colnames(APC) = colnames(score_mat)
	MI_APC = score_mat - APC

	diag(MI_APC) = 0
	results[, out.column] = NA
	results[, out.column] = MI_APC[as.matrix(results[,c('SFE_1','SFE_2')])]
	if(do.asinh.before)	results[,out.column] = asinh(results[,out.column])
	return(results)
}

#' APC correction. Extended version, aware of subdivision of the dataset.
#' 
#' @param results table of pairwise interactions
#' @param do.ASC whether calculating APC or ASC
#' @param do.fix_positivity Whether setting negative scores to 0
#' @param use.median Use median instead of average
#' @param column Which column to use as score
#' @return APC scores
APC.block <- function(results, divide.by=c('int_type', 'direction'), column='score', out.column='APC', verbose=FALSE, ...) {
	results[, out.column] = NA
	if(!is.null(divide.by) & length(divide.by)>0) {
		temp = results
		temp$group = apply(temp[,divide.by], 1, function(x) paste(x, collapse='||', sep='||'))
		APCs = lapply(unique(temp$group), function(x) {
			temp1 = results[which(temp$group == x),]
			# if(verbose) {
			# 	print(x)
			# 	print(dim(temp1))
			# }
			APC1 <- APC(temp1, column=column, out.column=out.column, ...)
		})
		for(i in 1:length(APCs)) {
			APC1 = APCs[[i]]
			common = intersect(rownames(APC1), rownames(results))
			results[common, out.column] = APC1[common, out.column]
		}
	} else {
		APC1 <- APC(results, column=column, out.column=out.column, ...)
		results[, out.column] = APC1[as.matrix(results[,c('SFE_1','SFE_2')])]
	}
	return(results)
}

#' Calculate ASC score
#' 
#' @param results table of pairwise interactions
#' @return a dataframe with ASC scores
#' 
add.ASC <- function(results, column, ...) {
	# Limit to non invariant
	results_redux = prune_invariants(results)

	# Run APC
	results_redux = APC.block(results=results_redux, column=column, use.median=TRUE, divide.by=c('int_type', 'direction'), do.asinh.before=TRUE, do.asinh.after=FALSE, do.ASC=TRUE, do.fix_positivity=FALSE, out.column='APC', ...)

	# Save results
	results$APC = NA
	results[rownames(results_redux),'APC'] = results_redux[rownames(results_redux),'APC']
	return(results)
}



# pa =  APC.block(results_redux, do.ASC=TRUE, do.fix_positivity=FALSE, divide.by=c('int_type', 'direction'), use.median=TRUE, column='stat_diff', verbose=TRUE)





# #------ #------ #------ #------ #------ #------ #------ #------ #------ #------ #------ #------ 

# results$N = MI$M.debug$corrected_N[cbind(results$SFE_1,results$SFE_2)]
# results$freq_min =  pmin(results$freq_1 , results$freq_2)
# results$freq_max =  pmax(results$freq_1 , results$freq_2)
# results$freq_prod =  (results$freq_1 * results$freq_2)
# results$freq_sum =  (results$freq_1 + results$freq_2)
# results$overlap_ratio = results$overlap / results$max_pos_overlap_block
# results$r_overlap_ratio = results$r_overlap / results$max_pos_overlap_block



# print('------> Calculating random scores...')
# # Get scores from randomized bam
# r_MI = sapply(MI$R.MI, function(x) x[cbind(results$SFE_1,results$SFE_2)])
# rownames(r_MI) = rownames(results)

# results$MI_ratio = results$MI / results$r_MI
# results$r_MI_sigmoid_max_block_correction_0.1_15_stat = as.vector(as.dist(Reduce('+', MI$R.MI)/length(MI$R.MI)))
# results$MI_sigmoid_max_block_correction_0.1_15_stat[results$MI_sigmoid_max_block_correction_0.1_15_stat < 1e-15] = 0  # Fix near to 0 MI becoming negative due 


















