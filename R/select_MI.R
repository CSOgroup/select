#' Binary Mutual Information. Works on single units, vectors or matrices
#' 
#' @param n_00
#' @param n_01
#' @param n_10
#' @param n_11
#' @param n_f1
#' @param n_f2
#' @return MI
#' 
binary.MI <- function(n_00, n_01, n_10, n_11, n_f1, n_f2, tolerance=0.00001) {
	# Normalize
	tot = n_00 + n_01 + n_10 + n_11
	tot[tot == 0] = NA # induce NA for those pairs with sample-size = 0
	f_00 = n_00 / tot
	f_01 = n_01 / tot
	f_10 = n_10 / tot
	f_11 = n_11 / tot
	f_f1 = n_f1 / tot
	f_f2 = n_f2 / tot

	# Joint entropy
	v_11 = f_11 * log2(f_11)
	v_10 = f_10 * log2(f_10)
	v_01 = f_01 * log2(f_01)
	v_00 = f_00 * log2(f_00)
	v_11[!is.finite(v_11)] = 0 # Fix infinity - happens when there are no samples for a given symbol.
	v_10[!is.finite(v_10)] = 0
	v_01[!is.finite(v_01)] = 0
	v_00[!is.finite(v_00)] = 0
	vt = v_11 + v_10 + v_01 + v_00

	# Marginal entropy
	mf1 = f_f1
	mf2 = f_f2
	mf1 = mf1 * log2(mf1) + ((1 - mf1) * log2(1 - mf1))
	mf2 = mf2 * log2(mf2) + ((1 - mf2) * log2(1 - mf2))
	mf1[!is.finite(mf1)] = 0
	mf2[!is.finite(mf2)] = 0

	# calculate MI
	MI = vt - mf1 - mf2
	MI[which((n_f1==0) | (n_f2==0))] = NA # fix cases with NO f1 or f2 alterations. Such cases do not make any sense!
	erroneous = MI[which(MI<0)]
	if(min(erroneous) < -tolerance) stop('Strong Negative MI detected. Something is really wrong with the input data.')
	MI[which(MI < 0)] = 0
	return(MI)
}

#' Debug version of Binary Mutual Information. Works on single units, vectors or matrices. 
#' 
#' @param n_00
#' @param n_01
#' @param n_10
#' @param n_11
#' @param n_f1
#' @param n_f2
#' @return debug structure containing MI and debug variables
#' 
binary.MI_debug <- function(n_00, n_01, n_10, n_11, n_f1, n_f2) {
	debug = list()
	debug$check_cons_1 = n_01 + n_11 - n_f1
	debug$check_cons_2 = n_10 + n_11 - n_f2

	# Normalize
	tot = n_00 + n_01 + n_10 + n_11
	tot[tot == 0] = NA # induce NA for those pairs with sample-size = 0
	f_00 = n_00 / tot
	f_01 = n_01 / tot
	f_10 = n_10 / tot
	f_11 = n_11 / tot
	f_f1 = n_f1 / tot
	f_f2 = n_f2 / tot
	debug$tot = tot
	debug$final_f = list('f_00'=f_00, 'f_01'=f_01, 'f_10'=f_10, 'f_11'=f_11, 'f_f1'=f_f1, 'f_f2'=f_f2)

	# Joint entropy
	v_11 = f_11 * log2(f_11)
	v_10 = f_10 * log2(f_10)
	v_01 = f_01 * log2(f_01)
	v_00 = f_00 * log2(f_00)
	v_11[!is.finite(v_11)] = 0 # Fix infinity - happens when there are no samples for a given symbol.
	v_10[!is.finite(v_10)] = 0
	v_01[!is.finite(v_01)] = 0
	v_00[!is.finite(v_00)] = 0
	vt = v_11 + v_10 + v_01 + v_00
	debug$joint_entropy = vt

	# Marginal entropy
	mf1 = f_f1
	mf2 = f_f2
	mf1 = mf1 * log2(mf1) + ((1 - mf1) * log2(1 - mf1))
	mf2 = mf2 * log2(mf2) + ((1 - mf2) * log2(1 - mf2))
	mf1[!is.finite(mf1)] = 0
	mf2[!is.finite(mf2)] = 0
	debug$single_entropy_f1 = mf1
	debug$single_entropy_f2 = mf2

	# calculate MI and MI_APC
	MI = vt - mf1 - mf2
	MI[which((n_f1==0) | (n_f2==0))] = NA # fix cases with NO f1 or f2 alterations. Such cases do not make any sense!
	debug$MI = MI
	
	return(debug)
}

#' Consistency enforcer for correction function for weight MI measure
#' 
#' @param 
#' @param 
#' @return Corrected N 
#' 
f_corr_fix_N_consistency <- function(N, marginal_f1, marginal_f2, oldN, pseudo_count = 1) {
	# which(N > oldN)
	# print(oldN)
	N = pmax(N, marginal_f1 + marginal_f2 + pseudo_count)
	N = pmin(N, oldN)
	N = round(N)
	return(N)
}

#' Correction function for weight MI measure
#' 
#' @param marginal_f1 Marginal frequency of event 1
#' @param marginal_f2 Marginal frequency of event 2
#' @param real.N Observed sample size
#' @return Corrected N
#' 
f_sigmoid_max_correction <- function(marginal_f1, marginal_f2, real.N, rho = 0.1, lambda = 15) {
	w1 = 1/(1+exp(-lambda*(marginal_f1/real.N - rho)))
	w2 = 1/(1+exp(-lambda*(marginal_f2/real.N - rho)))
	w1[which(marginal_f1==0)] = 0
	w2[which(marginal_f2==0)] = 0
	flag = (marginal_f1 == pmax(marginal_f1, marginal_f2))*1
	N = real.N * (w1*flag + w2*(1-flag))
	N = round(N)
	N = f_corr_fix_N_consistency(N, marginal_f1, marginal_f2, real.N)
	return(N)
}

#' Block aware correction function for weight MI measure
#' 
#' @param 
#' @param 
#' @return Corrected N 
#' 
f_sigmoid_max_block_correction <- function(marginal_f1, marginal_f2, real.N, rho = 0.1, lambda = 15, block_als) {
	als = block_als
	real.N = als$n.samples
	f_correction_algorithm = f_sigmoid_max_correction
	temp = lapply(als$sample.blocks, function(als) {
		overlap = as.matrix(als$alteration.pairwise$overlap)
		marginal_f1 = replicate(ncol(overlap), diag(overlap))
		colnames(marginal_f1) = rownames(marginal_f1)
		marginal_f2 = t(marginal_f1)
		# diag(overlap) = diag(marginal_f1)
		real.N = als$n.samples
		N = f_correction_algorithm(marginal_f1, marginal_f2, real.N = real.N, rho = rho, lambda = lambda)
		return(N)
	})
	temp = temp[which(!sapply(temp, is.null))]
	temp = Reduce('+', temp)
	N = f_corr_fix_N_consistency(temp, marginal_f1, marginal_f2, real.N)
	return(N)
}

#' Weighted MI measure
#' 
#' @param al The alteration landscape
#' @param alp Observed AL Overlap statistics
#' @return Pairwise weighted MI
#' 
corrected.binary.MI <- function(als, cache = NULL, f_correction_algorithm = f_sigmoid_max_block_correction, ...) {
	# debug = list()

	# Collect data
	overlap = as.matrix(als$alteration.pairwise$overlap)
	marginal_f1 = replicate(ncol(overlap), diag(overlap))
	colnames(marginal_f1) = rownames(marginal_f1)
	marginal_f2 = t(marginal_f1)
	# diag(overlap) = diag(marginal_f1)

	# Apply N correction
	N = cache$corrected_N
	if(is.null(N)) {
		real.N = als$n.samples
		N = f_correction_algorithm(marginal_f1, marginal_f2, real.N, ...)
	}

	# Derive corrected frequencies
	v_11 = overlap
	v_10 = marginal_f1 - v_11
	v_01 = marginal_f2 - v_11
	v_11 = v_11 / N
	v_10 = v_10 / N
	v_01 = v_01 / N
	v_00 = 1 - v_11 - v_01 - v_10
	mf1 = marginal_f1 / N
	mf2 = marginal_f2 / N

	# Calculate MI
	MI = binary.MI(v_00, v_01, v_10, v_11, mf1, mf2)
	# temp = binary.MI_debug(v_00, v_01, v_10, v_11, mf1, mf2)
	# debug = c(debug, temp)
	# return(debug)

	cache = list('corrected_N' = N)
	# debug = c(debug, temp)
	# debug$'cache' = cache
	# debug$'tot' = N
	# debug$corrected_N = N
	# debug$original_inc_matrix = list('inc_f1' = marginal_f1, 'inc_f2' = marginal_f2)

	return(list('measure'=MI, 'N'=N, 'cache' = cache)) #'debug' = debug, 
}


#' Generic interface to calculate pairwise scores between all the pairs of alterations. Used as front-end for MI measures
#' 
#' @param al The alteration landscape
#' @param alp Observed AL Overlap statistics
#' @return Pairwise ME P.vlaue
#' 
al.pairwise.alteration.custom_measure <- function(al, als, f = corrected.binary.MI, ...) {
	stats = list()
	temp_real <- f(als, ...)
	stats[['measure']] = temp_real$measure
	stats[['data']] = temp_real
	return(stats)
}

#' Compute custom score for the random ams
#' 
#' @param r.al Random Alteration landscape
#' @param r.als Random Alteration stats
#' @param cache constant data. Speed up calculations!
#' @return Custom Pairwise score for random ams
#' 
r.al.pairwise.alteration.custom_measure <- function(als, r.als, n.cores=1, verbose=TRUE, folder='./', ...) {
	stats = NULL

	f0 <- function(r.als, als, verbose=FALSE, f0.params) {
		if(verbose) print('Running single iteration of al.pairwise.alteration.stats...')
		f0.params[['als']] = list('n.samples' = als$n.samples, 'alteration.pairwise' = r.als)
		f0.params[['al']] = NULL
		stats <- do.call(al.pairwise.alteration.custom_measure, f0.params)
		return(stats)
	}

	if(n.cores>1) {
		# require(parallel)
		cl <- makeCluster(n.cores, outfile=paste(folder, "r.al.pairwise.alteration.custom_measure.log", sep='/'))
		clusterExport(cl, 'al.pairwise.alteration.custom_measure')
		clusterExport(cl, 'corrected.binary.MI')
		clusterExport(cl, 'binary.MI')
		# clusterExport(cl, 'f_sigmoid_max_block_correction')
		f0.params = list(...)
		environment(f0) <- .GlobalEnv
		stats <- clusterApply(cl=cl, r.als, f0, verbose = verbose, als = als, f0.params)
		stopCluster(cl)
		# } else {
		# 	stop('Error in registering a parallel cluster for randomization.')
		# }
	} else {
		stats = lapply(1:length(r.als), function(rid) {
			# if(verbose) print(rid)
			# r.M.pairwise = f0(r.als[[rid]], fake_al = fake_al, do.blocks=do.blocks, verbose=verbose)
			stats <- al.pairwise.alteration.custom_measure(al = NULL, als = list('n.samples' = als$n.samples, 'alteration.pairwise' = r.als[[rid]]), ...)
			return(stats)
		})
	}

	stats = lapply(stats, function(x) x$measure)
	return(stats)
}

#' Compute P vlaue for custom score
#' 
#' @param observed Random Alteration landscape
#' @param random Random Alteration stats
#' @return Custom Pairwise score for random ams
#' 
pairwise.p_MI <- function(observed, random) {
	observedValue = observed

	# faster & mem efficient
	temp.diff.pos.cum = Reduce(
		function(x,y) {
 			y = (observedValue - y) <= 0
			x + y
		}, random, init=0)
	temp.diff.pos.cum = temp.diff.pos.cum / length(random)

	# temp.diff.pos = lapply(random, function(x) (observedValue - x)<=0) #
	# temp.diff.pos.cum = Reduce('+', temp.diff.pos) / length(random)
	temp.pval =  temp.diff.pos.cum
	# p_MI = list('p.value' = temp.pval)
	return(temp.pval)
}
