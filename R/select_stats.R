# source('AL.R')

#' Create an Alteration Landscape Statistics (ALS) object
#' 
#' @param al Alteration landscape
#' @return An Alteration Landscape Statistics (ALS) object
#' 
new.ALS <- function(al) {
	if(is.null(al)) stop('Input al is NULL')
	als = list()
	# als = list('al' = al)
	class(als) <- "ALS";
	return(als)
}

#' Create an Alteration Matrix Statistics (AMS) object
#' 
#' @param al Alteration matrix
#' @return An Matrix Landscape Statistics (AMS) object
#' 
new.AMS <- function(am) {
	if(is.null(am)) stop('Input am is NULL')
	ams = list()
	class(ams) <- "AMS";
	return(ams)
}

#' Calculate simple statistics for an am
#' 
#' @param am Alteration matrix
#' @return Statistics (AMS object)
#' 
am.stats <- function(am) {
	ams = new.AMS(am)
	am = am * 1.0
	sample.alt.n = rowSums(am, na.rm=TRUE)
	feat.alt.n = colSums(am, na.rm=TRUE)
	numOfEdges <- sum(abs(am), na.rm=TRUE)
	temp = list('n.samples' = nrow(am), 'n.alterations' = ncol(am), 'n.occurrences' = numOfEdges, 'alterations.per.sample'=sample.alt.n, 'alteration.count'=feat.alt.n)
	ams[names(temp)] = temp
	return(ams)
}

#' Compute simple statistics for the AL, including accurate statistics on each sample/alteration block
#' 
#' @param al Alteration landscape
#' @return An Alteration Landscape Statistics (ALS) object
#' 
al.stats <- function(al) {
	als = new.ALS(al)
	temp = am.stats(al$am)
	ams = temp
	als[names(ams)] = ams

	# get stats for each single block
	blocks = get.blocks(al)
	
	als$sample.blocks = list()
	for(ib in 1:length(blocks$sample.blocks)) { fib = blocks$sample.blocks[[ib]]; if(length(fib)>0) {
		subM <- al$am[fib, , drop=FALSE]
		temp = am.stats(subM)
		als$sample.blocks[[ib]] = temp
	}}
	names(als$sample.blocks) = names(blocks$sample.blocks)

	als$alteration.blocks = list()
	for(ib in 1:length(blocks$alteration.blocks)) { fib = blocks$alteration.blocks[[ib]]; if(length(fib)>0) {
		subM <- al$am[, fib, drop=FALSE]
		temp = am.stats(subM)
		als$alteration.blocks[[ib]] = temp
	}}
	names(als$alteration.blocks) = names(blocks$alteration.blocks)

	# for(fib in blocks$sample.blocks) { if(length(fib)>0) {
	# 	for(sib in blocks$alteration.blocks) { if(length(sib)>0) {
	# 		subM <- al$am[fib, sib, drop=FALSE]
	# 		temp = am.stats(al$am)
	# 		als$blocks[[]] = als[names(temp)] = temp
return(als)
}


#' Compute pairwise overlap
#' 
#' @param am Alteration Matrix
#' @return The pairwise overlap between the alterations in 'dist' format
#' 
am.pairwise.alteration.overlap <- function(am) {
	A = am * 1
	overlap = as.dist(t(A) %*% A)
	return(overlap)
}

#' Compute pairwise coverage
#' 
#' @param am Alteration Matrix
#' @return The pairwise overlap between the alterations in 'dist' format
#' 
am.pairwise.alteration.coverage <- function(overlap_M, M.stats) {
	stats = list()
	overlap_M = as.matrix(overlap_M)
	# print(names(M.stats$alteration.count))
	# print(rownames(overlap_M))
	marginal_f1 = replicate(ncol(overlap_M),M.stats$alteration.count[rownames(overlap_M)])
	colnames(marginal_f1) = rownames(marginal_f1)
	marginal_f2 = t(marginal_f1)
	diag(overlap_M) = diag(marginal_f1)
	# me = marginal_f1 + marginal_f2 - 2*overlap_M # not used. Nov 2018
	# coverage = marginal_f2 + marginal_f1 - overlap_M # not used. Nov 2018
	stats[['overlap']] = Matrix:::Matrix(overlap_M)
	# stats[['coverage']] = coverage
	# stats[['me']] = me # not used. Nov 2018
	return(stats)
}

#' Compute overlap statistics for the AL
#' 
#' @param al Alteration landscape
#' @return Pairwise statistics
#' 
al.pairwise.alteration.stats <- function(al, als=NULL, do.blocks=FALSE) {
	if(is.null(als)) als = al.stats(al)
	M.overlap = am.pairwise.alteration.overlap(al$am)
	M.pairwise = am.pairwise.alteration.coverage(M.overlap, als)

	if(do.blocks) {
		blocks = get.blocks(al)
		pairwise.blocks = list()
		for(ib in 1:length(blocks$sample.blocks)) { fib = blocks$sample.blocks[[ib]]; if(length(fib)>0) {
			subM <- al$am[fib, , drop=FALSE]
			M.overlap.block = am.pairwise.alteration.overlap(subM)
			M.pairwise.block = am.pairwise.alteration.coverage(M.overlap.block, als$sample.blocks[[ib]])
			pairwise.blocks[[ib]] = M.pairwise.block
		}}
		names(pairwise.blocks) = names(blocks$sample.blocks)
		M.pairwise[['sample.blocks']] = pairwise.blocks
	}
	return(M.pairwise)
}


#' Compute overlap statistics for the random ams
#' 
#' @param al Alteration landscape
#' @return Pairwise statistics for random ams
#' 
r.al.pairwise.alteration.stats <- function(al, do.blocks=FALSE,  n.cores=1, verbose=TRUE, folder='./') {
	# require(foreach)
	r.M.pairwise = NULL
	if((is.null(al$r.am)) | (length(al$r.am)==0)) return(r.M.pairwise)
	fake_al = al
	fake_al$r.am = NULL

	f0 <- function(seq, fake_al, do.blocks=FALSE, verbose=FALSE) { #, pairwise_corr_MI_func, tosource) {
		# sapply(tosource, function(x) source(x, chdir=TRUE))
		require(select)
		if(verbose) print('Running single iteration of al.pairwise.alteration.stats...')
		tal = fake_al
		tal$am = as.matrix(seq)
		r.M.pairwise = select:::al.pairwise.alteration.stats(tal, NULL, do.blocks=do.blocks)
		return(r.M.pairwise)
	}

	if(n.cores>1) {
		# require(parallel)
		environment(f0) <- .GlobalEnv
		cl <- makeCluster(n.cores, outfile=paste(folder, "r.al.pairwise.alteration.stats.log", sep='/'))
		# clusterExport(cl, 'al.pairwise.alteration.stats')
		# clusterExport(cl, 'al.stats')
		# clusterExport(cl, 'am.stats')
		# clusterExport(cl, 'am.pairwise.alteration.overlap')
		# clusterExport(cl, 'am.pairwise.alteration.coverage')
		# clusterExport(cl, 'new.ALS')
		# clusterExport(cl, 'get.blocks')
		# clusterExport(cl, 'new.AMS')
		r.M.pairwise <- clusterApply(cl=cl, al$r.am, fun=f0, fake_al = fake_al, do.blocks=FALSE, verbose=verbose)
		stopCluster(cl)
		# } else {
		# 	stop('Error in registering a parallel cluster for randomization.')
		# }
	} else {
		r.M.pairwise = lapply(1:length(al$r.am), function(rid) {
			# if(verbose) print(rid)
			r.M.pairwise = f0(al$r.am[[rid]], fake_al = fake_al, do.blocks=do.blocks, verbose=verbose)
			return(r.M.pairwise)
		})
	}
	return(r.M.pairwise)
}


#' Compute the ME Pvalue starting from overlap statistics
#' 
#' @param alp Observed Overlap statistics
#' @param r.alp Random Overlap statistics
#' @return Pairwise ME P.vlaue
#' 
al.pairwise.alteration.ME_p <- function(alp, r.alp, enhanced.version = TRUE) {
	temp = lapply(r.alp, function(x) x$overlap) # extract coverage data
	temp1 = as.matrix(Reduce('+', temp)) / length(r.alp) # sum of coverage to decide whether performing left or right tail test
	observedValue = as.matrix(alp$overlap) # extract vector of observed values
	# Pos test: alternative: greater: tests whether the observed coverage is bigger than null distribution, that is, more mutual exclusive

	if(!enhanced.version) {
		temp.diff.pos = lapply(temp, function(x) (observedValue - as.matrix(x))>=0) # TRUE -> CO, FALSE -> ME
		temp.diff.pos.cum = Reduce('+', temp.diff.pos) / length(r.alp)
		rm(temp.diff.pos)
		temp.diff.neg = lapply(temp, function(x) (observedValue - as.matrix(x))<=0) # TRUE -> ME, FALSE -> CO
		temp.diff.neg.cum = Reduce('+', temp.diff.neg) / length(r.alp)
		rm(temp.diff.neg)
	} else {
		# efficient versions
		temp.diff.pos.cum = Reduce(
			function(x,y) {
	 			y = as.matrix(observedValue - y) >= 0
				x + y
			}, temp, init=0)
		temp.diff.pos.cum = temp.diff.pos.cum / length(r.alp)
		#
		temp.diff.neg.cum = Reduce(
			function(x,y) {
	 			y = as.matrix(observedValue - y) <= 0
				x + y
			}, temp, init=0)
		temp.diff.neg.cum = temp.diff.neg.cum / length(r.alp)
	}
	side = observedValue - temp1 < 0 # TRUE -> ME, FALSE -> CO
	temp.pval =  temp.diff.pos.cum*side + temp.diff.neg.cum*(!side)
	temp.sign = side * 2 - 1
	ME = list('p.value' = temp.pval, 'sign' = temp.sign)
	return(ME)
}

# Old version - based on coverage
#' Compute the ME Pvalue starting from overlap statistics
#' 
#' @param alp Observed Overlap statistics
#' @param r.alp Random Overlap statistics
#' @return Pairwise ME P.vlaue
#' 
# al.pairwise.alteration.ME_p <- function(alp, r.alp) {
# 	temp = lapply(r.alp, function(x) x$coverage) # extract coverage data
# 	temp1 = Reduce('+', temp) / length(r.alp) # sum of coverage to decide whether performing left or right tail test
# 	observedValue = alp$coverage # extract vector of observed values
# 	# Pos test: alternative: greater: tests whether the observed coverage is bigger than null distribution, that is, more mutual exclusive
# 	temp.diff.pos = lapply(temp, function(x) (observedValue - x)<=0) # TRUE -> CO, FALSE -> ME
# 	temp.diff.pos.cum = Reduce('+', temp.diff.pos) / length(r.alp)
# 	temp.diff.neg = lapply(temp, function(x) (observedValue - x)>=0) # TRUE -> ME, FALSE -> CO
# 	temp.diff.neg.cum = Reduce('+', temp.diff.neg) / length(r.alp)
# 	side = observedValue - temp1 > 0 # TRUE -> ME, FALSE -> CO
# 	temp.pval =  temp.diff.pos.cum*side + temp.diff.neg.cum*(!side)
# 	temp.sign = side * 2 - 1
# 	ME = list('p.value' = temp.pval, 'sign' = temp.sign)
# 	return(ME)
# }
