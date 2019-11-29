# source('AL.R')
#-------- #-------- #-------- #-------- #-------- 
# Randomization section

# .birewireEvents <- function(M, Q=100, events=NULL, verbose=FALSE) {
# 	uniqueTypes = NULL
# 	if(!is.null(events)) {
# 		eventTypes <- events$eventTypes[rownames(M)];
# 		# print(eventTypes)
# 		uniqueTypes <- unique(eventTypes);
# 		# print(uniqueTypes)
# 	}
# 	if(!is.null(uniqueTypes)) {
# 		for(eType in uniqueTypes) {
# 			subM <- M[which(eventTypes == eType), , drop=FALSE];
# 			numOfEdges <- sum(subM);
# 			# print(length(subM))
# 			subM <- birewire.rewire.bipartite(subM, max.iter=numOfEdges*Q, verbose=verbose);
# 			M[rownames(subM), ] <- subM;
# 		}
# 	} else {
# 		numOfEdges <- sum(M);
# 		M <- birewire.rewire.bipartite(M, max.iter=numOfEdges*Q, verbose=verbose);
# 	}
# 	return(M);
# }

#' Generate random binary matrix
#' 
#' @param am Input binary matrix (no row/col names!)
#' @param params parameter list
#' @return One random binary matrix
#' 
birewireBlocEvents <- function(M, Q=100, row.blocks=NULL, col.blocks=NULL, verbose=FALSE, switch_threshold = 15, seed=NULL) {
	# require(BiRewire)
	if(!is.null(seed)) set.seed(seed)
	if(is.null(row.blocks)) {
		row.blocks = list('all'= 1:nrow(M))
	}
	if(is.null(col.blocks)) {
		col.blocks = list('all'= 1:ncol(M))
	}
	missing_rows = which(! 1:nrow(M) %in% unlist(row.blocks))
	missing_cols = which(! 1:ncol(M) %in% unlist(col.blocks))
	if(length(missing_rows) > 0) { print("###### The following rows will not be reshuffled:"); print(missing_rows) }
	if(length(missing_cols) > 0) { print("###### The following cols will not be reshuffled:"); print(missing_cols) }
	for(fib in row.blocks) { if(length(fib)>0) {
			if(verbose) print(paste('\nConsidering rowBloc of', length(fib), 'rows.' ))
		for(sib in col.blocks) { if(length(sib)>0) {
			# if(verbose) print(paste('Permuting colBloc of', length(sib), 'cols.' ))
			subM <- M[fib, sib, drop=FALSE]
			numOfEdges <- sum(abs(subM*1.0))
			if(numOfEdges > 1) {
				if(verbose) print(paste('Permuting colBloc of', length(fib), 'rows and', length(sib),'cols,',numOfEdges,'edges...'))
				good_cols = length(which(apply(subM, 2, function(x) length(which(x>0)))>0))
				good_rows = length(which(apply(subM, 1, function(x) length(which(x>0)))>0))
				if(good_cols >= switch_threshold) {
					if(verbose) print(paste('Using birewire.rewire.bipartite as we have ', good_cols, ' >= ', switch_threshold, '...'))
					subM <- birewire.rewire.bipartite(subM, max.iter=numOfEdges*Q, verbose=verbose)
					M[fib, sib] <- subM
				} else {
					if(verbose) print(paste('Using simple row scrambling as we have ', good_cols, ' < ', switch_threshold, '...'))
					M1 = (apply(subM,2,function(x) sample(x))) # preserve alteration frequency but discard sample frequency
					M[fib, sib] <- M1
				}
			} else {
				if(verbose) print(paste('NOT Permuting colBloc of', length(fib), 'rows and', length(sib),'cols,',numOfEdges,'edges...'))
			}
		}}
	}}

	M = (M != 0) # * 1.0 # binarize to have occurrence counts. birewire acts strangely otherwise!!!
	M = Matrix:::Matrix(M) # convert Matrix back to Matrix object
	return(M)
}

#' Generate random alteration matrix. Enables parallel calculation
#' 
#' @param am Input alteration matrix
#' @param params parameter list
#' @return A collection of random Alteration matrices
#' 
gen.random.am <- function(M, N=1000, permuteFun=birewireBlocEvents, n.cores=1, folder='./', r.seed = 100, ...) {
	# require(foreach)
	set.seed(r.seed)
	seed.vector=runif(N, 0, 10^9)
	if(n.cores>1) {	
		# library(doParallel)
		cl <- makeCluster(n.cores, outfile=paste(folder, "gen.random.am.log", sep='/'))
		registerDoParallel(cl)  
		if( getDoParRegistered() ) {
			randomMs <- foreach(i=seed.vector) %dopar% permuteFun(M, seed=i, ...)
			stopCluster(cl)
		} else {
			stop('Error in registering a parallel cluster for randomization.')
		}
	} else {
		randomMs <- foreach(i=seed.vector) %do% permuteFun(M, seed=i, ...);
	}
	return(randomMs);
}

#' Generate random alteration landscapes
#' 
#' @param al The original alteration landscape
#' @param params parameter list
#' @return A collection of random Alteration Landscapes
#' 
gen.random.al <- function(al, N=1000, params=NULL, r.seed = 100, ...) {
	if(!class(al)=='AL') stop('An AL object is required.')
	blocks = get.blocks(al)
	al$r.seed = r.seed
	mdat <- gen.random.am(al$am, N=N, permuteFun=birewireBlocEvents, row.blocks=blocks$sample.blocks, col.blocks=blocks$alteration.blocks, r.seed = al$r.seed, ...)
	al$r.am = mdat
	return(al)
}

# save(dataset, file=paste(outfolder,'dataset_v2.RData',sep=''))
# save(dataset, file=paste(outfolder,'dataset_null_model_v2.RData',sep='')) # save the memo object.
