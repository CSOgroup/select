#' internal SELECT routine
#'
#' Run SELECT analysis on input al. All params and null model must be set and valid
#'
#' @import foreach
#' @import igraph
#' @import BiRewire
#' @import doParallel
#' @import parallel
#' @import Matrix
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom stats as.dist
#' @importFrom stats ecdf
#' @importFrom stats median
#' @importFrom stats runif
#' @importFrom utils combn
#' 
#' @return selection interactions between alterations
#'
#' #@examples
#' #select()
#'
 select.internal <- function(al, r.al,
	# sample.class,
	# alteration.class,
	folder = './',
	# r.seed = 110,
	n.cores = 40,
	vetos = NULL,
	# n.permut = 1000,
	# min.feature.support=5,
	# min.feature.freq=0.001,
	# remove.0.samples=TRUE,
	# remove.unknown.class.samples=TRUE,
	rho = 0.1,
	lambda = 15,
	save.intermediate.files = TRUE,
	# randomization.switch.threshold = 30,
	# max.memory.size=100,
	calculate_APC_threshold = TRUE,
	calculate_FDR = TRUE,
	verbose = TRUE,
	nbatches = 1,
	FDR.thresh = 0.1
	) {
	print(paste('-> Collecting event stats on observed GAM...'))
	# Get AL stats
	ptm <- proc.time()
	als = select:::al.stats(al) # Univariate stats
	alp = select:::al.pairwise.alteration.stats(al, als, do.blocks=TRUE) # Pairwise (alterations)
	als[['alteration.pairwise']] = alp
	als[['alteration.pairwise']]$sample.blocks = NULL
	for(i in names(als$sample.blocks)) als$sample.blocks[[i]][['alteration.pairwise']] = alp$sample.blocks[[i]]
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(als, file=paste(folder,'als.RData',sep=''))


	print(paste('-> Calculating wMI on observed GAM...'))
	# Get MI stats
	ptm <- proc.time()
	MI = select:::al.pairwise.alteration.custom_measure(al, als, f = corrected.binary.MI, block_als = als, f_correction_algorithm = f_sigmoid_max_block_correction, rho = rho, lambda = lambda)
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(MI, file=paste(folder,'MI.RData',sep=''))


	print(paste('-> Collecting event stats on null model...'))
	# Get AL stats for null model
	ptm <- proc.time()
	if(nbatches>1) {
		for(batch_n in 1:nbatches) {
			r.alp = list()
			rm(r.al)
			load(file=paste(folder,'r.al_batch.',batch_n,'.RData',sep='')) #r.al
			temp.r.alp = select:::r.al.pairwise.alteration.stats(r.al, do.blocks=FALSE,  n.cores=n.cores, verbose=FALSE, folder=folder)
			r.alp = c(r.alp, temp.r.alp)
		}
	} else {
		r.alp = select:::r.al.pairwise.alteration.stats(r.al, do.blocks=FALSE,  n.cores=n.cores, verbose=FALSE, folder=folder)	
	}
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(r.alp, file=paste(folder,'r.alp.RData',sep=''))

	# ME/CO naive P value
	print(paste('-> Comparing observed GAM vs null model...'))
	ptm <- proc.time()
	ME_p.value = select:::al.pairwise.alteration.ME_p(als$alteration.pairwise, r.alp)
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(ME_p.value, file=paste(folder,'ME_p.value.RData',sep=''))

	# Get MI stats for null model
	print(paste('-> Calculating wMI on null model...'))
	ptm <- proc.time()
	r.MI = select:::r.al.pairwise.alteration.custom_measure(als, r.alp, f = corrected.binary.MI, f_correction_algorithm = f_sigmoid_max_block_correction, rho = rho, lambda = lambda, cache=MI$data$cache, n.cores=1, verbose=TRUE, folder=folder)
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(r.MI, file=paste(folder,'r.MI.RData',sep=''))

	# MI-based P value
	print(paste('-> Calculating P-value based on wMI...'))
	ptm <- proc.time()
	MI_p.value = select:::pairwise.p_MI(MI$measure, r.MI)
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}
	if(save.intermediate.files) save(MI_p.value, file=paste(folder,'MI_p.value.RData',sep=''))

	MI_stats = list('stat'=MI$measure, 'p.value'=MI_p.value)
	ME_stats = list('p.value'=ME_p.value$p.value)

	# Background MI
	print(paste('-> Deriving average background wMI...'))
	ptm <- proc.time()
	# r.MI.elemvec = select:::listmat_to_elemvec(r.MI)
	# # r.MI.q.25.75 = sapply(r.MI.elemvec, function(x) sapply(x, function(y) max(quantile(y, c(0.25,0.5,0.75))))) # Better estimation of random MI by considering max within 25-75%
	# r.MI.q.25.75 = sapply(r.MI.elemvec, function(x) {
	# 	apply(matrixStats:::colQuantiles(x, probs = c(0.25,0.5,0.75)),1,max)
	# 	# sapply(x, function(y) max(quantile(y, c(0.25,0.5,0.75))))
	# }) # Better estimation of random MI by considering max within 25-75%
	# rm(r.MI.elemvec)
	r.MI.q.25.75 = select:::background_wMI(r.MI, probs = c(0.25,0.5,0.75)) # fast and efficient version
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}

	# Build interaction table
	print(paste('-> Building results table...'))
	alpi = select:::build.interaction.table(al, als, r.alp, patterns = list('wMI'=MI_stats, 'ME'=ME_stats, 'E.r.wMI'=list('stat'=r.MI.q.25.75)))
	rm(r.MI.q.25.75)

	# Prune and get APC
	print(paste('-> Calculating APC score...'))
	ptm <- proc.time()
	alpi = select:::prune_vetoed(alpi, vetos)
	alpi$MI_diff = alpi$wMI_stat - alpi$E.r.wMI_stat
	alpi <- select:::add.ASC(alpi, column='MI_diff')
	if(verbose) {
		print(proc.time() - ptm)
		print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
	}

	alpi$wMI_p.value_FDR = NA
	if(calculate_FDR) {
		print(paste('-> Performing wMI FDR analysis...'))
		ptm <- proc.time()
		FDR = select:::FDR_analysis(alpi, r.MI, column='wMI_p.value', FDR.thresh=FDR.thresh)
		alpi = FDR$alpi
		if(verbose) {
			print(proc.time() - ptm)
			print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
		}
	}

	alpi$select_score_good_cancer_cell_2017_criterion_1 = NA
	# Select optimal APC score threshold (removing 95% of FP)
	if(calculate_APC_threshold) {
		print(paste('-> Calculating APC threshold...'))
		TP.FP.thr = select:::establish_APC_threshold(alpi)
		threshold = TP.FP.thr$thresh['t']
		alpi$select_score_good_cancer_cell_2017_criterion_1 = FALSE
		alpi$select_score_good_cancer_cell_2017_criterion_1[which(alpi$APC >= threshold)] = TRUE
	}

	alpi$select_score = alpi$APC
	alpi$APC = NULL

	# Save pairwise results
	if(save.intermediate.files) save(alpi, file=paste(folder,'alpi.RData',sep=''))
	out.format = c('name','SFE_1','SFE_2','int_type', 'support_1', 'support_2', 'overlap', 'r_overlap', 'max_overlap', 'direction', 'freq_overlap', 'wMI_p.value', 'ME_p.value', 'wMI_p.value_FDR', 'select_score', 'select_score_good_cancer_cell_2017_criterion_1')
	if(save.intermediate.files)  write.table(alpi[,out.format], sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste(folder,'alpi.txt',sep=''))
	return(alpi)
}








#' SELECT analysis
#'
#' Run SELECT analysis on input gam
#'
#' @import foreach
#' @import igraph
#' @import BiRewire
#' @import doParallel
#' @import parallel
#' @import Matrix
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom stats as.dist
#' @importFrom stats ecdf
#' @importFrom stats median
#' @importFrom stats runif
#' @importFrom utils combn
#' 
#' @param M gam
#'
#' @return selection interactions between alterations
#'
#' #@examples
#' #select()
#'
#' @export
select <- function(
	M,
	sample.class,
	alteration.class,
	folder = './',
	r.seed = 110,
	n.cores = 40,
	vetos = NULL,
	n.permut = 1000,
	min.feature.support=5,
	min.feature.freq=0.001,
	remove.0.samples=TRUE,
	remove.unknown.class.samples=TRUE,
	rho = 0.1,
	lambda = 15,
	save.intermediate.files = TRUE,
	randomization.switch.threshold = 30,
	max.memory.size=100,
	calculate_APC_threshold = TRUE,
	calculate_FDR = TRUE,
	verbose = TRUE
	) {

	folder = paste(folder, '/', sep='')
	print(paste('Running SELECT... [', folder,']'))
	# Create analysis folder
	if(save.intermediate.files) dir.create(folder, recursive=TRUE, showWarnings = FALSE)
	# if(save.intermediate.files) dir.create(paste(folder,'log',sep='/'), recursive=TRUE, showWarnings = FALSE)
	
	print(paste('-> Parsing and Filtering GAM...'))
	# Create Alteration Landscape (AL) object
	al = new.AL(M)
	al$samples$sample.class = sample.class[rownames(al$am)]
	al$alterations$alteration.class = alteration.class[colnames(al$am)]

	# Filter AL
	al = select:::filter.al(al, list('min.feature.support'=min.feature.support, 'min.feature.freq'=min.feature.freq, 'remove.0.samples'=remove.0.samples, 'remove.unknown.class.samples'=remove.unknown.class.samples))
	if(save.intermediate.files) save(al, file=paste(folder,'al.RData',sep=''))


	print(paste('-> Generating null model...'))
	# Split up the random data in batches number of batches in order to avoid memory overflows
	estimated.r.size = as.numeric(object.size(al$am)) * n.permut
	if(max.memory.size < 1) max.memory.size = 1
	max.mem.size = max.memory.size * 2^30
	nbatches = 1
	n.permut.batch = n.permut
	if(estimated.r.size > max.mem.size) {
		nbatches = ceiling(estimated.r.size / max.mem.size)
		n.permut.batch = rep(ceiling(n.permut / nbatches), nbatches)
		excess = sum(n.permut.batch) - n.permut
		n.permut.batch[length(n.permut.batch)] = n.permut.batch[length(n.permut.batch)] - excess
		if(n.permut.batch[length(n.permut.batch)] < 2) n.permut.batch[length(n.permut.batch)] = 2
	}

	# Generate null model
	ptm <- proc.time()
	for(batch_n in 1:nbatches) {
		r.al = gen.random.al(al, n.cores=n.cores, N=n.permut.batch[batch_n], r.seed = r.seed - 1 + batch_n, switch_threshold = randomization.switch.threshold, verbose=FALSE, folder=folder)
		if(verbose) {
			print(proc.time() - ptm)
			print(paste('Memory usage:', gc()['Vcells',2], 'MB'))
		}
		if(save.intermediate.files & nbatches==1) save(r.al, file=paste(folder,'r.al.RData',sep=''))
		if(nbatches>1) save(r.al, file=paste(folder,'r.al_batch.',batch_n,'.RData',sep=''))
	}

	alpi = select.internal(al, r.al,
		folder = folder,
		n.cores = n.cores,
		vetos = vetos,
		rho = rho,
		lambda = lambda,
		save.intermediate.files = save.intermediate.files,
		calculate_APC_threshold = calculate_APC_threshold,
		calculate_FDR = calculate_FDR,
		verbose = verbose,
		nbatches = nbatches
	)

	return(alpi)
}
