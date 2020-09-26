#' SELECT analysis on an AL with null model
#'
#' Run SELECT analysis on input AL
#'
#' @import foreach
#' @import igraph
#' @import BiRewire
#' @import doParallel
#' @import parallel
#' @import matrixStats
#' @importFrom stats quantile
#' @importFrom utils write.table
#' @importFrom stats as.dist
#' @importFrom stats ecdf
#' @importFrom stats median
#' @importFrom stats runif
#' @importFrom utils combn
#' 
#' @param al Input AL with precalculated null model
#'
#' @return selection interactions between alterations
#'
#' #@examples
#' #select.on.AL()
#'
#' @export
select.on.AL <- function(
	al,
	# sample.class,
	# alteration.class,
	folder = './',
	r.seed = 110,
	n.cores = 40,
	vetos = NULL,
	# n.permut = 1000,
	min.feature.support=5,
	min.feature.freq=0.001,
	remove.0.samples=TRUE,
	remove.unknown.class.samples=TRUE,
	rho = 0.1,
	lambda = 15,
	save.intermediate.files = TRUE,
	# randomization.switch.threshold = 30,
	# max.memory.size=100,
	calculate_APC_threshold = TRUE,
	calculate_FDR = TRUE,
	verbose = TRUE
	) {

	folder = paste(folder, '/', sep='')
	if(save.intermediate.files) dir.create(folder, recursive=TRUE, showWarnings = FALSE)

	print(paste('Running SELECT... [', folder,']'))

	print(paste('-> Parsing and Filtering GAM...'))
	# Create Alteration Landscape (AL) object
	null_model = al
	input_al = al$al
	rm(al)
	al = select:::new.AL(as.matrix(input_al@am))
	sample.class = input_al@samples$sample.class
	alteration.class = input_al@alterations$alteration.class
	al$samples$sample.class = sample.class[rownames(input_al@am)]
	al$alterations$alteration.class = alteration.class[colnames(input_al@am)]

	# Filter AL
	al = select:::filter.al(al, list('min.feature.support'=min.feature.support, 'min.feature.freq'=min.feature.freq, 'remove.0.samples'=remove.0.samples, 'remove.unknown.class.samples' = remove.unknown.class.samples))
	if(save.intermediate.files) save(al, file=paste(folder,'al.RData',sep=''))

	# Adjust null model
	r.al = null_model
	r.al$al = al
	r.al$r.am = lapply(r.al$r.am, function(x) x[rownames(r.al$al$am), colnames(r.al$al$am)])

	nbatches = 1

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
