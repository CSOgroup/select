#' Assemble a table of pairwise interactions between alterations. Given N alterations, the table will contain N(N-1)/2 rows.
#' 
#' @param al Alteration Landscape
#' @param als AL statistics
#' @return a dataframe
#' 
interaction.table <- function(al, als) {
	features = colnames(al$am)
	num_features = length(features)
	if(! num_features == als$n.alterations) stop('Inconsistent al and als!')
	result_cols = c('SFE_1', 'SFE_2')
	results = data.frame(matrix(NA, nrow=num_features*(num_features-1)/2, ncol=length(result_cols)))
	colnames(results) = result_cols
	temp = rep(1:length(features), length(features))
	dim(temp) = c(length(features), length(features))
	diag(temp) = NA
	temp1 = t(temp)
	temp1[which(temp1 > temp)] = NA
	temp[which(temp < temp1)] = NA
	temp = cbind(as.vector(temp1), as.vector(temp))
	temp = temp[which(!is.na(temp[,1]) & !is.na(temp[,2])),,drop=FALSE]
 	results[, 'SFE_1'] = features[temp[,1]]
 	results[, 'SFE_2'] = features[temp[,2]]
	rownames(results) = paste(results[,'SFE_1'], results[,'SFE_2'], sep= ' - ')
	results[, 'name'] = rownames(results)

	results$type_1 = al$alterations$alteration.class[results$SFE_1]
	results$type_2 = al$alterations$alteration.class[results$SFE_2]
	temp = t(apply(results, 1, function(x) sort(c(x['type_1'], x['type_2']))))
	results$int_type = paste(temp[,1], temp[,2], sep= ' - ')

	results[,'support_1'] = als$alteration.count[results[,'SFE_1']]
	results[,'support_2'] = als$alteration.count[results[,'SFE_2']]
	results[,'freq_1'] = results[,'support_1']/als$n.samples
	results[,'freq_2'] = results[,'support_2']/als$n.samples

	results[,'overlap'] = as.vector(as.dist(als$alteration.pairwise$overlap))
	# results[,'coverage'] = as.vector(as.dist(als$alteration.pairwise$coverage))
	# results[,'freq_coverage'] = results[,'coverage'] / als$n.samples

	# Get max possible overlap
	if(!is.null(als$sample.blocks)) {
		A = lapply(als$sample.blocks, function(x) {
			a = rep(x$alteration.count, length(x$alteration.count))
			dim(a) = c(length(x$alteration.count),length(x$alteration.count))
			rownames(a) = names(x$alteration.count)
			colnames(a) = names(x$alteration.count)
			b = t(a)
			c = pmin(a,b)
		})
		A = Reduce('+', A)
		results$max_overlap = as.vector(as.dist(A))
	}

	results[,'freq_overlap'] = results[,'overlap'] / results[,'max_overlap']



	# results$max_naive_overlap = pmin(results$support_1, results$support_2)
	return(results)
}

#' Add some other statistics
#' 
#' @param results table of pairwise interactions
#' @param r.als random AL statistics
#' @return a dataframe with additional info
#' 
add.r.info <- function(results, r.als) {
	# Average random overlap
	temp = lapply(r.als, function(x) x$overlap)
	results[,'r_overlap'] = as.vector(as.dist(Reduce('+', temp) / length(r.als)))
	results[,'r_freq_overlap'] = results[,'r_overlap'] / results[,'max_overlap']

	# Difference between observed and expected overlap
	results[,'diff_overlap'] = results[,'overlap'] - results[,'r_overlap']
	results[,'abs_diff_overlap'] = abs(results[,'diff_overlap'])
	temp = sign(results[,'diff_overlap'])
	results[,'direction'] = 0
	results[temp==1,'direction'] = 'CO'
	results[temp==0,'direction'] = 'none'
	results[temp==-1,'direction'] = 'ME'

	return(results)
	# results[,'r_coverage'] = as.vector(as.dist(Reduce('+', lapply(dataset_props$R, function(x) x$pairwise$coverage)) / length(dataset_props$R)))
	# temp1 = listmat_to_elemvec(temp)
	# # prova = sapply(temp1, function(x) sapply(x, function(y) mad(y, na.rm=TRUE)))
	# prova = sapply(temp1, function(x) sapply(x, function(y) sd(y, na.rm=TRUE)))
	# results[,'r_overlap_mad'] = as.vector(as.dist(prova))
}


#' Add some other statistics
#' 
#' @param results table of pairwise interactions
#' @param r.als random AL statistics
#' @return a dataframe with additional info
#' 
add.patterns <- function(results, patterns, name) {
	if('stat' %in% names(patterns)) results[,paste(name,'_stat',sep='')] = as.vector(as.dist(patterns$stat))
	if('p.value' %in% names(patterns)) results[,paste(name,'_p.value',sep='')] = as.vector(as.dist(patterns$p.value))
	return(results)
	# for(i in names(MaxEnt)) {
	# 	if('empirical_P.value' %in% names(MaxEnt[[i]])) results[,paste(i,'_p.value',sep='')] = (MaxEnt[[i]]$empirical_P.value[rownames(results)])
	# 	if('fit_P.value' %in% names(MaxEnt[[i]])) results[,paste(i,'_p.value.fit',sep='')] = (MaxEnt[[i]]$fit_P.value[rownames(results)])
	# 	if('M' %in% names(MaxEnt[[i]])) results[,paste(i,'_stat',sep='')] = (MaxEnt[[i]]$M[rownames(results)])
	# 	if('R' %in% names(MaxEnt[[i]])) results[,paste(i,'_r_stat',sep='')] = apply(MaxEnt[[i]]$R[rownames(results), , drop=FALSE], 1, function(x) mean(x, na.rm=TRUE))
	# 	if('R' %in% names(MaxEnt[[i]])) results[,paste(i,'_r_mad',sep='')] = apply(MaxEnt[[i]]$R[rownames(results), , drop=FALSE], 1, function(x) mad(x, na.rm=TRUE))
	# }
}

#' Build interaction table
#' 
#' @param al Alteration Landscape
#' @param als AL statistics
#' @param r.alp random AL statistics
#' @param patterns patterns
#' @return a dataframe with all interaction info pooled together
#' 
build.interaction.table <- function(al, als, r.alp, patterns = NULL) {
	alpi = interaction.table(al, als)
	alpi = add.r.info(alpi, r.alp)
	for(i in names(patterns)) alpi = add.patterns(alpi, patterns = patterns[[i]], name=i)
	return(alpi)
}

#' Remove vetoed interactions (eg CNA between the same chromosome)
#' 
#' @param results table of pairwise interactions
#' @return a dataframe without vetoed interactions
#' 
prune_vetoed <- function(results, veto = NULL) {
	if(is.null(veto)) return(results)
	to_forget = as.vector(apply(veto, 1, function(x) c(paste(x[1],x[2], sep = ' - '),paste(x[2],x[1], sep = ' - '))))
	# intra_chr_results = results[which(rownames(results) %in% to_forget),]
	new_results =  results[which(!rownames(results) %in% to_forget),,drop=FALSE]
	return(new_results)
}

# results[,'descr_1'] = sapply(results[,'SFE_1'], function(x1) paste(dataset$mes$eventGraphIds[[x1]], collapse='/'))
# results[,'descr_2'] = sapply(results[,'SFE_2'], function(x1) paste(dataset$mes$eventGraphIds[[x1]], collapse='/'))

# f_correction_algorithms = config$MI$MI_f_correction_algorithms
# MI_c = lapply(f_correction_algorithms, function(x) {
# 	MI_name = paste(config$MI$prefix,x$name,sep='')
# 	MI = NULL
# 	if(file.exists(paste(outfolder,MI_name,'.RData',sep=''))) {
# 		load(file=paste(outfolder,MI_name,'.RData',sep=''))
# 		MI$name = MI_name
# 	}
# 	return(MI)
# })
# names(MI_c) = sapply(MI_c, function(x) x$name)
# patterns = list()
# temp = lapply(MI_c, function(x) list('p.value'=x$p.value, 'stat'=x$M.MI, 'N'=x$M.N, 'R.stat'=x$R.MI))
# patterns = c(patterns, temp)

