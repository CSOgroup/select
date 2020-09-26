require(select)

load('test_random_v1.RData')

sample.class = rep('SAMPLE', nrow(test_random_v1@am))
names(sample.class) = rownames(test_random_v1@am)

alteration.class = rep('ALT', ncol(test_random_v1@am))
names(alteration.class) = colnames(test_random_v1@am)

alpi <- select::select(
	M = test_random_v1@am,
	sample.class = sample.class,
	alteration.class = alteration.class,
	folder = './',
	r.seed = 110,
	n.cores = 40,
	vetos = NULL,
	n.permut = 100,
	min.feature.support=5,
	min.feature.freq=0.001,
	remove.0.samples=TRUE,
	remove.unknown.class.samples=TRUE,
	rho = 0.1,
	lambda = 15,
	save.intermediate.files = FALSE,
	randomization.switch.threshold = 30,
	max.memory.size=100,
	calculate_APC_threshold = TRUE,
	calculate_FDR = TRUE,
	verbose = TRUE
	)

require(dplyr)
options(width=200)

alpi %>% arrange(-select_score) %>% head(10)
alpi %>% arrange(wMI_p.value) %>% head(10)

alpi %>% count(select_score_good_cancer_cell_2017_criterion_1)
alpi %>% count(wMI_p.value_FDR)
alpi %>% nrow
