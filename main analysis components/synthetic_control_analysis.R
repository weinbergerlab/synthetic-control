#This is the analysis file. The functions used in this file are cointained in synthetic_control_functions.R
#There are two model variants: 
# *_full - Full synthetic control model with all covariates (excluding user-specified covariates).
# *_time - Trend adjustment using the specified variable (e.g., non-respiratory hospitalization or population size) as the denominator.

#############################
#                           #
#    System Preparations    #
#                           #
#############################

source('synthetic_control_functions.R', local = TRUE)

packages <- c('parallel', 'splines', 'lubridate', 'RcppRoll','pomp', 'BoomSpikeSlab', 'ggplot2', 'reshape','dummies')
packageHandler(packages, update_packages, install_packages)
sapply(packages, library, quietly = TRUE, character.only = TRUE)

#Detects number of available cores on computers. Used for parallel processing to speed up analysis.
n_cores <- detectCores()
set.seed(1)

###################################################
#                                                 #
# Directory setup and initialization of constants #
#                                                 #
###################################################

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
data_file <- paste(input_directory, file_name, sep = '')
prelog_data <- read.csv(data_file, check.names = FALSE)
groups <- as.character(unique(unlist(prelog_data[, group_name], use.names = FALSE)))
if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}

###############################################
#                                             #
# Data and covariate preparation for analysis #
#                                             #
###############################################


prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
#if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name,outcome_name)), groups)
time_points <- unique(ds[[1]][, date_name])


#Monthly dummies
if(n_seasons==4){x<-quarter(as.Date(time_points))}
if(n_seasons==12){x<-month(as.Date(time_points))}
if(n_seasons==3){
    x.m<-month(as.Date(time_points))
    x<-x.m
    x[x.m %in% c(1,2,3,4)]<-1
    x[x.m %in% c(5,6,7,8)]<-2
    x[x.m %in% c(9,10,11,12)]<-3
    }
season.dummies<-dummy(x)
season.dummies<-season.dummies[,-n_seasons]

ds <- lapply(ds, function(ds) {
	if (!(denom_name %in% colnames(ds))) {
		ds[denom_name] <- 0
	}
	return(ds)
})

sparse_groups <- sapply(ds, function(ds) {
	return(ncol(ds[!(colnames(ds) %in% c(date_name, group_name, denom_name, outcome_name, exclude_covar))]) == 0)
})
ds <- ds[!sparse_groups]
groups <- groups[!sparse_groups]

#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change,season.dummies=season.dummies,  intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {data[, outcome_name]})
outcome_plot=outcome
offset<- sapply(ds, FUN=function(data) exp(data[, denom_name]) )  #offset term on original scale; 1 column per age group


#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,       covars = covars_full), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome, covars = covars_time, trend=TRUE), groups)

###############################
#                             #
#        Main analysis        #
#                             #
###############################

#Start Cluster for CausalImpact (the main analysis function).
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(BoomSpikeSlab, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons'), environment())

impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)

stopCluster(cl)

#calculate WAIC
waic_full<-t(sapply(impact_full,waic_fun))
waic_time<-t(sapply(impact_time,waic_fun))


#Save the inclusion probabilities from each of the models.
inclusion_prob_full <- setNames(lapply(impact_full, inclusionProb), groups)
inclusion_prob_time <- setNames(lapply(impact_time, inclusionProb), groups)

#All model results combined
quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], denom_data = ds[[group]][, denom_name],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)}), groups)


#Model predicitons
pred_quantiles_full <- sapply(quantiles_full, getPred, simplify = 'array')
pred_quantiles_time <- sapply(quantiles_time, getPred, simplify = 'array')

#Pointwise RR and uncertainty for second stage meta analysis
log_rr_quantiles   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_quantiles}, simplify = 'array')
dimnames(log_rr_quantiles)[[1]] <- time_points
log_rr_sd   <- sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_sd}, simplify = 'array')
#log_rr_full_t_samples.prec<-sapply(quantiles_full,   FUN = function(quantiles) {quantiles$log_rr_full_t_samples.prec}, simplify = 'array')
saveRDS(log_rr_quantiles, file=paste0(output_directory, country, "_log_rr_quantiles.rds"))
saveRDS(log_rr_sd, file=paste0(output_directory, country, "_log_rr_sd.rds"))
#saveRDS(log_rr_full_t_samples.prec, file=paste0(output_directory, country, "_log_rr_full_t_samples.prec.rds"))


#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')

#Rate ratios for evaluation period.
rr_mean_full <- t(sapply(quantiles_full, getRR))
rr_mean_time <- t(sapply(quantiles_time, getRR))

rr_mean_full_intervals <- data.frame('Estimate (95% CI)'     = makeInterval(rr_mean_full[, 2], rr_mean_full[, 3], rr_mean_full[, 1]), check.names = FALSE, row.names = groups)
rr_mean_time_intervals <- data.frame('ITS Estimate (95% CI)' = makeInterval(rr_mean_time[, 2], rr_mean_time[, 3], rr_mean_time[, 1]), check.names = FALSE, row.names = groups)

colnames(rr_mean_time) <- paste('ITS', colnames(rr_mean_time))


cumsum_prevented <- sapply(groups, FUN = function(group, quantiles) {
	is_post_period <- which(time_points >= post_period[1])
	is_pre_period <- which(time_points < post_period[1])
	
	#Cumulative sum of prevented cases
	cases_prevented <- t(quantiles[[group]]$pred_samples) - outcome[, group]
	cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period, ], 2, cumsum)
	cumsum_cases_prevented_pre <- matrix(0, nrow = nrow(cases_prevented[is_pre_period, ]), ncol = ncol(cases_prevented[is_pre_period, ]))
	cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
	cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}, quantiles = quantiles_full, simplify = 'array')

################################
#                              #
#     Sensitivity Analyses     #
#                              #
################################

#Pred Sensitivity Analysis--tests effect of changing prior on Ncovars from 3 to 2 to 10
# cl <- makeCluster(n_cores)
# clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
# clusterExport(cl, c('doCausalImpact', 'predSensitivityAnalysis', 'inclusionProb', 'rrPredQuantiles', 'getPred', 'getRR', 'groups', 'ds', 'data_full', 'denom_name', 'outcome_mean', 'outcome_sd', 'intervention_date', 'eval_period', 'post_period', 'time_points', 'n_seasons'), environment())
# 
# sensitivity_analysis_pred_2  <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 2 ))), c('Lower CI', 'Point Estimate', 'Upper CI'))
# sensitivity_analysis_pred_10 <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 10))), c('Lower CI', 'Point Estimate', 'Upper CI'))
# 
# stopCluster(cl)
# 
# sensitivity_analysis_pred_2_intervals  <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_analysis_pred_2[, 2],  sensitivity_analysis_pred_2[, 3],  sensitivity_analysis_pred_2[, 1]),  row.names = groups, check.names = FALSE)
# sensitivity_analysis_pred_10_intervals <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_analysis_pred_10[, 2], sensitivity_analysis_pred_10[, 3], sensitivity_analysis_pred_10[, 1]), row.names = groups, check.names = FALSE)
# 
 bad_sensitivity_groups <- sapply(covars_full, function (covar) {ncol(covar) <= 3})
 sensitivity_covars_full <- covars_full[!bad_sensitivity_groups]
 sensitivity_ds <- ds[!bad_sensitivity_groups]
 sensitivity_impact_full <- impact_full[!bad_sensitivity_groups]
 sensitivity_groups <- groups[!bad_sensitivity_groups]

 #Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(BoomSpikeSlab, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
clusterExport(cl, c('sensitivity_ds', 'doCausalImpact',  'weightSensitivityAnalysis', 'rrPredQuantiles', 'sensitivity_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons',  'eval_period', 'post_period'), environment())

sensitivity_analysis_full <- setNames(parLapply(cl, sensitivity_groups, weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = time_points, intervention_date = intervention_date, n_seasons = n_seasons, outcome = outcome,  eval_period = eval_period, post_period = post_period), sensitivity_groups)

stopCluster(cl)

sensitivity_pred_quantiles  <- lapply(sensitivity_analysis_full, FUN = function(sensitivity_analysis) {
	pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
	for (sensitivity_index in 1:length(sensitivity_analysis)) {
		pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
	}
	return(pred_list)
})

#Table of rate ratios for each sensitivity analysis level
sensitivity_table <- t(sapply(sensitivity_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full, original_rr = rr_mean_full))
sensitivity_table_intervals <- data.frame('Estimate (95% CI)' = makeInterval(sensitivity_table[, 2],  sensitivity_table[, 3],  sensitivity_table[, 1]),
																					'Top Control 1' = sensitivity_table[, 'Top Control 1'],
																					'Inclusion Probability of Control 1' = sensitivity_table[, 'Inclusion Probability of Control 1'],
																					'Control 1 Estimate (95% CI)' = makeInterval(sensitivity_table[, 7],  sensitivity_table[, 8],  sensitivity_table[, 6]),
																					'Top Control 2' = sensitivity_table[, 'Top Control 2'],
																					'Inclusion Probability of Control 2' = sensitivity_table[, 'Inclusion Probability of Control 2'],
																					'Control 2 Estimate (95% CI)' = makeInterval(sensitivity_table[, 12],  sensitivity_table[, 13],  sensitivity_table[, 11]),
																					'Top Control 3' = sensitivity_table[, 'Top Control 3'],
																					'Inclusion Probability of Control 3' = sensitivity_table[, 'Inclusion Probability of Control 3'],
																					'Control 3 Estimate (95% CI)' = makeInterval(sensitivity_table[, 17],  sensitivity_table[, 18],  sensitivity_table[, 16]), check.names = FALSE)
rr_table <- cbind.data.frame(round(rr_mean_time[!bad_sensitivity_groups, ],2), sensitivity_table)
rr_table_intervals <- cbind('ITS Estimate (95% CI)' = rr_mean_time_intervals[!bad_sensitivity_groups, ], sensitivity_table_intervals)
