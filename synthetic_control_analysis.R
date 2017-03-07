#This is the analysis file. The functions used in this file are cointained in synthetic_control_functions.R
#There are three model variants: 
# *_full - Full synthetic control model with all covariates (excluding user-specified covariates).
# *_time - Trend adjustment using the specified variable (e.g., non-respiratory hospitalization or population size) as the denominator.

#############################
#                           #
#    System Preparations    #
#                           #
#############################

source('synthetic_control_functions.R', local = TRUE)

packages <- c('parallel', 'splines', 'zoo', 'lubridate', 'RcppRoll', 'CausalImpact', 'ggplot2', 'reshape')
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
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name)), groups)
time_points <- unique(ds[[1]][, date_name])

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
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change, intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(time_index = 1:nrow(covars)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {scale(data[, outcome_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, outcome_name])})
outcome_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset      <- sapply(ds, FUN = function(data) {data[, outcome_name] - data[, denom_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name] - data[, denom_name])})
outcome_offset      <- scale(outcome_offset)

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,        covars = covars_full, time_points = time_points), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome_offset, covars = covars_time, time_points = time_points), groups)

###############################
#                             #
#        Main analysis        #
#                             #
###############################

#Start Cluster for CausalImpact (the main analysis function).
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact', 'impactExtract', 'intervention_date', 'time_points', 'n_seasons'), environment())

impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)

stopCluster(cl)

#Save the inclusion probabilities from each of the models.
inclusion_prob_full <- setNames(lapply(impact_full, inclusionProb), groups)
inclusion_prob_time <- setNames(lapply(impact_time, inclusionProb), groups)

#All model results combined
quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], denom_data = ds[[group]][, denom_name], mean = outcome_mean[group],        sd = outcome_sd[group],        eval_period = eval_period, post_period = post_period)}), groups)
quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name], mean = outcome_offset_mean[group], sd = outcome_offset_sd[group], eval_period = eval_period, post_period = post_period, trend = TRUE)}), groups)

#Model predicitons
pred_quantiles_full <- sapply(quantiles_full, getPred, simplify = 'array')
pred_quantiles_time <- sapply(quantiles_time, getPred, simplify = 'array')

#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')

#Rate ratios for evaluation period.
rr_mean_full <- t(sapply(quantiles_full, getRR))
rr_mean_time <- t(sapply(quantiles_time, getRR))

colnames(rr_mean_time) <- paste('ITS', colnames(rr_mean_time))

cumsum_prevented <- sapply(groups, FUN = function(group, quantiles) {
	is_post_period <- which(time_points >= post_period[1])
	is_pre_period <- which(time_points < post_period[1])
	
	#Cumulative sum of prevented cases
	cases_prevented <- quantiles[[group]]$pred_samples - outcome_plot[, group]
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

#Pred Sensitivity Analysis
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact', 'impactExtract', 'predSensitivityAnalysis', 'inclusionProb', 'rrPredQuantiles', 'getPred', 'getRR', 'groups', 'ds', 'data_full', 'denom_name', 'outcome_mean', 'outcome_sd', 'intervention_date', 'eval_period', 'post_period', 'time_points', 'n_seasons'), environment())

sensitivity_analysis_pred_2  <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 2 ))), c('Lower CI', 'Point Estimate', 'Upper CI'))
sensitivity_analysis_pred_10 <- setNames(as.data.frame(t(parSapply(cl, groups, predSensitivityAnalysis, ds = ds, zoo_data = data_full, denom_name = denom_name, outcome_mean = outcome_mean, outcome_sd = outcome_sd, intervention_date = intervention_date, eval_period = eval_period, post_period = post_period, time_points = time_points, n_seasons = n_seasons, n_pred = 10))), c('Lower CI', 'Point Estimate', 'Upper CI'))

stopCluster(cl)

bad_sensitivity_groups <- sapply(covars_full, function (covar) {ncol(covar) <= 3})
sensitivity_covars_full <- covars_full[!bad_sensitivity_groups]
sensitivity_ds <- ds[!bad_sensitivity_groups]
sensitivity_impact_full <- impact_full[!bad_sensitivity_groups]
sensitivity_groups <- groups[!bad_sensitivity_groups]

#Weight Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE); library(RcppRoll, quietly = TRUE)})
clusterExport(cl, c('sensitivity_ds', 'doCausalImpact', 'impactExtract', 'weightSensitivityAnalysis', 'rrPredQuantiles', 'sensitivity_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons', 'outcome_mean', 'outcome_sd', 'eval_period', 'post_period'), environment())

sensitivity_analysis_full <- setNames(parLapply(cl, sensitivity_groups, weightSensitivityAnalysis, covars = sensitivity_covars_full, ds = sensitivity_ds, impact = sensitivity_impact_full, time_points = time_points, intervention_date = intervention_date, n_seasons = n_seasons, outcome = outcome, mean = outcome_mean, sd = outcome_sd, eval_period = eval_period, post_period = post_period), sensitivity_groups)

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
rr_table <- cbind(rr_mean_time[!bad_sensitivity_groups, ], sensitivity_table)