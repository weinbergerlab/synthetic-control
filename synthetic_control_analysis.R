#This is the analysis file. The functions are cointained in synthetic_control_functions.R
#There are six model variants: 
# *_offset - Simple pre-post analysis using the specified variable (e.g., non-respiratory hospitalization or population size) as an offset term (denominator). This is Model 1 from the paper.
# *_full - This is the full synthetic control model with all covariates.
# *_noj - This is the synthetic control model with bronchitis/bronchiolitis removed.
# *_ach - This is similar to _offset, but the specified variable is used as the sole covariate.
# *_none - No covariate or offsest.
# *_time - Trend adjustment with the same denominator as used in the _offset model.

rm(list = ls(all = TRUE)) #Clear workspace
gc()

#install.packages('devtools')
#library(devtools)
#devtools::install_github('google/CausalImpact')

library('parallel', quietly = TRUE)
library('splines', quietly = TRUE)
library('lubridate', quietly = TRUE)
library('CausalImpact', quietly = TRUE)
source('synthetic_control_functions.R')

#Detects number of available cores on computers. Used for parallel processing to speed up analysis.
n_cores <- detectCores()
set.seed(1)
par_defaults <- par(no.readonly = TRUE)

#############################
#                           #
#User-initialized constants.#
#                           #
#############################

country <- 'SC_Dec_Coverage groups'
factor_name <- 'age_group'
date_name <- 'date'
n_seasons <- NULL #12 for monthly, 4 for quarterly, 3 for trimester data.
do_weight_check <- FALSE
exclude_from_covars <- c('ACM-NoPCV', 'J00-06', 'J09-18', 'J20-22', 'J30-39',  'J40-47',  'J60-70', 'J80-84', 'J85_J86', 'J90-94', 'J95-99', 'A30-49', 'G00-09', 'H65-75', 'B95-98')

##################################################
#                                                #
#Directory setup and initialization of constants.#
#                                                #
##################################################

#Load file (.csv format). You can change the input and output directories to point to an appropriate spot on your computer.
input_directory <- paste('~/Documents/Synthetic Control Data/SC_Dec/', country, '/', sep = '')
output_directory <- paste(input_directory, 'results/', sep = '')
dir.create(output_directory, showWarnings = FALSE)
file_name <- paste('prelog', country, 'processed', 'data.csv', sep = '_')
pcv_file <- paste(input_directory, file_name, sep = '')
ds1a <- read.csv(pcv_file, check.names = FALSE)
age_groups <- paste('Age Group', unique(unlist(ds1a$age_group, use.names = FALSE)))

#Account for code-naming differences
#all_cause_pneu_name - Gives the outcome variable (e.g. pneumonia) 
#all_cause_name - Gives the name of denominator used for some of the analyses. Can be population size, non-respiratory hospitalization, etc.
SC_set <- c('SC1', 'SC2', 'SC3', 'SC_HDI', 'SC_HDI_Region', 'SC_National', 'SC_HDI_No_Pop', 'SC_HDI_Region_No_Pop', 'SC_National_No_Pop', 'SC_HDI_No_Pop_All', 'SC_HDI_Region_No_Pop_All', 'SC_National_No_Pop_All')
SC_subchapter_set <- c('SC_Subchapter_HDI', 'SC_Subchapter_HDI_Region', 'SC_Subchapter_National')
SC_Dec_set <- c('SC_Dec_Coverage groups', 'SC_Dec_HDI from Muni', 'SC_Dec_HDIxSuperRegion', 'SC_Dec_National', 'SC_Dec_Region')
if (country %in% SC_set || country %in% SC_subchapter_set) {
	all_cause_name <- 'ACM-NoJ'
	all_cause_pneu_name <- 'J12_18'
} else if (country %in% SC_Dec_set) {
	all_cause_name <- 'ACM-NoPCV'
	all_cause_pneu_name <- 'J12_18'
} else if (country == 'US') {
	all_cause_name <- 'ach_sid_noresp'
	all_cause_pneu_name <- 'm_ACP'
} else {
	all_cause_name <- 'ach_noj'
	all_cause_pneu_name <- 'J12_18'
}

#Date variables
#data_start_date gives training period start date.
#intervention_date gives training period end date.
#data_end_date gives end of data_set.
if (country == 'Brazil') {
	data_start_date <- as.Date('2004-01-01')
} else {
	data_start_date <- min(as.Date(ds1a[, date_name]))
}
if (country == 'Brazil' || country %in% SC_set) {
	intervention_date <- as.Date('2009-12-31') 
} else if (country %in% SC_subchapter_set || country %in% SC_Dec_set) {
	intervention_date <- as.Date('2009-04-30')
} else if (country == 'Chile') {
	intervention_date <- as.Date('2010-12-31') 
} else if (country == 'Ecuador') {
	intervention_date <- as.Date('2009-12-31') 
} else if (country == 'Mexico') {
	intervention_date <- as.Date('2005-12-31') 
} else if (country == 'US') {
	intervention_date <- as.Date('1999-12-31') 
}
data_end_date <- max(as.Date(ds1a[, date_name]))
pre_period <- c(data_start_date, intervention_date) #Define training period
post_period <- c(intervention_date + 1, data_end_date) #Define post-vaccine period.
#Define evaluation period.
if (country == 'Brazil') {
	eval_period <- c(as.Date(ds1a[, date_name][nrow(ds1a) - 23]), as.Date(ds1a[, date_name][nrow(ds1a)]))
} else if (country %in% SC_set || country %in% SC_subchapter_set || country %in% SC_Dec_set) {
	eval_period <- c(data_end_date - years(2), data_end_date)
} else if (country == 'Chile' || country == 'Ecuador') {
	eval_period <- c(as.Date(ds1a[, date_name][nrow(ds1a) - 11]), as.Date(ds1a[, date_name][nrow(ds1a)]))
} else if (country == 'Mexico') {
	eval_period <- c(as.Date('2010-01-01'), as.Date('2011-12-01'))
} else if (country == 'US') {
	eval_period <- c(as.Date('2003-01-01'), as.Date('2004-12-01'))
}

#Seasons
if (is.null(n_seasons)) {
	if (country %in% SC_set || country %in% SC_subchapter_set || country %in% SC_Dec_set) {
		n_seasons <- 3
	} else {
		n_seasons <- 12
	}	
}

#############################################
#                                           #
#Data and covariate preparation for analysis#
#                                           #
#############################################

ds1a[, date_name] <- as.Date(ds1a[, date_name])

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(unique(ds1a$age_group), FUN = logTransform, factor_name = factor_name, date_name = date_name, start_date = data_start_date, prelog_data = ds1a), age_groups)
data_start <- match(data_start_date, ds[[1]][, date_name])
time_points <- ds[[1]][, date_name][data_start:nrow(ds[[1]])]

ds <- lapply(ds, function(ds) {
	if (!(all_cause_name %in% colnames(ds))) {
		ds[all_cause_name] <- 0
	}
	return(ds)
})

#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars <- setNames(lapply(ds, FUN = function(ds_group) {
	if (country == 'Brazil') {
		#Eliminates effects from 2008 coding change
		covars <- ds_group[data_start:nrow(ds_group), 4:ncol(ds_group)]
		month_i <- as.factor(as.numeric(format(ds_group[, date_name][data_start:nrow(ds_group)], '%m')))
		spline <- setNames(as.data.frame(bs(1:nrow(covars), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
		year_2008 <- numeric(nrow(covars))
		year_2008[1:nrow(covars) >= match(as.Date('2008-01-01'), ds_group[, date_name])] <- 1
		data <- cbind.data.frame(year_2008, spline, month_i)
		trend <- lapply(covars, getTrend, data = data)
		covars <- covars - trend
	} else {
		covars <- ds_group[data_start:nrow(ds_group), 4:ncol(ds_group)]
	}
	if (intervention_date > as.Date('2009-09-01')) {
		covars$pandemic <- ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0))
	}
	covars <- as.data.frame(lapply(covars[, apply(covars, 2, var) != 0], scale), check.names = FALSE)
	if (!(all_cause_name %in% colnames(covars))) {
		covars[all_cause_name] <- 0
	}
	return(covars)
}), age_groups)
covars <- sapply(covars, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_from_covars)]})
covars_time <- setNames(lapply(covars, FUN = function(covars) {
	as.data.frame(list(time_index = 1:nrow(covars)))
}), age_groups)

#Standardize the outcome variables and save the original mean and SD for later analysis.
outcome <- sapply(ds, FUN = function(data) {scale(data[, all_cause_pneu_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, all_cause_pneu_name])})
outcome_sd <- sapply(ds, FUN = function(data) {sd(data[, all_cause_pneu_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset <- sapply(ds, FUN = function(data) {data[, all_cause_pneu_name] - data[, all_cause_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd <- sapply(ds, FUN = function(data) {sd(data[, all_cause_pneu_name] - data[, all_cause_name])})

#Combine the outcome, covariates, and time point information.
data_full   <- setNames(lapply(age_groups, makeTimeSeries, outcome = outcome,        covars = covars,      time_points = time_points, scale_outcome = FALSE), age_groups)
data_time   <- setNames(lapply(age_groups, makeTimeSeries, outcome = outcome_offset, covars = covars_time, time_points = time_points, scale_outcome = TRUE ), age_groups)

###############################
#                             #
#        Main analysis        #
#                             #
###############################

#Start Cluster for CausalImpact (the main analysis function).
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(CausalImpact, quietly = TRUE))
clusterExport(cl, c('ds', 'doCausalImpact', 'intervention_date', 'time_points', 'n_seasons'))

impact_full   <- setNames(parLapply(cl, data_full,   doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), age_groups)
impact_time   <- setNames(parLapply(cl, data_time,   doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE),  age_groups)

stopCluster(cl)

#Save the synthetic control response data.
dir.create(paste(output_directory, 'impact_responses/', sep = ''), showWarnings = FALSE)
impact_list <- list('impact_full' = impact_full, 'impact_time' = impact_time)
for (impact_set in names(impact_list)) {
	impact_series <- NULL
	for (group_factor in names(impact_list[[impact_set]])) {
		series <- data.frame('group_factor' = group_factor, 'date' = time(impact_list[[impact_set]][[group_factor]]$series), impact_list[[impact_set]][[group_factor]]$series, check.names = FALSE, row.names = NULL)
		impact_series <- rbind(impact_series, series)
	}
	write.csv(impact_series[, c('group_factor', 'date', 'point.pred')], paste(output_directory, 'impact_responses/', country, '_', impact_set, '_response.csv', sep = ''), row.names = FALSE)
}
rm(impact_list)

#Save the inclusion probabilities from each of the models.
inclusion_prob_full <- setNames(mapply(inclusionProb, age_groups, impact_full, SIMPLIFY = FALSE), age_groups)
inclusion_prob_time <- setNames(mapply(inclusionProb, age_groups, impact_time, SIMPLIFY = FALSE), age_groups)

#All model results combined
quantiles_full   <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(impact = impact_full[[age_group]],   all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period)}), age_groups)
quantiles_time   <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(impact = impact_time[[age_group]],   all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period)}), age_groups)

#Model predicitons
pred_quantiles_full   <- sapply(quantiles_full,   getPred, simplify = 'array')
pred_quantiles_time   <- sapply(quantiles_time,   getPred, simplify = 'array')

#Rolling rate ratios
rr_roll_full <- sapply(quantiles_full, FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array')
rr_roll_time <- sapply(quantiles_time, FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array')

#Rate ratios for evaluation period.
rr_mean_full   <- t(sapply(quantiles_full,   getRR))
rr_mean_time   <- t(sapply(quantiles_time,   getRR))

rr_col_names <- c('Lower CI', 'Point Estimate', 'Upper CI')

colnames(rr_mean_full)   <- rr_col_names
colnames(rr_mean_time)   <- paste('ITS', rr_col_names)

#Output the rate ratio estimates to a new file.
write.csv(rr_mean_full,   paste(output_directory, country, '_rr_full.csv', sep = ''))
write.csv(rr_mean_time,   paste(output_directory, country, '_rr_time_trend.csv', sep = ''))
write.csv(rr_roll_full,   paste(output_directory, country, '_rr_roll_full.csv', sep = ''))
write.csv(rr_roll_time,   paste(output_directory, country, '_rr_roll_time.csv', sep = ''))

plot_cumsum_prevented <- sapply(age_groups, FUN = function(age_group, quantiles) {
	pred_samples_post_full <- quantiles[[age_group]]$pred_samples_post_full
	
	post_period_start <- which(time_points == post_period[1]) 
	post_period_end <- which(time_points == post_period[2]) 
	is_post_period <- which(time_points >= post_period[1])
	is_pre_period <- which(time_points < post_period[1])
	
	#Cumulative sum of prevented cases
	cases_prevented <- pred_samples_post_full - outcome_plot[, age_group]
	cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period, ], 2, cumsum)
	cumsum_cases_prevented_pre <- apply(cases_prevented[is_pre_period, , drop = FALSE], 2, cumsum)
	cumsum_cases_prevented_pre[, ] <- 0
	cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
	plot_cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}, quantiles = quantiles_full, simplify = 'array')

#Full model predictions
pred_full <- matrix(data = 0, nrow = 0, ncol = 4, dimnames = list(NULL, c('pred_lcl', 'pred_median', 'pred_ucl', 'age_groups')))
for (age_group in age_groups) {
	pred_full <- rbind(pred_full, cbind(pred_quantiles_full[, ,age_group], rep(age_group, each = nrow(pred_quantiles_full[, ,age_group]))))
}
write.csv(pred_full, paste(output_directory, country, '_pred_mult_full.csv', sep = ''), row.names = FALSE)

##################################
#                                #
#    Visualization of results    #
#                                #
##################################

#Plot covars, models, and impacts
par(par_defaults)
invisible(lapply(age_groups, FUN = function(age_group) {
	#View scaled covariates
	matplot(covars[[age_group]], type = 'n')
	title(main = age_group, sub = 'Cumulative cases prevented')
	
	#View cumulative sums
	matplot(plot_cumsum_prevented[, , age_group], type = 'n')
	abline(h = 0)
	
	#Plot predictions
	par(mfrow = c(3, 1))
	plotPred(age_group, pred_quantiles_full[, , age_group])
	title(main = paste(age_group, 'Synthetic controls estimate'))

	plotPred(age_group, pred_quantiles_time[, , age_group])
	title(main = paste(age_group, 'Interupted time series estimate'))
	par(par_defaults)
	
	#Rolling rate ratio
	par(mfrow = c(2, 1))
	matplot(rr_roll_full[, , age_group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l')
	title(main = paste(age_group, 'Synthetic controls: rolling rate ratio'))
	abline(h = 1, lty = 2)
	
	matplot(rr_roll_time[, , age_group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l')
	title(main = paste(age_group, 'Interupted time series: rolling rate ratio'))
	abline(h = 1, lty = 2)
	par(par_defaults)
	
	#View cumulative sums
	matplot(plot_cumsum_prevented[, , age_group], type = 'l', col = 'black')
	abline(h = 0, lty = 2)
	title(main = age_group, sub = 'Cumulative cases prevented')
	return()
}))

log_rr_mean_full <- log(rr_mean_full)
log_rr_mean_full[is.na(log_rr_mean_full)] <- 1e6

min_max <- na.omit(log_rr_mean_full[log_rr_mean_full[1] < 100, ])

plotModel(log_rr_mean_full)

################################
#                              #
#     Sensitivity Analyses     #
#                              #
################################

#Start Cluster for Sensitivity Analysis
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(CausalImpact, quietly = TRUE))
clusterExport(cl, c('ds', 'doCausalImpact', 'sensitivityAnalysis', 'age_groups', 'intervention_date', 'outcome', 'time_points', 'n_seasons'))

#Sensitivity Analysis - top weighted variables are excluded and analysis is re-run.
sensitivity_analysis_full <- setNames(parLapply(cl, age_groups, sensitivityAnalysis, covars = covars, impact = impact_full, time_points = time_points, intervention_date = intervention_date, n_seasons = n_seasons), age_groups)

stopCluster(cl)

#Sensitivity Analysis 1 Model results
impact_sensitivity_analysis_1_full <- lapply(sensitivity_analysis_full, function(sensitivity_analysis) {sensitivity_analysis[[1]]$impact})
inclusion_prob_full_analysis_1 <- setNames(mapply(inclusionProb, age_groups, impact_sensitivity_analysis_1_full, SIMPLIFY = FALSE), age_groups)
quantiles_sensitivity_analysis_1_full <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(impact = impact_sensitivity_analysis_1_full[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period)}), age_groups)
rr_mean_sensitivity_analysis_1_full <- t(sapply(quantiles_sensitivity_analysis_1_full, getRR))
colnames(rr_mean_sensitivity_analysis_1_full) <- rr_col_names
write.csv(rr_mean_sensitivity_analysis_1_full, paste(output_directory, country, '_rr_sensitivity_analysis_1_full.csv', sep = ''))

#Table of rate ratios for each sensitivity analysis level
rr_table_full <- t(sapply(age_groups, rrTable, impact = impact_full, sensitivity_analysis = sensitivity_analysis_full, eval_period = eval_period))
write.csv(rr_table_full, paste(output_directory, country, '_rr_table_full.csv', sep = ''))
sensitivity_table <- t(sapply(age_groups, sensitivityTable, sensitivity_analysis = sensitivity_analysis_full))
write.csv(sensitivity_table, paste(output_directory, country, '_sensitivity_table.csv', sep = ''))
sensitivity_table <- read.csv(paste(output_directory, country, '_sensitivity_table.csv', sep = ''), check.names = FALSE, row.names = 1)
rr_table <- cbind(rr_mean_time, rr_table_full, sensitivity_table)
write.csv(rr_table, paste(output_directory, country, '_rr_table.csv', sep = ''))