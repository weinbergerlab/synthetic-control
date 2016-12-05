#Style Guides: https://www.r-bloggers.com/r-style-guide/, http://adv-r.had.co.nz/Style.html, https://google.github.io/styleguide/Rguide.xml, https://docs.google.com/document/d/1esDVxyWvH8AsX-VJa-8oqWaHLs4stGlIbk8kLc5VlII/edit

rm(list = ls(all = TRUE))
gc()

#Install relevant packages
#install.packages('devtools')
#library(devtools)
#devtools::install_github('google/CausalImpact')

#Include relevant packages
library(parallel, quietly = TRUE)
library(splines, quietly = TRUE)
library(CausalImpact, quietly = TRUE)

n_cores <- detectCores()

set.seed(1)

#Load file
country <- 'SC3'
input_directory <- paste('~/Desktop/PCV\ Project/Data/', country, '/', sep = '')
output_directory <- paste(input_directory, 'results/', sep = '')
dir.create(output_directory, showWarnings = FALSE)
file_name <- paste('prelog', country, 'processed', 'data.csv', sep = '_')
pcv_file <- paste(input_directory, file_name, sep = '')
ds1a <- read.csv(pcv_file)
age_groups <- paste('Age Group', unique(unlist(ds1a$age_group, use.names = FALSE)))

#Account for code-naming differences
if (country == 'SC1' || country == 'SC2' || country == 'SC3') {
	all_cause_name <- 'ACM_NoJ'
	all_cause_pneu_name <- 'J12_18'
	noj_name <- 'J20_22'
} else if (country == 'US') {
	all_cause_name <- 'ach_sid_noresp'
	all_cause_pneu_name <- 'm_ACP'
	noj_name <- 'm_466'
} else {
	all_cause_name <- 'ach_noj'
	all_cause_pneu_name <- 'J12_18'
	noj_name <- 'cJ20_J22'
}

#Date variables
if (country == 'Brazil') {
	data_start_date <- as.Date('2004-01-01')
} else {
	data_start_date <- min(as.Date(ds1a$date))
}
if (country == 'Brazil' || country == 'SC1' || country == 'SC2' || country == 'SC3') {
	data_intervention_date <- as.Date('2009-12-31') 
} else if (country == 'Chile') {
	data_intervention_date <- as.Date('2010-12-31') 
} else if (country == 'Ecuador') {
	data_intervention_date <- as.Date('2009-12-31') 
} else if (country == 'Mexico') {
	data_intervention_date <- as.Date('2005-12-31') 
} else if (country == 'US') {
	data_intervention_date <- as.Date('1999-12-31') 
}
data_end_date <- max(as.Date(ds1a$date))
pre_period <- c(data_start_date, data_intervention_date)
post_period <- c(data_intervention_date + 1, data_end_date)
if (country == 'Brazil' || country == 'SC1' || country == 'SC2' || country == 'SC3') {
	eval_period <- c(as.Date(ds1a$date[nrow(ds1a) - 23]), as.Date(ds1a$date[nrow(ds1a)]))
} else if (country == 'Chile' || country == 'Ecuador') {
	eval_period <- c(as.Date(ds1a$date[nrow(ds1a) - 11]), as.Date(ds1a$date[nrow(ds1a)]))
} else if (country == 'Mexico') {
	eval_period <- c(as.Date('2010-01-01'), as.Date('2011-12-01'))
} else if (country == 'US') {
	eval_period <- c(as.Date('2003-01-01'), as.Date('2004-12-01'))
}

#Seasons
if (country == 'SC1' || country == 'SC2' || country == 'SC3') {
	n_seasons <- 3
} else {
	n_seasons <- 12
}

logTransform <- function(age_group) {
	ds <- ds1a[ds1a$age_group == age_group,]
	ds <- ds[, colSums(is.na(ds)) == 0] #Delete columns with NAs
	ds <- ds[match(data_start_date, ds$date):nrow(ds),]
	ds[ds == 0] <- 0.5
	ds[, 3:ncol(ds)] <- log(ds[, 3:ncol(ds)])
	return(ds)
}

ds1a$date <- as.Date(ds1a$date)
ds <- lapply(unique(ds1a$age_group), FUN = logTransform)
names(ds) <- age_groups
data_start <- match(data_start_date, ds[[1]]$date)

time_points <- ds[[1]]$date[data_start:nrow(ds[[1]])]

getTrend <- function(covar_vector, data) {
	new_data <- data
	new_data[c('bs1', 'bs2', 'bs3', 'bs4')] <- 0
	new_data$month_i <- as.factor(1)
	covariate <- predict(glm(covar_vector~month_i + ., family = 'gaussian', data = data), type = 'response', newdata = new_data) #month_i is first to be the reference.
	return(covariate)
}

covars <- setNames(lapply(age_groups, FUN = function(age_group) {
	if (country == 'Brazil') {
		#Eliminates effects from 2008 coding change
		covars <- ds[[age_group]][data_start:nrow(ds[[age_group]]), 4:ncol(ds[[age_group]])]
		month_i <- as.factor(as.numeric(format(ds[[age_group]]$date[data_start:nrow(ds[[age_group]])], '%m')))
		spline <- setNames(as.data.frame(bs(1:nrow(covars), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
		year_2008 <- numeric(nrow(covars))
		year_2008[1:nrow(covars) >= match(as.Date('2008-01-01'), ds[[age_group]]$date)] <- 1
		data <- cbind.data.frame(year_2008, spline, month_i)
		trend <- lapply(covars, getTrend, data = data)
		covars <- covars - trend
	} else {
		covars <- ds[[age_group]][data_start:nrow(ds[[age_group]]), 4:ncol(ds[[age_group]])]
	}
	if (data_intervention_date > as.Date('2009-09-01')) {
		covars$pandemic <- ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0))
	}
	covars <- as.data.frame(lapply(covars[, apply(covars, 2, var) != 0], scale))
	return(covars)
}), age_groups)

#Handles case where noj_name variable did not have enough data to stay in analysis.
good_covars <- sapply(covars, FUN = function(age_group) {noj_name %in% colnames(age_group)})
covars <- covars[good_covars]
age_groups <- age_groups[good_covars]

covars_noj <- setNames(lapply(covars, FUN = function(covars) {
	covars[, -grep(noj_name, colnames(covars))] #Should remove all smoothness as well
}), age_groups)
covars_ach <- setNames(lapply(covars, FUN = function(covars) {
	covars[, all_cause_name, drop = FALSE]
}), age_groups)
covars_none <- setNames(lapply(covars, FUN = function(covars) {
	as.data.frame(list(constant = rep(1, times = nrow(covars))))
}), age_groups)
covars_time <- setNames(lapply(covars, FUN = function(covars) {
	as.data.frame(list(time_index = 1:nrow(covars)))
}), age_groups)

outcome <- sapply(ds, FUN = function(data) {scale(data[, all_cause_pneu_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, all_cause_pneu_name])})
outcome_sd <- sapply(ds, FUN = function(data) {sd(data[, all_cause_pneu_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset <- sapply(ds, FUN = function(data) {data[, all_cause_pneu_name] - data[, all_cause_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd <- sapply(ds, FUN = function(data) {sd(data[, all_cause_pneu_name] - data[, all_cause_name])})

doCausalImpact <- function(data, offset = FALSE) {
	n_iter = 2000
	
	y <- data[, 1]
	y[time_points >= as.Date(data_intervention_date)] <- NA
	sd_limit <- sd(y)
	sd <- sd(y, na.rm = TRUE)
	mean <- mean(y, na.rm = TRUE)
	
	post_period_response <- data[, 1]
	post_period_response <- as.vector(post_period_response[time_points >= as.Date(data_intervention_date)])
	
	sigma_prior_guess <- 1e-6
	prior_sample_size <- 1e6
	ss <- NA
	ss <- AddSeasonal(list(), y, nseasons = n_seasons, sigma.prior = SdPrior(sigma.guess = sigma_prior_guess, sample.size = prior_sample_size, upper.limit = sd_limit))
	ss <- AddLocalLevel(ss, y, sigma.prior = SdPrior(sigma.guess = sigma_prior_guess, sample.size = prior_sample_size, upper.limit = sd_limit), initial.state.prior = NormalPrior(mean, sd))
	
	if (offset) {
		bsts_model <- bsts(y, state.specification = ss, niter = n_iter, ping = 0, seed = 1)
	} else {
		x <- data[, -1] #Removes outcome column from dataset
		
		regression_prior_df <- 50
		n_pred <- max((ncol(x)/2), 1)
		exp_r2 <- 0.8
		
		bsts_model <- bsts(y~., data = x, state.specification = ss, niter = n_iter, expected.model.size = n_pred, prior.df = regression_prior_df, expected.r2 = exp_r2, ping = 0, seed = 1)	
	}
	CausalImpact(bsts.model = bsts_model, post.period.response = post_period_response)
}

data_full <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = outcome[, age_group], covars[[age_group]]), time_points)}), age_groups)
data_noj  <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = outcome[, age_group], covars_noj[[age_group]]), time_points)}), age_groups)
data_ach <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = outcome[, age_group], covars_ach[[age_group]]), time_points)}), age_groups)
data_none <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = outcome[, age_group], covars_none[[age_group]]), time_points)}), age_groups)
data_time <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = scale(outcome_offset[, age_group]), covars_time[[age_group]]), time_points)}), age_groups)
data_offset <- setNames(lapply(age_groups, FUN = function(age_group) {zoo(cbind(outcome = scale(outcome_offset[, age_group])), time_points)}), age_groups)

#Start Cluster for CausalImpact
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(CausalImpact, quietly = TRUE))
clusterExport(cl, c('ds', 'doCausalImpact', 'data_full', 'data_noj', 'data_ach', 'age_groups', 'data_intervention_date', 'time_points', 'n_seasons'))

impact_full <- setNames(parLapply(cl, data_full, doCausalImpact), age_groups)
impact_noj <- setNames(parLapply(cl, data_noj,   doCausalImpact), age_groups)
impact_ach <- setNames(parLapply(cl, data_ach,   doCausalImpact), age_groups)
impact_none <- setNames(parLapply(cl, data_none, doCausalImpact), age_groups)
impact_time <- setNames(parLapply(cl, data_time, doCausalImpact), age_groups)
impact_offset <- setNames(parLapply(cl, data_offset, doCausalImpact, offset = TRUE), age_groups)

stopCluster(cl)

inclusion_prob_full <- setNames(lapply(age_groups, FUN = function(age_group, impact) {setNames(as.data.frame(colMeans(impact[[age_group]]$model$bsts.model$coefficients != 0)), age_group)}, impact = impact_full), age_groups)
inclusion_prob_noj  <- setNames(lapply(age_groups, FUN = function(age_group, impact) {setNames(as.data.frame(colMeans(impact[[age_group]]$model$bsts.model$coefficients != 0)), age_group)}, impact = impact_noj),  age_groups)
inclusion_prob_ach  <- setNames(lapply(age_groups, FUN = function(age_group, impact) {setNames(as.data.frame(colMeans(impact[[age_group]]$model$bsts.model$coefficients != 0)), age_group)}, impact = impact_ach),  age_groups)
inclusion_prob_none <- setNames(lapply(age_groups, FUN = function(age_group, impact) {setNames(as.data.frame(colMeans(impact[[age_group]]$model$bsts.model$coefficients != 0)), age_group)}, impact = impact_none), age_groups)
inclusion_prob_time <- setNames(lapply(age_groups, FUN = function(age_group, impact) {setNames(as.data.frame(colMeans(impact[[age_group]]$model$bsts.model$coefficients != 0)), age_group)}, impact = impact_time), age_groups)

#Inclusion Probabilities
rrPredQuantiles <- function(data, all_cause_data, mean, sd, offset = FALSE) {
	burn <- SuggestBurn(0.1, data$model$bsts.model)
	#Posteriors
	state_samples <- rowSums(aperm(data$model$bsts.model$state.contributions[-(1:burn),,, drop = FALSE], c(1, 3, 2)), dims = 2)
	sigma_obs <- data$model$bsts.model$sigma.obs[-(1:burn)]
	#Sample from posterior predictive density over data
	obs_noise_samples <- matrix(rnorm(prod(dim(state_samples)), 0, sigma_obs), nrow = dim(state_samples)[1])
	y_samples <- state_samples + obs_noise_samples
	if (offset) {
		pred_samples_post <- exp(all_cause_data) * t(exp(y_samples * sd + mean))
	} else {
		pred_samples_post <- t(exp(y_samples * sd + mean))
	}
	plot_pred <- t(apply(pred_samples_post, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
	eval_indices <- match(eval_period[1], index(data$series$response)):match(eval_period[2], index(data$series$response))
	pred_eval_sum <- colSums(pred_samples_post[eval_indices,])
	if (offset) {
		eval_obs <- sum((exp(all_cause_data) * exp(data$series$response * sd + mean))[eval_indices])
	} else {
		eval_obs <- sum(exp((data$series$response[eval_indices] * sd + mean)))
	}
	eval_rr_sum <- eval_obs/pred_eval_sum
	rr <- quantile(eval_rr_sum, probs = c(0.025, 0.5, 0.975))
	mean_rate_ratio <- mean(eval_rr_sum)
	quantiles <- list(pred_samples_post_full = pred_samples_post, plot_pred = plot_pred, rr = rr, mean_rate_ratio = mean_rate_ratio)
	return(quantiles)
}

#Full model results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_full[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])}), age_groups)
pred_quantiles_full <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_full <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_full) <- c('rr_lcl', 'rr_median', 'rr_ucl')

#NoJ model results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_noj[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])}), age_groups)
pred_quantiles_noj <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_noj <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_noj) <- c('rr_lcl', 'rr_median', 'rr_ucl')

#Ach model results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_ach[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])}), age_groups)
pred_quantiles_ach <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_ach <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_ach) <- c('rr_lcl', 'rr_median', 'rr_ucl')

#No covariates model results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_none[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])}), age_groups)
pred_quantiles_none <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_none <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_none) <- c('rr_lcl', 'rr_median', 'rr_ucl')

#Time trend results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_time[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])}), age_groups)
pred_quantiles_time <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_time <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_time) <- c('rr_lcl', 'rr_median', 'rr_ucl')

#Ach offset results
quantiles <- setNames(lapply(age_groups, FUN = function(age_group) {rrPredQuantiles(data = impact_offset[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_offset_mean[age_group], sd = outcome_offset_sd[age_group], offset = TRUE)}), age_groups)
pred_quantiles_offset <- sapply(quantiles, FUN = function(quantiles) {quantiles$plot_pred}, simplify = 'array')
rr_mean_offset <- t(sapply(quantiles, FUN = function(quantiles) {quantiles$rr}))
colnames(rr_mean_offset) <- c('rr_lcl', 'rr_median', 'rr_ucl')

write.csv(rr_mean_full, paste(output_directory, country, '_rr_mult_full.csv', sep = ''))
write.csv(rr_mean_noj, paste(output_directory, country, '_rr_mult_noj.csv', sep = ''))
write.csv(rr_mean_ach, paste(output_directory, country, '_rr_mult_ach.csv', sep = ''))
write.csv(rr_mean_none, paste(output_directory, country, '_rr_mult_no_covars.csv', sep = ''))
write.csv(rr_mean_time, paste(output_directory, country, '_rr_mult_time_trend.csv', sep = ''))
write.csv(rr_mean_offset, paste(output_directory, country, '_rr_mult_offset.csv', sep = ''))

plot_cumsum_prevented <- sapply(age_groups, FUN = function(age_group) {
	pred_samples_post_full <- quantiles[[age_group]]$pred_samples_post_full
	
	post_period_start <- which(time_points == post_period[1]) 
	post_period_end <- which(time_points == post_period[2]) 
	is_post_period <- which(time_points >= post_period[1])
	is_pre_period <- which(time_points < post_period[1])
	
	#Cumulative sum of prevented cases
	cases_prevented <- pred_samples_post_full - outcome_plot[, age_group]
	cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period,], 2, cumsum)
	cumsum_cases_prevented_pre <- apply(cases_prevented[is_pre_period,, drop = FALSE], 2, cumsum)
	cumsum_cases_prevented_pre[,] <- 0
	cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
	plot_cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}, simplify = 'array')

#Full model predictions
pred_full <- matrix(data = 0, nrow = 0, ncol = 4, dimnames = list(NULL, c('pred_lcl', 'pred_median', 'pred_ucl', 'age_groups')))
for (age_group in age_groups) {
	pred_full <- rbind(pred_full, cbind(pred_quantiles_full[,,age_group], rep(age_group, each = nrow(pred_quantiles_full[,,age_group]))))
}
write.csv(pred_full, paste(output_directory, country, '_pred_mult_full.csv', sep = ''), row.names = FALSE)

plotPred <- function(age_group, data) {
	post_period_start <- which(time_points == post_period[1]) 
	post_period_end <- which(time_points == post_period[2])
	
	min_plot<-min(c(pred_quantiles_ach[,, age_group], pred_quantiles_full[,, age_group], outcome_plot[, age_group]))
	max_plot<-max(c(pred_quantiles_ach[,, age_group], pred_quantiles_full[,, age_group], outcome_plot[, age_group]))
	xx <- c(time_points[post_period_start:post_period_end], rev(time_points[post_period_start:post_period_end]))
	ci_poly <- c(data[post_period_start:post_period_end, 1], rev(data[post_period_start:post_period_end, 3]))
	plot(time_points, data[,2], type = 'l', col = 'darkgray', lwd = 1, bty = 'l', ylim = c(min_plot, max_plot))
	polygon(xx, ci_poly, lty = 0, col = 'lightgray')
	points(time_points, data[,2], type = 'l', col = 'white', lty = 2)
	points(time_points, outcome_plot[, age_group], type = 'l', col = 'black', lwd = 2)
	if (country == 'Brazil') {abline(v = as.Date('2008-01-01'), lty = 2)}
	abline(v = as.Date('2012-05-01'), lty = 2, col = 'lightgray') #pcv13
}

plotModel <- function(model) {
	plot(1:length(age_groups), model[1:length(age_groups), 'rr_median'], bty = 'l', ylim = c(min(min_max), max(min_max)))
	arrows(1:length(age_groups), model[1:length(age_groups), 'rr_lcl'], 1:length(age_groups), model[, 'rr_ucl'], code = 3, angle = 90, length = 0.0)
	abline(h = 0)
}

#Plot covars, models, and impacts
par_defaults <- par(no.readonly = TRUE)
invisible(lapply(age_groups, FUN = function(age_group) {
	#View scaled covariates
	matplot(covars[[age_group]], type = 'l')
	
	#View cumulative sums
	matplot(plot_cumsum_prevented[,, age_group], type='l')
	abline(h = 0)
	
	#Plot predictions
	par(mfrow = c(3, 1))
	#Ach
	plotPred(age_group, pred_quantiles_ach[,, age_group])
	#Full
	plotPred(age_group, pred_quantiles_full[,, age_group])
	par(par_defaults)
}))

#Plot results from offset model and full model
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1, 1)) #bltr
plot(rr_mean_offset[, 2], rr_mean_full[, 2], bty = 'l', xlim = c(0.5, 1.5), ylim = c(0.5, 1.5))
abline(v = 1, lty = 3)
abline(h = 1, lty = 3)
abline(a = 0, b = 1, lty = 2)

par(mfrow = c(3, 1))
par(mar = c(3, 2, 1, 1)) #bltr

log_rr_mean_offset <- log(rr_mean_offset)
log_rr_mean_offset[is.na(log_rr_mean_offset)] <- 1e6

log_rr_mean_full <- log(rr_mean_full)
log_rr_mean_full[is.na(log_rr_mean_full)] <- 1e6

min_max <- na.omit(cbind(log_rr_mean_offset[log_rr_mean_offset[1] < 100,], log_rr_mean_full[log_rr_mean_full[1] < 100,]))

#Offset Model
plotModel(log_rr_mean_offset)
#Full model
plotModel(log_rr_mean_full)

#Heatmaps
postProbHeatmap <- function(data) {
	post_prob <- Reduce(function(a, b) {
		ans <- merge(a, b, by = 'row.names', all = TRUE)
		row.names(ans) <- ans[, 'Row.names']
		ans[, !names(ans) %in% 'Row.names']
	}, data)
	post_prob <- post_prob[complete.cases(post_prob),]
	my_palette <- colorRampPalette(c('white','black'))(n = 299)
	heatmap(sqrt(as.matrix(post_prob[,1:length(age_groups)])), scale = 'none', col = my_palette)
	return(post_prob)
}

#All posterior probabilities for variables for all age groups
post_prob_full <- postProbHeatmap(inclusion_prob_full)
post_prob_noj <- postProbHeatmap(inclusion_prob_noj)
write.csv(post_prob_full, paste(output_directory, country, '_post_prob_full.csv', sep = ''))
write.csv(post_prob_noj, paste(output_directory, country, '_post_prob_noj.csv', sep = ''))

sensitivityAnalysis <- function(age_group, covars, impact) {
	par(mar = c(5, 4, 1, 2) + 0.1)
	covar_df <- covars[[age_group]]
	df <- ds[[age_group]]
	
	incl_prob <- plot(impact[[age_group]]$model$bsts.model, 'coefficients', cex.names = 0.5, main = age_group)$inclusion.prob 
	max_var <- names(incl_prob[length(incl_prob)])

	sensitivity_analysis <- vector('list', 3)
	
	for (i in 1:3) {
		df <- df[, names(df) != max_var]
		covar_df <- covar_df[, names(covar_df) != max_var]
		
		#Combine covars, outcome, date
		data <- zoo(cbind(outcome = outcome[, age_group], covar_df), time_points)
		
		sensitivity_analysis[[i]] <- doCausalImpact(data)
		
		incl_prob <- plot(sensitivity_analysis[[i]]$model$bsts.model, 'coefficients', cex.names = 0.5, main = paste(age_group, 'Analysis', i))$inclusion.prob
		max_var <- names(incl_prob[length(incl_prob)])
	}
	return(sensitivity_analysis)
}

#Start Cluster for Sensitivity Analysis
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(CausalImpact, quietly = TRUE))
clusterExport(cl, c('ds', 'doCausalImpact', 'sensitivityAnalysis', 'age_groups', 'data_intervention_date', 'outcome', 'time_points', 'n_seasons'))

#Sensitivity Analysis
sensitivity_analysis_full <- setNames(parLapply(cl, age_groups, sensitivityAnalysis, covars = covars, impact = impact_full), age_groups)
sensitivity_analysis_noj <- setNames(parLapply(cl, age_groups, sensitivityAnalysis, covars = covars_noj, impact = impact_noj), age_groups)

stopCluster(cl)

rr_table <- function(age_group, impact, sensitivity_analysis) {
	rr_pred_quantile <- rrPredQuantiles(data = impact[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])
	cred_int <- c('Initial' = round(rr_pred_quantile$mean_rate_ratio, 4), 'Initial .95' = paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
	cred_int_analyses <- lapply(1:length(sensitivity_analysis[[age_group]]), FUN = function(i) {
		rr_pred_quantile <- rrPredQuantiles(data = sensitivity_analysis[[age_group]][[i]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group])
		cred_int <- c(round(rr_pred_quantile$mean_rate_ratio, 4), paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
		names(cred_int) <- c(paste('Analysis', i), paste('Analysis', i, '.95'))
		return(cred_int)
	})
	c(cred_int, cred_int_analyses, recursive = TRUE)
}

#Table of Rate Ratios for each Analysis Level
rr_table_full <- t(sapply(age_groups, rr_table, impact = impact_full, sensitivity_analysis = sensitivity_analysis_full))
rr_table_noj <- t(sapply(age_groups, rr_table, impact = impact_noj, sensitivity_analysis = sensitivity_analysis_noj))
write.csv(rr_table_full, paste(output_directory, country, '_rr_table_full.csv', sep = ''))
write.csv(rr_table_noj, paste(output_directory, country, '_rr_table_noj.csv', sep = ''))
rr_table_full
rr_table_noj

#Weight Checker - counts how often a covariate appears with an inclusion probability above 5% of that of the covariate with the greatest inclusion probability.
all_vars <- unique(c(colnames(ds1a), unlist(sapply(covars, colnames))))
weight_checker <- rep(0, length = length(all_vars) - 2)
names(weight_checker) <- all_vars[3:length(all_vars)]
for (iteration in 1:10) {
	print(paste('Iteration', iteration, 'of', 10))
	invisible(lapply(age_groups, FUN = function(age_group) {
		par(mar = c(5, 4, 1, 2) + 0.1)
		covar_df <- covars[[age_group]]
		df <- ds[[age_group]]
		
		incl_prob <- plot(impact_full[[age_group]]$model$bsts.model, 'coefficients', cex.names = 0.5, main = age_group)$inclusion.prob 
		max_var <- names(incl_prob[length(incl_prob)])
		
		weight_checker <<- weight_checker + ifelse(is.na(incl_prob[names(weight_checker)]), 0, 
																							 ifelse(incl_prob[names(weight_checker)] >= 0.05*max(incl_prob), 1, 0)
		)
		
		for (i in 1:3) {
			df <- df[, names(df) != max_var]
			covar_df <- covar_df[, names(covar_df) != max_var]
			
			#Combine covars, outcome, date
			data <- zoo(cbind(outcome = outcome[, age_group], covar_df), time_points)
			impact <- doCausalImpact(data)
			
			incl_prob <- plot(impact$model$bsts.model, 'coefficients', cex.names = 0.5, main = paste(age_group, 'Analysis', i))$inclusion.prob
			max_var <- names(incl_prob[length(incl_prob)])
			
			#For weight checking: 
			weight_checker <<- weight_checker + ifelse(is.na(incl_prob[names(weight_checker)]), 0, 
																								 ifelse(incl_prob[names(weight_checker)] >= 0.05 * max(incl_prob), 1, 0)
			)
		}
	}))
}

barplot(weight_checker[order(weight_checker, decreasing = TRUE)], las = 2)
weight_check_table <- matrix(weight_checker[order(weight_checker, decreasing = TRUE)], dimnames = list(names(weight_checker)[order(weight_checker, decreasing = TRUE)], c('Counts')))
weight_check_table