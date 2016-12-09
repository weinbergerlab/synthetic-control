#This is the function file. It is called directly from the analysis file.

#Log-transform the data.
logTransform <- function(factor_name, factor_value, date_name, start_date, prelog_data) {
	ds <- prelog_data[prelog_data[, factor_name] == factor_value, ]
	ds <- ds[, colSums(is.na(ds)) == 0]
	ds <- ds[match(start_date, ds[, date_name]):nrow(ds), ]
	ds[ds == 0] <- 0.5
	ds[, !(colnames(ds) %in% c(factor_name, date_name))] <- log(ds[, !(colnames(ds) %in% c(factor_name, date_name))])
	return(ds)
}

#Used to adjust the Brazil data for a code shift in 2008.
getTrend <- function(covar_vector, data) {
	new_data <- data
	new_data[c('bs1', 'bs2', 'bs3', 'bs4')] <- 0
	new_data$month_i <- as.factor(1)
	trend <- predict(glm(covar_vector~month_i + ., family = 'gaussian', data = data), type = 'response', newdata = new_data) #month_i is first to be the reference.
	names(trend) <- NULL
	return(trend)
}

#Combine the outcome and covariates.
makeTimeSeries <- function(age_group, outcome, covars, time_points, scale_outcome) {
	if (scale_outcome) {
		ts <- zoo(cbind(outcome = scale(outcome[, age_group]), covars[[age_group]]), time_points)
	} else {
		ts <- zoo(cbind(outcome = outcome[, age_group], covars[[age_group]]), time_points)
	}
	return(ts)
}

#Main analysis function.
doCausalImpact <- function(zoo_data, intervention_date, time_points, n_seasons, trend = FALSE, offset = FALSE) {
	n_iter = 2000
	
	y <- zoo_data[, 1]
	y[time_points >= as.Date(intervention_date)] <- NA
	sd_limit <- sd(y)
	sd <- sd(y, na.rm = TRUE)
	mean <- mean(y, na.rm = TRUE)
	
	post_period_response <- zoo_data[, 1]
	post_period_response <- as.vector(post_period_response[time_points >= as.Date(intervention_date)])
	
	sigma_prior_guess <- 1e-6
	prior_sample_size <- 1e6
	ss <- NA
	ss <- AddSeasonal(list(), y, nseasons = n_seasons, sigma.prior = SdPrior(sigma.guess = sigma_prior_guess, sample.size = prior_sample_size, upper.limit = sd_limit))
	ss <- AddLocalLevel(ss, y, sigma.prior = SdPrior(sigma.guess = sigma_prior_guess, sample.size = prior_sample_size, upper.limit = sd_limit), initial.state.prior = NormalPrior(mean, sd))
	
	if (offset) {
		bsts_model <- bsts(y, state.specification = ss, niter = n_iter, ping = 0, seed = 1)
	} else if (trend){
		x <- zoo_data[, -1] #Removes outcome column from dataset
		bsts_model <- bsts(y~., data = x, state.specification = ss, prior.inclusion.probabilities = c(1.0,1.0), niter = n_iter,  ping = 0, seed = 1)	
	} else {
		x <- zoo_data[, -1] #Removes outcome column from dataset
		regression_prior_df <- 50
		exp_r2 <- 0.8
		n_pred <- max(ncol(x) / 2, 1)
		bsts_model <- bsts(y~., data = x, state.specification = ss, niter = n_iter, expected.model.size = n_pred, prior.df = regression_prior_df, expected.r2 = exp_r2, ping = 0, seed = 1)	
	}
	CausalImpact(bsts.model = bsts_model, post.period.response = post_period_response)
}

#Save inclusion probabilities.
inclusionProb <- function(age_group, impact) {
	return(setNames(as.data.frame(colMeans(impact$model$bsts.model$coefficients != 0)), age_group))
}

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.
rrPredQuantiles <- function(impact, all_cause_data, mean, sd, eval_period, offset = FALSE) {
	burn <- SuggestBurn(0.1, impact$model$bsts.model)
	#Posteriors
	state_samples <- rowSums(aperm(impact$model$bsts.model$state.contributions[-(1:burn),,, drop = FALSE], c(1, 3, 2)), dims = 2)
	sigma_obs <- impact$model$bsts.model$sigma.obs[-(1:burn)]
	#Sample from posterior predictive density over data
	obs_noise_samples <- matrix(rnorm(prod(dim(state_samples)), 0, sigma_obs), nrow = dim(state_samples)[1])
	y_samples <- state_samples + obs_noise_samples
	if (offset) {
		pred_samples_post <- exp(all_cause_data) * t(exp(y_samples * sd + mean))
	} else {
		pred_samples_post <- t(exp(y_samples * sd + mean))
	}
	plot_pred <- t(apply(pred_samples_post, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
	eval_indices <- match(eval_period[1], index(impact$series$response)):match(eval_period[2], index(impact$series$response))
	
	pred_eval_sum <- colSums(pred_samples_post[eval_indices, ])
	if (offset) {
		eval_obs <- sum((exp(all_cause_data) * exp(impact$series$response * sd + mean))[eval_indices])
	} else {
		eval_obs <- sum(exp((impact$series$response[eval_indices] * sd + mean)))
	}
	eval_rr_sum <- eval_obs/pred_eval_sum
	rr <- quantile(eval_rr_sum, probs = c(0.025, 0.5, 0.975))
	mean_rr <- mean(eval_rr_sum)
	
	plot_rr_date_start <- post_period %m-% months(24)  
	
	roll_rr_indices <- match(plot_rr_date_start[1], index(impact$series$response)):match(eval_period[2], index(impact$series$response))
	obs_full <- exp(impact$series$response * sd + mean)
	roll_sum_pred <- apply(pred_samples_post[roll_rr_indices, ], 2, rollsum, align = 'left', k = 12)
	roll_sum_obs <- rollsum(obs_full[roll_rr_indices], align = 'left', k = 12)
	roll_rr_est <- as.data.frame(sweep(1 / roll_sum_pred, 1, as.vector(roll_sum_obs), `*`))
	roll_rr <- t(apply(roll_rr_est, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))

	quantiles <- list(pred_samples_post_full = pred_samples_post, plot_pred = plot_pred, rr = rr, roll_rr = roll_rr, mean_rr = mean_rr)
	return(quantiles)
}

getPred <- function(quantiles) {
	return(quantiles$plot_pred)
}

getRR <- function(quantiles) {
	return(quantiles$rr)
}

#Plot predictions.
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

#Saves the posterior probabilities
postProbHeatmap <- function(data) {
	post_prob <- Reduce(function(a, b) {
		ans <- merge(a, b, by = 'row.names', all = TRUE)
		row.names(ans) <- ans[, 'Row.names']
		ans[, !names(ans) %in% 'Row.names']
	}, data)
	post_prob <- post_prob[complete.cases(post_prob), ]
	my_palette <- colorRampPalette(c('white','black'))(n = 299)
	heatmap(sqrt(as.matrix(post_prob[,1:length(age_groups)])), scale = 'none', col = my_palette)
	return(post_prob)
}

#Sensitivity analysis by dropping the top weighted covariates. 
sensitivityAnalysis <- function(age_group, covars, impact, time_points, intervention_date, n_seasons) {
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
		zoo_data <- zoo(cbind(outcome = outcome[, age_group], covar_df), time_points)
		
		sensitivity_analysis[[i]] <- list(removed_var = max_var, impact = doCausalImpact(zoo_data, intervention_date, time_points, n_seasons))
		
		incl_prob <- plot(sensitivity_analysis[[i]]$impact$model$bsts.model, 'coefficients', cex.names = 0.5, main = paste(age_group, 'Analysis', i))$inclusion.prob
		max_var <- names(incl_prob[length(incl_prob)])
	}
	return(sensitivity_analysis)
}

rrTable <- function(age_group, impact, sensitivity_analysis, eval_period) {
	rr_pred_quantile <- rrPredQuantiles(impact = impact[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period)
	cred_int <- c('Initial' = round(rr_pred_quantile$mean_rr, 4), 'Initial .95' = paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
	cred_int_analyses <- lapply(1:length(sensitivity_analysis[[age_group]]), FUN = function(i) {
		rr_pred_quantile <- rrPredQuantiles(impact = sensitivity_analysis[[age_group]][[i]]$impact, all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period)
		cred_int <- c(sensitivity_analysis[[age_group]][[i]]$removed_var, round(rr_pred_quantile$mean_rr, 4), paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
		names(cred_int) <- c(paste('Removed Variable', i), paste('Analysis', i), paste('Analysis', i, '.95'))
		return(cred_int)
	})
	c(cred_int, cred_int_analyses, recursive = TRUE)
}