#This is the function file. It is called directly from the analysis file.

packageHandler <- function(packages, update_packages = TRUE) {
	bad_packages <- list()
	for (package in packages) {
		tryCatch({
			find.package(package)
		}, error = function(e) {
			if (package == 'CausalImpact') {
				packageHandler('devtools')
				CI_desc <- do.call(rbind, strsplit(readLines('https://raw.githubusercontent.com/google/CausalImpact/master/DESCRIPTION'), split = ': '))
				rownames(CI_desc) <- CI_desc[, 1]
				CI_desc <- as.list(CI_desc[c('Package', 'Date', 'Copyright', 'Version', 'License', 'Imports', 'Depends', 'Suggests'), -1], rownames(CI_desc))
				CI_desc[c('Imports', 'Depends', 'Suggests')] <- lapply(CI_desc[c('Imports', 'Depends', 'Suggests')], FUN = function(string) {strsplit(string, split = ', ')[[1]]})
				sapply(c(CI_desc$Imports, CI_desc$Depends, CI_desc$Suggests), FUN = packageHandler, update_packages = update_packages)
				devtools::install_github('google/CausalImpact')
			} else {
				if (package %in% available.packages()) {
					install.packages(package)
				} else {
					bad_packages <<- append(bad_packages, package)
				}
			}
		}, warning = function(w) {
			paste(w, 'Shouldn\'t be here.')
		}, finally = {
			if (update_packages) {
				if (package == 'CausalImpact') {
					CI_desc <- do.call(rbind, strsplit(readLines('https://raw.githubusercontent.com/google/CausalImpact/master/DESCRIPTION'), split = ': '))
					rownames(CI_desc) <- CI_desc[, 1]
					CI_desc <- as.list(CI_desc[c('Package', 'Date', 'Copyright', 'Version', 'License', 'Imports', 'Depends', 'Suggests'), -1], rownames(CI_desc))
					if (packageVersion('CausalImpact') != CI_desc$Version) {
						devtools::install_github('google/CausalImpact')
					}
				} else {
					update.packages(package)
				}	
			}
		})
	}
	if (length(bad_packages) > 0) {
		if (length(bad_packages) == 1) {
			stop(paste('Package', paste('"', bad_packages, '"', sep = ''), 'is not available for', paste(version$version.string, '.', sep = '')))
		} else {
			stop(paste('Packages', paste(lapply(bad_packages, function(bad_package) {paste('"', bad_package, '"', sep = '')}), collapse = ', '), 'are not available for', paste(version$version.string, '.', sep = '')))
		}
	}
	return()
}

#Log-transform the data.
logTransform <- function(factor_name, factor_value, date_name, all_cause_name, all_cause_pneu_name, start_date, prelog_data) {
	ds <- prelog_data[prelog_data[, factor_name] == factor_value, ]
	ds <- ds[, colSums(is.na(ds)) == 0]
	ds <- ds[match(start_date, ds[, date_name]):nrow(ds), ]
	ds <- cbind(ds[, colnames(ds) %in% c(factor_name, date_name, all_cause_name, all_cause_pneu_name)], filterSparse(ds[, !(colnames(ds) %in% c(factor_name, date_name, all_cause_name, all_cause_pneu_name))]))
	ds[ds == 0] <- 0.5
	ds[, !(colnames(ds) %in% c(factor_name, date_name))] <- log(ds[, !(colnames(ds) %in% c(factor_name, date_name))])
	return(ds)
}

filterSparse <- function(dataset, threshold = 5) {
	return(dataset[, colMeans(dataset) > threshold, drop = FALSE])
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

impactExtract <- function(impact) {
	burn <- SuggestBurn(0.1, impact$model$bsts.model)
	#Posteriors
	state_samples <- rowSums(aperm(impact$model$bsts.model$state.contributions[-(1:burn), , , drop = FALSE], c(1, 3, 2)), dims = 2)
	sigma_obs <- impact$model$bsts.model$sigma.obs[-(1:burn)]
	#Sample from posterior predictive density over data
	obs_noise_samples <- matrix(rnorm(prod(dim(state_samples)), 0, sigma_obs), nrow = dim(state_samples)[1])
	y_samples <- state_samples + obs_noise_samples
	
	inclusion_probs <- sort(colMeans(impact$model$bsts.model$coefficients != 0))
	return(list(y_samples = y_samples, series = impact$series, inclusion_probs = inclusion_probs))
}

#Main analysis function.
doCausalImpact <- function(zoo_data, intervention_date, time_points, n_seasons = NULL, n_pred = 5, n_iter = 10000, trend = FALSE) {
	if (is.null(n_seasons) || is.na(n_seasons)) {
		n_seasons <- length(unique(month(time(zoo_data)))) #number of months
	}
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

	if (trend) {
		x <- zoo_data[, -1] #Removes outcome column from dataset
		bsts_model <- bsts(y~., data = x, state.specification = ss, prior.inclusion.probabilities = c(1.0, 1.0), niter = n_iter, ping = 0, seed = 1)	
	} else {
		x <- zoo_data[, -1] #Removes outcome column from dataset
		regression_prior_df <- 50
		exp_r2 <- 0.8
		bsts_model <- bsts(y~., data = x, state.specification = ss, niter = n_iter, expected.model.size = n_pred, prior.df = regression_prior_df, expected.r2 = exp_r2, ping = 0, seed = 1)
	}
	impact <- CausalImpact(bsts.model = bsts_model, post.period.response = post_period_response)
	colnames(impact$model$bsts.model$coefficients)[-1] <- names(zoo_data)[-1]
	impact_extract <- impactExtract(impact)
	return(impact_extract)
}

#Save inclusion probabilities.
inclusionProb <- function(age_group, impact) {
	return(setNames(as.data.frame(impact$inclusion_probs), age_group))
}

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.
rrPredQuantiles <- function(impact, all_cause_data = NULL, mean, sd, eval_period, post_period, trend = FALSE) {
	if (trend) {
		pred_samples_post <- exp(all_cause_data) * t(exp(impact$y_samples * sd + mean))
	} else {
		pred_samples_post <- t(exp(impact$y_samples * sd + mean))	
	}
	
	plot_pred <- t(apply(pred_samples_post, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
	eval_indices <- match(eval_period[1], index(impact$series$response)):match(eval_period[2], index(impact$series$response))
	
	pred_eval_sum <- colSums(pred_samples_post[eval_indices, ])
	
	if (trend) {
		eval_obs <- sum((exp(all_cause_data) * exp(impact$series$response * sd + mean))[eval_indices])
	} else {
		eval_obs <- sum(exp((impact$series$response[eval_indices] * sd + mean)))	
	}
	
	eval_rr_sum <- eval_obs/pred_eval_sum
	rr <- quantile(eval_rr_sum, probs = c(0.025, 0.5, 0.975))
	mean_rr <- mean(eval_rr_sum)
	
	plot_rr_date_start <- post_period %m-% months(24)
	roll_rr_indices <- match(plot_rr_date_start[1], index(impact$series$response)):match(eval_period[2], index(impact$series$response))
	if (trend) {
		obs_full <- exp(all_cause_data) * exp(impact$series$response * sd + mean)
	} else {
		obs_full <- exp(impact$series$response * sd + mean)
	}
	roll_sum_pred <- roll_sum(pred_samples_post[roll_rr_indices, ], 12)
	roll_sum_obs <- roll_sum(obs_full[roll_rr_indices], 12)
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
plotPred <- function(pred_quantiles, time_points, post_period, pred_quantiles_full, outcome_plot, country, title = NULL, sensitivity_pred_quantiles = NULL, sensitivity_title = 'Sensitivity Plots') {
	
	post_period_start <- which(time_points == post_period[1]) 
	post_period_end <- which(time_points == post_period[2])
	min_plot <- min(c(pred_quantiles_full, outcome_plot))
	max_plot <- max(c(pred_quantiles_full, outcome_plot))
	
	xx <- c(time_points[post_period_start:post_period_end], rev(time_points[post_period_start:post_period_end]))
	ci_poly <- c(pred_quantiles[post_period_start:post_period_end, 1], rev(pred_quantiles[post_period_start:post_period_end, 3]))
	plot(time_points, pred_quantiles[, 2], type = 'l', col = 'darkgray', lwd = 1, bty = 'l', ylim = c(min_plot, max_plot), main = title, xlab = 'Time', ylab = 'Number of Cases')
	polygon(xx, ci_poly, lty = 0, col = 'lightgray')
	points(time_points, pred_quantiles[, 2], type = 'l', col = 'white', lty = 2)
	points(time_points, outcome_plot, type = 'l', col = 'black', lwd = 2)
	if (country == 'Brazil') {abline(v = as.Date('2008-01-01'), lty = 2)}
	#abline(v = as.Date('2012-05-01'), lty = 2, col = 'lightgray') #pcv13
	
	if (!is.null(sensitivity_pred_quantiles) && length(sensitivity_pred_quantiles) == 3) {
		plot(time_points, pred_quantiles[, 2], type = 'l', col = 'darkgray', lwd = 1, bty = 'l', ylim = c(min_plot, max_plot), main = sensitivity_title, xlab = 'Time', ylab = 'Number of Cases')
		points(time_points, outcome_plot, type = 'l', col = 'black', lwd = 2)
		
		rgb_line <- c('red' = rgb(1, 0, 0, 1), 'green' = rgb(0, 1, 0, 1), 'blue' = rgb(0, 0, 1, 1))
		rgb_poly <- c('red' = rgb(1, 0, 0, 0.1), 'green' = rgb(0, 1, 0, 0.1), 'blue' = rgb(0, 0, 1, 0.1))
		
		#Due to plotting problems, a for loop can not be used here.
		points(time_points, sensitivity_pred_quantiles[[1]][, 2], type = 'l', col = rgb_line[1], lty = 1)
		points(time_points, sensitivity_pred_quantiles[[2]][, 2], type = 'l', col = rgb_line[2], lty = 1)
		points(time_points, sensitivity_pred_quantiles[[3]][, 2], type = 'l', col = rgb_line[3], lty = 1)
		
		ci_poly <- c(sensitivity_pred_quantiles[[1]][post_period_start:post_period_end, 1], rev(sensitivity_pred_quantiles[[1]][post_period_start:post_period_end, 3]))
		polygon(xx, ci_poly, lty = 0, col = rgb_poly[1])
		ci_poly <- c(sensitivity_pred_quantiles[[2]][post_period_start:post_period_end, 1], rev(sensitivity_pred_quantiles[[2]][post_period_start:post_period_end, 3]))
		polygon(xx, ci_poly, lty = 0, col = rgb_poly[2])
		ci_poly <- c(sensitivity_pred_quantiles[[3]][post_period_start:post_period_end, 1], rev(sensitivity_pred_quantiles[[3]][post_period_start:post_period_end, 3]))
		polygon(xx, ci_poly, lty = 0, col = rgb_poly[3])
		
		legend('topleft', legend = c('Outcome', 'Estimate', 'Sensitivity 1', 'Sensitivity 2', 'Sensitivity 3'), col = c('black', 'darkgray', rgb_line), lty = 1, cex = 0.8)
	}
}

plotModel <- function(model, age_groups, min_max) {
	plot(1:length(age_groups), model[1:length(age_groups), 'Point Estimate'], bty = 'l', ylim = c(min(min_max), max(min_max)))
	arrows(1:length(age_groups), model[1:length(age_groups), 'Lower CI'], 1:length(age_groups), model[, 'Upper CI'], code = 3, angle = 90, length = 0.0)
	abline(h = 0)
}

#Sensitivity analysis by dropping the top weighted covariates. 
weightSensitivityAnalysis <- function(age_group, covars, ds, impact, time_points, intervention_date, n_seasons, mean = NULL, sd = NULL, eval_period = NULL, post_period = NULL) {
	par(mar = c(5, 4, 1, 2) + 0.1)
	covar_df <- covars[[age_group]]
	df <- ds[[age_group]]
	
	incl_prob <- impact[[age_group]]$inclusion_probs
	max_var <- names(incl_prob)[length(incl_prob)]
	max_prob <- incl_prob[length(incl_prob)]
	sensitivity_analysis <- vector('list', 3)
	
	for (i in 1:3) {
		df <- df[, names(df) != max_var]
		covar_df <- covar_df[, names(covar_df) != max_var]
		#Combine covars, outcome, date
		zoo_data <- zoo(cbind(outcome = outcome[, age_group], covar_df), time_points)
		impact <- doCausalImpact(zoo_data, intervention_date, time_points, n_seasons)

		sensitivity_analysis[[i]] <- list(removed_var = max_var, removed_prob = max_prob, impact = impact)
		if (!is.null(mean) && !is.null(sd) && !is.null(eval_period) && !is.null(post_period)) {
			quantiles <- rrPredQuantiles(impact = sensitivity_analysis[[i]]$impact, mean = mean[age_group], sd = sd[age_group], eval_period = eval_period, post_period = post_period)
			sensitivity_analysis[[i]]$rr <- quantiles$rr
			sensitivity_analysis[[i]]$plot_pred <- quantiles$plot_pred
			names(sensitivity_analysis[[i]]$rr) <- c('Lower CI', 'Point Estimate', 'Upper CI')
			print(sensitivity_analysis[[i]]$rr)
		}
		
		incl_prob <- impact$inclusion_probs
		max_var <- names(incl_prob)[length(incl_prob)]
		max_prob <- incl_prob[length(incl_prob)]
	}
	return(sensitivity_analysis)
}

predSensitivityAnalysis <- function(group, ds, zoo_data, intervention_date, time_points, n_seasons, n_pred) {
	impact <- doCausalImpact(zoo_data[[group]], intervention_date, time_points, n_seasons, n_pred = n_pred)
	quantiles <- lapply(group, FUN = function(group) {rrPredQuantiles(impact = impact, all_cause_data = ds[[group]][, all_cause_name], mean = outcome_mean[group], sd = outcome_sd[group], eval_period = eval_period, post_period = post_period)})
	pred_quantiles_pred <- sapply(quantiles, getPred, simplify = 'array')
	
	rr_mean <- t(sapply(quantiles, getRR))
	rr_col_names <- c('Lower CI', 'Point Estimate', 'Upper CI')
	colnames(rr_mean) <- rr_col_names
	
	return(rr_mean)
}

rrTable <- function(age_group, impact, outcome_mean, outcome_sd, eval_period, post_period) {
	rr_pred_quantile <- rrPredQuantiles(impact = impact[[age_group]], all_cause_data = ds[[age_group]][, all_cause_name], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = eval_period, post_period = post_period)
	cred_int <- c(rr_pred_quantile$rr[1], rr_pred_quantile$rr[2], rr_pred_quantile$rr[3])
	rr_col_names <- c('Lower CI', 'Point Estimate', 'Upper CI')
	names(cred_int) <- paste('SC', rr_col_names)
	return(cred_int)
}

sensitivityTable <- function(age_group, sensitivity_analysis, original_rr = NULL) {
	top_controls <- lapply(1:length(sensitivity_analysis[[age_group]]), FUN = function(i) {
		top_control <- c(sensitivity_analysis[[age_group]][[i]]$removed_var, sensitivity_analysis[[age_group]][[i]]$removed_prob, sensitivity_analysis[[age_group]][[i]]$rr)
		names(top_control) <- c(paste('Top Control', i), paste('Inclusion Probability of Control', i), paste(names(sensitivity_analysis[[age_group]][[i]]$rr), i))
		return(top_control)
	})
	sensitivity_table <- c(original_rr[age_group, ], c(top_controls, recursive = TRUE))
	return(sensitivity_table)
}