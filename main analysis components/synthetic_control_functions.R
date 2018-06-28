#This is the function file. It is called directly from the analysis file.
packageHandler <- function(packages, update_packages = TRUE, install_packages = TRUE) {
	bad_packages <- list()
	for (package in packages) {
		if (install_packages) {
			tryCatch({
				find.package(package)
			}, error = function(e) {
									if (package %in% available.packages()) {
						install.packages(package, repos = 'http://cran.rstudio.com/')
					} else {
						bad_packages <<- append(bad_packages, package)
					}
			}, warning = function(w) {
				paste(w, 'Shouldn\'t be here.')
			}, finally = {
				if (update_packages) {
					
						update.packages(package, repos = 'http://cran.rstudio.com/')
					}	
			})
		}
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

#Rearrange date to YYYY-MM-DD format.
formatDate <- function(time_points) {
	time_points <- as_date(time_points)
	time_points <- as.Date(time_points, format = '%Y-%m-%d')
	return(time_points)
}

splitGroup <- function(ungrouped_data, group_name, group, date_name, start_date, end_date, no_filter = NULL) {
	ds <- ungrouped_data[ungrouped_data[, group_name] == group, ]
	ds <- ds[, colSums(is.na(ds)) == 0]
	ds <- ds[match(start_date, ds[, date_name]):match(end_date, ds[, date_name]), ]
	ds <- cbind(ds[, colnames(ds) %in% no_filter], filterSparse(ds[, !(colnames(ds) %in% no_filter)]))
	return(ds)
}

#Log-transform the covariate
logTransform <- function(prelog_data, no_log = NULL) {
  prelog_data[, !(colnames(prelog_data) %in% no_log)][prelog_data[, !(colnames(prelog_data) %in% no_log)] == 0] <- 0.5
  prelog_data[, !(colnames(prelog_data) %in% no_log)] <- log(prelog_data[, !(colnames(prelog_data) %in% no_log)])
  return(prelog_data)
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

makeCovars <- function(ds_group, code_change, intervention_date,season.dummies, time_points) {
	if (code_change) {
		#Eliminates effects from 2008 coding change
		covars <- ds_group[, 4:ncol(ds_group)]
		month_i <- as.factor(as.numeric(format(time_points, '%m')))
		spline <- setNames(as.data.frame(bs(1:nrow(covars), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
		year_2008 <- numeric(nrow(covars))
		year_2008[1:nrow(covars) >= match(as.Date('2008-01-01'), time_points)] <- 1
		data <- cbind.data.frame(year_2008, spline, month_i)
		trend <- lapply(covars, getTrend, data = data)
		covars <- covars - trend
	} else {
		covars <- ds_group[, 4:ncol(ds_group), drop = FALSE]
	}
	if (intervention_date > as.Date('2009-09-01')) {
		covars$pandemic <- ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0))
	}
	covars <- as.data.frame(lapply(covars[, apply(covars, 2, var) != 0, drop = FALSE], scale), check.names = FALSE)
	covars<-cbind(season.dummies,covars)
	return(covars)
}

#Combine the outcome and covariates.
makeTimeSeries <- function(group, outcome,  covars, trend=FALSE) {
  if(trend==FALSE){	return(cbind(outcome = outcome[, group],  covars[[group]]))}
  if(trend==TRUE){	return(cbind(outcome = outcome[, group],log.offset=log(offset[,group]+0.5),  covars[[group]]))}
}


#Main analysis function.
doCausalImpact <- function(zoo_data, intervention_date, time_points, n_seasons = NULL, n_pred = 5, n_iter = 10000, trend = FALSE) {
	if (is.null(n_seasons) || is.na(n_seasons)) {
		n_seasons <- length(unique(month(time(zoo_data)))) #number of months
	}
	y.pre <- zoo_data[time_points < as.Date(intervention_date), 1]
	y<-zoo_data[,1] #all y
	
	
	post_period_response <- zoo_data[, 1]
	post_period_response <- as.vector(post_period_response[time_points >= as.Date(intervention_date)])
	x <-as.matrix(zoo_data[, -1]) #Removes outcome column from dataset
	x.pre <-as.matrix(zoo_data[time_points < as.Date(intervention_date), -1]) #Removes outcome column from dataset
	
	post_period_response <- zoo_data[, 1]
	post_period_response <- as.vector(post_period_response[time_points >= as.Date(intervention_date)])
	regression_prior_df <- 50
	#instead try pogit packages?
	exp_r2 <- 0.8
   	if (trend) {
        #set inclusion prob=1 for all variab;es
    	  prior1=SpikeSlabPrior(cbind(1,x.pre), y=y.pre, prior.inclusion.probabilities =rep(1,14),prior.df = regression_prior_df, expected.r2 = exp_r2 )
    		bsts_model <- poisson.spike(y.pre~x.pre, prior=prior1 , niter = n_iter, ping = 0, seed = 1)	
    		covars<-x
		 }else {

			  denom <- ncol(x)-(n_seasons-1)
			  if(denom>0){
			  	prior.inclusion.probabilities = c( rep(1,n_seasons),  rep(n_pred/denom,denom) ) #force seasonality and intercept into model, repeat '1' 12 times, repeat inclusion prob by N of cnon-monthly covars
			  } else {	
			    prior.inclusion.probabilities = rep(1,12)   
			  }
  			  prior.inclusion.probabilities[prior.inclusion.probabilities>1] <- 1
  			  prior2=SpikeSlabPrior(cbind(1,x.pre),y=y.pre, prior.inclusion.probabilities = prior.inclusion.probabilities,prior.df = regression_prior_df, expected.r2 = exp_r2 )
  			  bsts_model <- poisson.spike(y.pre ~ x.pre,  niter = n_iter, prior=prior2 , ping = 0, seed = 1 )
  			  covars<-x
		  	}
	predict.bsts<-predict(bsts_model, newdata=cbind(rep(1, times=nrow(x)),x), burn=n_iter*0.1, mean.only=FALSE)
	beta.mat<-bsts_model$beta[-seq(from=1, to=n_iter*0.1, by=1),]

	coef.bsts<-SummarizeSpikeSlabCoefficients(bsts_model$beta, burn=n_iter*0.1, order=FALSE)
	inclusion_probs<- coef.bsts[,5]
	names(inclusion_probs)<-substring(names(inclusion_probs),2)
	impact <- list(covars,beta.mat,predict.bsts,inclusion_probs, post.period.response = post_period_response, observed.y=zoo_data[, 1])
	names(impact)<-c('covars','beta.mat','predict.bsts','inclusion_probs','post_period_response', 'observed.y' )
	return(impact)
}

#Save inclusion probabilities.
inclusionProb <- function(impact) {
	return(impact$inclusion_probs)
}

waic_fun<-function(impact,  eval_period, post_period, trend = FALSE) {
  covars.pre<-impact$covars[time_points < as.Date(intervention_date),]
  covars.pre<-cbind(rep(1,nrow(covars.pre)), covars.pre)
  y.pre<-impact$observed.y[time_points < as.Date(intervention_date)]
  reg.mean<-   exp(t(covars.pre %*% t(impact$beta.mat)))
  log.piece<-matrix(NA, nrow=nrow(reg.mean), ncol=ncol(reg.mean))
  for(j in 1:nrow(reg.mean)){
    log.piece[j,]<-dpois(y.pre, lambda=reg.mean[j,], log=TRUE)
  }
  llpd.part<-apply(log.piece,2, logmeanexp) #logmeanexp AVOIDS underflow issue
  llpd<-sum(  llpd.part)
  #PWAIC_1_piece1<-apply(log.piece,2, logmeanexp) #logmeanexp AVOIDS underflow issue
  PWAIC_1<-2*sum(llpd.part - colMeans(log.piece))
  WAIC_1<- -2*(llpd-PWAIC_1)

  temp<-rep(0,times=ncol(log.piece))
  for(j in 1:ncol(log.piece)){
    temp[j]<-var(log.piece[,j])
  }
  PWAIC_2<-sum(temp)
  WAIC_2<- -2*(llpd-PWAIC_2)
  #save output
  waics<-c(WAIC_1, WAIC_2)
  names(waics)<-c('waic_1','waic_2')
  return(waics)
}

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.

#Estimate the rate ratios during the evaluation period and return to the original scale of the data.
rrPredQuantiles <- function(impact, denom_data = NULL,  eval_period, post_period) {
  
    pred_samples <- impact$predict.bsts  
  
  pred <- t(apply(pred_samples, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  eval_indices <- match(which(time_points==eval_period[1]), (1:length(impact$observed.y))):match(which(time_points==eval_period[2]), (1:length(impact$observed.y)))
  
  pred_eval_sum <- colSums(pred_samples[eval_indices, ])
  
  
    eval_obs <- sum(impact$observed.y[eval_indices] )

  eval_rr_sum <- eval_obs/pred_eval_sum
  rr <- quantile(eval_rr_sum, probs = c(0.025, 0.5, 0.975))
  names(rr) <- c('Lower CI', 'Point Estimate', 'Upper CI')
  mean_rr <- mean(eval_rr_sum)
  
  plot_rr_start <- which(time_points==post_period[1]) - n_seasons
  roll_rr_indices <- match(plot_rr_start, (1:length(impact$observed.y))):match(which(time_points==eval_period[2]), (1:length(impact$observed.y)))
 
    obs_full <- impact$observed.y 
  
  roll_sum_pred <- roll_sum(pred_samples[roll_rr_indices, ], n_seasons)
  roll_sum_obs <- roll_sum(obs_full[roll_rr_indices], n_seasons)
  roll_rr_est <- roll_sum_obs/ roll_sum_pred
  roll_rr <- t(apply(roll_rr_est, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
	
     pred_samples_post<-pred_samples[eval_indices, ]
  
# obs_full[obs_full==0]<-0.5 #continuity correction for small sampls
   #  pred_quantiles<-t(apply(pred_samples, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
   #  matplot(impact$observed.y/pred_quantiles, type='l')
  log_rr_full_t_samples <- t(log(obs_full/pred_samples) )
  log_rr_full_t_quantiles<-t(apply(log_rr_full_t_samples, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  log_rr_full_t_sd<-t(apply(log_rr_full_t_samples, 2, sd, na.rm = TRUE))
# 	
# 	#Covariance matrix for pooled analysis
 	log_rr_full_t_samples.covar<-cov(log_rr_full_t_samples)
 	#post.indices<- which(time_points==post_period[1]):which(time_points==post_period[2])
 	#log_rr_full_t_samples.prec.post<-solve(log_rr_full_t_samples.covar[post.indices,post.indices]) #NOT INVERTIBLE?
# 	
  # quantiles <- list(pred_samples_post_full = pred_samples_post,roll_rr=roll_rr, log_rr_full_t_samples.prec=log_rr_full_t_samples.prec, log_rr_full_t_samples=log_rr_full_t_samples,log_rr_full_t_quantiles=log_rr_full_t_quantiles,log_rr_full_t_sd=log_rr_full_t_sd, plot_pred = plot_pred,log_plot_pred=log_plot_pred, log_plot_pred_SD=log_plot_pred_SD, rr = rr, mean_rate_ratio = mean_rate_ratio,rr.iter=rr.iter)
 # quantiles <- list(pred_samples = pred_samples, pred = pred, rr = rr, roll_rr = roll_rr, mean_rr = mean_rr)
   quantiles <- list(pred_samples = pred_samples, pred = pred, rr = rr, roll_rr = roll_rr, mean_rr = mean_rr, pred_samples_post_full = pred_samples_post,roll_rr=roll_rr, log_rr_full_t_quantiles=log_rr_full_t_quantiles,log_rr_full_t_sd=log_rr_full_t_sd, rr = rr)
   return(quantiles)
}

getPred <- function(quantiles) {
  return(quantiles$pred)
}

getRR <- function(quantiles) {
  return(quantiles$rr)
}

makeInterval <- function(point_estimate, upper_interval, lower_interval, digits = 2) {
  return(paste(round(as.numeric(point_estimate), digits), ' (', round(as.numeric(lower_interval), digits), ', ', round(as.numeric(upper_interval), digits), ')', sep = ''))
}

#Plot predictions.
plotPred <- function(pred_quantiles, time_points, post_period, ylim, outcome_plot, title = NULL, sensitivity_pred_quantiles = NULL, sensitivity_title = 'Sensitivity Plots', plot_sensitivity = FALSE) {
  
  post_period_start <- which(time_points == post_period[1]) 
  post_period_end <- which(time_points == post_period[2])
  post_dates <- c(time_points[post_period_start:post_period_end], rev(time_points[post_period_start:post_period_end]))
  
  if (!plot_sensitivity) {
    pred_plot <- ggplot() + 
      geom_polygon(data = data.frame(time = c(post_dates, rev(post_dates)), pred_bound = c(pred_quantiles[which(time_points %in% post_dates), 3], rev(pred_quantiles[which(time_points %in% post_dates), 1]))), aes_string(x = 'time', y = 'pred_bound'), alpha = 0.3) +
      geom_line(data = data.frame(time = time_points, outcome = outcome_plot), aes_string(x = 'time', y = 'outcome')) +
      geom_line(data = data.frame(time = time_points, pred_outcome = pred_quantiles[, 2]), aes_string(x = 'time', y = 'pred_outcome'), linetype = 'dashed', color = '#e41a1c') + 
      labs(x = 'Time', y = 'Number of Cases') + 
      ggtitle(title) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(pred_plot)
  } else if (!is.null(sensitivity_pred_quantiles)) {
    sensitivity_df <- data.frame('Outcome' = outcome_plot, 'Estimate' = pred_quantiles[, 2], 'Sensitivity 1' = sensitivity_pred_quantiles[[1]][, 2], 'Sensitivity 2' = sensitivity_pred_quantiles[[2]][, 2], 'Sensitivity 3' = sensitivity_pred_quantiles[[3]][, 2], check.names = FALSE)
    sensitivity_bound <- data.frame('Sensitivity 1' = c(sensitivity_pred_quantiles[[1]][which(time_points %in% post_dates), 3], rev(sensitivity_pred_quantiles[[1]][which(time_points %in% post_dates), 1])), 'Sensitivity 2' = c(sensitivity_pred_quantiles[[2]][which(time_points %in% post_dates), 3], rev(sensitivity_pred_quantiles[[2]][which(time_points %in% post_dates), 1])), 'Sensitivity 3' = c(sensitivity_pred_quantiles[[3]][which(time_points %in% post_dates), 3], rev(sensitivity_pred_quantiles[[3]][which(time_points %in% post_dates), 1])), check.names = FALSE)
    
    pred_plot <- ggplot() + 
      geom_polygon(data = melt(sensitivity_bound, id.vars = NULL), aes_string(x = rep(post_dates, ncol(sensitivity_bound)), y = 'value', fill = 'variable'), alpha = 0.3) +
      geom_line(data = melt(sensitivity_df, id.vars = NULL), aes_string(x = rep(time_points, ncol(sensitivity_df)), y = 'value', color = 'variable')) +
      scale_colour_manual(values = c('black', '#e41a1c', '#4daf4a', '#4daf4a', '#984ea3')) +
      scale_fill_hue(guide = 'none') +
      labs(x = 'Time', y = 'Number of Cases') + 
      ggtitle(sensitivity_title) + 
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(pred_plot)
  }
}

#Sensitivity analysis by dropping the top weighted covariates. 
weightSensitivityAnalysis <- function(group, covars, ds, impact, time_points, intervention_date, n_seasons, outcome,n_iter = 10000, eval_period = NULL, post_period = NULL) {
  par(mar = c(5, 4, 1, 2) + 0.1)
  covar_df <- as.matrix(covars[[group]])
  #colnames(covar_df)<-substring(colnames(covar_df), 2)
  
  incl_prob <- sort(impact[[group]]$inclusion_probs[-c(1:n_seasons)])
  max_var <- names(incl_prob)[length(incl_prob)]
  max_prob <- round(incl_prob[length(incl_prob)],2)
  sensitivity_analysis <- vector('list', 3)
  
  for (i in 1:3) {
    covar_df <- covar_df[, colnames(covar_df) != max_var]
  
    #Combine covars, outcome, date
    y <- outcome[, group]
    y.pre<-outcome[time_points < as.Date(intervention_date), group]
    covar_df.pre<-covar_df[time_points < as.Date(intervention_date),]
    post_period_response <- outcome[, group]
    post_period_response <- as.vector(post_period_response[time_points >= as.Date(intervention_date)])
    
    regression_prior_df <- 50
    exp_r2 <- 0.8
    n_pred<-3
    denom <- ncol(covar_df.pre)-(n_seasons-1)
    if(denom>0){
      prior.inclusion.probabilities = c( rep(1,n_seasons),  rep(n_pred/denom,denom) ) #force seasonality and intercept into model, repeat '1' 12 times, repeat inclusion prob by N of cnon-monthly covars
    } else {	
      prior.inclusion.probabilities = rep(1,12)   
    }
    prior.inclusion.probabilities[prior.inclusion.probabilities>1] <- 1
    prior2=SpikeSlabPrior(cbind(1,covar_df.pre),y=y.pre, prior.inclusion.probabilities = prior.inclusion.probabilities,prior.df = regression_prior_df, expected.r2 = exp_r2 )
    bsts_model <- poisson.spike(y.pre ~ covar_df.pre,  niter = n_iter, prior=prior2 , ping = 0, seed = 1 )
   
    predict.bsts<-predict(bsts_model, newdata=cbind(1,covar_df), burn=n_iter*0.1, mean.only=FALSE)
    coef.bsts<-SummarizeSpikeSlabCoefficients(bsts_model$beta, burn=n_iter*0.1, order=FALSE)
    inclusion_probs<-coef.bsts[,5]
    names(inclusion_probs)<-substring(names(inclusion_probs),9)
    
    impact_sens <- list(predict.bsts,inclusion_probs, post.period.response = post_period_response, observed.y=outcome[, group])
    names(impact_sens)<-c('predict.bsts','inclusion_probs','post_period_response', 'observed.y' )
    sensitivity_analysis[[i]] <- list(removed_var = max_var, removed_prob = max_prob)
   # if (!is.null(mean) && !is.null(sd) && !is.null(eval_period) && !is.null(post_period)) {
      quantiles <- rrPredQuantiles(impact = impact_sens,  eval_period = eval_period, post_period = post_period)
      sensitivity_analysis[[i]]$rr <- round(quantiles$rr,2)
      sensitivity_analysis[[i]]$pred <- quantiles$pred
  #  }
    
    incl_prob.sens <- sort(impact_sens$inclusion_probs[-c(1:n_seasons)])
    max_var <- names(incl_prob.sens)[length(incl_prob.sens)]
    max_prob <- round(incl_prob.sens[length(incl_prob.sens)],2)
  }
  return(sensitivity_analysis)
}

predSensitivityAnalysis <- function(group, ds, zoo_data, denom_name, outcome_mean, outcome_sd, intervention_date, eval_period, post_period, time_points, n_seasons , n_pred) {
  impact <- doCausalImpact(zoo_data[[group]], intervention_date, time_points, n_seasons, n_pred = n_pred)
  quantiles <- lapply(group, FUN = function(group) {rrPredQuantiles(impact = impact, denom_data = ds[[group]][, denom_name],  eval_period = eval_period, post_period = post_period)})
  rr_mean <- t(sapply(quantiles, getRR))
  return(rr_mean)
}

sensitivityTable <- function(group, sensitivity_analysis, original_rr = NULL) {
  top_controls <- lapply(1:length(sensitivity_analysis[[group]]), FUN = function(i) {
    top_control <- c(sensitivity_analysis[[group]][[i]]$removed_var, sensitivity_analysis[[group]][[i]]$removed_prob, sensitivity_analysis[[group]][[i]]$rr)
    names(top_control) <- c(paste('Top Control', i), paste('Inclusion Probability of Control', i), paste(names(sensitivity_analysis[[group]][[i]]$rr), i))
    return(top_control)
  })
  sensitivity_table <- c(original_rr[group, ], c(top_controls, recursive = TRUE))
  return(sensitivity_table)
}
