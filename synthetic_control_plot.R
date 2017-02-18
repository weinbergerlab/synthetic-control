#This is the plot file. Run it after the analysis file to visualize results if not using the rmarkdown report.

#Plot covars, models, and impacts
par(par_defaults)
invisible(lapply(groups, FUN = function(group) {
	#View scaled covariates
	matplot(covars_full[[group]], type = 'l', xlab = 'Log Covars')
	title(main = group)
	
	#Plot predictions
	plotPred(pred_quantiles_full[, , group], time_points, post_period, pred_quantiles_full[, , group], outcome_plot[, group], code_change, title = paste(group, 'Synthetic controls estimate'))
	plotPred(pred_quantiles_time[, , group], time_points, post_period, pred_quantiles_full[, , group], outcome_plot[, group], code_change, title = paste(group, 'Interupted time series estimate'))
	
	#Rolling rate ratio
	matplot(rr_roll_full[, , group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l')
	title(main = paste(group, 'Synthetic controls: rolling rate ratio'))
	abline(h = 1, lty = 2)
	
	matplot(rr_roll_time[, , group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l')
	title(main = paste(group, 'Interupted time series: rolling rate ratio'))
	abline(h = 1, lty = 2)
	
	#View cumulative sums
	matplot(cumsum_prevented[, , group], type = 'l', col = 'black')
	abline(h = 0, lty = 2)
	title(main = group, sub = 'Cumulative cases prevented')
	return()
}))

log_rr_mean_full <- log(rr_mean_full)
log_rr_mean_full[is.na(log_rr_mean_full)] <- 1e6

min_max <- na.omit(log_rr_mean_full[log_rr_mean_full[1] < 100, ])
plotModel(log_rr_mean_full, groups, min_max)