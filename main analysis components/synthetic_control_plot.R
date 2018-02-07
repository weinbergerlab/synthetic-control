#This is the plotting file. Run it after the analysis file to visualize results if not using the rmarkdown report.

plot_list <- setNames(vector(mode = 'list', length = length(groups)), groups)
for (group in groups) {
	plots <- list()
	
	#View scaled covariates
	covars.sub<-covars_full[[group]][,-c(1:12)]
	covar_plot <- ggplot(melt(covars.sub, id.vars = NULL), mapping = aes_string(x = rep(time_points, ncol(covars.sub)), y = 'value', group = 'variable', alpha = rep(inclusion_prob_full[[group]][colnames(covars.sub)], each = nrow(covars.sub)   ))) + 
		geom_line() + 
		labs(x = 'Time', y = 'Scaled Covariates') + 
		ggtitle(paste(group, 'Scaled Covariates Weighted by Inclusion Probability')) +
		theme_bw() +
		theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	#Plot predictions
	min_max <- c(min(c(pred_quantiles_full[, , group], outcome_plot[, group])), max(c(pred_quantiles_full[, , group], outcome_plot[, group])))
	pred_full_plot <-        plotPred(pred_quantiles_full[, , group], time_points, post_period, min_max, outcome_plot[, group], title = paste(group, 'Synthetic controls estimate'))
	pred_time_plot <-        plotPred(pred_quantiles_time[, , group], time_points, post_period, min_max, outcome_plot[, group], title = paste(group, 'Interupted time series estimate'))
	pred_sensitivity_plot <- plotPred(pred_quantiles_full[, , group], time_points, post_period, min_max, outcome_plot[, group], sensitivity_pred_quantiles = sensitivity_pred_quantiles[[group]], sensitivity_title = paste(group, 'Sensitivity Plots'), plot_sensitivity = TRUE)
	
	#matplot(pred_quantiles_full[, , 10], ylim=c(0,22000), type='l')	 ##Check
	#points(prelog_data[[10]]$J12_18)
	
	#Plot rolling rate ratio
	min_max <- c(min(rr_roll_full[, , group], rr_roll_time[, , group]), max(rr_roll_full[, , group], rr_roll_time[, , group]))
	rr_roll_full_plot <- ggplot(melt(as.data.frame(rr_roll_full[, , group]), id.vars = NULL), mapping = aes_string(x = rep(time_points[(length(time_points) - nrow(rr_roll_full[, , group]) + 1):length(time_points)], ncol(rr_roll_full[, , group])), y = 'value', linetype = 'variable')) + 
		geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Rolling Rate Ratio') + 
		ggtitle(paste(group, 'Synthetic Control Rolling Rate Ratio')) +
		coord_cartesian(ylim = min_max) +
		theme_bw() +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	rr_roll_time_plot <- ggplot(melt(as.data.frame(rr_roll_time[, , group]), id.vars = NULL), mapping = aes_string(x = rep(time_points[(length(time_points) - nrow(rr_roll_time[, , group]) + 1):length(time_points)], ncol(rr_roll_time[, , group])), y = 'value', linetype = 'variable')) + 
		geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Rolling Rate Ratio') + 
		ggtitle(paste(group, 'TT Rolling Rate Ratio')) +
		coord_cartesian(ylim = min_max) +
		theme_bw() +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	#Plot cumulative sums
	cumsum_prevented_plot <- ggplot(melt(as.data.frame(cumsum_prevented[, , group]), id.vars = NULL), mapping = aes_string(x = rep(time_points, ncol(cumsum_prevented[, , group])), y = 'value', linetype = 'variable')) + 
		geom_line() + geom_hline(yintercept = 1, linetype = 4) +
		labs(x = 'Time', y = 'Cumulative Sum Prevented') + 
		ggtitle(paste(group, 'Cumulative Number of Cases Prevented')) + 
		theme_bw() +
		theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(colour = NA, fill = 'transparent'), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	plot_list[[group]] <- list(covar_plot = covar_plot, 
															pred_full_plot = pred_full_plot, 
															pred_time_plot = pred_time_plot, 
															pred_sensitivity_plot = pred_sensitivity_plot, 
															rr_roll_full_plot = rr_roll_full_plot, 
															rr_roll_time_plot = rr_roll_time_plot, 
															cumsum_prevented_plot = cumsum_prevented_plot)
}

for (group in groups) {
	cat('###', group, '\n', sep = '')
	for (group_plot in plot_list[[group]]) {
		print(group_plot)
	}
	cat('\n\n')
}
