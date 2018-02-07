source('synthetic_control_functions_Shiny.R', local=FALSE)

packages <- c('parallel', 'shiny', 'splines',  'lubridate', 'RcppRoll', 'BoomSpikeSlab', 'ggplot2', 'reshape','dummies')
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

# library(shiny, quietly = TRUE)
# library(splines, quietly = TRUE)
# library(lubridate, quietly = TRUE)
# library(BoomSpikeSlab, quietly = TRUE)
# library(parallel, quietly = TRUE)
# library(RcppRoll, quietly = TRUE)
# library(reshape, quietly = TRUE)
# library(ggplot2, quietly = TRUE)
# library(dummies, quietly=TRUE)

#Set max file size
options(shiny.maxRequestSize = 100 * 1024 ^ 2) #100MB


shinyServer(function(input, output, session) {
	
	output$Sage <- renderImage(list(
		src = 'images/Sage Analytica Logo.png',
		contentType = 'image/png', 
		alt = 'Sage Analytica Logo', 
		height = 100
	), deleteFile = FALSE)
	output$Yale <- renderImage(list(
		src = 'images/Yale Logo.png',
		contentType = 'image/png', 
		alt = 'Yale Logo', 
		height = 100
	), deleteFile = FALSE)
	
	file_name <- reactive({
		data_file <- input$in_file
		if (is.null(data_file)) {return(NULL)}
		return(data_file$datapath)
	})
	input_data <- reactive({
		data_file <- input$in_file
		if (is.null(data_file)) {return(NULL)}
		data <- read.csv(data_file$datapath, check.names = FALSE, stringsAsFactors = FALSE)
		return(data)
	})
	observe({
		ds <- input_data()
		updateSelectInput(session = session, inputId = 'group', choices = colnames(ds))
	})
	observeEvent(input$group, {
		ds <- input_data()
		updateSelectInput(session = session, inputId = 'group_selected', choices = unique(ds[, input$group]), selected = unique(ds[, input$group]))
	})
	observeEvent(input$group, {
		ds <- input_data()
		updateSelectInput(session = session, inputId = 'date', choices = colnames(ds)[-c(match(c(input$group), colnames(ds)))])
	})
	observeEvent(input$date, {
		ds <- input_data()
		all_dates <- unique(ds[, input$date])
		updateSelectInput(session = session, inputId = 'variable', choices = colnames(ds)[-c(match(c(input$group, input$date), colnames(ds)))])
		updateDateRangeInput(session = session, inputId = 'training_range', start = as.character(all_dates[1]), end = as.character(all_dates[length(all_dates) / 2]), min = all_dates[1], max = all_dates[length(all_dates)])
		updateDateRangeInput(session = session, inputId = 'eval_range', start = as.character(all_dates[length(all_dates) / 2 + 1]), end = as.character(all_dates[length(all_dates)]), min = all_dates[1], max = all_dates[length(all_dates)])
	})
	observeEvent(input$variable, {
		ds <- input_data()
		choices <- colnames(ds)[-c(match(c(input$group, input$date, input$variable), colnames(ds)))]
		if (input$covariate_checkbox) {
			selected <- NULL
		} else {
			selected <- choices
		}
		updateSelectInput(session = session, inputId = 'covariate', choices = choices, selected = selected)
	})
	observeEvent(input$covariate, {
		ds <- input_data()
		choices <- colnames(ds)[-c(match(c(input$group, input$date, input$variable), colnames(ds)))]
		updateSelectInput(session = session, inputId = 'noj_denom', choices = choices)
	})
	
	input_vars <- eventReactive(input$run_CausalImpact, {
		input_data <- input_data()
		reformatted_date <- unique(formatDate(input_data[, input$date]))
		
		if (input$covariate_checkbox) {
			exclude_covar <- input$covariate
		} else {
			exclude_covar <- colnames(input_data)[-c(match(c(input$group, input$date, input$variable, input$covariate), colnames(input_data)))]
		}
		exclude_group <- unique(input_data[, input$group])[!(unique(input_data[, input$group]) %in% input$group_selected)]
		code_change <- input$adjust_covariates_checkbox
		
		input_directory  <- ''
		output_directory <- dirname(file_name())
		file_name        <- file_name()
		
		group_name   <- input$group
		date_name    <- input$date
		outcome_name <- input$variable
		denom_name   <- input$noj_denom
		
		groups <- as.character(unique(unlist(input_data[, group_name], use.names = FALSE)))
		if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}
		
		if (as.numeric(input$training_range[1]) < as.numeric(reformatted_date[1])) {
			start_date <- as.Date(reformatted_date[1])
		} else {
			start_date <- as.Date(reformatted_date[findInterval(as.numeric(input$training_range[1]), as.numeric(as.Date(unique(reformatted_date))))])
		}
		intervention_date <- input$training_range[2]
		if (as.numeric(input$eval_range[2]) > as.numeric(as.Date(reformatted_date[length(reformatted_date)]))) {
			end_date <- as.Date(reformatted_date[length(reformatted_date)])
		} else {
			end_date <- as.Date(reformatted_date[Position(function(date) {date >= input$eval_range[2]}, as.numeric(as.Date(reformatted_date)))])
		}
		time_points <- as.Date(as.character(reformatted_date[match(start_date, as.Date(reformatted_date)):match(end_date, as.Date(reformatted_date))]))
		pre_period <- input$training_range
		post_period <- c(as.Date(min(time_points[time_points >= intervention_date])), end_date)
		eval_period <- input$eval_range
		
		n_seasons <- length(unique(month(time_points)))
		#Monthly dummies
		if(n_seasons==4){x<-quarter(time_points)}
		if(n_seasons==12){x<-month(time_points)}
		if(n_seasons==3){
		  x.m<-month(time_points)
		  x<-x.m
		  x[x.m %in% c(1,2,3,4)]<-1
		  x[x.m %in% c(5,6,7,8)]<-2
		  x[x.m %in% c(9,10,11,12)]<-3
		}
		season.dummies<-dummy(x)
		season.dummies<-season.dummies[,-n_seasons]
		
		return(list(
			groups = groups,
			
			n_seasons     = n_seasons, 
			season.dummies     = season.dummies, 
			exclude_covar = exclude_covar, 
			exclude_group = exclude_group,
			code_change   = code_change,
			
			input_directory  = input_directory,
			output_directory = output_directory,
			file_name        = file_name,
			
			group_name   = group_name, 
			date_name    = date_name,
			outcome_name = outcome_name,
			denom_name   = denom_name,
			
			start_date        = start_date,
			intervention_date = intervention_date,
			end_date          = end_date,
			pre_period        = pre_period,
			post_period       = post_period,
			time_points       = time_points,
			eval_period       = eval_period
		))
	})
	
	prefiltered_data <- eventReactive(input_vars(), {
		groups <- input_vars()$groups
		
		group_name <- input_vars()$group_name
		date_name <- input_vars()$date_name
		outcome_name <- input_vars()$outcome_name
		denom_name <- input_vars()$denom_name
		
		start_date <- input_vars()$start_date
		end_date <- input_vars()$end_date
		
		prelog_data <- input_data()
		prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
		prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
		
		#Log-transform all variables, adding 0.5 to counts of 0.
		ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name)), groups)
		ds <- lapply(ds, function(ds) {
			if (!(denom_name %in% colnames(ds))) {
				ds[denom_name] <- 0
			}
			return(ds)
		})
		return(ds)
	})
	sparse_groups <- reactive({
		ds <- prefiltered_data()
		
		group_name   <- input_vars()$group_name
		date_name    <- input_vars()$date_name
		outcome_name <- input_vars()$outcome_name
		denom_name   <- input_vars()$denom_name
		
		exclude_covar <- input_vars()$exclude_covar
		
		sparse_groups <- sapply(ds, function(ds) {
			return(ncol(ds[!(colnames(ds) %in% c(date_name, group_name, denom_name, outcome_name, exclude_covar))]) == 0)
		})
		return(sparse_groups)
	})
	ds <- reactive(prefiltered_data()[!sparse_groups()])
	groups <- reactive(input_vars()$groups[!sparse_groups()])
	
	covars_full <- reactive({
		ds <- ds()
		groups <- groups()
		
		code_change <- input_vars()$code_change
		intervention_date <- input_vars()$intervention_date
		time_points <- input_vars()$time_points
		season.dummies<-input_vars()$season.dummies
		exclude_covar <- input_vars()$exclude_covar
		
		covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change, intervention_date = intervention_date, season.dummies=season.dummies,time_points = time_points), groups)
		covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
		return(covars_full)
	})
	covars_time <- reactive({
		covars_full <- covars_full()
		groups <- groups()
		season.dummies<-input_vars()$season.dummies
		setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(cbind(season.dummies,time_index = 1:nrow(covars))))}), groups)
	})
	
	outcome      <- reactive(sapply(ds(), FUN = function(data) {scale(data[, input_vars()$outcome_name])}))
	outcome_mean <- reactive(sapply(ds(), FUN = function(data) {mean( data[, input_vars()$outcome_name])}))
	outcome_sd   <- reactive(sapply(ds(), FUN = function(data) {sd(   data[, input_vars()$outcome_name])}))
	outcome_plot <- reactive(exp(t(t(outcome()) * outcome_sd() + outcome_mean())))
	outcome_offset      <- reactive(scale(   sapply(ds(), FUN = function(data) {   data[, input_vars()$outcome_name] - data[, input_vars()$denom_name]})))
	outcome_offset_mean <- reactive(colMeans(sapply(ds(), FUN = function(data) {   data[, input_vars()$outcome_name] - data[, input_vars()$denom_name]})))
	outcome_offset_sd   <- reactive(         sapply(ds(), FUN = function(data) {sd(data[, input_vars()$outcome_name] - data[, input_vars()$denom_name])}))
	
	data_full <- reactive(setNames(lapply(groups(), makeTimeSeries, outcome = outcome(),        covars = covars_full()), groups()))
	data_time <- reactive(setNames(lapply(groups(), makeTimeSeries, outcome = outcome_offset(), covars = covars_time()), groups()))
	
	impact_full <- reactive({
		data_full <- data_full()
		groups <- groups()
		
		intervention_date <- input_vars()$intervention_date
		time_points <- input_vars()$time_points
		n_seasons <- input_vars()$n_seasons
		season.dummies <- input_vars()$season.dummies
		
		# if (Sys.info()['sysname'] == 'Windows') {
		# 	cl <- makeCluster(n_cores)
		# 	clusterEvalQ(cl, {library(BoomSpikeSlab, quietly = TRUE); library(lubridate, quietly = TRUE)})
		# 	clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons','season.dummies'), environment())
		# } else {
		# 	cl <- makeForkCluster(n_cores)
		# }
		
		impact_full <- setNames(lapply( data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
		
		#stopCluster(cl)
		gc()
		return(impact_full)
	})
	impact_time <- reactive({
		data_time <- data_time()
		groups <- groups()
		
		intervention_date <- input_vars()$intervention_date
		time_points <- input_vars()$time_points
		n_seasons <- input_vars()$n_seasons
		
		# if (Sys.info()['sysname'] == 'Windows') {
		# 	cl <- makeCluster(n_cores)
		# 	clusterEvalQ(cl, {library(BoomSpikeSlab, quietly = TRUE); library(lubridate, quietly = TRUE)})
		# 	clusterExport(cl, c('doCausalImpact',  'intervention_date', 'time_points', 'n_seasons'), environment())
		# } else {
		# 	cl <- makeForkCluster(n_cores)
		# }
		
		impact_time <- setNames(lapply( data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)
		
	#	stopCluster(cl)
		gc()
		return(impact_time)
	})
	
	inclusion_prob_full <- reactive(setNames(lapply(impact_full(), inclusionProb), groups()))
	inclusion_prob_time <- reactive(setNames(lapply(impact_time(), inclusionProb), groups()))
	
	#All model results combined
	quantiles_full <- reactive({
		impact_full <- impact_full()
		groups <- groups()
		ds <- ds()
		outcome_mean <- outcome_mean()
		outcome_sd <- outcome_sd()
		
		eval_period <- input_vars()$eval_period
		post_period <- input_vars()$post_period
		denom_name <- input_vars()$denom_name
		time_points <- input_vars()$time_points
		n_seasons <- input_vars()$n_seasons
		
		quantiles_full <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_full[[group]], time_points=time_points, denom_data = ds[[group]][, denom_name],n_seasons=n_seasons, mean = outcome_mean[group], sd = outcome_sd[group], eval_period = eval_period, post_period = post_period)}), groups)
		return(quantiles_full)
	})
	quantiles_time <- reactive({
		impact_time <- impact_time()
		groups <- groups()
		ds <- ds()
		outcome_offset_mean <- outcome_offset_mean()
		outcome_offset_sd <- outcome_offset_sd()
		
		eval_period <- input_vars()$eval_period
		post_period <- input_vars()$post_period
		denom_name <- input_vars()$denom_name
		time_points <- input_vars()$time_points
		n_seasons <- input_vars()$n_seasons
		
		quantiles_time <- setNames(lapply(groups, FUN = function(group) {rrPredQuantiles(impact = impact_time[[group]], denom_data = ds[[group]][, denom_name], mean = outcome_offset_mean[group], time_points=time_points,n_seasons=n_seasons,sd = outcome_offset_sd[group], eval_period = eval_period, post_period = post_period, trend = TRUE)}), groups)
		return(quantiles_time)
	})
	
	#Model predicitons
	pred_quantiles_full <- reactive(sapply(quantiles_full(), getPred, simplify = 'array'))
	pred_quantiles_time <- reactive(sapply(quantiles_time(), getPred, simplify = 'array'))
	
	#Rolling rate ratios
	rr_roll_full <- reactive(sapply(quantiles_full(), FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array'))
	rr_roll_time <- reactive(sapply(quantiles_time(), FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array'))
	
	#Rate ratios for evaluation period.
	rr_mean_full <- reactive(t(sapply(quantiles_full(), getRR)))
	rr_mean_full_intervals <- reactive(data.frame('Estimate (95% CI)' = makeInterval(rr_mean_full()[, 2], rr_mean_full()[, 3], rr_mean_full()[, 1]), check.names = FALSE, row.names = groups()))
	rr_mean_time <- reactive({
		rr_mean_time <- t(sapply(quantiles_time(), getRR))
		colnames(rr_mean_time) <- paste('ITS', colnames(rr_mean_time))
		return(rr_mean_time)
	})
	rr_mean_time_intervals <- reactive(data.frame('ITS Estimate (95% CI)' = makeInterval(rr_mean_time()[, 2], rr_mean_time()[, 3], rr_mean_time()[, 1]), check.names = FALSE, row.names = groups()))
	
	cumsum_prevented <- reactive({
		groups <- groups()
		outcome_plot <- outcome_plot()
		quantiles_full <- quantiles_full()
		
		time_points <- input_vars()$time_points
		post_period <- input_vars()$post_period
		
		sapply(groups, FUN = function(group, quantiles) {
			is_post_period <- which(time_points >= post_period[1])
			is_pre_period <- which(time_points < post_period[1])
			
			#Cumulative sum of prevented cases
			cases_prevented <- t(quantiles[[group]]$pred_samples) - outcome_plot[, group]
			cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period, ], 2, cumsum)
			cumsum_cases_prevented_pre <- matrix(0, nrow = nrow(cases_prevented[is_pre_period, ]), ncol = ncol(cases_prevented[is_pre_period, ]))
			cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
			cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
		}, quantiles = quantiles_full, simplify = 'array')
	})
	
	
	# cumsum_prevented <- sapply(groups, FUN = function(group, quantiles) {
	#   is_post_period <- which(time_points >= post_period[1])
	#   is_pre_period <- which(time_points < post_period[1])
	# 
	#   #Cumulative sum of prevented cases
	#   cases_prevented <- t(quantiles[[group]]$pred_samples) - outcome_plot[, group]
	#   cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period, ], 2, cumsum)
	#   cumsum_cases_prevented_pre <- matrix(0, nrow = nrow(cases_prevented[is_pre_period, ]), ncol = ncol(cases_prevented[is_pre_period, ]))
	#   cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
	#   cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
	# }, quantiles = quantiles_full, simplify = 'array')

	
	bad_sensitivity_groups <- reactive(sapply(covars_full(), function (covar) {ncol(covar) <= 3}))
	sensitivity_covars_full <- reactive(covars_full()[!bad_sensitivity_groups()])
	sensitivity_ds <- reactive(ds()[!bad_sensitivity_groups()])
	sensitivity_impact_full <- reactive(impact_full()[!bad_sensitivity_groups()])
	sensitivity_groups <- reactive(groups()[!bad_sensitivity_groups()])
	#groups <- reactive(input_vars()$groups[!sparse_groups()])
	
	sensitivity_analysis_full <- reactive({
	  
		outcome_mean <- outcome_mean()
		outcome_sd <- outcome_sd()
		
		intervention_date <- input_vars()$intervention_date
		n_seasons <- input_vars()$n_seasons
		eval_period <- input_vars()$eval_period
		post_period <- input_vars()$post_period
		time_points <- input_vars()$time_points
		
		sensitivity_analysis_full <- setNames(lapply(sensitivity_groups(), weightSensitivityAnalysis, covars = sensitivity_covars_full(), outcome_sd=outcome_sd, outcome_mean=outcome_mean, n_seasons = n_seasons, time_points=time_points, ds = sensitivity_ds(), impact = sensitivity_impact_full(), intervention_date = intervention_date,  outcome = outcome(),  eval_period = eval_period, post_period = post_period), sensitivity_groups())
		#	impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
		
		#stopCluster(cl)
		gc()
		return(sensitivity_analysis_full)
		
	})
	
	sensitivity_pred_quantiles <- reactive({
	  time_points <- input_vars()$time_points
		sensitivity_pred_quantiles  <- lapply(sensitivity_analysis_full(), FUN = function(sensitivity_analysis) {
			pred_list <- vector(mode = 'list', length = length(sensitivity_analysis))
			for (sensitivity_index in 1:length(sensitivity_analysis)) {
				pred_list[[sensitivity_index]] <- getPred(sensitivity_analysis[[sensitivity_index]])
			}
			return(pred_list)
		})
	})
	sensitivity_table <- reactive(t(sapply(sensitivity_groups(), sensitivityTable, sensitivity_analysis = sensitivity_analysis_full(), original_rr = rr_mean_full())))
	sensitivity_table_intervals <- reactive({
	  time_points <- input_vars()$time_points
		sensitivity_table <- sensitivity_table()
		data.frame('Estimate (95% CI)' = makeInterval(sensitivity_table[, 2],  sensitivity_table[, 3],  sensitivity_table[, 1]),
							 'Top Control 1' = sensitivity_table[, 'Top Control 1'],
							 'Inclusion Probability of Control 1' = sensitivity_table[, 'Inclusion Probability of Control 1'],
							 'Control 1 Estimate (95% CI)' = makeInterval(sensitivity_table[, 7],  sensitivity_table[, 8],  sensitivity_table[, 6]),
							 'Top Control 2' = sensitivity_table[, 'Top Control 2'],
							 'Inclusion Probability of Control 2' = sensitivity_table[, 'Inclusion Probability of Control 2'],
							 'Control 2 Estimate (95% CI)' = makeInterval(sensitivity_table[, 12],  sensitivity_table[, 13],  sensitivity_table[, 11]),
							 'Top Control 3' = sensitivity_table[, 'Top Control 3'],
							 'Inclusion Probability of Control 3' = sensitivity_table[, 'Inclusion Probability of Control 3'],
							 'Control 3 Estimate (95% CI)' = makeInterval(sensitivity_table[, 17],  sensitivity_table[, 18],  sensitivity_table[, 16]), check.names = FALSE)
	})
	rr_table <- reactive(cbind(rr_mean_time()[!bad_sensitivity_groups(), ], sensitivity_table()))
	rr_table_intervals <- reactive({cbind('ITS Estimate (95% CI)' = rr_mean_time_intervals()[!bad_sensitivity_groups, ], sensitivity_table_intervals())})
	
	incl_probs <- reactive({
		groups <- groups()
		impact_full <- impact_full()
		
		n_seasons <- input_vars()$n_seasons
		incl_probs <- NULL
		for (group in groups) {
			incl_prob <- -1*sort(-impact_full[[group]]$inclusion_probs[-c(1:n_seasons)])[1:3]
			incl_prob <- data.frame('Group' = group, 'Greatest Inclusion Variable' = names(incl_prob)[1], 'Greatest Inclusion Probability' = incl_prob[1], 'Second Greatest Inclusion Variable' = names(incl_prob)[2], 'Second Greatest Inclusion Probability' = incl_prob[2], 'Third Greatest Inclusion Variable' = names(incl_prob)[3], 'Third Greatest Inclusion Probability' = incl_prob[3], check.names = FALSE)
			incl_probs <- rbind(incl_probs, incl_prob)
		}
		rownames(incl_probs) <- NULL
		return(incl_probs)
	})
	
	group_tabs <- reactive({
		progress <- Progress$new()
		on.exit(progress$close())
		progress$set(message = 'Setup', value = 0)
		groups <- groups()
		sparse_groups <- sparse_groups()
		time_points <- input_vars()$time_points
		outcome_plot <- outcome_plot()
		post_period  <- input_vars()$post_period
		progress$inc(amount = 0.25, message = 'Synthetic Control Analysis')
		rr_mean_full_intervals <- rr_mean_full_intervals()
		rr_mean_time_intervals <- rr_mean_time_intervals()
		rr_roll_full <- rr_roll_full()
		rr_roll_time <- rr_roll_time()
		inclusion_prob_full <- inclusion_prob_full()
		incl_probs <- incl_probs()
		progress$inc(amount = 0.25, message = 'Sensitivity Analysis')
		sensitivity_table_intervals <- sensitivity_table_intervals()
		#sensitivity_pred_2_intervals <- sensitivity_pred_2_intervals()
		#sensitivity_pred_10_intervals <- sensitivity_pred_10_intervals()
		sensitivity_pred_quantiles <- sensitivity_pred_quantiles()
		covars_full <- covars_full()
		pred_quantiles_full <- pred_quantiles_full()
		pred_quantiles_time <- pred_quantiles_time()
		cumsum_prevented <- cumsum_prevented()
		progress$inc(amount = 0.25, message = 'Plotting')
		source('synthetic_control_plot.R', local = TRUE)
		progress$inc(amount = 0.25, message = 'Finishing')
		tabs <- lapply(groups, FUN = function(group) {
			tab <- tabPanel(
				title = group, 
				renderPlot({print(plot_list[[group]]$covar_plot)}),
				renderPlot({print(plot_list[[group]]$pred_full_plot)}),
				renderPlot({print(plot_list[[group]]$pred_time_plot)}),
				renderPlot({print(plot_list[[group]]$pred_sensitivity_plot)}),
				renderPlot({print(plot_list[[group]]$rr_roll_full_plot)}),
				renderPlot({print(plot_list[[group]]$rr_roll_time_plot)}),
				renderPlot({print(plot_list[[group]]$cumsum_prevented_plot)})
			)
			return(tab)
		})
		
		summary_tab <- tabPanel(
			title = 'Summary',
			if (!is.null(names(sparse_groups[sparse_groups])) && length(names(sparse_groups[sparse_groups])) != 0) {
				renderTable(data.frame('Sparse Groups' = names(sparse_groups[sparse_groups]), check.names = FALSE), caption = 'Sparse Groups', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL))
			},
			renderTable(rr_mean_full_intervals,                 caption = 'Synthetic Control Quantiles',                                    caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL), rownames = TRUE),
			renderTable(rr_mean_time_intervals,                 caption = 'Interrupted Time Series Quantiles',                              caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL), rownames = TRUE),
			renderTable(incl_probs,                             caption = 'Inclusion Probabilities',                                        caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL)),
			renderTable(sensitivity_table_intervals,            caption = 'Weight Sensitivity Analysis',                                    caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL), rownames = TRUE)#,
			#renderTable(sensitivity_pred_2_intervals,  caption = 'Predictor Sensitivity Analysis where Number of Predictors = 2',  caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL), rownames = TRUE),
			#renderTable(sensitivity_pred_10_intervals, caption = 'Predictor Sensitivity Analysis where Number of Predictors = 10', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL), rownames = TRUE)
		)
		tabs <- append(list(summary_tab), tabs)
		return(tabs)
	})
	output$tab_view <- renderUI({
		groups <- groups()
		if (length(groups) == 0) {return(NULL)}
		tabs <- group_tabs()
		do.call(tabsetPanel, tabs)
	})
})
