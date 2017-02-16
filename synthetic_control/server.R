library(shiny, quietly = TRUE)
library(splines, quietly = TRUE)
library(CausalImpact, quietly = TRUE)
library(parallel, quietly = TRUE)
library(RcppRoll, quietly = TRUE)
source(file.path('R_scripts','processing_functions.R'))

#Set max file size
options(shiny.maxRequestSize = 100 * 1024 ^ 2) #100MB
n_cores <- detectCores()
progress <- NULL
value <- 0
counter <- 0

shinyServer(function(input, output, session) {
	
	output$download_preprocesser <- downloadHandler(
		filename = function() {file_name <- 'sample_preprocessing.R'},
		content = function(file) {file.copy('./downloads/sample_preprocessing.R', file)},
		contentType = 'application/R'
	)
	output$download_sample <- downloadHandler(
		filename <- function() {file_name <- 'example_data.csv'},
		content <- function(file) {file.copy('./downloads/example_data.csv', file)},
		contentType = 'text/csv'
	)
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
	
	progress <<- Progress$new()
	customIncProgress <- function(amount = 0, reset = FALSE) {
		value <<- amount + value
		progress$set(value = value)
		if (reset) {
			value <<- 0
		}
	}
	
	ds <- reactive({
		data_file <- input$in_file
		if (is.null(data_file)) {return(NULL)}
		data <- read.csv(data_file$datapath, check.names = FALSE)
		return(data)
	})
	observe({
		ds <- ds()
		updateSelectInput(session = session, inputId = 'age_group', choices = colnames(ds))
	})
	observeEvent(input$age_group, {
		ds <- ds()
		updateSelectInput(session = session, inputId = 'age_group_selected', choices = unique(ds[, input$age_group]), selected = unique(ds[, input$age_group]))
	})
	observeEvent(input$age_group, {
		ds <- ds()
		updateSelectInput(session = session, inputId = 'date', choices = colnames(ds)[-c(match(c(input$age_group), colnames(ds)))])
	})
	observeEvent(input$date, {
		ds <- ds()
		all_dates <- unique(ds[, input$date])
		updateSelectInput(session = session, inputId = 'variable', choices = colnames(ds)[-c(match(c(input$age_group, input$date), colnames(ds)))])
		updateDateRangeInput(session = session, inputId = 'training_range', start = as.character(all_dates[1]), end = as.character(all_dates[length(all_dates)/2]), min = all_dates[1], max = all_dates[length(all_dates)])
		updateDateRangeInput(session = session, inputId = 'eval_range', start = as.character(all_dates[length(all_dates)/2+1]), end = as.character(all_dates[length(all_dates)]), min = all_dates[1], max = all_dates[length(all_dates)])
	})
	observeEvent(input$variable, {
		ds <- ds()
		choices <- colnames(ds)[-c(match(c(input$age_group, input$date, input$variable), colnames(ds)))]
		if (input$covariate_checkbox) {
			selected <- NULL
		} else {
			selected <- choices
		}
		updateSelectInput(session = session, inputId = 'covariate', choices = choices, selected = selected)
	})
	observeEvent(input$covariate, {
		ds <- ds()
		choices <- colnames(ds)[-c(match(c(input$age_group, input$date, input$variable), colnames(ds)))]
		updateSelectInput(session = session, inputId = 'noj_denom', choices = choices)
	})
	
	observeEvent(input$run_CausalImpact, {
		if (is.null(progress)) {
			progress <<- Progress$new()
		}
		progress$set(value = 0, message = 'Starting...')
		customIncProgress(reset = TRUE)
	})
	locked_input <- eventReactive(input$run_CausalImpact, {return(input)})
	reformatted_date <- eventReactive(input$run_CausalImpact, {return(formatDate(ds()[, locked_input()$date]))})
	
	start_date <- eventReactive(input$run_CausalImpact, {
		reformatted_date <- reformatted_date()
		if (as.numeric(locked_input()$training_range[1]) < as.numeric(reformatted_date[1])) {
			start_date <- as.Date(reformatted_date[1])
		} else {
			start_date <- as.Date(reformatted_date[findInterval(as.numeric(locked_input()$training_range[1]), as.numeric(as.Date(unique(reformatted_date))))])
		}
	})
	end_date <- eventReactive(input$run_CausalImpact, {
		reformatted_date <- reformatted_date()
		if (as.numeric(locked_input()$eval_range[2]) > as.numeric(as.Date(reformatted_date[length(reformatted_date)]))) {
			end_date <- as.Date(reformatted_date[length(reformatted_date)])
		} else {
			end_date <- as.Date(reformatted_date[Position(function(date) {date >= locked_input()$eval_range[2]}, as.numeric(as.Date(reformatted_date)))])
		}
	})
	post_period <- eventReactive(time_points(), {
		time_points <- time_points()
		intervention_date <- locked_input()$training_range[2]
		c(time_points[which(abs(time_points - intervention_date) == min(abs(time_points[time_points >= intervention_date] - intervention_date)))], end_date())
	})
	time_points <- eventReactive(input$run_CausalImpact, {
		ds <- ds()
		reformatted_date <- reformatted_date()
		start_date <- start_date()
		end_date <- end_date()
		time_points <- as.Date(as.character(reformatted_date[match(start_date, as.Date(reformatted_date)):match(end_date, as.Date(reformatted_date))]))
	})
	age_groups <- eventReactive(input$run_CausalImpact, {paste(locked_input()$age_group, locked_input()$age_group_selected)})
	
	ds_log <- eventReactive(age_groups(), {
		ds <- ds()
		ds[, locked_input()$date] <- reformatted_date()
		age_groups <- locked_input()$age_group_selected
		ds <- setNames(lapply(age_groups, FUN = logTransform, factor_name = locked_input()$age_group, date_name = locked_input()$date, all_cause_name = locked_input()$noj_denom, all_cause_pneu_name = locked_input()$variable, start_date = start_date(), prelog_data = ds), age_groups())
		return(ds)
	})
	covar_list <- eventReactive(ds_log(), {
		ds <- ds_log()
		age_groups <- age_groups()
		covar_list <- sapply(age_groups, FUN = function(age_group) {
			setNames(ifelse(locked_input()$covariate %in% colnames(ds[[age_group]]), TRUE, FALSE), locked_input()$covariate)
		})
		if (length(locked_input()$covariate) <= 1) {
			covar_list <- t(matrix(covar_list, dimnames = list(age_groups, locked_input()$covariate)))
		}
		return(covar_list)
	})
	covars <- eventReactive(covar_list(), {
		age_groups <- age_groups()
		ds <- ds_log()
		covar_list <- covar_list()
		if (locked_input()$covariate_checkbox) {
			if (is.null(locked_input()$covariate)) {
				covariates <-  colnames(ds)[-c(match(c(locked_input()$age_group, locked_input()$date, locked_input()$variable), colnames(ds)))]
			} else {
				covariates <- colnames(ds)[-c(match(c(locked_input()$age_group, locked_input()$date, locked_input()$variable, locked_input()$covariate), colnames(ds)))]
			}
		} else {
			covariates <- locked_input()$covariate
		}
		if (length(covariates) == 0) {
			return(NULL)
		} else {
			covars <- setNames(lapply(age_groups, FUN = function(age_group) {
				good_covars <- rownames(covar_list[covar_list[, age_group] == TRUE, age_group, drop = FALSE])
				if (locked_input()$adjust_covariates_checkbox) {
					#Eliminates effects from 2008 coding change
					covars <- ds[[age_group]][, good_covars, drop = FALSE]
					month_i <- as.factor(as.numeric(format(as.Date(as.character(ds[[age_group]][, locked_input()$date])), '%m')))
					spline <- setNames(as.data.frame(bs(1:nrow(covars), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
					year_2008 <- numeric(nrow(covars))
					year_2008[1:nrow(covars) >= match(as.Date('2008-01-01'), as.Date(as.character(ds[[age_group]][, locked_input()$date])))] <- 1
					data <- cbind.data.frame(year_2008, spline, month_i)
					trend <- lapply(covars, getTrend, data = data)
					covars <- covars - trend
				} else {
					covars <- ds[[age_group]][, good_covars, drop = FALSE]
				}
				#if (locked_input()$loess_slidebar != 0) {
				#	span <- locked_input()$loess_slidebar
				#	covars <- as.data.frame(sapply(covars, FUN = function(covar) {predict(loess(covar~seq(nrow(covars)), span = span))}))
				#}
				covars <- as.data.frame(lapply(covars[, apply(covars, 2, var) != 0, drop = FALSE], scale))
				return(covars)
			}), age_groups)
			return(covars)
		}
	}) #Possibly add 2nd model and 3rd model for covars_ach and covars_noj
	#Add conditionals for model selection
	covars_time <- eventReactive(covars(), {
		setNames(lapply(covars(), FUN = function(covars) {
			as.data.frame(list(time_index = 1:nrow(covars)))
		}), age_groups())
	}) #Add conditionals for model selection
	run_sensitivity <- eventReactive(covars(), {
		results <- sapply(covars(), FUN = function(covars) {return(length(colnames(covars)) > 3)})
		return(locked_input()$run_sensitivity && results)
	})
	noj_denom <- eventReactive(ds_log(), {
		ds <- ds_log()
		setNames(sapply(ds_log(), FUN = function(ds) {
			if (locked_input()$noj_denom %in% colnames(ds)) {
				denom <- ds[, locked_input()$noj_denom]
			} else {
				denom <- ds[, 1]
				denom[] <- 0
				print(paste(locked_input()$noj_denom, 'is not in dataset.'))
			}
			if (locked_input()$no_noj_denom) {
				denom[] <- 0
			}
			return(denom)
		}), age_groups())})
	
	outcome <- eventReactive(ds_log(), {sapply(ds_log(), FUN = function(data) {scale(data[, locked_input()$variable])})})
	outcome_mean <- eventReactive(ds_log(), {sapply(ds_log(), FUN = function(data) {mean(data[, locked_input()$variable])})})
	outcome_sd <- eventReactive(ds_log(), {sapply(ds_log(), FUN = function(data) {sd(data[, locked_input()$variable])})})
	outcome_plot <- eventReactive(list(outcome(), outcome_sd(), outcome_mean()), {exp(t(t(outcome()) * outcome_sd() + outcome_mean()))})
	outcome_offset <- eventReactive(noj_denom(), {
		ds <- ds_log()
		noj_denom <- noj_denom()
		age_groups <- age_groups()
		null_outcome <- sapply(age_groups, FUN = function(age_group) {ds[[age_group]][, locked_input()$variable] - noj_denom[, age_group]})
		return(null_outcome)
	})
	outcome_offset_mean <- reactive({colMeans(outcome_offset())})
	outcome_offset_sd <- eventReactive(noj_denom(), {
		ds <- ds_log()
		noj_denom <- noj_denom()
		age_groups <- age_groups()
		sapply(age_groups, FUN = function(age_group) {sd(ds[[age_group]][, locked_input()$variable] - noj_denom[, age_group])})
	})
	
	#CausalImpact
	data_full <- eventReactive(list(outcome(), time_points()), {
		age_groups <- age_groups()
		setNames(lapply(age_groups, makeTimeSeries, outcome = outcome(), covars = covars(), time_points = time_points(), scale_outcome = FALSE), age_groups)
	})
	#Add conditional for model selection.
	data_time <- eventReactive(list(outcome(), time_points()), {
		age_groups <- age_groups()
		setNames(lapply(age_groups, makeTimeSeries, outcome = outcome(), covars = covars_time(), time_points = time_points(), scale_outcome = TRUE), age_groups)
	}) #Add conditional for model selection.
	data_offset <- eventReactive(list(outcome_offset(), time_points()), {
		age_groups <- age_groups()
		setNames(lapply(age_groups, makeTimeSeries, outcome = outcome_offset(), covars = NULL, time_points = time_points(), scale_outcome = TRUE), age_groups)
	})
	
	impact_full <- eventReactive(data_full(), {
		age_groups <- age_groups()
		zoo_data <- data_full()
		time_points <- time_points()
		intervention_date <- as.Date(as.character(locked_input()$training_range))[2]

		progress$set(message = 'Running Selected Model.')
		
		cl <- makeCluster(n_cores)
		clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
		clusterExport(cl, c('zoo_data', 'doCausalImpact', 'time_points', 'intervention_date', 'age_groups'), envir = environment())
		
		impact_full <- setNames(parLapply(cl, zoo_data, doCausalImpact, intervention_date = intervention_date, time_points = time_points), age_groups)
		
		stopCluster(cl)
		customIncProgress(amount = 0.2)
		
		return(impact_full)
	})
	impact_time <- eventReactive(data_time(), {
		age_groups <- age_groups()
		zoo_data <- data_time()
		time_points <- time_points()
		intervention_date <- as.Date(as.character(locked_input()$training_range))[2]
		
		progress$set(message = 'Running Time Trend Model.')
		
		cl <- makeCluster(n_cores)
		clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
		clusterExport(cl, c('zoo_data', 'doCausalImpact', 'time_points', 'intervention_date', 'age_groups'), envir = environment())
		
		impact_time <- setNames(parLapply(cl, zoo_data, doCausalImpact, intervention_date = intervention_date, time_points = time_points, trend = TRUE), age_groups)
		
		stopCluster(cl)
		
		customIncProgress(amount = 0.2)
		
		return(impact_time)
	})
	impact_offset <- eventReactive(outcome_offset(), {
		age_groups <- age_groups()
		zoo_data <- data_offset()
		time_points <- time_points()
		intervention_date <- as.Date(as.character(locked_input()$training_range))[2]

		progress$set(message = 'Running Offset Model.')
		
		cl <- makeCluster(n_cores)
		clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
		clusterExport(cl, c('zoo_data', 'doCausalImpact', 'time_points', 'intervention_date', 'age_groups'), envir = environment())
		impact_offset <- setNames(parLapply(cl, zoo_data, doCausalImpact, intervention_date = intervention_date, time_points = time_points, offset = TRUE), age_groups)
		stopCluster(cl)
		
		customIncProgress(amount = 0.2)
		return(impact_offset)
	})
	
	inclusion_prob_full <- eventReactive(impact_full(), {
		age_groups <- age_groups()
		impact <- impact_full()
		setNames(mapply(inclusionProb, age_groups, impact, SIMPLIFY = FALSE), age_groups)
	})
	inclusion_prob_time <- eventReactive(impact_time(), {
		age_groups <- age_groups()
		impact <- impact_time()
		setNames(mapply(inclusionProb, age_groups, impact, SIMPLIFY = FALSE), age_groups)
	})
	
	#Model results
	quantiles_full   <- eventReactive(impact_full(),   {setNames(lapply(age_groups(), FUN = function(age_group) {rrPredQuantiles(impact = impact_full()[[age_group]], mean = outcome_mean()[age_group], sd = outcome_sd()[age_group], eval_period = locked_input()$eval_range, post_period = post_period())}), age_groups())})
	quantiles_time   <- eventReactive(impact_time(),   {setNames(lapply(age_groups(), FUN = function(age_group) {rrPredQuantiles(impact = impact_time()[[age_group]], mean = outcome_mean()[age_group], sd = outcome_sd()[age_group], eval_period = locked_input()$eval_range, post_period = post_period())}), age_groups())})
	quantiles_offset <- eventReactive(impact_offset(), {setNames(lapply(age_groups(), FUN = function(age_group) {rrPredQuantiles(impact = impact_offset()[[age_group]], all_cause_data = noj_denom()[, age_group], mean = outcome_offset_mean()[age_group], sd = outcome_offset_sd()[age_group], eval_period = locked_input()$eval_range, post_period = post_period(), offset = TRUE)}), age_groups())})
	
	pred_quantiles_full   <- eventReactive(quantiles_full(),   {sapply(quantiles_full(),   getPred, simplify = 'array')})
	pred_quantiles_time   <- eventReactive(quantiles_time(),   {sapply(quantiles_time(),   getPred, simplify = 'array')})
	pred_quantiles_offset <- eventReactive(quantiles_offset(), {sapply(quantiles_offset(), getPred, simplify = 'array')})
	
	rr_roll_full <- eventReactive(quantiles_full(), sapply(quantiles_full(), FUN = function(quantiles_full) {quantiles_full$roll_rr}, simplify = 'array'))
	rr_roll_time <- eventReactive(quantiles_time(), sapply(quantiles_time(), FUN = function(quantiles_time) {quantiles_time$roll_rr}, simplify = 'array'))
	
	rr_mean_full   <- eventReactive(quantiles_full(),   {quantiles <- cbind(age_groups(), t(sapply(quantiles_full(),   getRR))); colnames(quantiles)[1] <- locked_input()$age_group; return(quantiles)})
	rr_mean_time   <- eventReactive(quantiles_time(),   {quantiles <- cbind(age_groups(), t(sapply(quantiles_time(),   getRR))); colnames(quantiles)[1] <- locked_input()$age_group; return(quantiles)})
	rr_mean_offset <- eventReactive(quantiles_offset(), {quantiles <- cbind(age_groups(), t(sapply(quantiles_offset(), getRR))); colnames(quantiles)[1] <- locked_input()$age_group; return(quantiles)})
	
	#Cumulative sum of prevented cases
	plot_cumsum_prevented <- eventReactive(quantiles_full(), {
		sapply(age_groups(), FUN = function(age_group, post_period, quantiles) {
			pred_samples_post_full <- quantiles[[age_group]]$pred_samples_post_full
			
			time_points <- time_points()
			post_period_start <- which(time_points == post_period[1]) 
			post_period_end <- which(time_points == post_period[2]) 
			is_post_period <- which(time_points >= post_period[1])
			is_pre_period <- which(time_points < post_period[1])
			
			cases_prevented <- pred_samples_post_full - outcome_plot()[, age_group]
			cumsum_cases_prevented_post <- apply(cases_prevented[is_post_period,], 2, cumsum)
			cumsum_cases_prevented_pre <- apply(cases_prevented[is_pre_period,, drop = FALSE], 2, cumsum)
			cumsum_cases_prevented_pre[,] <- 0
			cumsum_cases_prevented <- rbind(cumsum_cases_prevented_pre, cumsum_cases_prevented_post)
			plot_cumsum_prevented <- t(apply(cumsum_cases_prevented, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
		}, post_period = post_period(), quantiles = quantiles_full(), simplify = 'array')})
	
	sensitivity_analysis_full <- eventReactive(impact_full(), {
		if (!run_sensitivity()) {counter <<- counter + 1; return(counter)}
		age_groups <- age_groups()
		ds <- ds_log()
		covars <- covars()
		impact_full <- impact_full()
		time_points <- time_points()
		intervention_date <- as.Date(as.character(locked_input()$training_range))[2]
		outcome <- outcome()

		progress$set(message = 'Running Sensitivity Analysis.')
		
		cl <- makeCluster(n_cores)
		clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
		clusterExport(cl, c('ds', 'doCausalImpact', 'sensitivityAnalysis', 'age_groups', 'intervention_date', 'outcome', 'time_points'), envir = environment())
		
		#Sensitivity Analysis
		sensitivity_analysis <- setNames(parLapply(cl, age_groups, sensitivityAnalysis, ds = ds, covars = covars, impact = impact_full, intervention_date = intervention_date, outcome = outcome, time_points = time_points), age_groups)
		
		stopCluster(cl)
		
		customIncProgress(amount = 0.6)
		
		return(sensitivity_analysis)
	})
	
	#Table of Rate Ratios for each Analysis Level
	mean_rr_table <- eventReactive(sensitivity_analysis_full(), t({
		if (!run_sensitivity()) {counter <<- counter + 1; return(counter)}
		sensitivity_analysis <- sensitivity_analysis_full()
		age_groups <- isolate(age_groups())
		impact_full <- isolate(impact_full())
		ds <- isolate(ds_log())
		outcome_mean <- isolate(outcome_mean())
		outcome_sd <- isolate(outcome_sd())
		sapply(age_groups, FUN = function(age_group) {
			rr_pred_quantile <- rrPredQuantiles(impact = impact_full[[age_group]], all_cause_data = noj_denom()[, age_group], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = locked_input()$eval_range, post_period = post_period())
			cred_int <- c(age_group, 'Initial' = round(rr_pred_quantile$mean_rr, 4), 'Initial .95' = paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
			names(cred_int)[1] <- locked_input()$age_group
			cred_int_analyses <- lapply(1:length(sensitivity_analysis[[age_group]]), FUN = function(i) {
				rr_pred_quantile <- rrPredQuantiles(impact = sensitivity_analysis[[age_group]][[i]], all_cause_data = noj_denom()[, age_group], mean = outcome_mean[age_group], sd = outcome_sd[age_group], eval_period = locked_input()$eval_range, post_period = post_period())
				cred_int <- c(round(rr_pred_quantile$mean_rr, 4), paste('(', round(rr_pred_quantile$rr[1], 4), ',', round(rr_pred_quantile$rr[3], 4), ')', sep = ''))
				names(cred_int) <- c(paste('Analysis', i), paste('Analysis', i, '.95'))
				return(cred_int)
			})
			c(cred_int, cred_int_analyses, recursive = TRUE)
		})
	}))
	
	covar_warning <- reactive({
		setNames(lapply(age_groups(), FUN = function(age_group) {
			covar_list <- covar_list()
			bad_covars <- rownames(covar_list[covar_list[, age_group] == FALSE, age_group, drop = FALSE])
			if (is.null(bad_covars) || length(bad_covars) == 0) {
				return('')
			} else if (length(bad_covars) == 1) {
				return(paste('Note that', bad_covars, 'was removed due to insufficient data.'))
			} else {
				return(paste('Note that', paste(bad_covars, collapse = ", "), 'were removed from the analysis due to insufficient data.'))
			}
		}), age_groups())
	})
	age_group_tabs <- eventReactive(mean_rr_table(), {
		ds <- ds_log()
		covars <- covars()
		age_groups <- age_groups()
		pred_quantiles_full <- pred_quantiles_full()
		pred_quantiles_time <- pred_quantiles_time()
		pred_quantiles_offset <- pred_quantiles_offset()
		inclusion_prob <- inclusion_prob_full()
		tabs <- lapply(age_groups, FUN = function(age_group) {
			dates <- as.Date(as.character(ds[[age_group]][, locked_input()$date]))
			tab <- tabPanel(
				title = age_group, 
				renderPlot(plotPred(data = pred_quantiles_full()[, , age_group], time_points = time_points(), post_period = post_period(), pred_quantiles_full = pred_quantiles_full()[, , age_group], pred_quantiles_offset = pred_quantiles_offset()[, , age_group], outcome_plot = outcome_plot()[, age_group], fix_2008 = locked_input()$covariate_checkbox, title = 'Full Synthetic Control Plot')),
				renderPlot(plotPred(data = pred_quantiles_time()[, , age_group], time_points = time_points(), post_period = post_period(), pred_quantiles_full = pred_quantiles_full()[, , age_group], pred_quantiles_offset = pred_quantiles_offset()[, , age_group], outcome_plot = outcome_plot()[, age_group], fix_2008 = locked_input()$covariate_checkbox, title = 'Interrupted Time Series Plot')),
				renderPlot(plotPred(data = pred_quantiles_offset()[,, age_group], time_points = time_points(), post_period = post_period(), pred_quantiles_full = pred_quantiles_full()[, , age_group], pred_quantiles_offset = pred_quantiles_offset()[, , age_group], outcome_plot = outcome_plot()[, age_group], fix_2008 = locked_input()$covariate_checkbox, offset = TRUE, title = 'Simple Pre-post Comparison Plot')),
				renderPlot({matplot(rr_roll_full()[, , age_group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l', main = paste(age_group, 'Synthetic controls: rolling rate ratio'), ylab = 'Rate Ratio'); abline(h = 1, lty = 2)}),
				renderPlot({matplot(rr_roll_time()[, , age_group], type = 'l', xlim = c(12, 72), ylim = c(0.3, 1.7), col = 'black', bty = 'l', main = paste(age_group, 'Interupted time series: rolling rate ratio'), ylab = 'Rate Ratio'); abline(h = 1, lty = 2)}),
				renderPlot(plot(impact_full()[[age_group]]$model$bsts.model, 'coefficients', main = age_group)),
				renderPlot(plot(zoo(covars[[age_group]][, rownames(inclusion_prob[[age_group]])[2:length(rownames(inclusion_prob[[age_group]]))]], dates), plot.type = 'single', type = 'l', main = 'Weight Assigned to Different Covariates', xlab = 'Time', ylab = 'Scaled Values', col = rgb(0, 0, 0, alpha = unlist(inclusion_prob[[age_group]][2:length(rownames(inclusion_prob[[age_group]])),])))), #col = rainbow(ncol(covars[[age_group]]))
				renderText(covar_warning()[[age_group]])
			)
			return(tab)
		})
		if (!run_sensitivity()) {mean_rr_table <- NULL; sensitivity_header <- NULL} else {mean_rr_table <- mean_rr_table(); sensitivity_header <- "Rate ratios and credible intervals from sensitivity analyses."}
		summary_tab <- tabPanel(
			title = 'Summary',
			h3("Rate ratios and credible intervals from main analyses."),
			renderTable(rr_mean_full(), caption = 'Full Synthetic Control Model Quantiles', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL)),
			renderTable(rr_mean_time(), caption = 'Time Trend Model Quantiles', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL)),
			renderTable(rr_mean_offset(), caption = 'Offset Model Quantiles', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL)),
			h3(sensitivity_header),
			renderTable(mean_rr_table, caption = 'Rate Ratios and Credibility Intervals', caption.placement = getOption('xtable.caption.placement', 'top'), caption.width = getOption('xtable.caption.width', NULL))
		)
		tabs <- append(list(summary_tab), tabs)
		return(tabs)
	})
	output$tab_view <- renderUI({
		age_groups <- age_groups()
		#if (length(age_groups) > 20) {
		#	print('Too many factors - likely chose the wrong column.')
		#	return(NULL)
		#} else if (length(age_groups) == 0) {return(NULL)}
		if (length(age_groups) == 0) {return(NULL)}
		tabs <- age_group_tabs()
		if (!is.null(progress)) {
			progress$set(value = 1.0)
			progress$close()
			progress <<- NULL
		}
		do.call(tabsetPanel, tabs)
	})
})