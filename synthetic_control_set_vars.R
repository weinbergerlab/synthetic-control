#Used when running analysis to programmatically set constants and load data for analysis.

#Clear workspace
rm(list = ls(all = TRUE))
gc()

library(lubridate)

countries <- c('Brazil', 'Ecuador', 'Chile', 'Mexico', 'US',
							 'SC1', 'SC2', 'SC3', 'SC_HDI', 'SC_HDI_Region', 'SC_National', 'SC_HDI_No_Pop', 'SC_HDI_Region_No_Pop', 'SC_National_No_Pop', 'SC_HDI_No_Pop_All', 'SC_HDI_Region_No_Pop_All', 'SC_National_No_Pop_All',
							 'SC_Subchapter_HDI', 'SC_Subchapter_HDI_Region', 'SC_Subchapter_National',
							 'SC_Dec/Coverage groups', 'SC_Dec/HDI from Muni', 'SC_Dec/HDIxSuperRegion', 'SC_Dec/National', 'SC_Dec/Region')

for (country in countries) {
print(country)

#############################
#                           #
#User-initialized constants.#
#                           #
#############################

group_name <- 'age_group'
date_name <- 'date'
n_seasons <- NULL #12 for monthly, 4 for quarterly, 3 for trimester data.
exclude_covar <- c('ACM-NoPCV', 'J00-06', 'J09-18', 'J20-22', 'J30-39', 'J40-47', 'J60-70', 'J80-84', 'J85-J86', 'J90-94', 'J95-99', 'A30-49', 'G00-09', 'H65-75', 'B95-98')
exclude_covar <- c(exclude_covar, gsub('-', '_', exclude_covar))
update_packages <- TRUE
install_packages <- TRUE
install_pandoc <- TRUE

##################################################
#                                                #
#Directory setup and initialization of constants.#
#                                                #
##################################################

#Load file (.csv format). You can change the input and output directories to point to an appropriate spot on your computer.
input_directory <- paste('~/Documents/Synthetic Control Data/', country, '/', sep = '')
country <- gsub('/', '_', country)
output_directory <- paste(input_directory, 'results/', sep = '')
output_directory <- paste(output_directory, format(Sys.time(), '%Y-%m-%d-%H%M%S'), '/', sep = '')
file_name <- paste('prelog', country, 'processed', 'data.csv', sep = '_')
pcv_file <- paste(input_directory, file_name, sep = '')
ds1a <- read.csv(pcv_file, check.names = FALSE)
ds1a[, date_name] <- as.Date(ds1a[, date_name])

#Account for code-naming differences
#outcome_name - Gives the outcome variable (e.g. pneumonia) 
#denom_name - Gives the name of denominator used for some of the analyses. Can be population size, non-respiratory hospitalization, etc.
SC_set <- c('SC1', 'SC2', 'SC3', 'SC_HDI', 'SC_HDI_Region', 'SC_National', 'SC_HDI_No_Pop', 'SC_HDI_Region_No_Pop', 'SC_National_No_Pop', 'SC_HDI_No_Pop_All', 'SC_HDI_Region_No_Pop_All', 'SC_National_No_Pop_All')
SC_subchapter_set <- c('SC_Subchapter_HDI', 'SC_Subchapter_HDI_Region', 'SC_Subchapter_National')
SC_Dec_set <- c('SC_Dec_Coverage groups', 'SC_Dec_HDI from Muni', 'SC_Dec_HDIxSuperRegion', 'SC_Dec_National', 'SC_Dec_Region', 'SC_Dec/Coverage groups', 'SC_Dec/HDI from Muni', 'SC_Dec/HDIxSuperRegion', 'SC_Dec/National', 'SC_Dec/Region')
if (country %in% SC_set || country %in% SC_subchapter_set) {
	denom_name <- 'ACM_NoJ'
	outcome_name <- 'J12_18'
} else if (country %in% SC_Dec_set) {
	denom_name <- 'ACM_NoPCV'
	outcome_name <- 'J12_18'
} else if (country == 'US') {
	denom_name <- 'ach_sid_noresp'
	outcome_name <- 'm_ACP'
} else {
	denom_name <- 'ach_noj'
	outcome_name <- 'J12_18'
}

#Date variables
#start_date gives training period start date.
#intervention_date gives training period end date.
#end_date gives end of data_set.
if (country == 'Brazil') {
	start_date <- as.Date('2004-01-01')
} else {
	start_date <- min(as.Date(ds1a[, date_name]))
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
end_date <- max(as.Date(ds1a[, date_name]))
pre_period <- c(start_date, intervention_date) #Define training period
post_period <- c(unique(ds1a[, date_name])[which(abs(unique(ds1a[, date_name]) - intervention_date) == min(abs(unique(ds1a[, date_name])[unique(ds1a[, date_name]) >= intervention_date] - intervention_date)))], end_date) #Define post-vaccine period.
#Define evaluation period.
if (country == 'Brazil') {
	eval_period <- c(as.Date(ds1a[, date_name][nrow(ds1a) - 23]), as.Date(ds1a[, date_name][nrow(ds1a)]))
} else if (country %in% SC_set || country %in% SC_subchapter_set || country %in% SC_Dec_set) {
	eval_period <- c(end_date - years(2), end_date)
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

if (country == 'Brazil') {
	code_change <- TRUE
} else {
	code_change <- FALSE
}

#Run analysis and generate report
source('synthetic_control_report.R', local = TRUE)
}