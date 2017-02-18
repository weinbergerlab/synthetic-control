#This is the file used to set variables to be used in analysis, as well as to run the analysis.
#Make sure *_analysis.R, *_report.R, *_report.Rmd, *_functions.R, and *_plot.R are all in the same folder as this file.
#Model the setup shown in this file, then run this file from the console using source('This file's directory/This file's name'). 

#Clear the workspace
rm(list = ls(all = TRUE))
gc()

#Set the working directory
#This step only works as written if this file is run using the source() command. Otherwise, skip this step and set manually.
tryCatch({
	setwd(dirname(sys.frame(1)$ofile))
}, error = function(e) {
	cat('Warning: \n  Could not programmatically set working directory. \n  Try running this file by using "source(\'This file\'s directory/This file\'s name\')" from the command line, \n   or set the working directory manually using "setwd()".')
}, warning = function(w) {}, finally = {})

#Used to check for relevant packages and update them if out of date or install them if not installed.
update_packages  <- TRUE #Whether to update outdated packages.
install_packages <- TRUE #Whether to install missing packages.
install_pandoc   <- TRUE #Whether to install pandoc, which requires an external installer, and rmarkdown, a package that depends on pandoc's successful installation.

#Assign variable values
country     <- 'Brazil'  #Country or region name.
n_seasons   <- 3         #Number of months (seasons) per year. 12 for monthly, 4 for quarterly, 3 for trimester data.
exclude     <- c('ACM-NoPCV', 'J00-06', 'J09-18', 'J20-22', 'J30-39', 'J40-47', 'J60-70', 'J80-84', 'J85-J86', 'J90-94', 'J95-99', 'A30-49', 'G00-09', 'H65-75', 'B95-98') #User-defined list of covariate columns to exclude from all analyses.
exclude     <- c(exclude, gsub('-', '_', exclude))
code_change <- FALSE     #Used for Brazil data. Set to TRUE to adjust for year 2008 coding changes; otherwise, set to FALSE.

input_directory  <- 'https://raw.githubusercontent.com/weinbergerlab/synthetic-control/master/Datasets%20for%20PNAS/' #Directory (or URL) containing input data file.
output_directory <- paste('~/Desktop/Results/', sep = '')                                                             #Directory where results will be saved.
output_directory <- paste(output_directory, format(Sys.time(), '%Y-%m-%d-%H%M%S'), '/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
file_name        <- 'Dataset%20S1%20Brazil.csv'                                                                       #Name of file containing data for analysis. Must be a .csv file.

group_name   <- 'age_group' #Name of column containing group labels.
date_name    <- 'date'      #Name of column containing dates.
outcome_name <- 'J12_18'    #Name of column containing outcome.
denom_name   <- 'ACM_NoPCV' #Name of column containing denominator to be used in offset.

start_date        <- as.Date('2003-01-01') #Indicates the date of the first data point.
intervention_date <- as.Date('2009-04-30') #Indicates the date of intervention in the data.
end_date          <- as.Date('2013-12-01') #Indicates the date of the last data point.
pre_period        <- as.Date(c('2003-01-01', '2009-04-30')) #Range over which the data is trained for the CausalImpact model.
post_period       <- as.Date(c('2009-05-01', '2013-12-01')) #Range from the intervention date to the end date.
eval_period       <- as.Date(c('2012-01-01', '2013-12-01')) #Range over which rate ratio calculation will be performed.

#Run analysis, but don't generate HTML report
#source('synthetic_control_analysis.R')
#source('synthetic_control_plot.R')

#Run analysis and generate HTML report
source('synthetic_control_report.R')