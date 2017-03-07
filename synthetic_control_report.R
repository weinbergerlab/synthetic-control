#Call this file after defining the relevant variables to generate the HTML report.
source('synthetic_control_functions.R', local = TRUE)
packages <- c('curl', 'evaluate', 'digest', 'formatR', 'highr', 'markdown', 'stringr', 'yaml', 'Rcpp', 'htmltools', 'caTools', 'bitops', 'knitr', 'jsonlite', 'base64enc', 'rprojroot', 'rmarkdown')
packageHandler(packages, update_packages, install_packages)
if (install_pandoc) {
	if (Sys.info()['sysname'] == 'Windows') {
		packageHandler(c('installr', 'rmarkdown'), update_packages, install_packages)
		if (!rmarkdown::pandoc_available()) {
			installr::install.pandoc()
		}
	} else {
		if (!(rmarkdown::pandoc_available())) {
			warning('This system cannot programmatically install/update Pandoc. To install/update Pandoc, visit "https://pandoc.org/installing.html".')
		}
	}
}
if (!exists('exclude_group')) {exclude_group <- c()}

sapply(packages, library, quietly = TRUE, character.only = TRUE)
param_list <- list(update_packages = update_packages,
							install_packages = install_packages,
							
							country       = country, 
							n_seasons     = n_seasons, 
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
							eval_period       = eval_period)
run_pandoc <- rmarkdown::pandoc_available()
if (run_pandoc) {
	rmarkdown::render('synthetic_control_report.Rmd', output_file = 'Synthetic Control Report.html', output_dir = output_directory, params = param_list, envir = environment())	
} else {
	knitr::knit('synthetic_control_report.Rmd', envir = environment())
	markdown::markdownToHTML('synthetic_control_report.md', output = paste(output_directory, 'synthetic_control_report.html', sep = ''))
	file.remove('synthetic_control_report.md')
	unlink('figure/', recursive = TRUE)
}
source('synthetic_control_write_results.R', local = TRUE)