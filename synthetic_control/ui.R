library(shiny, quietly = TRUE)

shinyUI(fluidPage(
	
	# Application title
	titlePanel('Synthetic Control'),
	
	sidebarLayout(sidebarPanel(
		fileInput(
			inputId = 'in_file',
			label = 'Choose a file* to upload.', 
			accept = c(
				'.csv'
			)
		), 
		tags$hr(),
		selectInput(
			inputId = 'age_group', 
			label = 'Select column for subgroups (e.g., age category).',
			choices = NULL, 
			selected = NULL,
			selectize = TRUE
		), 
		selectInput(
			inputId = 'age_group_selected', 
			label = 'Select subgroup factors.',
			choices = NULL, 
			selected = NULL,
			selectize = TRUE, 
			multiple = TRUE
		),
		selectInput(
			inputId = 'date', 
			label = 'Select date column.',
			choices = NULL, 
			selected = NULL,
			selectize = TRUE
		), 
		selectInput(
			inputId = 'variable', 
			label = 'Select outcome variable.',
			choices = NULL, 
			selected = NULL,
			selectize = TRUE
		), 
		checkboxInput(
			inputId = 'covariate_checkbox',
			label = 'List of synthetic controls to exclude (if checked).',
			value = FALSE
		),
		selectInput(
			inputId = 'covariate',
			label = 'Select synthetic control columns.',
			choices = NULL,
			selected = NULL,
			selectize = TRUE, 
			multiple = TRUE
		),
		#sliderInput(
		#	inputId = 'loess_slidebar',
		#	label = 'Loess smoothing (0 for none, 1 for full).',
		#	min = 0, 
		#	max = 1, 
		#	value = 0, 
		#	step = 0.01
		#),
		checkboxInput(
			inputId = 'adjust_covariates_checkbox',
			label = 'Adjust for year 2008 coding change (Brazil only).',
			value = FALSE
		),
		conditionalPanel(
			condition = '!input.no_noj_denom',
			selectInput(
				inputId = 'noj_denom',
				label = 'Select denominator for offset model.', 
				choices = NULL, 
				selectize = TRUE
			)
		),
		checkboxInput( #TODO: Implement
			inputId = 'no_noj_denom',
			label = 'No denominator (will use 1 as the denominator).',
			value = FALSE
		),
		tags$hr(),
		dateRangeInput(
			inputId = 'training_range',
			label = 'Training Range (YYYY-MM-01)',
			start = NULL, end = NULL
		),
		dateRangeInput(
			inputId = 'eval_range',
			label = 'Evaluation Range (YYYY-MM-01)',
			start = NULL, end = NULL
		),
		tags$hr(),
		conditionalPanel(condition = '(input.covariate != null && input.covariate.length > 3)', 
										 checkboxInput(
										 	inputId = 'run_sensitivity',
										 	label = 'Run sensitivity analysis.',
										 	value = TRUE
										 )),
		actionButton(
			inputId = 'run_CausalImpact', 
			label = 'Run Causal Impact'
		)
	),
	
	mainPanel(align = 'center',
						uiOutput('tab_view')
	)), 
	
	fluidRow(
		column(12, align = 'justify',
		hr(),
		span('*To use this applet properly, the file should be a .csv and the date must be in the YYYY-MM-DD format. We also suggest preprocessing using a method similar to that found in '),
		downloadLink(
			outputId = 'download_preprocesser',
			label = 'this R script,'
		),
		span('which you are free to download and modify. We provide '),
		downloadLink(
			outputId = 'download_sample',
			label = 'an example dataset'
		), 
		span('for testing and comparison purposes.'),
		span('For additional information and answers to commonly asked questions, visit the applet\'s'),
		a(target = '_blank', 'help page.', href = 'https://apps.google.com/products/sites/'),
		
		hr(), 
		span('This applet uses the Causal Impact R package developed by Google. For more information about the package, visit the'),
		a(target = '_blank', 'Causal Impact GitHub page.', href = 'https://google.github.io/CausalImpact/'), 
		
		hr(), 
		span('This app was developed by VIPR (Vaccine ImPact Research), a collaborative network of investigators based at Yale University (Weinberger lab), George Washington University (Lone Simonsen), and Sage Analytica focused on evaluating the impacts of vaccines around the world. VIPR works on developing and using cutting-edge data analysis approaches to estimate the impacts of vaccines from imperfect data sources. The project is supported by a grant from the Bill and Melinda Gates Foundation to Yale University (PI: Weinberger). The synthetic control analyses are based on the model and code developed by '),
		a(target = '_blank', 'Brodersen et al.', href = 'https://research.google.com/pubs/pub41854.html') 
	)),
	fluidRow(
		column(6, align = 'center',
			imageOutput('Sage')
		),
		column(6, align = 'center',
			imageOutput('Yale')
		)
	)
))