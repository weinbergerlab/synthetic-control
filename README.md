# Synthetic Control

##Analysis
To run the analyses, use the code found in "synthetic_control_analysis.R". Save this and the the "synthetic_control_functions.R" code in the same directory, and change the input directory to point to the folder containing the data.

To generate reports in HTML format, knit the RMarkdown file "synthetic_control_output.Rmd". Make sure both "synthetic_control_functions.R" and "synthetic_control_analysis_function.R" are in the same folder. To set the input directory, in the call to sytheticControls of the .Rmd file set the argument `use_defaults` to FALSE. You can then manually set the input directory (as well as the other arguments necessary for analysis). 
To change the output directory of the .Rmd, change the assignment of the variable `out_dir` found in the knit portion of the YAML header.

##Sample Data

The sample data provided are CSV files properly formatted for use with the synthetic control analysis. There are five files, each containing data from different countries. The data are subsetted by age group, which can be translated using the list below.

|Numeric Representation | Age Group |
| --- | ------------ |
|  0  | < 3 months |
|  1  | 3-<12 months |
|  2  | 12-24 months |
|  3  | 24-59 months |
|  4  | 5-<18 years |
|  5  | 18-<40 years |
|  6  | 40-<65 years |
|  7  | 65-<80 years |
|  8  | 80+ years |
|  9  | <12 months |

For the US, cells with fewer than 10 counts are replaced with 9999 due to privacy considerations.

##Web Applet

The synthetic control applet is a resource-intensive process. To handle concurrent users, there are multiple instances, listed below.

https://weinbergerlab.shinyapps.io/synthetic_control/

https://weinbergerlab.shinyapps.io/synthetic_control_1/

https://weinbergerlab.shinyapps.io/synthetic_control_2/
