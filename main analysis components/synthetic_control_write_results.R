#Use this script after running the analysis script to write important results to file.

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

#Output the rate ratio estimates to a new file.
write.csv(rr_mean_full, paste(output_directory, country, '_rr_full.csv', sep = ''))
write.csv(rr_mean_time, paste(output_directory, country, '_rr_time_trend.csv', sep = ''))
write.csv(rr_mean_full_intervals, paste(output_directory, country, '_rr_full_intervals.csv', sep = ''))
write.csv(rr_mean_time_intervals, paste(output_directory, country, '_rr_time_trend_intervals.csv', sep = ''))
write.csv(rr_roll_full, paste(output_directory, country, '_rr_roll_full.csv', sep = ''))
write.csv(rr_roll_time, paste(output_directory, country, '_rr_roll_time.csv', sep = ''))

#WAICs
write.csv(waic_full, paste(output_directory, country, '_waic_full.csv', sep = ''))
write.csv(waic_time, paste(output_directory, country, '_waic_time.csv', sep = ''))

#Output the sensitivity analysis rate ratio estimates to a new file.
# write.csv(sensitivity_analysis_pred_2,  paste(output_directory, country, '_sensitivity_analysis_pred_2.csv',  sep = ''))
# write.csv(sensitivity_analysis_pred_10, paste(output_directory, country, '_sensitivity_analysis_pred_10.csv', sep = ''))
# write.csv(sensitivity_analysis_pred_2_intervals,  paste(output_directory, country, '_sensitivity_analysis_pred_2_intervals.csv',  sep = ''))
# write.csv(sensitivity_analysis_pred_10_intervals, paste(output_directory, country, '_sensitivity_analysis_pred_10_intervals.csv', sep = ''))

#Tables for rate ratios.
write.csv(sensitivity_table, paste(output_directory, country, '_sensitivity_table.csv', sep = ''))
write.csv(sensitivity_table_intervals, paste(output_directory, country, '_sensitivity_table_intervals.csv', sep = ''))
write.csv(rr_table, paste(output_directory, country, '_rr_table.csv', sep = ''))
write.csv(rr_table_intervals, paste(output_directory, country, '_rr_table_intervals.csv', sep = ''))
