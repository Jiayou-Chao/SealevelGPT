source("ComparePlots.R")

GPT_predict <- function(location, year, measurement.type = "sa", add_VLM = F, VLM_na_interp = T, na.rm = F, equal_time = FALSE, ...) {
    #' Returns a description of the prediction for a given year and save it to a local file.
    #' @param location The name of the station.
    #' @param year The year to predict.

    prediction_table <- RSL_predict(location, measurement.type = measurement.type, add_VLM = add_VLM, VLM_na_interp = VLM_na_interp, na.rm = na.rm, equal_time = equal_time, ...)
    prediction_table <- prediction_table[[1]]
    # keep only 2 decimal places
    prediction_table[, 2:ncol(prediction_table)] <- round(prediction_table[, 2:ncol(prediction_table)], 2)
    # return(prediction_table)
    predictions <- prediction_table[prediction_table$Year == year, ]
    #   The column names of the CSV table are "Year", "mean_no", "up95_no", "low95_no", "up99_no", "low99_no", "mean_cop", "up95_cop", "low95_cop", "up99_cop", "low99_cop", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585". "mean_no" means the mean of the prediction under no restriction scenario. "up95_no" means the upper bound of the 95% confidence interval of the prediction under no restriction scenario. "low95_no" means the lower bound of the 95% confidence interval of the prediction under no restriction scenario. "up99_no" means the upper bound of the 99% confidence interval of the prediction under no restriction scenario. "low99_no" means the lower bound of the 99% confidence interval of the prediction under no restriction scenario. "mean_cop" means the mean of the prediction under the restriction scenario. "up95_cop" means the upper bound of the 95% confidence interval of the prediction under the restriction scenario. "low95_cop" means the lower bound of the 95% confidence interval of the prediction under the restriction scenario. "up99_cop" means the upper bound of the 99% confidence interval of the prediction under the restriction scenario. "low99_cop" means the lower bound of the 99% confidence interval of the prediction under the restriction scenario. "ssp119" means the prediction under the SSP1-1.9 scenario. "ssp126" means the prediction under the SSP1-2.6 scenario. "ssp245" means the prediction under the SSP2-4.5 scenario. "ssp370" means the prediction under the SSP3-7.0 scenario. "ssp460" means the prediction under the SSP4-6.0 scenario. "ssp585" means the prediction under the SSP5-8.5 scenario.
    result_string <- paste0(
        "The predicted sea level in ", location, " in ", year, " is ", predictions$mean_no, " m under no restriction scenario and ", predictions$mean_cop, "m under restricted scenario. The 95% confidence interval of the prediction under no restriction scenario is [", predictions$low95_no, ", ", predictions$up95_no, "]m. The 95% confidence interval of the prediction under restricted scenario is [", predictions$low95_cop, ", ", predictions$up95_cop, "]m. The 99% confidence interval of the prediction under no restriction scenario is [", predictions$low99_no, ", ", predictions$up99_no, "]m. The 99% confidence interval of the prediction under restriction scenario is [", predictions$low99_cop, ", ", predictions$up99_cop, "]m.",
        ifelse(!is.na(predictions$ssp119), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP1-1.9 scenario is ", predictions$ssp119, "m."), ""),
        ifelse(!is.na(predictions$ssp126), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP1-2.6 scenario is ", predictions$ssp126, "m."), ""),
        ifelse(!is.na(predictions$ssp245), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP2-4.5 scenario is ", predictions$ssp245, "m."), ""),
        ifelse(!is.na(predictions$ssp370), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP3-7.0 scenario is ", predictions$ssp370, "m."), ""),
        ifelse(!is.na(predictions$ssp460), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP4-6.0 scenario is ", predictions$ssp460, "m."), ""),
        ifelse(!is.na(predictions$ssp585), paste0(" The predicted sea level in ", location, " in ", year, " under the SSP5-8.5 scenario is ", predictions$ssp585, "m."), "")
    )
    # Write the result to a local file called "prediction.txt".
    write(result_string, file = "temporary_prediction.txt")
    return(result_string)
}
