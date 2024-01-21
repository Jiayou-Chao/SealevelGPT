library(magrittr)
library(tidyverse)
library(corrplot)
DATA_DIR <- "data/plots"
station.names=list.dirs(DATA_DIR, full.names = F, recursive = F)

##########
#### New function for converting daily VLM to monthly (11/6/23)
convert_daily_VLM_to_monthly <- function(data, date_col = 'month', value_col = 'VLM', na.rm = F) {
  #' Converts daily VLM data to monthly VLM data. The output is YYYY-MM-01 format.
  #' Deprecated. Use convert_daily_to_monthly() instead.
  #' 
  #' 
  warning("Deprecated. Use convert_daily_to_monthly() instead.")
  if (value_col == 'VLM') {
    data1 <- data %>%
              mutate(month1 = month(get(date_col)), year = year(get(date_col))) %>%
              group_by(year, month1) %>%
              summarise_all(mean)
    data1$month <- as.Date(paste0(data1$year,'-',data1$month1,'-01'))
    # return(select(data1, month, VLM))
    data1 %<>% 
      ungroup() %>% 
      select(month, VLM)
    return(data1)
    
  }
  data <- data %>%
    mutate(month1 = month(get(date_col)), year = year(get(date_col))) %>%
    group_by(year, month1) %>%
    summarise(value = mean(get(value_col), na.rm = na.rm))
  data %>% select(-month1)
}

convert_daily_to_monthly <- function(data, month_col = 'month', na.rm = F) {
  #' Converts daily VLM data to monthly VLM data. The output is YYYY-MM-01 format.
  #' Columns with numeric values will be averaged. Columns with non-numeric values will be dropped.

  data_num <- data %>% 
    select(where(is.numeric), month)

  data_num <- data_num %>%
    mutate(month1 = month(get(month_col)), year = year(get(month_col))) %>%
    group_by(year, month1) %>%
    summarise_all(mean, na.rm = na.rm)
  data_num$month <- as.Date(paste0(data_num$year,'-',data_num$month1,'-01'))
  data_num %>% 
    ungroup() %>% 
    select(-month1, -year) %>% 
    select(month, everything())

}


##########

# convert_daily_VLM_to_monthly <- function(daily_VLM) {
#   #' Convert daily VLM data to monthly VLM data
#   #' @param daily_VLM a data frame with daily VLM data
#   #' @return a data frame with monthly VLM data
#   #' @export
#   #' @examples
#   #' convert_daily_VLM_to_monthly(daily_VLM)
#   daily_VLM <- daily_VLM %>% 
#     mutate(month = format(Date, '%Y-%m')) %>% 
#     group_by(month) %>% 
#     summarise(VLM = mean(VLM))
#   daily_VLM$month <- as.Date(daily_VLM$month)
#   return(daily_VLM)
# }

# merge_monthly_data <- function(..., by = 'month') {
#   #' Merge monthly VLM data
#   #' @param ... a list of data frames with monthly VLM data
#   #' @param by the column name to merge by
#   #' @return a data frame with merged monthly VLM data
#   #' @export
#   #' @examples
#   #' merge_monthly_data(monthly_VLM_1, monthly_VLM_2, by = 'month')
#   Reduce(function(x, y) merge(x, y, by = by), list(...))
# }

merge_monthly_data <- function(datasets, by = 'month') {
  #' Merge monthly VLM data
  #' @param datasets a list of data frames with monthly VLM data
  #' @param by the column name to merge by
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' monthly_VLM_1 <- get_VLM_by_station(station.names[1], method = 1)
  #' monthly_VLM_2 <- get_VLM_by_station(station.names[1], method = 2)
  #' merge_monthly_data(c(monthly_VLM_1, monthly_VLM_2), by = 'month')
  #' 
  #' 
  for (dataset in datasets) {
    dataset[by] <- as.Date(dataset[[by]])
  }
  Reduce(function(x, y) merge(x, y, by = by, all = TRUE), datasets)
}

beautiful_corplot <- function(cor_data) {
  #' Make a beautiful correlation plot
  #' @param cor_data a data frame with correlation data, usually the output of cor().
  #' @return a correlation plot
  #' @export
  #' @examples
  #' data <- merge_VLM_station("south_beach")
  #' cor_data <- data %>% select(-month) %>% cor(use = "complete.obs")
  #' beautiful_corplot(cor_data)
  corrplot(cor_data,
         method = 'circle',
         type = 'lower',
         sig.level = c(.001, .01, .05), 
         tl.pos="lt", 
         tl.col="black", tl.cex=1.3, 
         tl.offset=0.2,
         cl.pos="r",
         insig = "label_sig",
         pch.cex = 1.3,
         pch.col="red",
         cl.cex = 1.3)
  corrplot(cor_data,  type="upper", method="number",
          col="coral4",  tl.pos="n", cl.pos="n", number.cex = 1.2, add=T,diag=F)
  recordPlot()
}
