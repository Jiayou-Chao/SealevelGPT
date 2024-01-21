library(maps)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(zoo)
library(tidyverse)
library(docstring)
source("utils.R")
library(geosphere)

vlm_rename_dict <- function(type = "method") {
  #' Get the rename dictionary for VLM data
  #' @param type character, "method" or "rename". The names of the methods are "1", "2", "3", "4", which correspond to "ULR", "NGL", "JPL", "GFZ".
  #' @return a dictionary
  #' @export
  #' @examples
  #' vlm_rename_dict("method")
  method_dict <- c("1" = "ULR", "2" = "NGL", "3" = "JPL", "4" = "GFZ")
  rename_dict <- c("ULR" = "1", "NGL" = "2", "JPL" = "3", "GFZ" = "4")
  if (type == "method") {
    return(method_dict)
  }
  else if (type == "rename") {
    return(rename_dict)
  }
  else {
    stop("Invalid type. Please use 'method' or 'rename'.")
  }
}

get_VLM_old<-function(lan_input,lon_input,verbose=FALSE){
  location<-read.table('data/VLM/JPL_location.txt',head=TRUE)
  location<-location[location$X2=="POS",]
  distance=(location$N-lan_input)^2+(location$E-lon_input)^2
  VLM_station_name=location[which.min(distance),]$Name
  if(verbose){
    print('The info of the nearest VLM station: ')
    print(location[which.min(distance),])
  }

  VLM_ts=read.table(paste('https://sideshow.jpl.nasa.gov/pub/JPL_GPS_Timeseries/repro2018a/post/point/',VLM_station_name,'.series',sep='')
                    ,head=FALSE)
  
  
  
  
  VLM_ts=data.frame(as.Date(date_decimal(VLM_ts$V1)),VLM_ts$V4)
  names(VLM_ts)=c('Date','VLM')
  VLM_ts$month <- format(VLM_ts$Date, '%Y-%m')
  VLM_ts_monthly<- data.frame(VLM = tapply(VLM_ts$VLM,VLM_ts$month, mean))
  
  
  all_month=paste(rownames(VLM_ts_monthly),'-01',sep='')
  
  VLM_ts_monthly=data.frame(month=all_month,
                            VLM=VLM_ts_monthly$VLM)
  VLM_ts_monthly <- VLM_ts_monthly %>% mutate(month = as.Date(month))
  
  return(VLM_ts_monthly)
  
}

read_neu <- function(neu_path, convert_decimal_year = T) {
  # Test file: 'data/VLM/ULR7a_neu/d0ARK_99999M018_ULR7.neu'
  # list files: list.files('data/VLM/ULR7a_neu', pattern = '*.neu', full.names = TRUE)
  neu <- read.table(neu_path, header = FALSE, comment.char = "#")
  # Set column names to #  Year      DN (m)      DE (m)      DU (m)   SDN     SDE     SDU
  names(neu) <- c("Year", "DN", "DE", "DU", "SDN", "SDE", "SDU")
  if (convert_decimal_year) {
    neu$Year <- date_decimal(neu$Year)
  }
  return(neu)
}

read_station_table <- function(table_path) {
  station_table <- read.table(table_path, header = FALSE, skip = 13)
  names(station_table) <- c("Site", "DOMES", "Lon", "Lat", "T_GPS", "Data", "V_GPS", "S_GPS", "MODEL")
  return(station_table)
}

read_ULR_table <- function(ULR_table_path = 'data/VLM/ulr7_vertical_velocities.txt') {
  ULR_table <- read_station_table(ULR_table_path)
  return(ULR_table)
}

read_NGL_table <- function(NGL_table_path = 'data/VLM/vertical_velocities_table_ngl14.txt') {
  NGL_table <- read_station_table(NGL_table_path)
  return(NGL_table)
}

read_JPL_table <- function(JPL_table_path = 'data/VLM/vertical_velocities_table_jpl14.txt') {
  JPL_table <- read_station_table(JPL_table_path)
  return(JPL_table)
}

read_GFZ_table <- function(GFZ_table_path = 'data/VLM/vertical_velocities_table_gt3.txt') {
  GFZ_table <- read_station_table(GFZ_table_path)
  return(GFZ_table)
}

.get_VLM <- function(lat, lon, method = 'default', choose_which = 1,
                     max_distance = 100000,
                     ULR_data_dir = 'data/VLM/ULR7a_neu/',
                     NGL_data_dir = 'data/VLM/NGL14/',
                     JPL_data_dir = 'data/VLM/JPL14/',
                     GFZ_data_dir = 'data/VLM/GT3_neu/',
                     return_VLM_only = T,
                     convert_daily_to_monthly = T,
                     choose_multiple = FALSE,
                     verbose = F,
                     return_VLM_station_info_only = F) {
  #' Get VLM data from the nearest station
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int. 'default', 'ULR', 'NGL', 'JPL', 'GFZ'. Alternatively, you can use 0, 1, 2, 3, 4.
  #' The `default` method uses the code from Guanchao's original code to keep the backward compatibility.
  #' @param choose_which integer, choose the nth nearest station
  #' @param max_distance numeric, the maximum distance to search for the nearest station.
  #' If there is no station within the max_distance, an empty data frame will be returned.
  #' @param return_VLM_only logical, if TRUE, only return VLM data, if FALSE, return more details including distance and speed in 3 directions.
  #' @param convert_daily_to_monthly logical, if TRUE, convert daily VLM data to monthly VLM data.
  #' @param choose_multiple logical, if TRUE, choose multiple stations from 1 to choose_which. Otherwise, choose the choose_whichth nearest station.
  #' @param verbose logical, if TRUE, print the info of the station.
  #' @param return_VLM_station_info_only logical, if TRUE, only return the info of the station.
  #' @return VLM data.
  #' @details Distance is calculated using the distHaversine function from the geosphere package. The distHaversine function uses the Haversine formula, which is a good general-purpose formula for distances on the Earth's surface. 
  #' The alternative distVincentySphere function is a bit more accurate for longer distances.
  #' @export
  #' @examples
  #' get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' get_VLM_by_station(station.names[3], method = 'ULR', choose_which = 5, choose_multiple = T, return_VLM_station_info_only = T, max_distance = 50)
  #' 
  if (method == 'default' | method == 0) {
    return(get_VLM_old(lat, lon, verbose))
  }
  else if (method == 'ULR' | method == 1) {
    station_table <- read_ULR_table()
    file_suffix <- paste0('_ULR7.neu')
    data_dir <- ULR_data_dir
  }
  else if (method == 'NGL' | method == 2) {
    station_table <- read_NGL_table()
    file_suffix <- paste0('_NGL14.neu')
    data_dir <- NGL_data_dir
  }
  else if (method == 'JPL' | method == 3) {
    station_table <- read_JPL_table()
    file_suffix <- paste0('_JPL14.neu')
    data_dir <- JPL_data_dir
  }
  else if (method == 'GFZ' | method == 4) {
    station_table <- read_GFZ_table()
    file_suffix <- paste0('_GT3.neu')
    data_dir <- GFZ_data_dir
  }
  else {
    stop(paste0('Invalid method. Please use "default", "ULR", "NGL", "JPL", "GFZ", or 0, 1, 2, 3, 4. You used: ', method, '.'))
  }
  
  # station_table$distance <- (station_table$Lat - lat)^2 + (station_table$Lon - lon)^2
  # Calculating the distance with a newer method. The result might be different. The unit is km.
  # The distHaversine function uses the Haversine formula, which is a good general-purpose formula for distances on the Earth's surface. The distVincentySphere function is a bit more accurate for longer distances.
  station_table$distance <- distHaversine(c(lon, lat), station_table[, c('Lon', 'Lat')]) / 1000
  station_table <- station_table[station_table$distance < max_distance, ] 
  empty_data_frame <- data.frame(month = character(), VLM = numeric())
  # Check if there is any station within the max_distance
  if (nrow(station_table) == 0) {
    # warning(paste0('There is no station within ', max_distance, ' km.'))
    return(empty_data_frame)
  }
  if (choose_which > nrow(station_table)) {
    warning(paste0('There are only ', nrow(station_table), ' stations within ', max_distance, ' km.'))
    if (choose_multiple) {
      choose_which <- nrow(station_table)
    }
    else {
      return(empty_data_frame)
    }
  }
  station_table <- station_table[order(station_table$distance), ]
  get_station_path <- function(choose_which) {
    station_name <- paste0('d', station_table$Site[choose_which], '_', station_table$DOMES[choose_which])
    station_path <- paste0(data_dir, station_name, file_suffix)
    return(station_path)
  }

  check_file_exist <- function(choose_which) {
    station_path <- get_station_path(choose_which)
    return(file.exists(station_path))
  }

  while((!check_file_exist(choose_which)) && (choose_which > 1)) {
    choose_which <- choose_which - 1
  }

  station_path <- get_station_path(choose_which)
  is_file_exist <- file.exists(station_path)
  while ((!is_file_exist) && (choose_which < nrow(station_table))) {
    choose_which <- choose_which + 1
    # station_name <- paste0('d', station_table$Site[choose_which], '_', station_table$DOMES[choose_which])
    # station_path <- paste0(data_dir, station_name, file_suffix)
    # is_file_exist <- file.exists(station_path)
    is_file_exist <- check_file_exist(choose_which)
    if (verbose) {
      print(paste0('The ', choose_which, 'th nearest station is not available. Try the ', choose_which + 1, 'th nearest station.'))
    }
  }


  station_path <- get_station_path(choose_which)

  if (verbose) {
    print(paste0('The requested location: ', lat, ', ', lon))
    print(paste0('The info of the nearest station: '))
    print(station_table[choose_which, ])
  }
  if (return_VLM_station_info_only) {
    if (choose_multiple) {
      station_table <- station_table[1:choose_which, ]
      station_table <- station_table %>% 
        mutate(site = paste(method, 1:choose_which, sep = '_'))
      return(station_table)
    }
    station_table <- station_table[choose_which, ]
    station_table$site <- paste(method, choose_which, sep = '_')
    return(station_table)
  }

  neu <- read_neu(station_path)

  if (return_VLM_only) {
    neu  <- neu %>% 
      rename(VLM = DU, month = Year) %>% 
      select(month, VLM)

    if (convert_daily_to_monthly | choose_multiple) {
      neu <- convert_daily_to_monthly(neu)
    }
  }

  if (choose_multiple) {
    # neu  <- neu %>% 
    #   rename(!!paste("VLM", choose_which, sep = "_") := DU, month = Year) %>% 
    #   select(month, !!paste("VLM", choose_which, sep = "_"))
    # neu <- convert_daily_to_monthly(neu)

    for (i in 1:(choose_which - 1)) {
      # station_name <- paste0('d', station_table$Site[i], '_', station_table$DOMES[i])
      # station_path <- paste0(data_dir, station_name, file_suffix)
      if (check_file_exist(i)) {
        station_path <- get_station_path(i)
        neu1 <- read_neu(station_path)
        neu1  <- neu1 %>% 
          rename(!!(paste("VLM", i, sep = "_")) := DU, month = Year) %>% 
          select(month, !!paste("VLM", i, sep = "_"))
        if (convert_daily_to_monthly) {
          neu1 <- convert_daily_to_monthly(neu1)
        }
        neu <- merge(neu, neu1, by = 'month', all = TRUE)
      }
    }
  }

  return(neu)
}

get_VLM_longest <- function(lat, lon, 
                method = 'default',
                choose_which = 3,
                convert_daily_to_monthly = TRUE,
                return_VLM_station_info_only = FALSE,
                ...) {
  #' Get VLM data with the longest time series from the nearby stations
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param choose_which integer, the maximum number of stations to choose from. When choose_which = 1, it is the same as get_VLM.
  #' @param convert_daily_to_monthly logical, it has to be TRUE for get_VLM_longest.
  #' @param ... other parameters passed to .get_VLM. Run `docstring(.get_VLM)` or `?.get_VLM` for details.
  #' @return VLM data
  #' @export 
  #' @examples
  #' get_VLM_longest(lat = 40, lon = 40, method = 'ULR', choose_which = 3, convert_daily_to_monthly = F)
  #' 
  #' 

  # convert_daily_to_monthly has to be TRUE for get_VLM_longest
  # ...$convert_daily_to_monthly <- T



  longest_station <- 0
  if (choose_which == 1) {
    return(.get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which, ...))
  }
  for (.choose_which in 1:choose_which) {
    vlm_data <- .get_VLM(lat, lon, method = method, choose_which = .choose_which, 
                        # convert_daily_to_monthly = T,
                        return_VLM_station_info_only = FALSE,
                        ...)
    # print(paste("nrow(vlm_data):", nrow(vlm_data)))
    if (nrow(vlm_data) > longest_station) {
      longest_station <- nrow(vlm_data)
      longest_station_data <- vlm_data
      if (return_VLM_station_info_only) {
        longest_station_data <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = .choose_which, 
                                        return_VLM_station_info_only = TRUE,
                                        ...)
        longest_station_data$site <- paste(method, .choose_which, sep = '_')
        if (.choose_which > 1) {

          closest_station <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = 1, 
                                    return_VLM_station_info_only = TRUE,
                                    ...)
          closest_station$site <- paste(method, 1, sep = '_')
          longest_station_data <- rbind(closest_station, longest_station_data)
        }
        
      }
    }
  }
  return(longest_station_data)
}

get_VLM_mean <- function(lat, lon,
            methods = c('ULR', 'NGL', 'JPL', 'GFZ'),
            choose_which = 1,
            return_VLM_only = TRUE,
            convert_daily_to_monthly = TRUE,
            return_VLM_station_info_only = FALSE,
            ...) {
  #' Get VLM data Using the four methods and take the mean.
  #' If one method is not available, impute it with the mean of the other methods.
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param methods character vector, the methods to use.
  #' @param choose_which integer, the maximum number of stations to choose from
  #' @param convert_daily_to_monthly logical, it has to be TRUE for get_VLM_mean. It doesn't really have any effect.
  #' @param ... other parameters passed to get_VLM
  #' @return VLM data
  #' @examples
  #' get_VLM_mean(lat = 40, lon = 40, methods = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #' get_VLM_mean(lat = 30, lon = 30, methods = c('ULR', 'NGL', 'JPL', 'GFZ'), return_VLM_station_info_only = T)
  #' 

  for (method in methods) {
    vlm_data <- .get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which, convert_daily_to_monthly = T, return_VLM_only = TRUE, return_VLM_station_info_only = return_VLM_station_info_only, ...)
    if (method == methods[1]) {
      if (return_VLM_station_info_only) {
        vlm_data_all <- vlm_data
        vlm_data_all$site <- paste(method, choose_which, sep = '_')
      }
      else {
        vlm_data_all <- vlm_data  %>% 
          rename(!!(method) := VLM)
      }
    }
    else {
      if (return_VLM_station_info_only) {
        vlm_data$site <- paste(method, 1, sep = '_')
        vlm_data_all <- rbind(vlm_data_all, vlm_data)
      }
      else {
        vlm_data_all <- merge(vlm_data_all, vlm_data, by = 'month', all = TRUE) %>% 
          rename(!!(method) := VLM)
      }
    }
  }

  if (return_VLM_station_info_only) {
    return(vlm_data_all)
  }

  vlm_data_all <- vlm_data_all %>%
    mutate_at(vars(methods), ~ifelse(is.na(.), mean(., na.rm = T), .))
  vlm_data_all <- vlm_data_all %>%
    mutate(VLM = rowMeans(.[methods], na.rm = T)) %>% 
    select(month, VLM)

  return(vlm_data_all)

}

get_VLM_circlemean <- function(lat, lon,
                               method = c('ULR', 'NGL', 'JPL', 'GFZ'),
                               max_distance = 1000,
                               choose_multiple = TRUE,
                               convert_daily_to_monthly = TRUE,
                               na.rm = TRUE,
                               choose_which = 100,
                               return_VLM_station_info_only = FALSE,
                               verbose = FALSE,
                               ...){
  #' Get VLM data Using the mean of the stations within a circle.
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method character or character vector, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param max_distance numeric, the maximum distance to search for the nearest station. The unit is km.
  #' @param ... other parameters passed to .get_VLM
  #' @return data frame with VLM data
  #' @examples
  #' get_VLM_circlemean(lat = 30, lon = 30, method = 'ULR', max_distance = 50)
  #' get_VLM_circlemean(lat = 30, lon = 30, method = c('ULR', 'NGL', 'JPL', 'GFZ'), max_distance = 50)
  #' 
  #' 
  datalist <- list()
  for (.method in method) {
    data <- .get_VLM(lat=lat, lon=lon, method = .method, max_distance = max_distance,
                     choose_which = choose_which, choose_multiple = TRUE, 
                     return_VLM_station_info_only = return_VLM_station_info_only, ...)
    if (return_VLM_station_info_only) {
      data <- data %>% 
        mutate(site = paste0(.method, '_', row_number()))
    }
    else {
      data <- data %>% 
        rename_with(~ paste(.method, .x, sep = '_'), -month)
    }
    # data <- data %>% 
    #   rename_with(~ paste(.method, .x, sep = '_'), -month)
    # datalist <- c(datalist, data)
    datalist[[length(datalist) + 1]] <- data
  }
  # if (verbose) {
  #   print(paste('The number of stations within', max_distance, 'km:', length(datalist)))
  # } # The code is wrong.
  if (return_VLM_station_info_only) {
    data <- Reduce(function(x, y) rbind(x, y), datalist)
    return(data)
  }
  data <- Reduce(function(x, y) merge(x, y, by = 'month', all = TRUE), datalist)
  data <- data %>% 
    # mutate(VLM = rowMeans(.[method], na.rm = na.rm)) %>%  
    rowwise() %>%
    mutate(VLM = mean(c_across(ends_with("VLM")), na.rm = na.rm)) %>%
    select(month, VLM)
  return(data)
}

get_VLM <- function(lat, lon,
                    method = 'default',
                    VLM_retreive_method = 'deterministic',
                    ...) {
  #' Get VLM data from the nearest station
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character, int or vector. 'default', 'ULR', 'NGL', 'JPL', 'GFZ'. Alternatively, you can use 0, 1, 2, 3, 4.
  #' When `VLM_retreive_method` is `mean` or `circlemean`, `method` can be a vector of methods that are used to calculate the mean.
  #' The `default` method uses the code from Guanchao's original code to keep the backward compatibility.
  #' @param VLM_retreive_method: character, 'deterministic', 'longest', 'mean' or 'circlemean'. 'deterministic' uses the exact specified station. 
  #' 'longest' uses the station with the longest time series constrained by `choose_which`. Alternatively, you can `get_VLM_longest()`
  #' 'mean' uses the mean of the VLM data from different stations constrained by `choose_which`. Alternatively, you can `get_VLM_mean()`.
  #' 'circlemean' uses the mean of the VLM data from different stations within a circle constrained by `max_distance`. Alternatively, you can `get_VLM_circlemean()`.
  #' @param choose_which integer, choose the nth nearest station. When `VLM_retreive_method` is `longest`, it is the maximum number of stations to choose from.
  #' @param return_VLM_only logical, if TRUE, only return VLM data, if FALSE, return more details including distance and speed in 3 directions.
  #' @param convert_daily_to_monthly logical, if TRUE, convert daily VLM data to monthly VLM data.
  #' @param verbose logical, if TRUE, print the info of the station.
  #' @param return_VLM_station_info_only logical, if TRUE, only return the info of the station.
  #' @param ... other parameters passed to `.get_VLM`. Run `docstring(.get_VLM)` or `?.get_VLM` for details.
  #' @return VLM data
  #' @export
  #' @examples
  #' get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' get_VLM(lat = 40, lon = 40, method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'mean')
  #' get_VLM_mean(lat = 40, lon = 40, methods = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #' get_VLM(lat = 40, lon = 40, method = 'NGL', VLM_retreive_method = 'longest', choose_which = 3, verbose = T)
  #' get_VLM_longest(lat = 40, lon = 40, method = 'NGL', choose_which = 3, verbose = T)
  #' get_VLM(lat = 30, lon = 30, method = 'ULR', VLM_retreive_method = 'circlemean', max_distance = 500)
  #' get_VLM(lat = 30, lon = 30, method = c('ULR', 'NGL', 'JPL', 'GFZ'), VLM_retreive_method = 'circlemean', max_distance = 500
  #' 
  if (VLM_retreive_method == 'deterministic') {
    return(.get_VLM(lat=lat, lon=lon, method = method, ...))
  }
  else if (VLM_retreive_method == 'longest') {
    return(get_VLM_longest(lat=lat, lon=lon, method = method, ...))
  }
  else if (VLM_retreive_method == 'mean') {
    return(get_VLM_mean(lat=lat, lon=lon, methods = method, ...))
  }
  else if (VLM_retreive_method == 'circlemean') {
    return(get_VLM_circlemean(lat=lat, lon=lon, method = method, ...))
  }
  else {
    stop(paste0('Invalid VLM_retreive_method. Please use "deterministic", "longest", or "mean". You used: ', VLM_retreive_method, '.'))
  }
}

station2lonlat <- function(station_name, data_dir=DATA_DIR){
  #' returns a list with the longitude and latitude of the station
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @export
  #' @return a list with the longitude and latitude of the station. The list has two elements: lon and lat.
  
  possible_station_names <- dir(data_dir)
  if(!(station_name %in% possible_station_names)){
    print("possible station names:")
    print(possible_station_names)
    stop("station name not found")
  }
  station_info <- read_csv(file.path(data_dir,station_name,"GMR.csv"), ,show_col_types = FALSE)
  return(list("lon"=station_info$lon,"lat"=station_info$lat))
}

get_VLM_by_station <- function(station_name, data_dir=DATA_DIR, ...){
  #' get VLM data by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data. Usually just use the default.
  #' @param ... other parameters passed to get_VLM. See `docstring(get_VLM)` or `?get_VLM` for details.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' get_VLM_by_station("astoria")
  #' get_VLM_by_station("astoria",method="JPL", choose_which=2)
  station_info <- station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  # print(paste("getting VLM for station",station_name,"at lon",lon,"lat",lat))
  return(get_VLM(lat=lat,lon=lon, ...))
}

vlm_plot <- function(lat, lon, method = 1, choose_which = 1) {
  #' plot VLM data based on the location
  #' @param lat latitude of the location
  #' @param lon longitude of the location
  #' @param method: character or int, 'default', 'ULR', 'NGL', 'JPL', 'GFZ', or 0, 1, 2, 3, 4.
  #' @param choose_which integer, choose the nth nearest station
  #' @return a ggplot object
  vlm_data <- get_VLM(lat=lat, lon=lon, method = method, choose_which = choose_which)
  p <- ggplot(vlm_data, aes(x = month, y = VLM)) +
    geom_line() +
    theme_bw() +
    labs(x = "Month", y = "VLM (m)")
  return(p)
}

vlm_plot_station <- function(station_name, data_dir=DATA_DIR, method = "ULR", ...) {
  #' get VLM data by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @param ... other parameters passed to get_VLM_by_station and vlm_source_plot.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' vlm_source_plot_by_station("astoria")
  #' vlm_source_plot_by_station("astoria",method="JPL", choose_which=2)
  station_info <- station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  p <- vlm_plot(lat, lon, method = method, ...)
  method_dict <- vlm_rename_dict("method")
  if (method %in% names(method_dict)) {
    method <- method_dict[[method]]
  }
  p <- p + labs(title = paste(station_name, 'VLM (m) from', method))
  return(p)
}

vlm_sources_plots <- function(lat, lon, methods = c(1, 2, 3, 4), choose_which = 1) {
  vlm_data <- get_VLM(lat=lat, lon=lon, method = methods[1], choose_which = choose_which)
  vlm_data <- vlm_data %>% 
    rename(!!as.character(methods[1]) := VLM)
  for (method in methods[-1]) {
    vlm_data <- merge(vlm_data, get_VLM(lat, lon, method = method, choose_which = choose_which), by = "month", all = TRUE)
    vlm_data <- vlm_data %>% 
      rename(!!as.character(method) := VLM)
  }
  method_dict <- c("1" = "ULR", "2" = "NGL", "3" = "JPL", "4" = "GFZ")
  rename_dict <- c("ULR" = "1", "NGL" = "2", "JPL" = "3", "GFZ" = "4")
  vlm_data <- vlm_data %>% 
    # rename(ULR = `1`, NGL = `2`, JPL = `3`, GFZ = `4`) %>% 
    rename(any_of(rename_dict)) %>% 
    select(month, any_of(c('ULR', 'NGL', 'JPL', 'GFZ')))
  # p  <- vlm_data %>% 
  #   gather(key = "source", value = "vlm", -month) %>% 
  #   ggplot(aes(x = month, y = vlm, color = source)) +
  #   geom_line() +
  #   theme_bw() +
  #   theme(legend.position = "bottom") +
  #   labs(x = "Month", y = "VLM (m)")
  p <- ggplot(vlm_data, aes(x = month))
  
  for (method in methods) {
    if (as.character(method) %in% names(method_dict)) {
      method <- method_dict[[method]]
    }
    p <- p + geom_line(data = na.omit(vlm_data[, c("month", method)]),
      aes(x = month, y = .data[[method]], color = !!method))
  }
  p <- p + theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Month", y = "VLM (m)")
  return(p)
}

vlm_sources_plots_by_station <- function(station_name, data_dir=DATA_DIR, ...) {
  #' Compare VLM data from different sources by station.
  #' @param station_name character, the name of the station
  #' @param data_dir character, the directory of the data
  #' @param ... other parameters passed to get_VLM_by_station and vlm_sources_plots.
  #' @return a data frame with the VLM data for the station
  #' @examples
  #' @export
  #' vlm_sources_plots_by_station("astoria")
  #' vlm_sources_plots_by_station("astoria",method="JPL", choose_which=2)
  station_info <- station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  p <- vlm_sources_plots(lat=lat, lon=lon, ...)
  return(p)
}

.station2color <- function(station) {
  # if station starts with ULR
  if (startsWith(station, 'ULR')) {
    return("#1b9e77")
  }
  else if (startsWith(station, 'NGL')) {
    return("#d95f02")
  }
  else if (startsWith(station, 'JPL')) {
    return("#7570b3")
  }
  else if (startsWith(station, 'GFZ')) {
    return("#e7298a")
  }
  else {
    return("#000000")
  }
}

station2color <- function(station) {
  #' get the color of the station
  #' @param station vector of characters, the name of the station
  #' @return a vector of colors
  #' 
  #' 
  if (is.character(station)) {
    return(sapply(station, .station2color))
  }
  else {
    return(.station2color(station))
  }
}

get_VLM_coords <- function(lat,lon, methods = c('ULR','NGL','JPL','GFZ'), ...){
  coords_df = data.frame(site = NULL,lat = NULL,lon=NULL,distance=NULL)
  i = 1
  for (meth in methods){
    station_table = get_VLM(lat=lat,lon=lon,method=meth,return_VLM_station_info_only=TRUE, ...)
    if ('site' %in% names(station_table)){
      # tempdf = data.frame('site'=station_table$site,lat=station_table$Lat,lon=station_table$Lon,distance = (station_table$Lat - lat)^2 + (station_table$Lon - lon)^2)
      tempdf = data.frame('site'=station_table$site,lat=station_table$Lat,lon=station_table$Lon,distance = station_table$distance)
    } else {
      # tempdf = data.frame('site'=meth,lat=station_table$Lat,lon=station_table$Lon,distance = (station_table$Lat - lat)^2 + (station_table$Lon - lon)^2)
      
    }
    coords_df = rbind(coords_df,tempdf)
    i=i+1
  }
  coords_df <- coords_df  %>% 
    mutate(color = station2color(site))

  return(coords_df)
}

get_loc_info<- function(station_name,data_dir=DATA_DIR,save=TRUE, ...){
  #' get VLM station info by station.
  #' @param station_name character, the name of the station
  #' @param save logical, if TRUE, save the data to a csv file
  #' @param ... other parameters passed to get_VLM_by_station and get_VLM_coords.
  #' @return a data frame with the VLM station information for the station. An example:
  #' > get_VLM_coords(3,3)
  #'    site    lat    lon distance   color
  #' 1 ULR_1 6.3847 2.4500 381.6908 #1b9e77
  #' 2 NGL_1 6.3847 2.4500 381.6908 #d95f02
  #' 3 JPL_1 0.3539 9.6721 798.6556 #7570b3
  #' 4 GFZ_1 6.3847 2.4500 381.6908 #e7298a
  #' @examples
  #' get_VLM_coords_by_station(station.names[6], method=1,choose_which = 11, VLM_retreive_method = "longest")
  station_info = station2lonlat(station_name, data_dir=data_dir)
  lon <- station_info$lon
  lat <- station_info$lat
  df = data.frame(site = station_name,lat = lat,lon=lon,distance=0, color="#000000")
  df = rbind(df,get_VLM_coords(lat=lat,lon=lon, ...))
  
  if (save){
    write.csv(df,file.path(data_dir,station_name,"VLM_distances.csv"))
  }
  
  return(df)
}

get_VLM_coords_by_station <- function(...) get_loc_info(...)

get_snapshot<-function(df,centlat,centlon){
  #' get a snapshot of the VLM stations around a location
  #' @param df data frame, the output of `get_VLM_coords` or `get_loc_info`.
  #' @param centlat numeric, the latitude of the center of the snapshot
  #' @param centlon numeric, the longitude of the center of the snapshot
  #' @return a data frame with the boundaries of the snapshot
  #' @examples
  #' df = get_VLM_coords(30,30)
  #' get_snapshot(df,30,30)
  radius_lat = 1.1*max(abs(df$lat-centlat))
  radius_lon = 1.1*max(abs(df$lon-centlon))
  r = max(radius_lon,radius_lat)
  
  df = data.frame(uplat = centlat+r,lowlat = centlat-r,uplon = centlon+r,lowlon=centlon-r)
  df
}

get_snapshot_station <- function(station_name, data_dir=DATA_DIR){
  #' get VLM station info by station.
  #' @param station_name character, the name of the station
  df <- get_loc_info(station_name, data_dir=data_dir)
  centlat <- df$lat[1]
  centlon <- df$lon[1]
  get_snapshot(df,centlat,centlon)
}

vlm_map_df <- function(coords){
  #' plots a map of the VLM stations around a location.
  #' @param coords data frame, usually the output of `get_VLM_coords` or `get_loc_info`.
  #' @return a ggplot object
  #' @examples
  #' get_VLM_coords_by_station(station.names[6], method=1,choose_which = 11, VLM_retreive_method = "longest") 
  #'    %>% vlm_map_df()
  #' get_VLM_coords_by_station(station.names[6],choose_which = 11, VLM_retreive_method = "mean") %>% vlm_map_df()
  #initialize map
  us_map = map_data("state")
  
  #get snapshot
  snap = get_snapshot(coords,coords$lat[1],coords$lon[1])
  print(coords$color)
  #Add map to base plot
  base_map <- ggplot() + geom_polygon(data=us_map, aes(x=long, y=lat, group=group), 
                                      color="black", fill="white") +
    coord_cartesian(xlim = c(snap$lowlon,snap$uplon),ylim = c(snap$lowlat,snap$uplat))
  print(unique(coords$color))
  map_w_points <- base_map + geom_point(data=coords,aes(x=lon, y=lat, color = color),show.legend = FALSE) + scale_color_manual(values = unique(coords$color))
  map_w_points+ geom_text_repel(data = coords, max.overlaps = Inf, aes(x = lon, y = lat, label = site),
                                                                                                 box.padding = 0.25, point.padding = 0.25) 
}

vlm_map <- function(station_name, data_dir = DATA_DIR,savedistance=FALSE,saveplot=FALSE){
  coords = get_loc_info(station_name,data_dir = data_dir,save=savedistance)

  map_w_points <- vlm_map_df(coords)

  if(saveplot){
    ggsave(file.path(data_dir,station_name,"VLM_Locations.jpg"))
  }
  
  map_w_points
  
}

vlm_map_station <- function(...) vlm_map(...)

merge_VLM_df <- function(datasets, col_names = c('ULR', 'NGL', 'JPL', 'GFZ')) {
  #' Merge monthly VLM data
  #' @param datasets a list of data frames with monthly VLM data. The datasets should come from `get_VLM` or `get_VLM_by_station`.
  #' The column names should be `month` and `VLM`.
  #' @param col_names character vector, the col_names of the column names. It should have the same length as `datasets`.
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' monthly_VLM_1 <- get_VLM(lat = 40, lon = 40, method = 'ULR')
  #' monthly_VLM_2 <- get_VLM(lat = 40, lon = 40, method = 'NGL')
  #' merge_VLM(list(monthly_VLM_1, monthly_VLM_2), col_names = c('ULR', 'NGL', 'JPL', 'GFZ'))
  #' 
  #' 
  for (i in 1:length(datasets)) {
    datasets[[i]] <- datasets[[i]] %>% 
      rename(!!(col_names[i]) := VLM) %>% 
      mutate(month = as.Date(month))
  }
  Reduce(function(x, y) merge(x, y, by = 'month', all = TRUE), datasets)
}

merge_VLM_station <- function(station, data_dir = DATA_DIR, methods = c('ULR', 'NGL', 'JPL', 'GFZ')) {
  #' Merge monthly VLM data by station
  #' @param station character, the name of the station
  #' @param data_dir character, the directory of the data. Usually just use the default.
  #' @param methods character vector, the methods to use. Default is c('ULR', 'NGL', 'JPL', 'GFZ').
  #' @return a data frame with merged monthly VLM data
  #' @export
  #' @examples
  #' merge_VLM_station(station.names[6])
  #' 
  #' 
  datasets <- list()
  for (method in methods) {
    datasets[[method]] <- get_VLM_by_station(station, data_dir = data_dir, method = method)
  }
  merge_VLM_df(datasets, col_names = methods)
}

if (sys.nframe() == 0){
  # test codes
  test_data <- get_VLM(lan_input=40,lon_input=40)
  test_data
}
