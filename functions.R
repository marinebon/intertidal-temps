# libraries ----
if (!require(librarian)){
  install.packages("librarian")
  library(librarian)
}
librarian::shelf(
  # time-series
  caTools, tools, dygraphs, xts,
  # spatial
  sf, leaflet,
  # tidyverse / other
  fs, glue, here, lubridate, stringr, tidyverse, purrr, yaml)

# paths & variables ----
user <- Sys.info()[["user"]]

# set dir_gdata as filepath for robomussels data on Google Drive 
dir_gdata <- case_when(
  #user == "bbest"       ~ "/Users/bbest/Downloads/robomusseldata20201030",
  user == "bbest"       ~ "/Volumes/GoogleDrive/My Drive/projects/mbon-p2p/data/rocky/MARINe/robomusseldata20201030",
  user == "cdobbelaere" ~ "/Users/cdobbelaere/Documents/robomussels/robomusseldata20201030")
dir_avg <- file.path(dirname(dir_gdata), "robomusseldata20201030_avg")
stopifnot(any(dir.exists(c(dir_gdata, dir_avg))))


# define functions ----

# convert sf to st_point; add sf geometry list column & coord ref sys
xy2pt <- function(x, y){
  st_point(c(x, y)) %>% 
    st_sfc(crs = 4326)
}


# find where metadata ends for gdata files
find_skip <- function(file, pattern, n = 20) { 
  min(grep(pattern, read_lines(file, n_max = n)))
}

# function to read temp csv files for gdata
read_tempcsv <- function(path) {
  data <- tibble(
    read_csv(
      path,
      skip = (find_skip(file = path, pattern = "^time,") - 1))) %>% 
    drop_na() 
  if ("temp" %in% names(data)) {
    data <- data %>% rename(Temp_C = temp)
  }
  if (TRUE %in% grepl("-", data$time)) {
    data <- data %>% mutate(time = parse_date_time(time,"y-m-d H:M:S"))
  } else if (TRUE %in% grepl("/", data$time)) {
    data <- data %>% mutate(time = parse_date_time(time, "m/d/y H:M"))
  }
  data <- data %>% 
    mutate(path = path)
  
  if ("site" %in% colnames(data)) {
    data$site <- data$site
  } else if ("sensor" %in% colnames(data)) {
    data$site <- data$sensor
  } else data$site <- as.character(NA)

  if (TRUE %in% grepl("sun", data$path)) {
    data$zone <- "sun"
  } else if (TRUE %in% grepl("shade", data$path)) {
    data$zone <- "shade"
  } else if (TRUE %in% grepl("exposed", data$path)) {
    data$zone <- "exposed"
  } else data$zone <- NA
  
  if ("sensor" %in% colnames(data)) {data <- data %>% select(-sensor)}
  return(data)
}


# read and clean metadata for gdata
read_metadata <- function(path) { 
  
  n_start <- find_skip(file = path, pattern = "^time,")
  
  # if metadata present:
  if (n_start != Inf) {
    
    metadata <- tibble(
      raw_data  = read_lines(
        path,
        n_max = n_start - 1)) 
    
    metadata <- metadata %>% 
      mutate_if(
        is.character, 
        function(x) {Encoding(x) <- "latin1"; return(x)}) %>% 
      filter(str_detect(raw_data, ":")) %>% 
      separate(raw_data, c("key", "value"), ": ", convert = T) %>% 
      pivot_wider(names_from = "key", values_from = "value")
    
    if ("custom name" %in% names(metadata)) {
      metadata <- metadata %>% select(-"custom name")
    }
    
    if ("coords" %in% names(metadata)) {
      if (TRUE %in% grepl(",", metadata$coords)) {
        metadata <- metadata %>%
          separate(coords, c("lat", "lon"), ", ")
      } else {
        metadata <- metadata %>%
          separate(coords, c("lat", "lon"), " ")
      }
    }
    
    # fix coords
    if (TRUE %in% grepl("S", metadata$lat)) {
      metadata$lat <- gsub("S", "", metadata$lat)
      metadata <- metadata %>% mutate(lat = glue("-{metadata$lat}"))
    } else if (TRUE %in% grepl("N", metadata$lat)) {
      metadata$lat <- metadata$lat %>%
        gsub("N", "", .) 
    }
    
    if (TRUE %in% grepl("E", metadata$lon)) {
      metadata$lon <- gsub("E", "", metadata$lon)
      metadata <- metadata %>% mutate(lon = glue("-{metadata$lon}"))
    } else if (TRUE %in% grepl("W", metadata$lon)) {
      metadata$lon <- metadata$lon %>% 
        gsub("W", "", .) 
    }
    
    metadata <- metadata %>%
      mutate_all(list(~gsub(",", "", .))) %>% 
      select(-"accuracy (m)")
    
    if ("X1" %in% colnames(metadata)) {
      metadata <- metadata %>% select(-X1)
    }
    
    if (TRUE %in% grepl("ample", names(metadata))) {
      metadata <- metadata %>%
        rename(number_of_samples = grep("ample", names(metadata), value = T))
    }
    
    # remove non-numeric characters
    metadata <- metadata %>% mutate(across(!`serial number`, as.numeric)) 
    
    # convert back to character for filling NAs; add path & file names
    metadata <- metadata %>% 
      mutate(across(everything(), as.character)) %>% 
      mutate(path = as.character(path)) 
  }
  
  if (TRUE %in% grepl("sun", metadata$path)) {
    metadata$zone <- "sun"
  } else if (TRUE %in% grepl("shade", metadata$path)) {
    metadata$zone <- "shade"
  } else if (TRUE %in% grepl("exposed", metadata$path)) {
    metadata$zone <- "exposed"
  } else metadata$zone <- metadata$`serial number`
  
  return(metadata)
}


# combine data files for MARINe sites
dataCombiner <- function(data_site_zone) {
  # message("reading in sites")
  
  # site and zone names
  if ("id" %in% colnames(data_site_zone)) {
    site <- unique(data_site_zone$id) 
  }
  if ("zone" %in% colnames(data_site_zone)) {
    zone <- unique(data_site_zone$zone)
  }
  
  # read temp file corresponding to each path name 
  if (file_ext(data_site_zone$path) == "csv") {
    temp_data <- bind_rows(lapply(data_site_zone$path, read_csv, col_types = cols())) %>% 
      rename(Temp_C = sst) %>% 
      mutate(time   = parse_date_time(date, "y-m-d")) %>% 
      select(-date)
    
  } else {
    temp_data <- bind_rows(lapply(data_site_zone$path, read_tsv)) %>% 
      mutate(time = parse_date_time(Time_GMT, "m/d/y H:M")) %>%
      select(-Time_GMT)
  }
  
  temp_data <- temp_data %>% 
    drop_na() %>% 
    group_by(time) %>% 
    summarize(Temp_C = mean(Temp_C)) %>% 
    mutate(
      site = if("id" %in% colnames(data_site_zone)) site else NA,
      zone = if("zone" %in% colnames(data_site_zone)) zone else NA) %>% 
    relocate(site)
  
  return(temp_data)
}

dailyQuantilesData <- function(data) {
  data <- data %>% 
    mutate(
      day = floor_date(time, unit = "day")) %>%
    group_by(day) %>%
    distinct(day, .keep_all = T) %>% 
    mutate(
      temp_c_q10 = quantile(Temp_C, 0.1),
      temp_c_q90 = quantile(Temp_C, 0.9),
      temp_c_avg = mean(Temp_C),
      temp_c_min = min(Temp_C),
      temp_c_max = max(Temp_C)) %>% 
    select(-time, -Temp_C) 
  if (!("path" %in% colnames(data))) {
    data <- data %>% 
      gather("metric", "Temp_C", c(-1, -2, -3)) %>% 
      select(-zone, zone)
  } else if ("path" %in% colnames(data)) {
    data <- data %>% 
      gather("metric", "Temp_C", c(-1, -2, -3, -4)) %>% 
      select(-zone, zone)
  }
    
  return(data)
}

dailyQuantiles_gdata <- function(data) {
  data <- data %>% 
    group_by(zone) %>% 
    mutate(Temp_C = mean(Temp_C)) %>% 
    ungroup() %>% 
    mutate(
      day = floor_date(time, unit = "day")) %>%
    group_by(day) %>%
    distinct(day, .keep_all = T) %>% 
    mutate(
      temp_c_q10 = quantile(Temp_C, 0.1),
      temp_c_q90 = quantile(Temp_C, 0.9),
      temp_c_avg = mean(Temp_C),
      temp_c_min = min(Temp_C),
      temp_c_max = max(Temp_C)) %>% 
    select(-time, -Temp_C) 
  if (!("path" %in% colnames(data))) {
    data <- data %>% 
      gather("metric", "Temp_C", c(-1, -2, -3)) %>% 
      select(-zone, zone)
  } else if ("path" %in% colnames(data)) {
    data <- data %>% 
      gather("metric", "Temp_C", c(-1, -2, -3, -4)) %>% 
      select(-zone, zone)
  }
  
  return(data)
}




# read in smoothed csv and convert to xts for dygraphs
get_xts <- function(path) {
    d_smoothed <- read_csv(path) %>%
      mutate(
        day = parse_date_time(day, "ymd")) %>%
      mutate_at(vars(-1), funs(as.numeric))
    x_smoothed <- xts(select(d_smoothed, -day), order.by = d_smoothed$day)
    return(x_smoothed)
}
  
get_dygraph <- function(xts) {
  
  # map color palette to zone names
  pal  <- c("#3D2C9A", "#3E98C5", "#4A9A78", "#F7BD33", "#D74B00")
  zone_colors <- setNames(pal, names(xts)) 
  
  # plot
  dygraph <- dygraph(xts, main = "Daily Temperature") %>%
    dyHighlight(
      highlightCircleSize = 5, 
      highlightSeriesBackgroundAlpha = 0.2,
      hideOnMouseOut = TRUE) %>%
    dyOptions(
      # use only colors corresponding to zones that 
      # exist in the site's xts data
      colors = as.character(
        zone_colors[names(xts)]),
      connectSeparatedPoints = FALSE) %>% 
    dyOptions(
      fillGraph = FALSE, fillAlpha = 0.4) %>%
    dyRangeSelector()
  
  return(dygraph)
}
  