# libraries ----
if (!require(librarian)){
  install.packages("librarian")
  library(librarian)
}
librarian::shelf(
  # time-series
  caTools, tools, dygraphs, xts,
  #spatial
  sf, leaflet,
  # tidyverse
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

# combine data files for MARINe sites
dataCombiner <- function(data_site_zone) {
  # message("reading in sites")
  
  # store site and zone names
  if ("site" %in% colnames(data_site_zone)) {
    site <- unique(data_site_zone$site) 
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
      site = if("site" %in% colnames(data_site_zone)) site else NA,
      zone = if("zone" %in% colnames(data_site_zone)) zone else NA) %>% 
    relocate(site)
}

dailyQuantilesData <- function(data) {
  data %>% 
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
    select(-time, -Temp_C) %>% 
    gather("metric", "Temp_C", c(-1, -2, -3)) %>% 
    select(-zone, zone)
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
  