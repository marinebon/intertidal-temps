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
  fs, glue, here, lubridate, stringr, tidyverse, purrr)



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
  zone_colors <- setNames(pal, zones) 
  
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
  
  
# paths & variables ----
user <- Sys.info()[["user"]]

# set dir_gdata as filepath for robomussels data on Google Drive 
dir_gdata <- case_when(
  #user == "bbest"       ~ "/Users/bbest/Downloads/robomusseldata20201030",
  user == "bbest"       ~ "/Volumes/GoogleDrive/My Drive/projects/mbon-p2p/data/rocky/MARINe/robomusseldata20201030",
  user == "cdobbelaere" ~ "/Users/cdobbelaere/Documents/robomussels/robomusseldata20201030")
dir_avg <- file.path(dirname(dir_gdata), "robomusseldata20201030_avg")
stopifnot(any(dir.exists(c(dir_gdata, dir_avg))))



# sites ----

## read in site data ----

# get [Pole to Pole](https://marinebon.org/p2p) sites & convert to sf 
d_psites <- read_csv(here::here("data/p2p_sites.csv"))

d_psites_sf <- d_psites %>% 
  st_as_sf(
    coords = c("lon", "lat"), remove = F,
    crs    = 4326) # geographic coordinate ref system


# get MARINe sites, prep for combining temp data txts, & convert to sf
# from robo metadata
d_msites <- read_csv("Robomussel metadata.csv") %>% 
  rename(
    microsite_id   = `microsite id`,
    logger_type    = `logger type`,
    tidal_height_m = `tidal height (m)`,
    wave_exposure  = `wave exposure`,
    start_date     = `start date`,
    end_date       = `end date`)

d_mfiles  <- tibble(
  path            = list.files(dir_gdata, ".*\\.txt", full.names = T),
  file            = basename(path),
  microsite_year  = file %>% path_ext_remove()) %>% 
  separate(microsite_year, c("msite", "year"), "_", convert = T)

# join file names with their associated metadata
d_msites <- d_mfiles %>% 
  left_join(
    d_msites, by = c("msite" = "microsite_id")) %>% 
  drop_na() # sum(is.na(d_msites$site)): n = 2

# convert to sf for joining with p2p
d_msites_sf <- d_msites %>% 
  st_as_sf(
    coords = c("longitude", "latitude"), remove = F,
    crs    = 4326) # geographic coordinate ref system (WGS1984)
 

## join MARINe and p2p sites by nearest feature ----

# join 4 nearby sites by nearest feature
# (sites that are both MARINe and p2p that we have data for)
sites_joined <- st_join(
  d_msites_sf, 
  d_psites_sf %>% 
    select(-country) %>% 
    rename(pgeometry = geometry), 
  join = st_nearest_feature)

# add MARINe & p2p geom columns,
# then use to calc distance between sites in km
sites_joined <- sites_joined %>% 
  st_drop_geometry() %>% 
  mutate(
    geom_marine = map2(longitude, latitude, xy2pt),
    geom_p2p    = map2(lon      , lat     , xy2pt),
    dist_km     = map2_dbl(geom_marine, geom_p2p, st_distance) / 1000) %>% 
  st_as_sf( 
    coords = c("lon", "lat"), remove = F,
    crs    = 4326)



## smooth MARINe sites ----

# split by each unique site & zone combination

d_filtered <- split(sites_joined, list(sites_joined$site, sites_joined$zone), drop = T, sep = "_")

# loop through all site & zone combinations and smooth 
# store dailyq for each in a list
dailyq <- list()

for (i in 1:length(d_filtered)){ # i = 
  
  # combine data for all zones of each site
  d_site_zone <- dataCombiner(data_site_zone = d_filtered[[i]])
  
  # get site and zone names
  site <- unique(d_filtered[[i]][["site"]])
  if ("zone" %in% colnames(d_filtered[[i]])) {
    zone <- unique(d_filtered[[i]][["zone"]])
    site_zone <-  paste0(site, "_", zone)
  } else {
    zone <- ""
    site_zone <- site
  }
  
  # smooth data
  d_dailyq <- dailyQuantilesData(d_site_zone)
  
  # assign global names to local objects
  stringname <- paste0(site_zone, "_dailyq")
  
  # populate dailyq list
  dailyq[[stringname]] <- d_dailyq
  
}

# create vector of ordered zones for when we convert zones to factor
zones <- c("Low", "Lower-Mid", "Mid", "Upper-Mid", "High")

# bind temps for all sites & zones
d <- dailyq %>% 
  bind_rows %>% 
  mutate(
    site = as.factor(site),
    zone = factor(zone, levels = zones, ordered = T))
 

# write smoothed avg temp data for each site to a unique csv
for (i in 1:length(levels(d$site))) { # for each site i
  
  site <- levels(d$site)[i]
  
  d_site <- d %>% 
    filter(site == !!site) %>%
    ungroup()
  
  # Filter out avgs
  d_site_avg <- d_site %>% 
    filter(metric == "temp_c_avg") %>% 
    select(-site, -metric) %>% 
    mutate(
      zone = factor(zone, zones, ordered = T)) %>% 
    arrange(zone, day) %>% 
    pivot_wider(day, names_from = zone, values_from = Temp_C) %>% 
    arrange(day) 
  
  write_csv(
    d_site_avg,
    paste0(getwd(), "/data_smoothed/", site, ".csv"))
  
}



## organize report output ----
# get start & end dates and other metadata
sites_joined_summarized <- sites_joined %>% 
  mutate(num_sites = length(unique(site))) %>% 
  group_by(site) %>% 
  mutate(num_loggers = length(unique(msite))) %>% 
  # mutate(
  #   start_date = min(start_date),
  #   end_date   = max(end_date)) %>% 
  # group_by(site, location) %>% 
  distinct(site, .keep_all = T) %>% 
  select(site, location, region, country, num_loggers, num_sites) %>% 
  rename(site_id = site, site_name = location)

sites_smoothed <- 
  tibble(
    path        = list.files(glue("{getwd()}/data_smoothed"), full.names = T),
    file        = basename(path),
    site_id     = file %>% path_ext_remove()) %>% 
    # ,xts        = glue("x_{site_id}"),) 
  left_join(
    sites_joined_summarized,
    by = c("site_id" = "site_id"), copy = F, keep = F) %>% 
  select(site_id, site_name, everything())

