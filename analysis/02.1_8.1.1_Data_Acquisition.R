# load libraries
library(tidyverse)
library(geodata)

# download SRTM data
srtm_path <- file.path("data-raw", "elevation", "srtm_38_03.tif")
if (!file.exists(srtm_path)) {
  geodata::elevation_3s(lat = 46.6756, lon = 7.85480, path = "data-raw")
}
dem <- terra::rast(srtm_path)

# load libraries
library(MODISTools)

# download and save phenology data
phenology_path <- file.path("data-raw", "phenology.rds")
if (!file.exists(phenology_path)) {
  phenology <- MODISTools::mt_subset(
    product = "MCD12Q2", lat = 46.6756, lon = 7.85480,
    band = "Greenup.Num_Modes_01", start = "2012-01-01", end = "2012-12-31",
    km_lr = 100, km_ab = 100, site_name = "swiss", internal = TRUE, progress = FALSE
  )
  saveRDS(phenology, phenology_path)
} else {
  phenology <- readRDS(phenology_path)
}

# download and save land cover data
land_cover_path <- file.path("data-raw", "land_cover.rds")
if (!file.exists(land_cover_path)) {
  land_cover <- MODISTools::mt_subset(
    product = "MCD12Q1", lat = 46.6756, lon = 7.85480,
    band = "LC_Type1", start = "2012-01-01", end = "2012-12-31",
    km_lr = 100, km_ab = 100, site_name = "swiss", internal = TRUE, progress = FALSE
  )
  saveRDS(land_cover, land_cover_path)
} else {
  land_cover <- readRDS(land_cover_path)
}