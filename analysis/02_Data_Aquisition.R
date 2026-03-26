# load libraries
library(tidyverse)
library(geodata)

# download SRTM data
# This stores file srtm_38_03.tif in 
# subfolder elevation of tempdir()
geodata::elevation_3s(
  lat = 46.6756,
  lon = 7.85480,
  path = tempdir()
)

# read the downloaded data
# use file.path() to combine
# a directory path with a filename
dem <- terra::rast(
  file.path(
    tempdir(),
    "elevation",
    "srtm_38_03.tif"
  )
)
# load libraries
library(MODISTools)

# download and save phenology data
phenology <- MODISTools::mt_subset(
  product = "MCD12Q2",
  lat = 46.6756,
  lon = 7.85480,
  band = "Greenup.Num_Modes_01",
  start = "2012-01-01",
  end = "2012-12-31",
  km_lr = 100,
  km_ab = 100,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)
# screening of data
phenology <- phenology |>
  mutate(
    value = ifelse(value > 32656, NA, value),
    value = as.numeric(format(as.Date("1970-01-01") + value, "%j")),
    value = ifelse (value < 200, value, NA)
  )
phenology_raster <- MODISTools::mt_to_terra(
  phenology,
  reproject = TRUE
)
# crop the dem
dem <- terra::crop(
  x = dem,
  y = phenology_raster
)
# resample the dem using
# the mean DEM value in a
# MODIS pixel
dem <- terra::resample(
  x = dem,
  y = phenology_raster,
  method = "average"
)

# mask the locations which
# have no data
dem <- terra::mask(
  dem,
  is.na(phenology_raster),
  maskvalues = TRUE
)
# download and save land cover data
land_cover <- MODISTools::mt_subset(
  product = "MCD12Q1",
  lat = 46.6756,
  lon = 7.85480,
  band = "LC_Type1",
  start = "2012-01-01",
  end = "2012-12-31",
  km_lr = 100,
  km_ab = 100,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)
land_cover_raster <- MODISTools::mt_to_terra(
  land_cover,
  reproject = TRUE
)