library(MODISTools)
lat <- 43.5
lon <- -74.5
years <- 2001:2010

## MODIS Phenology-Daten laden
phenology <- mt_subset(
  product = "MCD12Q2",
  lat = lat,
  lon = lon,
  band = c("Greenup.Num_Modes_01","Maturity.Num_Modes_01"),
  start = "2001-01-01",
  end   = "2010-12-31",
  km_lr = 100, km_ab = 100,
  site_name = "adirondacks",
  internal = TRUE, progress = FALSE
)
write.csv(phenology, "data-raw/phenology.csv", row.names = FALSE)

## Land Cover-Daten laden
land_cover <- mt_subset(
  product = "MCD12Q1",
  lat = lat, lon = lon,
  band = "LC_Type1",
  start = "2010-01-01", end = "2010-12-31",
  km_lr = 100, km_ab = 100,
  site_name = "adirondacks_lc",
  internal = TRUE
)
write.csv(land_cover, "data-raw/land_cover.csv", row.names = FALSE)

dem_us <- geodata::elevation_30s(country="USA", path="data-raw/")
