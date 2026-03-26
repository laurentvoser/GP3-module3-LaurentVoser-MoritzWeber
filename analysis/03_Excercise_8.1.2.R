library(MODISTools)
library(terra)
library(tidyverse)

# daten laden
bands <- mt_bands("MCD12Q2")
print(bands)

lat <- 43.5
lon <- -74.5

phenology <- mt_subset(
  product = "MCD12Q2",
  lat = lat,
  lon = lon,
  band = c(
    "Greenup.Num_Modes_01",
    "Maturity.Num_Modes_01"
  ),
  start = "2001-01-01",
  end = "2010-12-31",
  km_lr = 100,
  km_ab = 100,
  site_name = "adirondacks",
  internal = TRUE,
  progress = FALSE
)

# Daten aufbereiten & Raster erstellen
phenology <- phenology |>
  mutate(
    value = ifelse(value > 30000, NA, value),
    doy = as.numeric(format(as.Date("1970-01-01") + value, "%j")),
    year = as.numeric(substr(calendar_date, 1, 4))
  )
greenup <- phenology |>
  filter(band == "Greenup.Num_Modes_01")

greenup_raster <- mt_to_terra(greenup, reproject = TRUE)

maturity <- phenology |>
  filter(band == "Maturity.Num_Modes_01")

maturity_raster <- mt_to_terra(maturity, reproject = TRUE)

# Landcover filtern (2010)
land_cover <- mt_subset(
  product = "MCD12Q1",
  lat = lat,
  lon = lon,
  band = "LC_Type1",
  start = "2010-01-01",
  end = "2010-12-31",
  km_lr = 100,
  km_ab = 100,
  site_name = "adirondacks_lc",
  internal = TRUE
)

lc_raster <- mt_to_terra(land_cover, reproject = TRUE)

# nur Broadleaf + Mixed Forest (IGBP 4+5)
mask_forest <- lc_raster %in% c(4,5)

greenup_raster <- mask(greenup_raster, mask_forest, maskvalues = FALSE)
maturity_raster <- mask(maturity_raster, mask_forest, maskvalues = FALSE)


# LTM und STDEV 2001 - 2009
# Filter für Jahre 2001-2009
greenup_hist <- greenup_raster[[1:9]]  # Layer 1 bis 9 entsprechen 2001-2009
maturity_hist <- maturity_raster[[1:9]]

# LTM und SD berechnen
greenup_ltm <- mean(greenup_hist, na.rm = TRUE)
greenup_sd <- terra::app(greenup_hist, sd, na.rm = TRUE)

maturity_ltm <- mean(maturity_hist, na.rm = TRUE)
maturity_sd <- terra::app(maturity_hist, sd, na.rm = TRUE)


# Abweichung 2010
# Layer für 2010
greenup_2010 <- greenup_raster[[10]]
maturity_2010 <- maturity_raster[[10]]

# Frühgrün (< LTM - 1 SD)
early_greenup <- greenup_2010 < (greenup_ltm - greenup_sd)

# Spätmaturity (> LTM + 1 SD)
late_maturity <- maturity_2010 > (maturity_ltm + maturity_sd)

# Plot zur Kontrolle
par(mfrow=c(1,2))
terra::plot(early_greenup, main="Early Greenup 2010")
terra::plot(late_maturity, main="Late Maturity 2010")

# Digital Elevation
library(geodata)
library(terra)

# DEM der USA (30 arcsec) herunterladen
dem_us <- geodata::elevation_30s(country = "USA", path = tempdir())

# DEM auf dasselbe CRS wie die Greenup-Raster reprojizieren
dem_us <- terra::project(dem_us, crs(greenup_raster))

# DEM auf dasselbe Grid wie Greenup 2010 bringen (Resampling)
dem_resampled <- terra::resample(dem_us, greenup_raster[[10]], method = "bilinear")

# Maskieren der relevanten Pixel
dem_early <- terra::mask(dem_resampled, early_greenup, maskvalues = FALSE)
dem_late  <- terra::mask(dem_resampled, late_maturity, maskvalues = FALSE)

# Werte extrahieren & zu numerischen Vektoren machen
vals_early <- as.numeric(terra::values(dem_early, na.rm = TRUE))
vals_late  <- as.numeric(terra::values(dem_late, na.rm = TRUE))
vals_rest  <- as.numeric(terra::values(dem_resampled, na.rm = TRUE))

# Nur Pixel, die nicht in early oder late sind
vals_rest  <- vals_rest[!(vals_rest %in% c(vals_early, vals_late))]

# Prüfen, dass Vektoren nicht leer sind
vals_early <- if(length(vals_early)==0) NA else vals_early
vals_late  <- if(length(vals_late)==0) NA else vals_late
vals_rest  <- if(length(vals_rest)==0) NA else vals_rest

# Boxplot erstellen
boxplot(vals_early, vals_late, vals_rest,
        names = c("Early Greenup", "Late Maturity", "Other"),
        main = "Elevation vs Phenology Patterns",
        ylab = "Elevation (m)")