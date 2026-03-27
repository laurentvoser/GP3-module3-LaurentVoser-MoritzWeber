library(MODISTools)
library(terra)
library(tidyverse)
library(geodata)
library(tmap)

# --- Parameter ---
lat <- 43.5
lon <- -74.5
years <- 2001:2010

# --- MODIS Phenology-Daten laden ---
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

# Aufbereiten

phenology <- phenology %>%
  mutate(
    value = ifelse(value > 30000, NA, value),
    doy   = as.numeric(format(as.Date("1970-01-01") + value, "%j")),
    year  = as.numeric(substr(calendar_date,1,4))
  )

greenup <- phenology %>% filter(band=="Greenup.Num_Modes_01")
maturity <- phenology %>% filter(band=="Maturity.Num_Modes_01")

greenup_raster <- mt_to_terra(greenup, reproject=TRUE)
maturity_raster <- mt_to_terra(maturity, reproject=TRUE)

# --- Landcover filtern (2010) ---
land_cover <- mt_subset(
  product = "MCD12Q1",
  lat = lat, lon = lon,
  band = "LC_Type1",
  start = "2010-01-01", end = "2010-12-31",
  km_lr = 100, km_ab = 100,
  site_name = "adirondacks_lc",
  internal = TRUE
)
lc_raster <- mt_to_terra(land_cover, reproject = TRUE)
mask_forest <- lc_raster %in% c(4,5)


# --- DOY-Raster erstellen ---
# Werte in DOY umwandeln

# Use raw (unmasked) rasters as source for DOY conversion
greenup_raw  <- mt_to_terra(greenup,   reproject = TRUE)
maturity_raw <- mt_to_terra(maturity,  reproject = TRUE)

# Convert the full values matrix (preserves spatial order)
to_doy_matrix <- function(r) {
  v <- values(r)                         # N-cells × L-layers matrix
  v[v > 30000] <- NA
  matrix(
    as.numeric(format(as.Date("1970-01-01") + as.vector(v), "%j")),
    nrow = nrow(v), ncol = ncol(v)
  )
}

greenup_raster_doy  <- greenup_raw
maturity_raster_doy <- maturity_raw
values(greenup_raster_doy)  <- to_doy_matrix(greenup_raw)
values(maturity_raster_doy) <- to_doy_matrix(maturity_raw)

# Re-align mask grid (MCD12Q1 and MCD12Q2 may not snap perfectly after reproject)
mask_forest_aligned <- resample(mask_forest, greenup_raster_doy, method = "near")

# Apply mask AFTER DOY conversion
greenup_raster_doy  <- mask(greenup_raster_doy,  mask_forest_aligned, maskvalues = FALSE)
maturity_raster_doy <- mask(maturity_raster_doy, mask_forest_aligned, maskvalues = FALSE)
# --- LTM und SD 2001-2009 ---
greenup_hist <- greenup_raster_doy[[1:9]]
maturity_hist <- maturity_raster_doy[[1:9]]

greenup_ltm <- mean(greenup_hist, na.rm=TRUE)
greenup_sd  <- terra::app(greenup_hist, sd, na.rm=TRUE)
maturity_ltm <- mean(maturity_hist, na.rm=TRUE)
maturity_sd  <- terra::app(maturity_hist, sd, na.rm=TRUE)

# --- Abweichungen 2010 ---
greenup_2010 <- greenup_raster_doy[[10]]
maturity_2010 <- maturity_raster_doy[[10]]

early_greenup <- greenup_2010 < (greenup_ltm - greenup_sd)
late_maturity  <- maturity_2010 > (maturity_ltm + maturity_sd)

# --- Kontrolle Plot ---
par(mfrow=c(1,2))
terra::plot(early_greenup, main="Early Greenup 2010")
terra::plot(late_maturity, main="Late Maturity 2010")

# --- DEM / Höhe ---
dem_us <- geodata::elevation_30s(country="USA", path=tempdir())
dem_us <- terra::project(dem_us, crs(greenup_raster))
dem_resampled <- terra::resample(dem_us, greenup_raster_doy[[10]], method="bilinear")

dem_early <- terra::mask(dem_resampled, early_greenup, maskvalues=FALSE)
dem_late  <- terra::mask(dem_resampled, late_maturity, maskvalues=FALSE)

vals_early <- as.numeric(terra::values(dem_early, na.rm=TRUE))
vals_late  <- as.numeric(terra::values(dem_late, na.rm=TRUE))
vals_rest  <- as.numeric(terra::values(dem_resampled, na.rm=TRUE))
neither <- mask_forest_aligned & !early_greenup & !late_maturity
dem_rest <- terra::mask(dem_resampled, neither, maskvalues = FALSE)
vals_rest <- as.numeric(terra::values(dem_rest, na.rm = TRUE))

vals_early <- if(length(vals_early)==0) NA else vals_early
vals_late  <- if(length(vals_late)==0) NA else vals_late
vals_rest  <- if(length(vals_rest)==0) NA else vals_rest

boxplot(vals_early, vals_late, vals_rest,
        names=c("Early Greenup","Late Maturity","Other"),
        main="Elevation vs Phenology Patterns",
        ylab="Elevation (m)")

# --- Karte mit OpenStreetMap Hintergrund ---
tmap_mode("plot")

safe_tm_shape <- function(rast, col, title) {
  if(all(is.na(values(rast)))) return(NULL)
  tm_shape(rast) + tm_raster(palette=col, alpha=0.6, title=title)
}

tm <- tm_tiles("OpenStreetMap") +
  tm_shape(mask_forest) + tm_raster(palette="lightgreen", alpha=0.3, title="Forest Mask")

tm_early <- safe_tm_shape(early_greenup, "yellow", "Early Greenup")
tm_late  <- safe_tm_shape(late_maturity, "red", "Late Maturity")

if(!is.null(tm_early)) tm <- tm + tm_early
if(!is.null(tm_late))  tm <- tm + tm_late

tm <- tm + tm_layout(legend.outside=TRUE) + tm_title("Phenology Patterns 2010")
tm