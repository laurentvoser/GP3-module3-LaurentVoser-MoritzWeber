# load libraries
library(tidyverse)
library(geodata)

# download SRTM data
srtm_path <- file.path("data-raw", "elevation", "srtm_38_03.tif")

dem <- terra::rast(srtm_path)

# load libraries
library(MODISTools)

# download and save phenology data
phenology_path <- file.path("data-raw", "phenology.rds")

phenology <- readRDS(phenology_path)

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
land_cover_path <- file.path("data-raw", "land_cover.rds")

land_cover <- readRDS(land_cover_path)

land_cover_raster <- MODISTools::mt_to_terra(
  land_cover,
  reproject = TRUE
)
p <- ggplot() +
  tidyterra::geom_spatraster(data = dem) +
  scale_fill_viridis_c(
    na.value = NA,
    name = "altitude (m)"
  ) +
  theme_bw()

p2 <- ggplot() +
  tidyterra::geom_spatraster(data = phenology_raster) +
  scale_fill_viridis_c(
    na.value = NA,
    name = "DOY"
  ) +
  theme_bw()

dir.create("graphs", showWarnings = FALSE)

# compositing with patchwork package
library(patchwork)
combined_plot <- p + p2 +
  plot_layout(ncol = 1) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

ggsave(file.path("graphs", "dem_phenology.png"), combined_plot)

# convert to data frame and merge
dem_df <- as.vector(dem)
phenology_df <- as.vector(phenology_raster)
sct_df <- data.frame(
  altitude = dem_df,
  doy = phenology_df
)

scatter_plot <- ggplot(data = sct_df, aes(altitude, doy)) +
  geom_hex() +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE, colour = "white", lty = 2) +
  labs(x = "altitude (m)", y = "MODIS vegetation greenup (DOY)") +
  theme_bw()

ggsave(file.path("graphs", "scatter_altitude_doy.png"), scatter_plot)

# fit a linear regression to the data of the figure above
# (for the pre-processing see the collapsed code of the figure)
fit <- lm(doy ~ altitude, data = sct_df)
print(summary(fit))