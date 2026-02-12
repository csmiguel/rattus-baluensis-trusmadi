###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# Oct 2025
###.............................................................................
#GOAL: compute area per elevation and plot it, as requested by reviewer 2
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(raster)
library(terra)

# input file
p_tm_raster <- "map/n05_n06_e116_1arc_v3.tif"

# read raster
tm_raster <- rast(p_tm_raster)

# define extent: xmin, xmax, ymin, ymax
e <- ext(116.45, 116.62, 5.52, 5.67)

# crop raster
tm_crop <- crop(tm_raster, e)
tm_utm <- project(tm_crop, "EPSG:32650")  # WGS84 / UTM zone 50N

pdf("output/map-2000-contour-line.pdf")
# plot cropped raster
plot(tm_utm)

# Create contour at 2000 m
contour_2000 <- as.contour(tm_utm, levels = 2000)

# Add contour to plot
lines(contour_2000, col = "red", lwd = 2)

# Add scale bar (length in km)
scalebar(d = 10000, xy = c(447000, 611000), type = "bar", below = "10 km")
dev.off()

breaks <- seq(2000, 2700, 100)
areas <- cellSize(tm_utm, unit = "m")
# Calculate area per bin
area_per_bin <- sapply(seq_along(breaks[-length(breaks)]), function(i) {
  # mask for pixels in this band
  band_mask <- tm_utm >= breaks[i] & tm_utm < breaks[i+1]
  # sum of areas (m²)
  sum(values(areas * band_mask), na.rm = TRUE)
})
# Assemble result as data.frame (convert to km²)
result <- data.frame(
  elev_min = breaks[-length(breaks)],
  elev_max = breaks[-1],
  area_km2 = area_per_bin / 1e6
)
pdf("output/binned-elevation-trusmadi.pdf")
graphics::barplot(height = result$area_km2,
                  names.arg = paste0(breaks, "-", breaks[-1])[1:(length(breaks)-1)],
                  ylab = "km^2",
                  xlab = "Elevation (m), binned every 100 m")
dev.off()
