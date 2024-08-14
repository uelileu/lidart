# setup
setwd(here::here())

library(sf)
library(lidR)
library(ggplot2)
library(ggthemes)
library(tidyverse)

if (!dir.exists("temp")) {
  dir.create("temp")
}
if (!dir.exists("img")) {
  dir.create("img")
}


# dl laz file -------------------------------------------------------------
url.laz <- "https://maps.zh.ch/download/hoehen/2022/lidar/2683000_1247000.laz"
if (!file.exists("temp/2683000_1247000.laz")) {
  download.file(url = url.laz,
                destfile = "temp/2683000_1247000.laz")
}

# load shape
outline <- read_sf("geom/cut_geom.gpkg")

# load point cloud
laz <- readLAS("temp/2683000_1247000.laz")

# clip and tranform point cloud
clipTransform <- function(laz, outline) {
  subs <- clip_roi(las = laz, geometry = outline)

  # subs <- readLAS("luz.las")
  subset_df <- data.frame(
    subs$X,
    subs$Y,
    subs$Z
  )
  
  # Ensure points is a numeric matrix
  points <- as.matrix(subset_df)
  # Calculate centroid for X and Y only (keep Z fixed)
  centroid_xy <- colMeans(points[, 1:2])
  
  # Center the points in the XY plane
  points_centered_xy <- sweep(points[, 1:2], 2, centroid_xy)
  # Covariance matrix for XY plane
  cov_matrix_xy <- cov(points_centered_xy)
  # Eigen decomposition for XY plane
  eigen_decomp_xy <- eigen(cov_matrix_xy)
  # Eigenvectors for rotation matrix in XY plane
  eigenvectors_xy <- eigen_decomp_xy$vectors
  # Rotation matrix for XY plane (Z remains unchanged)
  rotation_matrix_xy <- eigenvectors_xy
  # Rotate the points in the XY plane
  rotated_xy <- points_centered_xy %*% rotation_matrix_xy

  # Combine rotated XY with the original Z values
  rotated_points <- cbind(rotated_xy, points[, 3])
  # Translate back in the XY plane if needed
  # rotated_points[, 1:2] <- sweep(rotated_points[, 1:2], 2, centroid_xy, FUN = "+")
  # begin X and Y by zero
  rotated_points[, 1] <- rotated_points[, 1] - min(rotated_points[, 1])
  rotated_points[, 2] <- rotated_points[, 2] - min(rotated_points[, 2])

  rot_df <- as.data.frame(rotated_points)
  names(rot_df) <- c("X", "Y", "Z")
  rot_df$distance = abs(rot_df$Y - mean(rot_df$Y)) / 10
  
  return(rot_df)
}

cutting <- clipTransform(laz, outline)


# plot --------------------------------------------------------------------

ggplot(cutting) + 
  geom_point(aes(X, Z, alpha = distance), colour="white", size = 1.1) +
  coord_equal(ratio = 1) +
  ggthemes::theme_tufte() +
  theme(panel.background = element_rect(fill = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  )  


ggsave("img/profile.png", height = 120, width = 120, units = "cm", dpi = 200,
       limitsize = FALSE)
knitr::plot_crop("img/profile.png")
