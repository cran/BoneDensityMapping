## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BoneDensityMapping)
library(rgl)


## -----------------------------------------------------------------------------
url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url, scan_filepath, mode = "wb")
CTscan <- import_scan(scan_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url2, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

landmark_path <- system.file("extdata", "test_femur.mrk.json",
                             package = "BoneDensityMapping")
landmarks <- import_lmks(landmark_path)

## -----------------------------------------------------------------------------
landmark_check(surface_mesh, landmarks, threshold = 1.0)

## -----------------------------------------------------------------------------
bone_scan_check(surface_mesh, CTscan)

## -----------------------------------------------------------------------------
vertices <- t(surface_mesh$vb)[, c(1:3)]
surface_density <- surface_normal_intersect(surface_mesh, vertices, normal_dist = 3.0, CTscan)

## -----------------------------------------------------------------------------
surface_colors <- color_mapping(surface_density, maxi = 2000, mini = 0)

## -----------------------------------------------------------------------------
plot_mesh(surface_mesh, surface_colors, title = 'Bone surface', legend_maxi = 2000, legend_mini = 0)

## -----------------------------------------------------------------------------
internal_fill <- fill_bone_points(surface_mesh, 1)

## -----------------------------------------------------------------------------
internal_density <- voxel_point_intersect(internal_fill, CTscan)

## -----------------------------------------------------------------------------
internal_colors <- color_mapping(internal_density, maxi = 2000, mini = 0)

## -----------------------------------------------------------------------------
plot_cross_section_bone(surface_mesh, surface_colors, IncludeSurface = FALSE, internal_fill, internal_colors, slice_axis = 'x', slice_val = 0.5, legend_maxi = 2000, legend_mini = 0)

## -----------------------------------------------------------------------------
patients <- list()

base_url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/"
patients <- list()

for (i in 1:6) {
  id <- paste0("SCAP00", i)
  
  mesh_url <- paste0(base_url, id, ".stl")
  ct_url <- paste0(base_url, id, "_resampled.nii.gz")
  
  mesh_path <- tempfile(fileext = ".stl")
  ct_path <- tempfile(fileext = ".nii.gz")
  
  download.file(mesh_url, mesh_path, mode = "wb")
  download.file(ct_url, ct_path, mode = "wb")
  
  landmark_path <- system.file("extdata", paste0(id, "_landmarks.fcsv"), package = "BoneDensityMapping")
  
  patients[[id]] <- list(
    mesh = import_mesh(mesh_path),
    landmarks = import_lmks(landmark_path),
    ct = import_scan(ct_path)
  )
}

## -----------------------------------------------------------------------------
patients$SCAP001$surface_points <- surface_points_template(patients$SCAP001$mesh, patients$SCAP001$landmarks, 5000)  

## -----------------------------------------------------------------------------
patient_sides <- data.frame(
  ID = sprintf("SCAP%03d", 1:6),
  Side = c("LEFT", "RIGHT", "LEFT", "LEFT", "RIGHT", "RIGHT"),
  stringsAsFactors = FALSE
)

for (i in 2:6) {
  id <- sprintf("SCAP%03d", i)

  mirror_side <- if (patient_sides$Side[i] == "RIGHT") "x" else FALSE

  patients[[id]]$surface_points <- surface_points_new(
    patients[[id]]$mesh,
    patients[[id]]$landmarks,
    patients$SCAP001$surface_points,
    mirror = mirror_side,
    plot_check = FALSE
  ) 
}

## -----------------------------------------------------------------------------
for (i in 1:6) {
  id <- sprintf("SCAP%03d", i)

  patients[[id]]$surface_density <- surface_normal_intersect(
    patients[[id]]$mesh,
    patients[[id]]$surface_points,
    normal_dist = 1.0,
    patients[[id]]$ct
  )
}

## -----------------------------------------------------------------------------
for (i in 1:6) {
  id <- sprintf("SCAP%03d", i)
  
  patients[[id]]$colored_mesh <- color_mesh(
    patients$SCAP001$mesh,
    patients$SCAP001$surface_points,
    patients[[id]]$surface_density,
    maxi = 2100, 
    mini = 0
  )
}

plot_mesh(patients$SCAP001$colored_mesh, title = 'SCAP001')
plot_mesh(patients$SCAP002$colored_mesh, title = 'SCAP002')
plot_mesh(patients$SCAP003$colored_mesh, title = 'SCAP003')
plot_mesh(patients$SCAP004$colored_mesh, title = 'SCAP004')
plot_mesh(patients$SCAP005$colored_mesh, title = 'SCAP005')
plot_mesh(patients$SCAP006$colored_mesh, title = 'SCAP006')


## -----------------------------------------------------------------------------
# Collect all surface density vectors into a matrix
density_matrix <- do.call(cbind, lapply(1:6, function(i) {
  id <- sprintf("SCAP%03d", i)
  patients[[id]]$surface_density
}))

# Average across columns (i.e. across patients), row by row
average_density <- rowMeans(density_matrix)

average_colored_mesh <- color_mesh(
    patients$SCAP001$mesh,
    patients$SCAP001$surface_points,
    average_density
  )

plot_mesh(average_colored_mesh, title = 'avg')

## -----------------------------------------------------------------------------
# Define patient IDs
young_ids <- c(1, 4, 6)
old_ids <- c(2, 3, 5)

# Helper function to get surface densities
get_densities <- function(indices) {
  do.call(cbind, lapply(indices, function(i) {
    id <- sprintf("SCAP%03d", i)
    patients[[id]]$surface_density
  }))
}

average_young_density <- rowMeans(get_densities(young_ids))
average_old_density <- rowMeans(get_densities(old_ids))

# Create colored meshes for each group
young_colored_mesh <- color_mesh(patients$SCAP001$mesh, patients$SCAP001$surface_points, average_young_density)
old_colored_mesh <- color_mesh(patients$SCAP001$mesh, patients$SCAP001$surface_points, average_old_density)

close3d()
plot_mesh(young_colored_mesh, title = "Young Group")
plot_mesh(old_colored_mesh, title = "Old Group")


density_diff <- average_young_density - average_old_density
max_abs <- max(abs(density_diff))

# Color scheme: red = young > old, blue = old > young
my_colors <- c("blue", "white", "red")

# Create colored mesh of differences
diff_colored_mesh <- color_mesh(
  patients$SCAP001$mesh,
  patients$SCAP001$surface_points,
  density_diff,
  maxi = max_abs,
  mini = -max_abs,
  color_sel = my_colors
)

# Plot it
plot_mesh(diff_colored_mesh, title = "Young - Old Density Difference", legend_maxi = max_abs, legend_mini = -max_abs, legend_color_sel = my_colors)

