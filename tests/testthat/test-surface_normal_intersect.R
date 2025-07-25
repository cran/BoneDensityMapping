url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url2, scan_filepath, mode = "wb")
nifti <- import_scan(scan_filepath)

landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
landmarks <- import_lmks(landmark_path)

test_that("surface_normal_intersect returns expected output", {
  # Generate surface points
  mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 10)

  # Run the function
  mat_peak <- surface_normal_intersect(
    surface_mesh, mapped_coords, normal_dist = 3.0,
    nifti
  )

  # Check output type and length
  expect_length(mat_peak, nrow(mapped_coords))

  # Check for NAs (optional depending on what you expect)
  expect_false(any(is.na(mat_peak)))
})
