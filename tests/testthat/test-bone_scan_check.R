url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url2, scan_filepath, mode = "wb")
nifti <- import_scan(scan_filepath)

test_that("bone_scan_check works with mesh input", {
  df <- bone_scan_check(surface_mesh, nifti, return_limits = TRUE)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Axis", "Mesh_Min", "Mesh_Max", "Scan_Min", "Scan_Max") %in% colnames(df)))
})

test_that("bone_scan_check works with matrix input", {
  coords <- t(surface_mesh$vb)[, 1:3]
  df <- bone_scan_check(coords, nifti, return_limits = TRUE)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Axis", "Mesh_Min", "Mesh_Max", "Scan_Min", "Scan_Max") %in% colnames(df)))
})

test_that("bone_scan_check works with data.frame input", {
  coords_df <- as.data.frame(t(surface_mesh$vb)[, 1:3])

  result <- bone_scan_check(coords_df, nifti)
  expect_null(result)  # returns invisible(NULL) when return_limits = FALSE
})

test_that("bone_scan_check errors on unsupported surface_mesh types", {
  expect_error(bone_scan_check(list(1, 2, 3), nifti, return_limits = TRUE),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
  expect_error(bone_scan_check(NULL, nifti, return_limits = TRUE),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
})

test_that("bone_scan_check errors when mesh extends outside scan volume", {
  mesh_outside <- surface_mesh
  mesh_outside$vb[1, ] <- mesh_outside$vb[1, ] + 10000

  expect_error(
    bone_scan_check(mesh_outside, nifti, return_limits = TRUE),
    "Mesh not within scan volume."
  )
})

