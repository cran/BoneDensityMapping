url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url2, scan_filepath, mode = "wb")
nifti <- import_scan(scan_filepath)

test_that("import_mesh correctly imports a valid .stl mesh", {
  expect_true(inherits(surface_mesh, "mesh3d"))
  expect_true(!is.null(surface_mesh$vb))  # vertices
  expect_true(!is.null(surface_mesh$it))  # triangle indices
})

test_that("import_mesh errors on unsupported file format", {
  fake_file <- tempfile(fileext = ".txt")
  writeLines("not a mesh", fake_file)

  expect_error(import_mesh(fake_file))
})

