url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url2, scan_filepath, mode = "wb")
nifti <- import_scan(scan_filepath)

test_that("import_scan reads valid .nii file", {

  expect_s4_class(nifti, "nifti")
  expect_true(length(scan) > 0)
})

test_that("import_scan errors on unsupported file type", {
  fake_path <- tempfile(fileext = ".txt")
  file.create(fake_path)

  expect_error(import_scan(fake_path), "Unsupported file type")
})

test_that("import_scan errors when file does not exist", {
  missing_file <- testthat::test_path("testdata/does_not_exist.nii")

  expect_error(import_scan(missing_file))
})
