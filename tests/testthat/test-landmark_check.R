url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

lmk_path <- system.file("extdata", "test_femur.fcsv", package = "BoneDensityMapping")
landmarks <- import_lmks(lmk_path)

test_that("landmark_check confirms all landmarks within threshold", {
  expect_message(
    landmark_check(surface_mesh, landmarks, threshold = 2.0),
    "All landmarks are on bone surface."
  )
})


