test_that("import_lmks correctly imports .json file", {
  json_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
  df <- import_lmks(json_path)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("lmk_id", "x", "y", "z"))
  expect_true(nrow(df) > 0)
  expect_type(df$x, "double")
  expect_type(df$lmk_id, "character")
})

test_that("import_lmks returns correct data for .fcsv file", {
  fcsv_path <- test_path("testdata/test_landmarks.fcsv")
  df <- import_lmks(fcsv_path)

  expect_s3_class(df, "data.frame")
  expect_named(df, c("lmk_id", "x", "y", "z"))
})

test_that("import_lmks errors on unsupported file type", {
  expect_error(import_lmks("bad_file_type.txt"), "Unsupported file type")
})

