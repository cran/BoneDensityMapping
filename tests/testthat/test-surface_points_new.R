url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP001.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
mesh1 <- import_mesh(bone_filepath)

landmark_path <- system.file("extdata", "SCAP001_landmarks.fcsv",
                            package = "BoneDensityMapping")
lmks1 <- import_lmks(landmark_path)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP001.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url2, bone_filepath, mode = "wb")
mesh2 <- import_mesh(bone_filepath)

landmark_path <- system.file("extdata", "SCAP002_landmarks.fcsv",
                             package = "BoneDensityMapping")
lmks2 <- import_lmks(landmark_path)

test_that("surface_points_new returns valid remapped coordinates", {

  # Generate template and remap
  template <- surface_points_template(mesh1, lmks1, 100)
  result <- surface_points_new(mesh2, lmks2, template)

  # Basic checks
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(template))
  expect_true(all(is.finite(result)))
  expect_type(result[1, 1], "double")
  mesh_coords <- t(mesh2$vb)[, 1:3]
  for (i in 1:3) {
    expect_true(all(result[, i] >= min(mesh_coords[, i])))
    expect_true(all(result[, i] <= max(mesh_coords[, i])))
  }
})
