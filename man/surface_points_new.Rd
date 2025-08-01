% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoneDensityFunctions.R
\name{surface_points_new}
\alias{surface_points_new}
\title{New mapped surface points from template}
\usage{
surface_points_new(
  surface_mesh,
  landmarks,
  template,
  mirror = FALSE,
  plot_check = FALSE
)
}
\arguments{
\item{surface_mesh}{List. Mesh data imported via ply_import function}

\item{landmarks}{Data frame. Contains 3D coords of landmarks}

\item{template}{Data frame. 3D coords of remapped surface points}

\item{mirror}{Logical or character. Set to "x", "y", or "z" to mirror the
mesh and landmarks across that axis before remapping.}

\item{plot_check}{Logical. If TRUE, generates a 3D plot showing the mirrored
mesh, mirrored landmarks, remapped surface points, and original template
points to visually verify correct orientation and laterality.}
}
\value{
Data frame. 3D coords of remapped surface points
}
\description{
New mapped surface points from template
}
\examples{
\donttest{
  url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP001.stl"
  bone_filepath <- tempfile(fileext = ".stl")
  download.file(url, bone_filepath, mode = "wb")
  scap_001_mesh <- import_mesh(bone_filepath)
  landmark_path <- system.file("extdata", "SCAP001_landmarks.fcsv",
                               package = "BoneDensityMapping")
  scap_001_lmk <- import_lmks(landmark_path)
  template_coords <- surface_points_template(scap_001_mesh, scap_001_lmk,
                                             1000)
  url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP002.stl"
  bone_filepath <- tempfile(fileext = ".stl")
  download.file(url2, bone_filepath, mode = "wb")
  scap_002_mesh <- import_mesh(bone_filepath)
  landmark_path <- system.file("extdata", "SCAP002_landmarks.fcsv",
                               package = "BoneDensityMapping")
  scap_002_lmk <- import_lmks(landmark_path)
  scap_002_remapped <- surface_points_new(scap_002_mesh, scap_002_lmk,
                                          template_coords, mirror = "x",
                                          plot_check = FALSE)
}
}
\author{
Scott Telfer \email{scott.telfer@gmail.com} Adapted from geomorph
}
