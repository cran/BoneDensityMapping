% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoneDensityFunctions.R
\name{import_lmks}
\alias{import_lmks}
\title{import landmark coordinates}
\usage{
import_lmks(landmark_path, x = 1, y = 1, z = 1)
}
\arguments{
\item{landmark_path}{String. File path to landmark data. .json or .fcsv
format}

\item{x}{Integer Value to apply to convert mesh i.e. -1 will mirror x coords}

\item{y}{Integer Value to apply to convert mesh i.e. -1 will mirror y coords}

\item{z}{Integer Value to apply to convert mesh i.e. -1 will mirror z coords}
}
\value{
dataframe. Columns are landmark name, x, y, and z coordinates
}
\description{
import landmark coordinates
}
\examples{
landmark_path <- system.file("extdata", "test_femur.mrk.json",
                             package = "BoneDensityMapping")
landmarks <- import_lmks(landmark_path)
}
\author{
Scott Telfer \email{scott.telfer@gmail.com}
}
