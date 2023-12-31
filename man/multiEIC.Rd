% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiEIC.R
\name{multiEIC}
\alias{multiEIC}
\title{Plot stacked EICs for multiple files}
\usage{
multiEIC(files, mz, mzmin, mzmax, ppm = 100)
}
\arguments{
\item{files}{Vector of mzML files}

\item{mz}{peak mz}

\item{mzmin}{min mz range (for EIC sum)}

\item{mzmax}{max mz range (for EIC sum)}

\item{ppm}{mz expansion in rt-mz-map}
}
\value{
plot (base R) with viridis color for intensity
}
\description{
Uses mzR to plot summed intensities (within mzmin-mzmax range) vs rt
as well as intensities in rt-mz-map (intensity as color)
The mz-axis of the rt-mz-map is expanded +/- ppm
}
