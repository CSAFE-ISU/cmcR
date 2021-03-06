% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcess.R
\name{preProcess_crop}
\alias{preProcess_crop}
\title{Remove observations from the exterior of interior of a breech face scan}
\usage{
preProcess_crop(
  x3p,
  region = "exterior",
  radiusOffset = 0,
  croppingThresh = 1,
  agg_function = median,
  scheme = 3,
  high_connectivity = FALSE,
  tolerance = 0
)
}
\arguments{
\item{x3p}{an x3p object containing the surface matrix of a cartridge case
scan}

\item{region}{dictates whether the observations on the "exterior" or
"interior" of the scan are removed}

\item{radiusOffset}{number of pixels to add to estimated breech face radius.
This is commonly a negative value (e.g., -30 for region = "exterior") to
trim the cartridge case primer roll-off from the returned, cropped surface
matrix or a positive value (e.g., 200 for region = "interior") to remove
observations around the firing pin impression hole.}

\item{croppingThresh}{minimum number of non-NA pixels that need to be in a
row/column for it to not be cropped out of the breech face scan exterior}

\item{agg_function}{the breech face radius estimation procedure returns a
number of radius estimates. This argument dictates the function used to
aggregate these into a final estimate.}

\item{scheme}{argument for imager::imgradient}

\item{high_connectivity}{argument for imager::label}

\item{tolerance}{argument for imager::label}
}
\value{
An x3p object containing the surface matrix of a breech face
impression scan where the observations on the exterior/interior of the
breech face scan surface.
}
\description{
Remove observations from the exterior of interior of a breech face scan
}
\note{
The radius estimation procedure tends to over-estimate the desired
radius values. As such, a lot of the breech face impression "roll-off" is
included in the final scan. Excessive roll-off can bias the calculation of
the CCF. As such, we can manually shrink the radius estimate (-30 or -30
seems to work well for the Fadul cartridge cases) so that little to no
roll-off is included in the final processed scan.

The radius estimation procedure is effective at estimating the radius of
the firing pin hole. Unfortunately, it is often desired that more than just
observations in firing pin hole are removed. In particular, the plateaued
region surrounding the firing pin impression hole does not come into contact
with the breech face of a firearm and is thus unwanted in the final,
processed scan. The radiusOffset argument must be tuned (around 200 seems to
work well for the Fadul cartridge cases) to remove these unwanted
observations.
}
\examples{

#Process fadul1.1 "from scratch" (takes > 5 seconds to run)
\dontrun{
nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"

fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))

fadul1.1_extCropped <- preProcess_crop(x3p = fadul1.1,
                                       radiusOffset = -30,
                                       region = "exterior")

fadul1.1_extIntCropped <- preProcess_crop(x3p = fadul1.1_extCropped,
                                          radiusOffset = 200,
                                          region = "interior")

x3pListPlot(list("Original" = fadul1.1,
                 "Exterior Cropped" = fadul1.1_extCropped,
                 "Exterior & Interior Cropped" = fadul1.1_extIntCropped ))
}
}
