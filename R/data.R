#' Downsampled cartridge case scans "Fadul 1-1" and "Fadul 1-2" from Fadul et
#' al. (2011)
#'
#' x3p (XML 3D Surface Profile) objects containing downsampled versions of the
#' "Fadul 1-1" and "Fadyl 1-2" cartridge case from Fadul et al. (2011). The scans
#' have been downsampled by a factor of 8 using the sample_x3p function from the
#' x3ptools package due to file size constraints.
#'
#'
#' @format An x3p object containing a surface matrix and metainformation
#'   concerning the conditions under which the scan was taken \describe{
#'   \item{header.info}{size and resolution of scan}
#'
#'   \item{surface.matrix}{spatially-ordered matrix of elements representing the
#'   height values of the cartridge case surface at particular locations}
#'
#'   \item{feature.info}{provides structure for storing surface data}
#'
#'   \item{general.info}{information concerning the author of the scan and
#'   capturing device}
#'
#'   \item{matrix.info}{provides link to surface measurements in binary format}
#'
#'   }
#'
#' @seealso T. Fadul, G. Hernandez, S. Stoiloff, and G. Sneh. An Empirical Study
#'   to Improve the Scientific Foundation of Forensic Firearm and Tool Mark
#'   Identification Utilizing 10 Consecutively Manufactured Slides, 2011.
#'
#' @seealso \url{https://github.com/heike/x3ptools}
#'
#' @source
#' \url{https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/Details/2d9cc51f-6f66-40a0-973a-a9292dbee36d}
#'
#' @name fadulData_raw
#' @keywords datasets
"fadul1.1_raw"

#' @rdname fadulData_raw
"fadul1.2_raw"

#' Processed versions of the fadul1.1_raw and fadul1.2_raw datasets using
#' functions from the cmcR package
#'
#' An x3p (XML 3D Surface Profile) object containing a processed version of the
#' "Fadul 1-1" and "Fadul 1-2" cartridge cases from Fadul et al. (2011). The
#' scans have been downsampled by a factor of 8 and processed using functions
#' from the cmcR package.
#'
#'
#' @format An x3p object containing a surface matrix and metainformation
#'   concerning the conditions under which the scan was taken \describe{
#'   \item{header.info}{size and resolution of scan}
#'
#'   \item{surface.matrix}{spatially-ordered matrix of elements representing the
#'   height values of the processed cartridge case surface at particular
#'   locations}
#'
#'   \item{feature.info}{provides structure for storing surface data}
#'
#'   \item{general.info}{information concerning the author of the scan and
#'   capturing device}
#'
#'   \item{matrix.info}{provides link to surface measurements in binary format}
#'
#'   }
#'
#' @seealso T. Fadul, G. Hernandez, S. Stoiloff, and G. Sneh. An Empirical Study
#'   to Improve the Scientific Foundation of Forensic Firearm and Tool Mark
#'   Identification Utilizing 10 Consecutively Manufactured Slides, 2011.
#'
#' @seealso \url{https://github.com/heike/x3ptools}
#'
#' @source
#' \url{https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/Details/2d9cc51f-6f66-40a0-973a-a9292dbee36d}
#'
#' @name fadulData_processed
#' @keywords datasets
"fadul1.1_processed"

#' @rdname fadulData_processed
"fadul1.2_processed"
