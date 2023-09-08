#' Identify Regions of Interest (ROIs) for trace ions bleeding through large RT regions
#'
#' Using findmzROI from XCMC, author: Johannes Rainer
#'
#' @param file mzML to scout for ROIs
#' @param prefilter Criterion for ROI detection: number of scans with intensity over a limit (defaults to c(400, 2000))
#' @param noise Noise level (defaults to 500)
#' @param ppm Between-scan ppm allowance (defaults to 40)
#'
#' @return Data frame mz, scans, etc
#' @export
#'
#' @examples
#' getROI(mzMLfile)
getROI <- function(file, prefilter = c(400, 2000), noise = 500, ppm = 40) {

  cat('Reading file', file, '\n')

  rawdata <- readMSData(file, msLevel = 1, verbose = FALSE)

  cat('Finished reading.\n')

  mzs <- mz(rawdata)
  mz <- unlist(mzs)
  ints <- intensity(rawdata)
  int <- unlist(ints)
  scantime <- rtime(rawdata)
  # plot(scantime)
  valsPerSpect <- lengths(mzs)
  scanindex <- xcms:::valueCount2ScanIndex(valsPerSpect)
  scanrange <- c(1, length(scantime))
  minCentroids <- 4 # This is just a guess - something to try to get the ROI detection to run :)

  cat('Detecting ROIs\n')

  roiList <- .Call("findmzROI",
                   mz, int, scanindex,
                   as.double(c(0.0, 0.0)),
                   as.integer(scanrange),
                   as.integer(length(scantime)),
                   as.double(ppm * 1e-6),
                   as.integer(minCentroids),
                   as.integer(prefilter),
                   as.integer(noise),
                   PACKAGE ='xcms')

  cat('Finished ROI detection.\n\n\n')

  if(length(roiList) > 0) {

    for (i in 1:length(roiList)) roiList[[i]] <- as.data.frame(roiList[[i]])
    roiList <- do.call(rbind, roiList)

    roiList <- roiList[order(roiList$mz), ]

    return(roiList)

  } else return(NULL)

}
