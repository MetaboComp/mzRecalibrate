#' Identify Regions of Interest (ROIs) for trace ions bleeding through large RT regions
#'
#' Finds and merges ROIs from multiple files
#'
#' @param files character vector with pathnames to mzML files
#' @param dmz mz allowance to merge ROIs from multiple files
#' @param ... Pass on arguments to getROI()
#'
#' @return Data frame with candidates
#' @export
#'
#' @examples
#' candidates <- findCandidates(mzMLfiles)
findCandidates <- function(files, dmz = 0.002, ...) {

  # Support function to identify and merge very close features
  mergeROIs <- function(ROIs, dmz) {
    # Calculate mz distances between hits
    roiDist <- dist(ROIs$mz) %>% as.matrix()
    roiDist[roiDist == 0] <- NA
    roiDist[lower.tri(roiDist)] <- NA

    # Are any ROIs really close in mz?
    if(any(roiDist < dmz, na.rm = TRUE)) {
      whichDist <- which(roiDist < dmz, arr.ind = TRUE)
      removes <- numeric()

      # Merge hits
      for (d in 1:nrow(whichDist)) {
        merge1 <- whichDist[d, 1]
        merge2 <- whichDist[d, 2]
        # ROIs[c(merge1, merge2),]
        # Should only be merged if they don't overlap
        if (!any(table(c(ROIs[merge1, ]$scmin:ROIs[merge1, ]$scmax, ROIs[merge2, ]$scmin:ROIs[merge2, ]$scmax)) == 2)) {
          ROIs[merge1, ]$mzmin <- min(c(ROIs[merge1, ]$mzmin, ROIs[merge2, ]$mzmin))
          ROIs[merge1, ]$mzmax <- max(c(ROIs[merge1, ]$mzmax, ROIs[merge2, ]$mzmax))
          ROIs[merge1, ]$mz <- ((ROIs[merge1, ]$mz * ROIs[merge1, ]$length) + (ROIs[merge2, ]$mz * ROIs[merge2, ]$length)) / (ROIs[merge1, ]$length + ROIs[merge2, ]$length)
          ROIs[merge1, ]$scmin <- min(c(ROIs[merge1, ]$scmin, ROIs[merge2, ]$scmin))
          ROIs[merge1, ]$scmax <- max(c(ROIs[merge1, ]$scmax, ROIs[merge2, ]$scmax))
          ROIs[merge1, ]$length <- ROIs[merge1, ]$scmax - ROIs[merge1, ]$scmin + 1
          removes <- c(removes, merge2)
        }
      }
      if (length(removes) > 0) ROIs <- ROIs[-removes, ]
    }
    return(ROIs)
  }

  # Loop thru all files to stack up ROIs
  for (i in 1:length(files)) {

    # Read in first file -> Start of roiList
    if (i == 1) {

      roiList <- getROI(files[i], ...) # Find ROIs
      roiList <- subset(roiList, select = -intensity) # Intensity is hard to merge -> drop it
      roiList$fileCount <- 1 # Counter for number of times an ROI appears
      roiList <- mergeROIs(roiList, dmz = dmz) # Merge signals that appear with multiple ROIs (e.g. early & late)

    } else {

      # As for the 1st file - See above for comments
      roiNew <- getROI(files[i], ...)
      roiNew <- subset(roiNew, select = -intensity)
      roiNew$fileCount <- 1
      roiNew <- mergeROIs(roiNew, dmz = dmz)

      # Merge, add or discard new ROIs
      for (r in 1:nrow(roiNew)) {

        # Potential matches to previous hits (by mz)
        potentialROI <- roiNew[r,]$mz > roiList$mzmin & roiNew[r,]$mz < roiList$mzmax

        # Merge if matching 1 ROI exclusively
        if (sum(potentialROI) == 1) {
          merge1 <- which(potentialROI)
          roiList[merge1, ]$mzmin <- min(c(roiList[merge1, ]$mzmin, roiNew[r,]$mzmin))
          roiList[merge1, ]$mzmax <- max(c(roiList[merge1, ]$mzmax, roiNew[r,]$mzmax))
          roiList[merge1, ]$mz <- ((roiList[merge1, ]$mz * roiList[merge1, ]$length) + (roiNew[r,]$mz * roiNew[r,]$length)) / (roiList[merge1, ]$length + roiNew[r,]$length)
          roiList[merge1, ]$scmin <- min(c(roiList[merge1, ]$scmin, roiNew[r,]$scmin))
          roiList[merge1, ]$scmax <- max(c(roiList[merge1, ]$scmax, roiNew[r,]$scmax))
          roiList[merge1, ]$length <- roiList[merge1, ]$scmax - roiList[merge1, ]$scmin + 1
          roiList[merge1, ]$fileCount <- roiList[merge1, ]$fileCount + 1
        }

        # Add if not matching any previous ROI
        if (sum(potentialROI) == 0) roiList[nrow(roiList) + 1, ] <- roiNew[r,]
      }
    }
  }

  roiList <- roiList[order(roiList$mz), ]

  return(roiList)
}
