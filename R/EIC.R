#' Plot EICs for one or multiple mzML files
#'
#' Uses mzR to produce 2 plots (per file):
#' - Summed intensities (within mzmin-mzmax range) vs rt
#' - Intensities in rt-mz-map (intensity as color)
#'
#' The mz-axis of the rt-mz-map is expanded +/- ppmPlot for visualization.
#' Multiplot default is 2 x 1 for single file and nFile x 2 for multiple files.
#' The multiplot format can be manually set using mfrow as an argument
#'
#' @param files mzML file(s)
#' @param mz peak mz
#' @param mzmin min mz range (for EIC sum; defaults to (1 - ppmEIC * 1e-6) * mz)
#' @param mzmax max mz range (for EIC sum; defaults to (1 + ppmEIC * 1e-6) * mz)
#' @param ppmEIC Defines mz bin width for EIC (used if mzmin and mzmax are not specified)
#' @param ppmPlot Defines plot range in rt-mz view
#' @param rtmin min rt range (defaults to first rt)
#' @param rtmax max rt range (defaults to last rt)
#' @param mfrow
#'
#' @return plot (base R) with viridis color for intensity
#' @export
EIC <- function(files, mz, ppmEIC = 20, mzmin, mzmax, ppmPlot = 100, rtmin, rtmax, mfrow) {

  nFile <- length(files)

  # Make automated mfrow (if not specified)
  if(missing(mfrow)) {
    if (nFile == 1) {
      mfrow = c(2, 1)
    } else {
      mfrow = c(nFile, 2)
    }
  }

  # Calculate mzmin and mzmax (if not specified) from ppmEIC
  if(missing(mzmin)) mzmin <- (1 - ppmEIC * 1e-6) * mz
  if(missing(mzmax)) mzmax <- (1 + ppmEIC * 1e-6) * mz

  # Plot set-up
  par(mfrow = mfrow, mar = c(2, 4, 0, 0) + .5, oma = c(0, 0, 2, 0))

  for (i in 1:length(files)) {
    file <- files[i]

    # Open connection to mzML file
    ms <- openMSfile(file)
    hd <- header(ms)

    # Extract actual rts (instead of scans)
    RT <- hd$retentionTime

    # Extract RT ranges (if rtmin and/or rt max not specified)
    if (missing(rtmin)) rtmin <- min(RT)
    if (missing(rtmax)) rtmax <- max(RT)

    # Convert RTs to scans
    scanMin <- which.min(abs(RT - rtmin))
    scanMax <- which.min(abs(RT - rtmax))
    nScan <- scanMax - scanMin + 1

    # Allocate summed intensities for EIC
    eicInt <- numeric(nScan)

    # Allocate vectors for showing rt-mz
    mzPlot <- numeric() # mz
    intPlot <- numeric() # intensity
    sPlot <- numeric() # scan number

    # for rt-mz visualization - expand with ppm
    dmz <- ppmPlot * 1e-6 * mz

    # loop thru scans to extract points and EIC
    for (s in 1:nScan) {
      scan <- s + scanMin - 1
      spect <- peaks(ms, scan)
      mzs <- spect[,1]
      int <- spect[,2]
      # plot(mzs, int, pch = '.')

      # Extract EIC summed intensities within mzmin-mzmax range
      eicInt[s] <- sum(int[mzs >= mzmin & mzs <= mzmax])

      # for rt-mz visualization - expand with ppm
      mzPlotNew <- mzs[mzs >= (mz - dmz) & mzs <= (mz + dmz)]

      # If any mz within said range store them for later plotting
      if (length(mzPlotNew) > 0) {
        mzPlot <- c(mzPlot, mzPlotNew) # Add the new mzs
        intPlot <- c(intPlot, int[mzs >= (mz - dmz) & mzs <= (mz + dmz)]) # Add corresponding intensities
        sPlot <- c(sPlot, rep(scan, length(mzPlotNew))) # Add at what scan
      }
    }

    # Plot EIC
    colEIC <- viridis::viridis(10, alpha = 0.75)[as.numeric(cut(eicInt,breaks = 10))]
    plot(RT[scanMin:scanMax], eicInt, type = 'l', col = 'grey40', las = 1, ylab = '', main = file, cex.main = .85)
    abline(h = 1000, col = 'grey80')
    points(RT[scanMin:scanMax], eicInt, pch = 15, cex = .5, col = colEIC)
    legend('topleft', legend = c('sumInt (mzmin-mazmax)', '1000'), lty = rep(1, 2), col = c('grey40', 'grey80'), bty = 'n', cex = 0.75)

    # Plot rt-mz map
    if(length(sPlot) > 0) {

    }
    plot(RT[sPlot], mzPlot, xlim = c(rtmin, rtmax), pch = 16, cex = .5, type = 'n', ylim = c(mz - dmz, mz + dmz), las = 1, ylab = '', main = file, cex.main = .85)
    abline(h = c(mz, mzmin, mzmax), col = c('grey40', 'grey80', 'grey80'))
    if (length(intPlot) > 0) {
      colSpect <- viridis::viridis(10)[as.numeric(cut(intPlot,breaks = 10))]
      points(RT[sPlot], mzPlot, pch = 16, cex = .5, col = colSpect)
    }
    legend('topright', legend = c('m/z peak', 'm/z min', 'm/z max'), lty = rep(1, 3), col = c('grey40', 'grey80', 'grey80'), bty = 'n', cex = 0.75)

  }
  mtext(text = paste0('m/z ', round(mz, 5),' (',round(mzmin, 5),'-',round(mzmax, 5),')'), outer = TRUE, line = 1)
}
