#' Plot EICs
#'
#' Uses mzR to plot summed intensities (within mzmin-mzmax range) vs rt
#' as well as intensities in rt-mz-map (intensity as color)
#' The mz-axis of the rt-mz-map is expanded +/- ppm
#'
#' @param file mzML file
#' @param mz peak mz
#' @param mzmin min mz range (for EIC sum)
#' @param mzmax max mz range (for EIC sum)
#' @param ppm mz expansion in rt-mz-map
#' @param rtmin min rt range (defaults to first rt)
#' @param rtmax max rt range (defaults to last rt)
#'
#' @return plot (base R) with viridis color for intensity
#' @export
EIC <- function(file, mz, mzmin, mzmax, rtmin, rtmax, ppm = 100) {
  ms <- openMSfile(file)
  hd <- header(ms)

  if (missing(rtmin)) rtmin <- min(hd$retentionTime)
  if (missing(rtmax)) rtmax <- max(hd$retentionTime)

  scanMin <- which.min(abs(hd$retentionTime - rtmin))
  scanMax <- which.min(abs(hd$retentionTime - rtmax))
  nScan <- scanMax - scanMin + 1

  # Allocate summed intensities for EIC
  eicInt <- numeric(nScan)

  # Allocate vectors for showing rt-mz
  mzPlot <- numeric() # mz
  intPlot <- numeric() # intensity
  sPlot <- numeric() # scan number

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
    dmz <- ppm * 1e-6 * mz
    mzPlotNew <- mzs[mzs >= (mz - dmz) & mzs <= (mz + dmz)]

    # If any mz within said range store them for later plotting
    if (length(mzPlotNew) > 0) {
      mzPlot <- c(mzPlot, mzPlotNew) # Add the new mzs
      intPlot <- c(intPlot, int[mzs >= (mz - dmz) & mzs <= (mz + dmz)]) # Add corresponding intensities
      sPlot <- c(sPlot, rep(scan, length(mzPlotNew))) # Add at what scan
    }
  }

  # Extract actual rts (instead of scans)
  RT <- hd$retentionTime

  # Plot set-up
  par(mfrow = c(2, 1), mar = c(2, 4, 0, 0) + .5, oma = c(0, 0, 2, 0))

  # Plot EIC
  colEIC <- viridis::viridis(10, alpha = 0.75)[as.numeric(cut(eicInt,breaks = 10))]
  plot(RT[scanMin:scanMax], eicInt, type = 'l', col = 'grey40', las = 1, ylab = '')
  abline(h = 1000, col = 'grey80')
  points(RT[scanMin:scanMax], eicInt, pch = 15, cex = .5, col = colEIC)
  legend('topleft', legend = c('sumInt (mzmin-mazmax)', '1000'), lty = rep(1, 2), col = c('grey40', 'grey80'), bty = 'n', cex = 0.75)

  # Plot rt-mz map
  plot(RT[sPlot], mzPlot, xlim = c(rtmin, rtmax), pch = 16, cex = .5, type = 'n', ylim = c(mz - dmz, mz + dmz), las = 1, ylab = '')
  abline(h = c(mz, mzmin, mzmax), col = c('grey40', 'grey80', 'grey80'))
  colSpect <- viridis::viridis(10)[as.numeric(cut(intPlot,breaks = 10))]
  points(RT[sPlot], mzPlot, pch = 16, cex = .5, col = colSpect)
  legend('topright', legend = c('m/z peak', 'm/z min', 'm/z max'), lty = rep(1, 3), col = c('grey40', 'grey80', 'grey80'), bty = 'n', cex = 0.75)
  mtext(text = paste0('m/z ', round(mz, 5),' (',round(mzmin, 5),'-',round(mzmax, 5),')'), outer = TRUE, line = 1)
  mtext(text = file, outer = TRUE, line = 0, cex = 0.75)

}
