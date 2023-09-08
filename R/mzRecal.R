#' Recalibrate mass axis for MS1 (mzML format) - One file only
#'
#' Please see mzRecalibrate for standard interface for multiple files and parallelization
#'
#' @param file mzML file
#' @param mzCandidates data frame
#' @param ppmFind ppm limit to discover candidate from data
#' @param ppmQuantile ppm quantile limit to exclude potential candidates from calibration (deviation from observed ppms)
#' @param ppmLimit ppm quantile limit to exclude potential candidates from calibration (mean absolute error)
#' @param intensityLimit minimum signal intensity (to avoid "frayed rope")
#' @param calibLimit Lower limit, above which there must be a candidate in order to perform calibration
#' @param twoStage If TRUE, performs a rough calibration in-between pre-calibration (using last scan settings) and fine calibration
#' @param preCalibrate If TRUE, performs a pre-calibration of the mass axis using the calibration line from the last scan
#' @param ppmFindRough ppm limit to discover candidate from data in rough calibration (if twoStage = TRUE)
#' @param ppmQuantileRough ppm quantile limit to exclude potential candidates from rough calibration
#' @param plot If TRUE plots calibration metrics (mz and number of calibrants; intercept and slope)
#' @param jpg If TRUE and plot = TRUE plots to jpg file in "mzRecal_log/" folder
#' @param save If TRUE saves recalibrated file in "mzRecal/" folder
#' @param verbose If TRUE reveals additional verbose information (if save = TRUE, sinks to txt file in "mzRecal_log/" folder)
#' @param ...
#'
#' @return saves recalibrated file
#' @export
#'
#' @examples
#' # To be added
mzRecal <- function(file,
                    mzCandidates,
                    ppmFind = 50,
                    ppmQuantile = 0.2,
                    ppmLimit = 5,
                    intensityLimit = 1000,
                    calibLimit = 300,
                    preCalibrate = TRUE,
                    twoStage = TRUE,
                    ppmFindRough = 100,
                    ppmQuantileRough = 0.1,
                    plot = TRUE,
                    jpg = TRUE,
                    save = TRUE,
                    verbose = FALSE,
                    ...) {

  # Sink log
  if (save & verbose) {
    if(!dir.exists(paste0(dirname(file), '/mzRecal_log'))) dir.create(paste0(dirname(file), '/mzRecal_log'))
    sink(file = paste0(dirname(file), '/mzRecal_log/', basename(file), '.txt'))
  }

  # Read raw data and extract header info
  MS <- readMSData(file, msLevel = 1, verbose = FALSE)
  nScan <- length(MS)
  nCandidates <- length(mzCandidates)
  assayDataEnvironment <- new.env()

  # Allocate storage objects
  mzCalibrators <- matrix(nrow = length(mzCandidates), ncol = nScan, dimnames = list(mzCandidates,featureNames(MS)))
  intCalibrators <- matrix(nrow = length(mzCandidates), ncol = nScan, dimnames = list(mzCandidates,featureNames(MS)))
  nCalibrators <- numeric(nScan)
  ppmCalibrators <- list()
  intercept <- numeric(nScan)
  slope <- numeric(nScan)
  coefsLastScan <- c(0, 1)

  ####################################
  # Loop thru & calibrate all scans
  ####################################

  for (s in 1:nScan) {

    # Criterion for whether to abort calibration and calibrae using last scan
    calibrateUsingLast = FALSE

    # Extract spectral data
    spectrum <- MS[[s]]
    mzs <- mz(spectrum)
    ints <- intensity(spectrum)
    intBool <- ints > intensityLimit
    # plot(mzs, ints, type = 'l')

    ####################################
    # Pre-calibration
    ####################################

    if (!preCalibrate | s == 1) {
      ####################################
      # Force pre-calibration to 1:1 line
      ####################################

      coefs <- c(0, 1) # intercept = 0, slope = 1 by default
      coefsLast <- coefs

    }

    # Predict mz axis values from pre-calibration
    mzPreCal <- coefs[1] + coefs[2] * mzs

    if (twoStage) {

      ####################################
      # Rough calibration
      ####################################

      #########################################
      # Find candidates for rough calibration
      #########################################

      # Allocate vector for candidates
      mzFind <- rep(NA, length(mzCandidates)) # In recorded data
      intFind <- rep(NA, length(mzCandidates)) # In recorded data
      mzFindRough <- rep(NA, length(mzCandidates)) # In rough pre-calibrated data

      # Find candidate matches in the scan
      for (c in 1:nCandidates) {
        mzRef <- mzCandidates[c]
        dmz <- ppmFindRough * 1e-6 * mzRef
        # How many potential hits for each candidate in the ppm range?
        boolHit <- (mzPreCal >= (mzRef - dmz)) & (mzPreCal <= (mzRef + dmz)) & intBool
        # Limit to those that have one hit only
        if (sum(boolHit) == 1) {
          mzFind[c] <- mzs[boolHit]
          mzFindRough[c] <- mzPreCal[boolHit]
          intFind[c] <- ints[boolHit]
        }
      }

      # Exclude hits with too deviating ppm difference
      ppm <- ((mzFindRough - mzCandidates) / mzCandidates) * 1e6
      ppmLimits <- quantile(ppm, probs = c(ppmQuantileRough, 1 - ppmQuantileRough), na.rm = TRUE)
      removes <- which((ppm < ppmLimits[1]) | (ppm > ppmLimits[2]))
      ppm[removes] <- NA
      mzFind[removes] <- NA
      intFind[removes] <- NA

      # View(cbind(mzFind, mzFindRough, mzCandidates, ppm))

      #########################################
      # Perform rough calibration
      #########################################

      # At least 2 points needed for lm & at least one calibrant with mz > calibLimit
      if(((nCandidates - sum(is.na(mzFind))) >= 2) & any(mzFind > calibLimit, na.rm = TRUE)) {
        coefs <- coef(.lm.fit(cbind(rep(1, sum(!is.na(mzFind))), mzFind[!is.na(mzFind)]), mzCandidates[!is.na(mzFind)]))
      } else {
        # If not - carry forward previous calibration
        coefs <- coefsLast
      }

      # Store calibration as last
      coefsLast <- coefs

      # Use rough mass calibration to get rough pre-calibrations
      mzRough <- coefs[1] + coefs[2] * mzs

    } else mzRough <- mzPreCal

    #########################################
    # Find candidates for fine calibration
    #########################################

    # Allocate vector for candidates
    mzFind <- rep(NA, length(mzCandidates)) # In recorded data
    intFind <- rep(NA, length(mzCandidates)) # In recorded data
    mzFindFine <- rep(NA, length(mzCandidates)) # In rough pre-calibrated data

    # Find candidate matches in the scan
    for (c in 1:nCandidates) {
      mzRef <- mzCandidates[c]
      dmz <- ppmFind * 1e-6 * mzRef
      # How many potential hits for each candidate in the ppm range?
      boolHit <- (mzRough >= (mzRef - dmz)) & (mzRough <= (mzRef + dmz)) & intBool
      # Limit to those that have one hit only
      if (sum(boolHit) == 1) {
        mzFind[c] <- mzs[boolHit]
        mzFindFine[c] <- mzRough[boolHit]
        intFind[c] <- ints[boolHit]
      }
    }
    #
    #     # diagnostic plots
    #     par(mfrow = c(2, 2), mar = c(2, 2, 0, 0) + .5)
    #     # par(mfrow = c(1, 1), mar = c(2, 2, 0, 0) + .5)
    #     plot(mzCandidates, mzCandidates - mzFind)
    #     abline(h = c(0, 0.002, -0.002), col = c('grey40', 'grey80', 'grey80'), lty = c(1, 2, 2))
    #     plot(mzCandidates, mzCandidates - mzFindFine)
    #     abline(h = c(0, 0.002, -0.002), col = c('grey40', 'grey80', 'grey80'), lty = c(1, 2, 2))
    #     plot(mzCandidates, (mzCandidates - mzFind)/mzCandidates)
    #     abline(h = c(0, 10e-6, -10e-6), col = c('grey40', 'grey80', 'grey80'), lty = c(1, 2, 2))
    #     plot(mzCandidates, (mzCandidates - mzFindFine)/mzCandidates)
    #     abline(h = c(0, 10e-6, -10e-6), col = c('grey40', 'grey80', 'grey80'), lty = c(1, 2, 2))

    # Exclude hits with too deviating ppm difference
    ppm <- ((mzFindFine - mzCandidates) / mzCandidates) * 1e6
    ppmLimits <- quantile(ppm, probs = c(ppmQuantile, 1 - ppmQuantile), na.rm = TRUE)
    removes <- which((ppm < ppmLimits[1]) | (ppm > ppmLimits[2]))
    ppm[removes] <- NA
    mzFind[removes] <- NA
    intFind[removes] <- NA

    # View(cbind(mzFind, mzFindRough, mzFindFine, mzCandidates, ppm))

    #########################################
    # Perform fine calibration
    #########################################

    # At least 2 points needed for lm & at least one calibrant with mz > calibLimit
    if(((nCandidates - sum(is.na(mzFind))) >= 2) & any(mzFind > calibLimit, na.rm = TRUE)) {
      # calib <- lm(post ~ pre, data = mzFindDF)
      coefs <- coef(.lm.fit(cbind(rep(1, sum(!is.na(mzFind))), mzFind[!is.na(mzFind)]), mzCandidates[!is.na(mzFind)]))

    } else {
      # If not - carry forward last scan's calibration
      if(verbose) cat('Scan', s, '- Lack of calibrants - Using calibration from last scan.\n')
      coefs <- coefsLastScan
    }

    # Calculate calibrated mass axis
    mzFine <- coefs[1] + coefs[2] * mzs

    # Allocate vector for calibrated candidates
    mzCalibrated <- rep(NA, length(mzCandidates)) # In rough pre-calibrated data

    # Find calibrated candidate matches in the scan
    for (c in 1:nCandidates) {
      mzRef <- mzCandidates[c]
      dmz <- ppmFind * 1e-6 * mzRef
      # How many potential hits for each candidate in the ppm range?
      boolHit <- (mzFine >= (mzRef - dmz)) & (mzFine <= (mzRef + dmz)) & intBool
      # Limit to those that have one hit only
      if (sum(boolHit) == 1) {
        mzCalibrated[c] <- mzFine[boolHit]
      }
    }

    mzCalibrated[is.na(mzFind)] <- NA
    ppm <- ((mzCandidates - mzCalibrated) / mzCandidates) * 1e6
    # View(cbind(mzFind, mzFindRough, mzFindFine, mzCalibrated, mzCandidates, ppm))

    ppm <- ppm[complete.cases(ppm)]

    # Find if ppms are off
    if(length(ppm) == 0) {
      if(verbose) {
        cat('Scan', s, '- No calibrant found after calibration - ')
      }
      calibrateUsingLast = TRUE
    } else {
      if(mean(abs(ppm)) > ppmLimit) {
        if(verbose) {
          cat('Scan', s, '- Average calibrant ppm > ppmLimit of',ppmLimit,'- ')
        }
        calibrateUsingLast = TRUE
      }
      if(length(ppm) >= 2 && !is.na(t.test(ppm, mu = 0)$p.value)) {
        if(t.test(ppm, mu = 0)$p.value < 0.05) {
          if(verbose) {
            cat('Scan', s, '- Calibrant ppm different from zero at p < 0.05 - ')
          }
          calibrateUsingLast = TRUE
        }
      }
    }
    # if ppms are off use last scan calibration
    if(calibrateUsingLast) {
      if(verbose) cat('Using calibration from last scan.\n')
      mzFind <- rep(NA, length(mzFind))
      intFind <- rep(NA, length(intFind))
      mzFine <- coefsLastScan[1] + coefsLastScan[2] * mzs
    } else {
      if (length(ppm) > 0) ppm <- data.frame(scan = s, ppm = ppm)
    }

    # Store calibration data
    mzCalibrators[,s] <- mzFind
    intCalibrators[,s] <- intFind
    nCalibrators[s] <- sum(!is.na(mzFind))
    if(nCalibrators[s] > 0) ppmCalibrators[[s]] <- ppm
    intercept[s] <- coefs[1]
    slope[s] <- coefs[2]

    # Store re-calibrated mass axis
    spectrum@mz <- mzFine
    assign(featureNames(MS)[s], spectrum, envir = assayDataEnvironment)

    # Store calibration as last
    coefsLast <- coefs
    coefsLastScan <- coefs

  }

  # convert ppmCalibrators to DF
  ppmCalibrators <- do.call(rbind, ppmCalibrators)

  # Open jpg connection
  if(jpg & plot) {
    if(!dir.exists(paste0(dirname(file), '/mzRecal_log'))) dir.create(paste0(dirname(file), '/mzRecal_log'))
    jpeg(filename = paste0(dirname(file), '/mzRecal_log/', basename(file),'.jpg'), width = 1000, height = 2000, pointsize = 30)
  }

  # Plot
  if(plot) {
    par(mfrow = c(5, 1), mar = c(2, 6, 0, 0) + .5, oma = c(0, 0, 2, 0))
    matplot(rtime(MS), t(mzCalibrators), type = 'l', lty = 1, ylab = 'm/z', las = 1)
    plot(rtime(MS), nCalibrators, type = 'l', ylab = 'nCalib', las = 1)
    plot(rtime(MS)[ppmCalibrators$scan], ppmCalibrators$ppm, pch = '.', ylab = 'ppmCalib', las = 1)
    plot(rtime(MS), intercept, type = 'l', ylab = 'intercept', las = 1)
    plot(rtime(MS), slope, type = 'l', ylab = 'slope', las = 1)
    mtext(text = file, outer = TRUE, cex = 0.75)
  }

  # Close jpg
  if(jpg & plot) dev.off()

  # Store spectra in assayData slot
  lockEnvironment(assayDataEnvironment, bindings = TRUE)
  MS@assayData <- assayDataEnvironment

  # Write the new, re-calibrated mzML file
  if (save) {
    if(!dir.exists(paste0(dirname(file), '/mzRecal'))) dir.create(paste0(dirname(file), '/mzRecal'))
    writeMSData(MS, file = paste0(dirname(file), '/mzRecal/', basename(file)), copy = TRUE)
  }

  # Close sink
  if (save & verbose) sink()

  # Return mz re-calibration metrics
  return(invisible(list(mzCalibrators = mzCalibrators,
                        intCalibrators = intCalibrators,
                        nCalibrators = nCalibrators,
                        intercept = intercept,
                        slope = slope,
                        rtime = rtime(MS))))
}
