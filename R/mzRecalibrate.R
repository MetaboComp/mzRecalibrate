#' Recalibrate mzML files
#'
#' Distributes recalibration to multiple threads if a parallel backend is registered
#'
#' @param files mzML file names (character vector)
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
#' @param parallel If TRUE and backend is registered (doParallel), will perform calculations in parallel
#' @param ... Additional arguments
#'
#' @return
#' @export
#'
#' @examples
#' # To be added
mzRecalibrate <- function(files,
                          mzCandidates,
                          ppmFind = 50,
                          ppmQuantile = 0.2,
                          ppmLimit = 5,
                          intensityLimit = 1000,
                          calibLimit = 300,
                          twoStage = TRUE,
                          preCalibrate = TRUE,
                          ppmFindRough = 100,
                          ppmQuantileRough = 0.1,
                          plot = TRUE,
                          jpg = TRUE,
                          save = TRUE,
                          verbose = TRUE,
                          parallel = TRUE,
                          ...) {

  library(doParallel)

  if (parallel) "%doVersion%" <- get("%dopar%") else "%doVersion%" <- get("%do%") # Parallel vs serial

  pkg <- c('MSnbase', 'mzRecalibrate')

  iteration <- foreach(file = files, .packages = pkg) %doVersion% {
    mzRecal(file,
            mzCandidates = mzCandidates,
            ppmFind = ppmFind,
            ppmQuantile = ppmQuantile,
            ppmLimit = ppmLimit,
            intensityLimit = intensityLimit,
            calibLimit = calibLimit,
            twoStage = twoStage,
            preCalibrate = preCalibrate,
            ppmFindRough = ppmFindRough,
            ppmQuantileRough = ppmQuantileRough,
            plot = plot,
            jpg = jpg,
            save = save,
            verbose = verbose,
            ...)
  }
}
