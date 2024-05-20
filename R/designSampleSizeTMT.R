designSampleSizeTMT = function(
    data, desiredFC, FDR = 0.05, numSample = TRUE, power = 0.9,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL
) {

  # MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
  #                                     log_file_path, "MSstats_sampleSize_log_")
  # getOption("MSstatsLog")("INFO", "** MSstats - designSampleSize function")
  # getOption("MSstatsLog")("INFO", paste0("Desired fold change = ", 
  #                                        paste(desiredFC, collapse=" - ")))
  # getOption("MSstatsLog")("INFO", paste0("FDR = ", FDR))
  # getOption("MSstatsLog")("INFO", paste0("Power = ", power))
  
  var_component = .getVarComponentTMT(data)

  median_sigma_error = median(var_component[["Error"]], na.rm = TRUE)
  median_sigma_subject = .getMedianSigmaSubject(var_component)
  median_sigma_run = .getMedianSigmaRun(var_component)
  getOption("MSstatsLog")("INFO", "Calculated variance component. - okay")
  
  
  
  ## power calculation
  if (isTRUE(power)) {
    delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
    desiredFC = 2 ^ delta
    power_output = .calculatePower(desiredFC, FDR, delta, median_sigma_error, 
                                   median_sigma_subject, median_sigma_run, 
                                   numSample)        
    var_comp = (median_sigma_error / numSample + median_sigma_subject / numSample + median_sigma_run / numSample)
    CV = round( (2 * var_comp) / desiredFC, 3)
    getOption("MSstatsLog")("INFO", "Power is calculated. - okay")
    sample_size = data.frame(desiredFC, numSample, FDR, 
                             power = power_output, CV)
  }	
  
  if (is.numeric(power)) {
    delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
    desiredFC = 2 ^ delta
    ## Large portion of proteins are not changing
    m0_m1 = 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
    alpha = power * FDR / (1 + (1 - FDR) * m0_m1)
    if (isTRUE(numSample)) {
      numSample = .getNumSample(desiredFC, power, alpha, delta,
                                median_sigma_error, median_sigma_subject,
                                median_sigma_run)
      var_comp = (median_sigma_error / numSample + median_sigma_subject / numSample + median_sigma_run / numSample)
      CV = round(2 * var_comp / desiredFC, 3)
      getOption("MSstatsLog")("INFO", "The number of sample is calculated. - okay")
      sample_size = data.frame(desiredFC, numSample, FDR, power, CV)
    }
  } 
  sample_size
}


#' Get variances from models fitted by the groupComparison function
#' @param fitted_models FittedModels element of groupComparison output
#' @keywords internal
.getVarComponentTMT = function(fitted_models) {
  Protein = NULL
  
  result = data.table::rbindlist(
    lapply(fitted_models, function(fit) {
      if (!is.null(fit)) {
        if (!is(fit, "lmerMod")) {
          error = summary(fit)$sigma^2
          subject <- NA
          group_subject <- NA
          run <- NA
          mix_run <- NA
        } else {
          stddev = c(sapply(lme4::VarCorr(fit), function(el) attr(el, "stddev")), 
                     attr(lme4::VarCorr(fit), "sc"))
          error = stddev[names(stddev) == ""]^2
          if (any(names(stddev) %in% "SUBJECT.(Intercept)")) {
            subject = stddev["SUBJECT.(Intercept)"]^2
          } else {
            subject = NA
          }
          if (any(names(stddev) %in% "SUBJECT:GROUP.(Intercept)")) {
            group_subject = stddev["SUBJECT:GROUP.(Intercept)"]^2
          } else {
            group_subject = NA
          }
          if (any(names(stddev) %in% "Run.(Intercept)")) {
            run = stddev["Run.(Intercept)"]^2
          } else {
            run = NA
          }
          if (any(names(stddev) %in% "Mixture:TechRepMixture.(Intercept)")) {
            mix_run = stddev["Mixture:TechRepMixture.(Intercept)"]^2
          } else {
            mix_run = NA
          }
          
        }
        list(Error = error,
             Subject = subject,
             GroupBySubject = group_subject,
             Run = run,
             RunByMixture=mix_run)
      } else {
        NULL
      }
    })
  )
  result[, Protein := 1:.N]
  result
}


#' Get median per subject or group by subject
#' @param var_component data.frame, output of .getVarComponent
#' @importFrom stats median
#' @keywords internal
.getMedianSigmaRun = function(var_component) {
  if (sum(!is.na(var_component[, "RunByMixture"])) > 0) {
    median_sigma_run = median(var_component[["RunByMixture"]], na.rm=TRUE)
  } else {
    if (sum(!is.na(var_component[, "Run"])) > 0) {
      median_sigma_run = median(var_component[["Run"]], na.rm=TRUE)
    } else {
      median_sigma_run = 0
    }
  }
  median_sigma_run
}

#' Get median per run or run by mix
#' @param var_component data.frame, output of .getVarComponent
#' @importFrom stats median
#' @keywords internal
.getMedianSigmaSubject = function(var_component) {
  if (sum(!is.na(var_component[, "GroupBySubject"])) > 0) {
    median_sigma_subject = median(var_component[["GroupBySubject"]], na.rm=TRUE)
  } else {
    if (sum(!is.na(var_component[, "Subject"])) > 0) {
      median_sigma_subject = median(var_component[["Subject"]], na.rm=TRUE)
    } else {
      median_sigma_subject = 0
    }
  }
  median_sigma_subject
}


#' Power calculation
#' @inheritParams designSampleSize
#' @param delta difference between means (?)
#' @param median_sigma_error median of error standard deviation
#' @param median_sigma_subject median standard deviation per subject
#' @importFrom stats qnorm
#' @keywords internal
.calculatePower = function(desiredFC, FDR, delta, median_sigma_error, 
                           median_sigma_subject, median_sigma_run, numSample) {
  m0_m1 = 99
  t = delta / sqrt(2 * (median_sigma_error / numSample + median_sigma_subject / numSample + median_sigma_run / numSample))
  powerTemp = seq(0, 1, 0.01)
  power = numeric(length(t))
  for (i in seq_along(t)) {
    diff = qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
    min(abs(diff), na.rm = TRUE)
    power[i] = powerTemp[order(abs(diff))][1]
  }
  power
}


#' Get sample size
#' @inheritParams designSampleSize
#' @inheritParams .calculatePower
#' @param alpha significance level
#' @param delta difference between means (?)
#' @importFrom stats qnorm
#' @keywords internal
.getNumSample = function(desiredFC, power, alpha, delta, median_sigma_error, 
                         median_sigma_subject, median_sigma_run){
  z_alpha = qnorm(1 - alpha / 2)
  z_beta = qnorm(power)
  aa = (delta / (z_alpha + z_beta)) ^ 2
  numSample = round(2 * (median_sigma_error + median_sigma_subject + median_sigma_run) / aa, 0)
  numSample
}