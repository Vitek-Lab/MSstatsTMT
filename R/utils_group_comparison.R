#' @keywords internal
.checkGroupComparisonInput = function(input) {
  required_cols = c("Protein", "BioReplicate", "Abundance", "Run", 
                    "Channel", "Condition", "TechRepMixture", "Mixture")
  if (!all(required_cols %in% colnames(input))) {
    missing_cols = !(required_cols %in% colnames(input))
    missing_msg = paste(required_cols[missing_cols],
                        collapse = ", ")
    if (sum(missing_cols) == 1) {
      stop(paste("Please check the required input. ** columns :",
                 missing_msg, "is missed."))
    } else {
      stop(paste("Please check the required input. ** columns :",
                 missing_msg, ", are missed."))
    }
  }
  if (data.table::uniqueN(input$Condition) < 2){
    stop(paste("Please check the Condition column in annotation file.", 
               "There must be at least two conditions!"))
  }
  input
}

#' check whether pairwise comparison. If pairwise, generate a contrast matrix.
#' @importFrom stats coef
#' @keywords internal
#' @return a contrast matrix
.checkContrastMatrix = function(contrast_matrix) { 
  # TODO: use MSstatsdev::MSstatsContrastMatrix
  # TODO: add checking data.frame/validity in MSstatsContrastMatrix
  
  groups <- NULL
  ## check whether contrast.matrix is pairwise or matrix
  if (!(is.matrix(contrast.matrix) | is.data.frame(contrast.matrix))) {
    if(length(contrast.matrix) == 1){ # contrast.matrix is a character
      if(contrast.matrix != "pairwise"){ 
        stop("contrast.matrix must be 'pairwise' or a contrast matrix.")
      } else{ # create constrast matrix for pairwise comparison
        contrast.matrix = .makeContrast(groups)
      }
    } else{
      stop("contrast.matrix must be 'pairwise' or a contrast matrix.")
    }
  } else{ # contrast.matrix is a matrix or data frame
    contrast.matrix = as.matrix(contrast.matrix)
    
    if (!all(colnames(contrast.matrix) %in% groups)) {
      stop("Please check the contrast.matrix. Column names of contrast.matrix must be matched with conditions!")
    }
    
    if(!is.numeric(contrast.matrix)){
      stop("Please check the contrast.matrix. The elements of the contrast matrix must be all numeric!")
    }
  }
  ncomp = nrow(contrast.matrix)
}

#' @keywords internal
.fitModelTMT = function(single_protein, has_single_subject, has_techreps, 
                        has_biomixtures, has_single_run, has_Repeated_Measures) {
  if (has_single_subject) { # no biological variation
    if (has_techreps & has_biomixtures) { # multiple mixtures and tech MS runs
      fit = fit_Mix_TechRep_Group_model(single_protein)
      if (inherits(fit, "try-error")) { # full model is not applicable 
        fit = fit_Run_Group_model(single_protein) # fit the reduced model with only run effect
      }
      if (inherits(fit, "try-error")) { # second model not applicable - fit one-way anova model
        fit = fit_Group_model(single_protein) 
      }
    } else {
      if (has_techreps | has_biomixtures) { # multiple runs
        # fit the reduced model with only run effect
        fit = fit_Run_Group_model(single_protein)
        if (inherits(fit, "try-error")) { # run model not applicable - fit one-way anova model
          fit = fit_Group_model(single_protein) 
        }
      } else { # single run
        fit = fit_Group_model(single_protein) 
      }
    }
  } else { # the data has biological variation 
    if(has_Repeated_Measures){ # time course design
      if (has_biomixtures) { # multiple mixtures
        if (has_techreps) { # multiple tech MS runs
          
          ####### TODO: need to double check the full model ######
          fit = fit_Mix_TechRep_Group_Sub_model(single_protein) # fit the full model with mixture, techrep, subject effects
          if (!inherits(fit, "try-error")) { # full model is not applicable, fit the model with run and subject effects
            fit = fit_Run_Group_Sub_model(single_protein) 
          }
          if (inherits(fit, "try-error")) { # second model not applicable - fit model with only subject effect
            fit = fit_Group_Sub_model(single_protein) 
          }
          if (inherits(fit, "try-error")) { # subject model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        } else { # single MS run per mixture
          fit = fit_Group_Sub_model(single_protein) # fit the model with subject effect
          if (inherits(fit, "try-error")) { # subject model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        }
      } else { # single mixture
        if (has_techreps) { # multiple tech MS runs
          fit = fit_Run_Group_Sub_model(single_protein) # fit the model with run and subject effects
          if (inherits(fit, "try-error")) { # full model not applicable - fit model with only run effect
            fit = fit_Group_Sub_model(single_protein) 
          }
          if (inherits(fit, "try-error")) { # model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        } else { # single MS run
          fit = fit_Group_Sub_model(single_protein) 
          if (inherits(fit, "try-error")) { # model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        }
      }
    } else{  # group comparison design
      if (has_biomixtures) { # multiple mixtures
        if (has_techreps) { # multiple tech MS runs
          fit = fit_Mix_TechRep_Group_Sub_model(single_protein) # fit the full model with mixture, techrep, subject effects
          if (!inherits(fit, "try-error")) { # full model is not applicable
            fit = fit_Run_Group_Sub_model(single_protein) # fit the reduced model with run and subject effects
          }
          if (inherits(fit, "try-error")) { # full model not applicable - fit model with only run effect
            fit = fit_Run_Group_model(single_protein) 
          }
          if (inherits(fit, "try-error")) { # run model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        } else {
          fit = fit_Run_Group_model(single_protein) # fit the reduced model with only subject effect
          if (inherits(fit, "try-error")) { # subject model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        }
      } else { # single mixture
        if (has_techreps) { # multiple tech MS runs
          fit = fit_Run_Group_Sub_model(single_protein) # fit the reduced model with run and subject effects
          if (inherits(fit, "try-error")) { # model not applicable - fit one-way anova model
            fit = fit_Group_model(single_protein) 
          }
        } else { # single MS run per mixture
          fit = fit_Group_model(single_protein) 
        }
      }
    }
  }
  fit
}

#' @keywords internal
.addVarianceInformation = function(fitted_model, protein_name) {
  if (!inherits(fitted_model, "try-error")) {
    if (inherits(fitted_model, "lm")) { # single run case 
      av = anova(fitted_model)
      coefs = coef(fitted_model)
      s2_df = av["Residuals", "Df"]
      if(s2_df == 0) {
        s2 = 0
      } else {
        s2 = av["Residuals", "Mean Sq"] # use error variance for testing
      }
    } else { 
      av = anova(fitted_model)
      coefs = lme4::fixef(fitted_model)
      s2_df = av$DenDF
      s2 = av$'Mean Sq'/av$'F value'
    }
  } else {
    s2 = NA
    s2_df = NA
    coefs = NA
  }
  list(fitted_model = list(fitted_model),
       protein = protein_name,
       variance = s2,
       variance_df = s2_df,
       coefs = list(coefs)
  )
}


###############################################################
## make contrast matrix for pairwise comparisons
###############################################################
#' @keywords internal
.makeContrast = function(groups) {
  
  ncomp = length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast.matrix = matrix(rep(0, length(groups) * ncomp), ncol = length(groups))
  colnames(contrast.matrix) = groups
  
  count = 0
  contrast.matrix.rownames = NULL
  for(j in seq_len(length(groups)-1)){
    for(k in (j+1):length(groups)){
      
      count = count + 1
      # save row name
      contrast.matrix.rownames = c(contrast.matrix.rownames, paste(groups[j], groups[k], sep = "-"))
      # set constrast value
      contrast.matrix[count, groups[j]] = 1
      contrast.matrix[count, groups[k]] = -1
    }
  }
  rownames(contrast.matrix) = contrast.matrix.rownames
  
  return(contrast.matrix)
}

## check whether single subject per mixture and group
#' @keywords internal
.checkSingleSubject = function(annotation) {
  Subject <- NULL
  count_subjects = annotation[, .(NumSubjects = uniqueN(Subject)),
                              by = c("Group")]
  all(count_subjects$NumSubjects <= 1)
}


## check .checkTechReplicate
#' @keywords internal
.checkTechReplicate = function(annotation) {
  Run <- NULL
  count_runs = annotation[, .(NumRuns = uniqueN(Run)),
                              by = c("Mixture")]
  any(count_runs$NumRuns > 1)
}


## check whether there are multiple biological mixtures
#' @keywords internal
.checkMulBioMixture = function(annotation) {
  uniqueN(annotation$Mixture) > 1
}


# check whether there is only single run
#' @keywords internal
.checkSingleRun = function(annotation) {
  uniqueN(annotation$Run) == 1
}

# check whether the data has repeated measures
#' @keywords internal
.checkRepeatedMeasures = function(annotation) {
  Group <- NULL
  count_groups = annotation[, .(NumGroups = uniqueN(Group)),
                              by = c("Subject")]
  any(count_groups$NumGroups > 1)
}

## fit the full model with mixture, techrep and subject effects
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' multiple mixtures, multiple technical replicate runs per mixture and biological variation
fit_Mix_TechRep_Group_Sub_model = function(data) {
  suppressMessages(try(
    lmerTest::lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) +  # whole plot
                     Group + #subplot
                     (1|Subject), data = data), TRUE))
}

## fit the reduced model with run and subject effects
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' single mixture with multiple technical replicate runs
fit_Run_Group_Sub_model = function(data) {
  suppressMessages(try(
    lmerTest::lmer(Abundance ~ 1 + (1|Run) +  # whole plot
                     Group + #subplot
                     (1|Subject), data = data), TRUE))
}

## fit the reduced model with mixture and techrep effects
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures with multiple technical replicate runs
fit_Mix_TechRep_Group_model = function(data) {
  suppressMessages(try(
    lmerTest::lmer(
      Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group,
      data = data
    ), TRUE
  ))
}

## fit the reduced with only run effect
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures or multiple technical replicate runs
#' or if the data has multiple mixtures but single technical replicate MS run
fit_Run_Group_model = function(data) {
  suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1|Run) + Group,
                                      data = data), TRUE))
}

## fit one-way anova model
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has single run
fit_Group_model = function(data) {
  suppressMessages(try(lm(Abundance ~ 1 + Group, data = data), TRUE))
}

## fit subject model for repeated measures design
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit subject model if the data uses repeated measures design
fit_Group_Sub_model = function(data) {
  suppressMessages(try(
    lmerTest::lmer(
      Abundance ~ 1 + Group + (1|Subject),
      data = data
    ), TRUE
  ))
}


#' perform statistical inference for single protein and single contrast
#' @keywords internal
.handleSingleContrastTMT = function(contrast, fit, single_protein, coefs, 
                                    protein, groups, s2_posterior, rho, vss,
                                    df_prior, s2_df) {
  positive.groups = colnames(contrast)[contrast > 0]
  negative.groups = colnames(contrast)[contrast < 0]
  
  if (any(positive.groups %in% groups) & 
      any(negative.groups %in% groups)) {
    cm = .makeContrastSingleTMT(fit, contrast, single_protein, coefs)
    FC = (cm %*% coefs)[, 1]
    if (inherits(fit, "lm")) {
      se2.post = diag(t(cm) %*% summary(fit)$cov.unscaled %*% cm) * s2_posterior
      df.post = s2_df + df_prior
    } else {
      # Acknowlege: Tyler Bradshawthis contributed to this part of implementation
      vcov = fit@vcov_beta
      se2 = as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
      ## calculate posterior variance
      vcov.post = fit@pp$unsc() * s2_posterior
      se2.post = as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
      ## calculate posterior df
      g = .mygrad(function(x)  vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
      denom = try(t(g) %*% fit@vcov_varpar %*% g, silent=TRUE)
      if (inherits(denom, "try-error")) {
        df.post = s2_df + df_prior
      } else{
        df.post = 2*(se2)^2/denom + df_prior
      }
    }
    
    t = FC / sqrt(se2.post)[1, 1]
    p = 2*pt(-abs(t), df = df.post)
    se = sqrt(se2.post)[1, 1]
    DF = df.post[1, 1]
    if (s2_df == 0) {
      issue = "SingleMeasurePerCondition"
    } else{
      issue = NA
    }
  } else {
    result = .handleMissingConditionTMT(single_protein, contrast[1, ])
    FC = result[["logFC"]]
    p = NA
    se = NA
    DF = NA
    issue = result[["issue"]]
  }
  
  list(Protein = protein, comparison = row.names(contrast),
       log2FC = FC, pvalue = p, SE = se, DF = DF, issue = issue)
}


## check the reason for results with NA
#' @keywords internal
#' check the possible reason for untestable comparison
.handleMissingConditionTMT = function(input, contrast){
  groups = as.character(unique(input$Group)) 
  positive = names(contrast)[contrast > 0]
  negative = names(contrast)[contrast < 0]
  
  if(is.null(positive) | is.null(negative)){
    stop("Please check the contrast.matrix. 
         Each row must have both positive and negative values,
         and their sum must be 1!")
  } # TODO: do this check earlier
  
  if (any(positive %in% groups) & any(negative %in% groups)) {
    logFC = NA
    issue = "unfittableModel"
  } else {
    if (all(!positive %in% groups) & any(negative %in% groups)) {
      logFC = -Inf
      issue = "oneConditionMissing"
    } else {
      if (any(positive %in% groups) & all(!negative %in% groups)) {
        logFC = Inf
        issue = "oneConditionMissing"
      } else {
        logFC = NA
        issue = "completeMissing"
      }
    }
  }
  list(logFC = logFC, issue = issue)
}

#' Make a contrast
#' @keywords internal
#' @return a contrast vector
.makeContrastSingleTMT = function(fit, contrast, single_protein, coefs) {
  sub_groups = as.character(unique(single_protein$Group))
  positive = names(contrast)[contrast > 0]
  negative = names(contrast)[contrast < 0]
  
  if (!(all(positive %in% sub_groups) & all(negative %in% sub_groups))) {
    contrast_updated = unname(contrast[, sub_groups])
    contrast_updated = .normalizeContrastGroups(contrast_updated, "negative")
    contrast_updated = .normalizeContrastGroups(contrast_updated, "positive")
    contrast[] = 0
    contrast[1, sub_groups] = contrast_updated
  }
  
  if (inherits(fit, "lm")) {
    coef_name = names(coefs)
  } else {
    coef_name = names(coefs)
  }
  
  temp = coef_name[grep("Intercept", coef_name)]
  if (length(temp) == 0) {
    intercept_c = NULL
  } else {
    intercept_c = rep(0, length(temp))
    names(intercept_c) = temp
  }
  
  temp = coef_name[grep("Group", coef_name)]
  if (length(temp) == 0) {
    group_c = NULL
  } else {
    tempcontrast = contrast[1, sub_groups]
    group_c = tempcontrast[gsub("^Group", "", temp)] 
    names(group_c) = temp
  }
  
  new_contrast = c(intercept_c, group_c)
  if (inherits(fit, "lm")) {
    new_contrast = new_contrast[!is.na(coefs)]
  } else {
    new_contrast = new_contrast[!is.na(coefs)]
  }
  new_contrast
}

#' @keywords internal
.normalizeContrastGroups = function(contrast, type) {
  if (type == "negative") {
    condition = contrast < 0
  } else {
    condition = contrast > 0
  }
  new_values = contrast[condition]
  new_values = new_values * abs(1 / sum(new_values, na.rm = TRUE))
  contrast[condition] = new_values
  contrast
}