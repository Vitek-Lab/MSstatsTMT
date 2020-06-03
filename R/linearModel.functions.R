###############################################################
## make contrast matrix for pairwise comparisons
###############################################################
#' @keywords internal
.makeContrast <- function(groups) {
  
  ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast.matrix <- matrix(rep(0, length(groups) * ncomp), ncol = length(groups))
  colnames(contrast.matrix) <- groups
  
  count <- 0
  contrast.matrix.rownames <- NULL
  for(j in seq_len(length(groups)-1)){
    for(k in (j+1):length(groups)){
      
      count <- count + 1
      # save row name
      contrast.matrix.rownames <- c(contrast.matrix.rownames, paste(groups[j], groups[k], sep = "-"))
      # set constrast value
      contrast.matrix[count, groups[j]] <- 1
      contrast.matrix[count, groups[k]] <- -1
    }
  }
  rownames(contrast.matrix) <- contrast.matrix.rownames
  
  return(contrast.matrix)
}

###############################################################
## check single subject within each condition in each mixture
###############################################################
#' @keywords internal
.checkSingleSubject <- function(annotation) {

    temp <- unique(annotation[, c("Mixture", "Group", "Subject")])
    temp$Group <- factor(temp$Group)
    temp$Mixture <- factor(temp$Mixture)
    temp1 <- xtabs(~ Mixture+Group, data=temp)
    singleSubject <- all(temp1 <= "1")

    return(singleSubject)
}

#############################################
## check .checkTechReplicate
#############################################
#' @keywords internal
.checkTechReplicate <-  function(annotation) {

    temp <- unique(annotation[, c("Mixture", "Run")])
    temp$Mixture <- factor(temp$Mixture)
    temp1 <- xtabs(~ Mixture, data=temp)
    TechReplicate <- all(temp1 != "1")

    return(TechReplicate)
}

#############################################
## check whether there are multiple biological mixtures
#############################################
#' @keywords internal
.checkMulBioMixture <-  function(annotation) {

    temp <- unique(annotation[, "Mixture"])
    temp <- as.vector(as.matrix(temp))

    return(length(temp)>1)
}

#############################################
## check whether there is only single run
#############################################
#' @keywords internal
.checkSingleRun <-  function(annotation) {

    temp <- unique(annotation[, "Run"])
    temp <- as.vector(as.matrix(temp))

    return(length(temp)==1)
}

#############################################
## fit the full model with mixture, techrep and subject effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' multiple mixtures, multiple technical replicate runs per mixture and biological variation
fit_full_model <- function(data) {
  
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) +  # whole plot
                                           Group + #subplot
                                           (1|Subject:Group:Mixture), data = data), TRUE))
  
  if(!inherits(fit, "try-error")){
    return(fit)
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
}

#############################################
## fit the reduced model with run and subject effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' single mixture with multiple technical replicate runs
fit_reduced_model_techrep <- function(data) {
  
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1|Run) +  # whole plot
                                           Group + #subplot
                                           (1|Subject:Group), data = data), TRUE))
  
  if(!inherits(fit, "try-error")){
    return(fit)
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
}

#############################################
## fit the reduced model with mixture and techrep effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures with multiple technical replicate runs
fit_full_model_spikedin <- function(data) {
  
  fit  <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) 
                                         + Group, data = data), TRUE))
  
  if(!inherits(fit, "try-error")){
    return(fit)
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
}

#############################################
## fit the reduced with only run effect
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures or multiple technical replicate runs
#' or if the data has multiple mixtures but single technical replicate MS run
fit_reduced_model_mulrun <- function(data) {
  
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1|Run) + Group, data = data), TRUE))
  
  if(!inherits(fit, "try-error")){
    return(fit)
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
}

#############################################
## fit one-way anova model
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has single run
fit_reduced_model_onerun <- function(data) {
  
  fit <- suppressMessages(try(lm(Abundance ~ 1 + Group, data = data), TRUE))
  
  if(!inherits(fit, "try-error")){
    return(fit)
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
  
}

#############################################
## fit the proper linear model for each protein
#############################################
#' @importFrom lme4 fixef
#' @import lmerTest
#' @importFrom stats vcov
#' @importFrom dplyr filter
#' @keywords internal
#' fit the proper linear model for each protein
.linear.model.fitting <- function(data){
  
  Abundance <- Group <- Protein <- NULL
  
  data$Protein <- as.character(data$Protein) ## make sure protein names are character
  proteins <- as.character(unique(data$Protein)) ## proteins
  num.protein <- length(proteins)
  linear.models <- list() # linear models
  s2.all <- NULL # sigma^2
  s2_df.all <- NULL # degree freedom of sigma^2
  pro.all <- NULL # testable proteins
  coeff.all <- list() # coefficients
  ## do inference for each protein individually
  for(i in seq_along(proteins)) {
    
    message(paste("Model fitting for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    # sub_groups <- as.character(unique(sub_data$Group))
    # if(length(sub_groups) == 1){
    #   stop("Only one condition!")
    # }
    ## Record the annotation information
    sub_annot <- unique(sub_data[, c('Run', 'Channel', 'Subject',
                                     'Group', 'Mixture', 'TechRepMixture')])
    
    ## check the experimental design
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_bioMixture <- .checkMulBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)
    
    if(sub_singleSubject){ # no biological variation within each condition and mixture
      if(sub_TechReplicate & sub_bioMixture){ # multiple mixtures and technical replicates
        # fit the full model with mixture and techrep effects for spiked-in data
        fit <- fit_full_model_spikedin(sub_data)
        
        if(is.null(fit)){ # full model is not applicable 
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data)
          
        } 
        
        if(is.null(fit)){ # the second model is not applicable
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data) 
          
        }
      } else{
        if(sub_TechReplicate | sub_bioMixture){ # multiple mixtures or multiple technical replicates
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data)
          
          if(is.null(fit)){ # the second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data) 
            
          }
        } else{ # single run case
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data) 
          
        }
      }
    } else{ # biological variation exists within each condition and mixture
      if (sub_bioMixture) {  # multiple biological mixtures
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the full model with mixture, techrep, subject effects
          fit <- fit_full_model(sub_data) 
          
          if(is.null(fit)){ # full model is not applicable
            # fit the reduced model with run and subject effects
            fit <- fit_reduced_model_techrep(sub_data) 
          }
          
          if(is.null(fit)){ # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data) 
          }
          
        } else { # single technical replicate MS run
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data) 
          
          if(is.null(fit)){ # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data) 
          }
          
        }
      } else { # single biological mixture
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the reduced model with run and subject effects
          fit <- fit_reduced_model_techrep(sub_data)
          
          if(is.null(fit)){ # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data) 
          }
          
        } else { # single run
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data) 
          
        } # single technical replicate MS run
      } # single biological mixture
    } # biological variation
    
    ## estimate variance and df from linear models
    if(!is.null(fit)){ # the model is fittable
      if(inherits(fit, "lm")){# single run case 
        ## Estimate the coeff from fixed model
        av <- anova(fit)
        coeff <- coef(fit)
        
        s2_df <- av["Residuals", "Df"]
        
        if(s2_df == 0){
          s2 <- 0
          
        } else{
          # use error variance for testing
          s2 <- av["Residuals", "Mean Sq"]
          
        }

        linear.models[[proteins[i]]] <- list(model = fit)
        
      } else{ 
        ## Estimate the coeff from lmerTest model
        rho <- list() ## environment containing info about model
        rho <- .rhoInit(rho, fit, TRUE) ## save lmer outcome in rho envir variable
        rho$A <- .calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma
        
        av <- anova(rho$model)
        coeff <- lme4::fixef(rho$model)
        s2_df <- av$DenDF
        s2 <- av$'Mean Sq'/av$'F value'
        
        linear.models[[proteins[i]]] <- rho 
      }

      pro.all <-  c(pro.all, proteins[i])
      s2.all <- c(s2.all, s2)
      s2_df.all <- c(s2_df.all, s2_df)
      coeff.all[[proteins[i]]] <- coeff
      
    } else{ # the model is not fittble
      # message(proteins[i], " is untestable due to no enough measurements.")
      linear.models[[proteins[i]]] <- "unfittable" 
      pro.all <-  c(pro.all, proteins[i])
      s2.all <- c(s2.all, NA)
      s2_df.all <- c(s2_df.all, NA)
      coeff.all[[proteins[i]]] <- NA
      
    }
  } # for each protein
  names(s2.all) <- proteins
  names(s2_df.all) <- proteins
  
  return(list(protein = pro.all, 
              model = linear.models, 
              s2 = s2.all, 
              s2_df = s2_df.all, 
              coeff = coeff.all))
}

#############################################
## check the reason for results with NA
#############################################
#' @keywords internal
#' check the possible reason for untestable comparison
.issue.checking <- function(data, 
                            contrast.matrix){
  
  ## choose each comparison
  contrast.matrix.sub <- contrast.matrix
  
  # groups in the sub data
  sub_groups <- as.character(unique(data$Group)) 
  # groups with positive coefficients
  positive.groups <- names(contrast.matrix.sub)[contrast.matrix.sub>0]
  # groups with negative coefficients
  negative.groups <- names(contrast.matrix.sub)[contrast.matrix.sub<0]
  
  if(is.null(positive.groups) | is.null(negative.groups)){
    stop("Please check the contrast.matrix. 
         Each row must have both positive and negative values,
         and their sum must be 1!")
  }
  
  if(any(positive.groups %in% sub_groups) & 
     any(negative.groups %in% sub_groups)){
    logFC = NA
    issue = "unfittableModel"
    
  } else{
    # more than one condition
    if(all(!positive.groups %in% sub_groups) & 
       any(negative.groups %in% sub_groups)){
      logFC = (-Inf)
      issue = "oneConditionMissing"
      
    } else{
      if(any(positive.groups %in% sub_groups) & 
         all(!negative.groups %in% sub_groups)){
        logFC = Inf
        issue = "oneConditionMissing"
        
      } else{
        logFC = NA
        issue = "completeMissing"
        
      }
    }
  }
  
  return(list(logFC = logFC, issue = issue))
}

#############################################
## make constrast
#############################################
#	MSstats
#' @importFrom stats coef
#' @importFrom lme4 fixef
#' @keywords internal
.make.contrast.single <- function(fit, contrast, sub_data) {

  ## when there are some groups which are all missing
  sub_groups <- as.character(levels(sub_data[, c("Group")]))
  
  # groups with positive coefficients
  positive.groups <- names(contrast)[contrast>0]
  # groups with negative coefficients
  negative.groups <- names(contrast)[contrast<0]
  
  # if some groups not exist in the protein data
  if(!(all(positive.groups %in% sub_groups) & 
       all(negative.groups %in% sub_groups))){
    
    contrast.single <- contrast[sub_groups]
    
    ## tune the coefficients of positive groups so that their summation is 1
    temp <- contrast.single[contrast.single > 0]
    temp <- temp*(1/sum(temp, na.rm = TRUE))
    contrast.single[contrast.single > 0] <- temp
    
    ## tune the coefficients of positive groups so that their summation is 1
    temp2 <- contrast.single[contrast.single < 0]
    temp2 <- temp2*abs(1/sum(temp2, na.rm = TRUE))
    contrast.single[contrast.single < 0] <- temp2
    
    ## set the coefficients of non-existing groups to zero
    contrast[] <- 0
    contrast[sub_groups] <- contrast.single
  }

  if (inherits(fit, "lm")) {
    coef_name <- names(stats::coef(fit))
  } else {
    coef_name <- names(lme4::fixef(fit))
  }
  
  ## intercept
  temp <- coef_name[grep("Intercept", coef_name)]
  intercept_c <- rep(0, length(temp))
  names(intercept_c) <- temp
  if (length(temp) == 0) {
    intercept_c <- NULL
  }
  
  ## group
  temp <- coef_name[grep("Group", coef_name)]
  tempcontrast <- contrast[sub_groups]
  group_c <- tempcontrast[gsub("Group", "", temp)] 
  names(group_c) <- temp
  if (length(temp) == 0) {
    group_c<-NULL
  }
  
  ## combine all
  newcontrast <- c(intercept_c, group_c)
  if(inherits(fit, "lm")) {
    contrast1 <- newcontrast[!is.na(stats::coef(fit))]
  } else {
    contrast1 <- newcontrast[!is.na(lme4::fixef(fit))]
  }
  
  return(contrast1)
}


# retired fuction (2020.04.13)
# #############################################
# ## get the unscaled covariance matrix
# #############################################
# #	statOmics, MSqRob hurdle model
# #	Created 2020
# .getVcovUnscaled <- function(model){
# 
#   if(inherits(fixed.model, "lm")){
#     vcov <- summary(model)$cov.unscaled
# 
#   } else{
#     p <- ncol(lme4::getME(model,"X"))
#     q <- nrow(lme4::getME(model,"Zt"))
#     Ct <- rbind2(t(lme4::getME(model,"X")),lme4::getME(model,"Zt"))
#     Ginv <- Matrix::solve(Matrix::tcrossprod(lme4::getME(model,"Lambda"))+Matrix::Diagonal(q,1e-18))
#     vcovInv <- Matrix::tcrossprod(Ct)
#     vcovInv[((p+1):(q+p)),((p+1):(q+p))] <- vcovInv[((p+1):(q+p)),((p+1):(q+p))]+Ginv
# 
#     #remove rows with only zeros, making it uninvertible
#     defined <- rowSums(as.matrix(vcovInv==0))!=ncol(vcovInv)
#     defined[is.na(defined)] <- TRUE
#     vcovInv <- vcovInv[defined, defined, drop=FALSE]
# 
#     #Estimated variance-covariance matrix vcov:
#     vcov <- tryCatch(as.matrix(Matrix::solve(vcovInv)), error=function(e){
#       return(vcovInv*NA)
#     })
# 
#     rownames(vcov) <- colnames(vcovInv)
#     colnames(vcov) <- rownames(vcovInv)
#   }
# 
#   return(vcov)
# }