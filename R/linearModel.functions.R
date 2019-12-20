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

.checkSingleRun <-  function(annotation) {

    temp <- unique(annotation[, "Run"])
    temp <- as.vector(as.matrix(temp))

    return(length(temp)==1)
}

#############################################
## fit the full model with mixture, techrep and subject effects
#############################################
#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' multiple mixtures, multiple technical replicate runs per mixture and biological variation
fit_full_model <- function(data) {
  
  fit.mixed <- suppressMessages(try(lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) +  # whole plot
                                           Group + #subplot
                                           (1|Subject:Group:Mixture), data = data), TRUE))
  
  fit.fixed <- suppressMessages(try(lm(Abundance ~ 1 + Mixture + Mixture:TechRepMixture +  # whole plot
                                         Group +
                                         Subject:Group:Mixture, data = data), TRUE))
  
  if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
    return(list(fixed = fit.fixed, mixed = fit.mixed, subject = "Subject:Group:Mixture"))
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
}

#############################################
## fit the reduced model with run and subject effects
#############################################
#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has 
#' single mixture with multiple technical replicate runs
fit_reduced_model_techrep <- function(data) {
  
  fit.mixed <- suppressMessages(try(lmer(Abundance ~ 1 + (1|Run) +  # whole plot
                                           Group + #subplot
                                           (1|Subject:Group), data = data), TRUE))
  
  fit.fixed <- suppressMessages(try(lm(Abundance ~ 1 + Run +  # whole plot
                                         Group +
                                         Subject:Group, data = data), TRUE))
  
  if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
    return(list(fixed = fit.fixed, mixed = fit.mixed, subject = "Subject:Group"))
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
  
}

#############################################
## fit the reduced model with mixture and techrep effects
#############################################
#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures with multiple technical replicate runs
fit_full_model_spikedin <- function(data) {
  
  fit.mixed <- suppressMessages(try(lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) 
                                         + Group, data = data), TRUE))
  fit.fixed <- suppressMessages(try(lm(Abundance ~ 1 + Mixture + Mixture:TechRepMixture 
                                       + Group, data = data), TRUE))
  
  if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
    return(list(fixed = fit.fixed, mixed = fit.mixed, subject = "None"))
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
  
}

#############################################
## fit the reduced with only run effect
#############################################
#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures or multiple technical replicate runs
#' or if the data has multiple mixtures but single technical replicate MS run
fit_reduced_model_mulrun <- function(data) {
  
  fit.mixed <- suppressMessages(try(lmer(Abundance ~ 1 + (1|Run) + Group, data = data), TRUE))
  fit.fixed <- suppressMessages(try(lm(Abundance ~ 1 + Run + Group, data = data), TRUE))
  
  if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
    return(list(fixed = fit.fixed, mixed = fit.mixed, subject = "None"))
  } else{ # if the parameters are not estimable, return null
    return(NULL)
  }
  
}

#############################################
## fit one-way anova model
#############################################
#' @import lme4
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
#' @import lme4
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
  df.all <- NULL # degree freedom
  pro.all <- NULL # testable proteins
  ## do inference for each protein individually
  count = 0 # count the testable proteins
  # issue <- NA
  for(i in seq_along(proteins)) {
    
    message(paste("Model fitting for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
    TEST <-  FALSE
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    
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
      if(class(fit) == "lm"){# single run case 
        ## Estimate the group variance from fixed model
        av <- anova(fit)
        # use error variance for testing
        MSE <- av["Residuals", "Mean Sq"]
        df <- av["Residuals", "Df"]
        TEST <-  TRUE
        
      } else{ ## fit linear mixed model
        if(fit$subject=="None"){ # no technical replicates
          # Estimate the group variance and df
          varcomp <- as.data.frame(VarCorr(fit$mixed))
          # use error variance for testing
          MSE <- varcomp[varcomp$grp == "Residual", "vcov"] 
          av <- anova(fit$fixed)
          df <- av["Residuals", "Df"] # degree of freedom
          TEST <-  TRUE
          
        } else{ # multiple technical replicates
          if(fit$subject=="Subject:Group:Mixture"){ # multiple biological mixtures
            # Estimate the group variance and df
            varcomp <- as.data.frame(VarCorr(fit$mixed))
            av <- anova(fit$fixed)
            if(any(grepl("Subject",  rownames(av)))){
              # use subject variance for testing
              MSE <- varcomp[varcomp$grp == "Subject:Group:Mixture", "vcov"]
              df <- av[rownames(av)[grepl("Subject",  rownames(av))], "Df"] # degree of freedom
              TEST <-  TRUE
              
            }
          } else{ # single biological mixture
            # Estimate the group variance and df
            varcomp <- as.data.frame(VarCorr(fit$mixed))
            av <- anova(fit$fixed)
            if(any(grepl("Subject",  rownames(av)))){
              # use subject variance for testing
              MSE <- varcomp[varcomp$grp == "Subject:Group", "vcov"]
              df <- av[rownames(av)[grepl("Subject",  rownames(av))], "Df"] # degree of freedom
              TEST <-  TRUE
              
            }
          }
        }
      }
      
      if(TEST){ # make sure all the parameters are estimable
        count <- count + 1
        linear.models[[count]] <- fit 
        pro.all <-  c(pro.all, proteins[i])
        s2.all <- c(s2.all, MSE)
        df.all <- c(df.all, df)
        
      } else{ # the standard error is not estimable
        linear.models[[proteins[i]]] <- NA 
        pro.all <-  c(pro.all, proteins[i])
        s2.all <- c(s2.all, NA)
        df.all <- c(df.all, NA)
        
      }
    } else{ # the model is not fittble
      linear.models[[proteins[i]]] <- "unfittable" 
      pro.all <-  c(pro.all, proteins[i])
      s2.all <- c(s2.all, NA)
      df.all <- c(df.all, NA)
      
    }
  } # for each protein
  
  return(list(protein = pro.all, model = linear.models, s2 = s2.all, df = df.all))
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