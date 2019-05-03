#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has multiple mixtures, multiple technical replicate runs per mixture and biological variation
fit_full_model <- function(data) {

    fit.mixed <- try(lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) +  # whole plot
                              Group + (1|Group:Mixture) + (1|Group:Mixture:TechRepMixture) + #subplot
                              (1|Subject:Group:Mixture), data = data), TRUE)

    fit.fixed <- try(lm(Abundance ~ 1 + Mixture + Mixture:TechRepMixture +  # whole plot
                            Group + Group:Mixture + Group:Mixture:TechRepMixture +
                            Subject:Group:Mixture, data = data), TRUE)

    if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
        return(list(fixed = fit.fixed, mixed = fit.mixed))
    } else{ # if the parameters are not estimable, return null
        return(NULL)
    }
}


#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has multiple mixtures without technical replicate runs
fit_reduced_model_biomixture <- function(data) {

    fit.mixed <- try(lmer(Abundance ~ 1 + (1|Run) + Group, data = data), TRUE)
    fit.fixed <- try(lm(Abundance ~ 1 + Run + Group, data = data), TRUE)

    if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
        return(list(fixed = fit.fixed, mixed = fit.mixed))
    } else{ # if the parameters are not estimable, return null
        return(NULL)
    }

}

#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has single mixture with multiple technical replicate runs
fit_reduced_model_techrep <- function(data) {

    fit.mixed <- try(lmer(Abundance ~ 1 + (1|Run) +  # whole plot
                              Group + (1|Subject:Condition), data = data), TRUE)

    fit.fixed <- try(lm(Abundance ~ 1 + Run +  # whole plot
                              Group + Subject:Condition, data = data), TRUE)

    if((!inherits(fit.mixed, "try-error")) & (!inherits(fit.fixed, "try-error"))){
        return(list(fixed = fit.fixed, mixed = fit.mixed))
    } else{ # if the parameters are not estimable, return null
        return(NULL)
    }

}

#' @import lme4
#' @keywords internal
#' fit the whole plot and subplot model if the data has single run
fit_reduced_model_onerun <- function(data) {

    fit <- lm(Abundance ~ 1 + Group, data = data)

    if(!inherits(fit, "try-error")){
        return(fit)
    } else{ # if the parameters are not estimable, return null
        return(NULL)
    }

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
    singleSubject <- all(temp1 == "1")

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


#' @import limma
#' @import statmod
#' @keywords internal
#' estimate the prior variance and df
.estimate.prior.var <- function(data,
                                bioMixture,
                                TechReplicate,
                                singleSubject,
                                singleRun){

    Subject <- Abundance <- Protein <- NULL

    proteins <- as.character(unique(data$Protein)) ## proteins
    SIGMA_all <- NULL
    DF_all <- NULL
    for(i in 1:length(proteins)) {
        sub_data <- data %>% dplyr::filter(Protein == proteins[i])
        ## fit the linear model based on experimental design
        if (bioMixture) {  ## multiple biological mixtures
            if (!TechReplicate | singleSubject) {
                fit <- fit_reduced_model_biomixture(sub_data) # fit linear model
                varcomp <- as.data.frame(VarCorr(fit$mixed))
                sigma <- varcomp[varcomp$grp == "Residual", "vcov"] # residual variance
                av <- anova(fit$fixed)
                df <- av["Residuals", "Df"] # degree of freedom
                SIGMA_all <- c(SIGMA_all, sigma)
                DF_all <- c(DF_all, df)

            } else {
                fit <- fit_full_model(sub_data) # fit linear model
                varcomp <- as.data.frame(VarCorr(fit$mixed))
                sigma <- varcomp[varcomp$grp == "Residual", "vcov"] # residual variance
                av <- anova(fit$fixed)
                df <- av["Residuals", "Df"] # degree of freedom
                SIGMA_all <- c(SIGMA_all, sigma)
                DF_all <- c(DF_all, df)

            }
        } else { ## single biological mixture
            if (singleSubject) { # each condition has one subject
                fit <- fit_reduced_model_biomixture(sub_data) # fit linear model
                varcomp <- as.data.frame(VarCorr(fit$mixed))
                sigma <- varcomp[varcomp$grp == "Residual", "vcov"] # residual variance
                av <- anova(fit$fixed)
                df <- av["Residuals", "Df"] # degree of freedom
                SIGMA_all <- c(SIGMA_all, sigma)
                DF_all <- c(DF_all, df)

            } else { ## no single subject
                if (!TechReplicate) { # no technical replicates
                    fit <- fit_reduced_model_onerun(sub_data) # fit linear model
                    av <- anova(fit)
                    sigma <- av["Residuals", "Mean Sq"] # residual variance
                    df <- av["Residuals", "Df"] # degree of freedom
                    SIGMA_all <- c(SIGMA_all, sigma)
                    DF_all <- c(DF_all, df)

                } else { # with technical replicates
                    fit <- fit_reduced_model_techrep(sub_data) # fit linear model
                    varcomp <- as.data.frame(VarCorr(fit$mixed))
                    sigma <- varcomp[varcomp$grp == "Residual", "vcov"] # residual variance
                    av <- anova(fit$fixed)
                    df <- av["Residuals", "Df"] # degree of freedom
                    SIGMA_all <- c(SIGMA_all, sigma)
                    DF_all <- c(DF_all, df)
                }
            }
        } ## single biological mixture
    } ## for all the proteins

    ## Squeeze the residual variances together by computing
    ## empirical Bayes posterior means.
    ## require pacakage statmod
    eb_fit <- limma::squeezeVar(SIGMA_all, DF_all)

    if(is.infinite(eb_fit$df.prior)){
        return(list(df.prior = 0,
                    s2.prior = 0))
    } else{
        return(list(df.prior = eb_fit$df.prior,
                    s2.prior = eb_fit$var.prior))
    }
}
