#' @import statmod
#' @importFrom limma squeezeVar
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom dplyr %>% group_by filter
#' @importFrom stats aggregate anova coef lm median medpolish model.matrix p.adjust pt t.test xtabs
#' @keywords internal

.proposed.model <- function(data,
                            moderated = TRUE,
                            contrast.matrix = "pairwise",
                            adj.method = "BH") {

    Abundance <- Group <- Protein <- NULL

    groups <- as.character(unique(data$Group)) # groups
    if(length(groups) < 2){
        stop("Please check the Condition column in annotation file. There must be at least two conditions!")
    }
    
    ## contrast matrix can be matrix or character vector.
    if(is.matrix(contrast.matrix)){
      # comparison come from contrast.matrix
      if (!all(colnames(contrast.matrix) %in% groups)) {
        stop("Please check the contrast.matrix. Column names of contrast.matrix must be matched with conditions!")
      }
    } else{  
      # create constrast matrix for pairwise comparison
      contrast.matrix <- .makeContrast(groups)
    }
    ncomp <- nrow(contrast.matrix)
   
    ###################### start fitting linear model ######################
    ## fit the linear model for each protein
    fitted.models <- .linear.model.fitting(data)
    
    ## perform empirical bayes moderation
    if(moderated){ ## moderated t statistic
      ## Estimate the prior variance and degree freedom
      
      eb_input_s2 <- fitted.models$s2[fitted.models$s2_df!=0]
      eb_input_df <- fitted.models$s2_df[fitted.models$s2_df!=0]
      
      eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
        
      if(is.infinite(eb_fit$df.prior)){
        df.prior = 0
        s2.prior = 0
      } else{
        df.prior = eb_fit$df.prior
        s2.prior = eb_fit$var.prior
      }
    } else { ## ordinary t statistic
      s2.prior <- 0
      df.prior <- 0
    }

    # extract the linear model fitting results
    proteins <- fitted.models$protein # proteins
    s2.all <- fitted.models$s2  # group variance
    s2_df.all <- fitted.models$s2_df  # degree freedom of s2
    lms <- fitted.models$model # linear models
    coeff.all <- fitted.models$coeff # coefficients
    
    num.protein <- length(proteins)
    res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) ## store the inference results
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF", "issue")
    data$Group <- as.factor(data$Group) # make sure group is factor
    data$Run <- as.factor(data$Run)
    nrun <- length(unique(data$Run)) # check the number of MS runs in the data
    count <- 0
    
    message(paste0("Testing for ", num.protein , " proteins:"))
    pb <- txtProgressBar(max=num.protein, style = 3)
    
    for(i in seq_along(proteins)){
      #message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
      
      ## get the data for protein i
      sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
      ## record the contrast matrix for each protein
      sub.contrast.matrix <- contrast.matrix
      
      sub_groups <- as.character(unique(sub_data$Group))
      sub_groups <- sort(sub_groups) # sort the groups based on alphabetic order

      ## get the linear model for proteins[i]
      fit <- lms[[proteins[i]]]
      s2 <- s2.all[proteins[i]]
      s2_df <- s2_df.all[proteins[i]]
      coeff <- coeff.all[[proteins[i]]]
      
      if(!inherits(fit, "try-error")){ ## check the model is fittable 
        
        s2.post <- (s2.prior * df.prior + s2 * s2_df)/(df.prior + s2_df)
        
        if(!inherits(fit, "lm")){
          ## prepare for df calculation
          rho <- list() ## environment containing info about model
          rho <- suppressMessages(.rhoInit(rho, fit, TRUE)) ## save lmer outcome in rho envir variable
          # rho$A <- .calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma
          vss <- .vcovLThetaL(fit)
        }
        
        ## Compare one specific contrast
        # perform testing for required contrasts
        for(j in seq_len(nrow(sub.contrast.matrix))){
          count <- count + 1
          res[count, "Protein"] <- proteins[i] ## protein names
          res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison

          # groups with positive coefficients
          positive.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j,]>0]
          # groups with negative coefficients
          negative.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j,]<0]
          # make sure at least one group from each side of the contrast exist
          if(any(positive.groups %in% sub_groups) & 
             any(negative.groups %in% sub_groups)){
            
            contrast.matrix.single <- as.vector(sub.contrast.matrix[j,])
            names(contrast.matrix.single) <- colnames(sub.contrast.matrix)
            
            cm <- .make.contrast.single(fit, contrast.matrix.single, sub_data)
            
            ## logFC
            FC <- (cm%*%coeff)[,1]
            
            ## variance and df
            if(inherits(fit, "lm")){
 
              se2.post <- diag(t(cm) %*% summary(fit)$cov.unscaled %*% cm)*s2.post
              df.post <- s2_df + df.prior
              
            } else{
              
              # Acknowlege: Tyler Bradshawthis contributed to this part of implementation
              vcov <- fit@vcov_beta
              se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
              
              ## calculate posterior variance
              vcov.post <- fit@pp$unsc() * s2.post
              se2.post <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
              
              ## calculate posterior df
              g <- .mygrad(function(x)  vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
              denom <- try(t(g) %*% fit@vcov_varpar %*% g, silent=TRUE)
              if(inherits(denom, "try-error")) {
                
                df.post <- s2_df + df.prior
              } else{
                
                df.post <- 2*(se2)^2/denom + df.prior
              }
            }
            
            ## calculate the t statistic
            t <- FC/sqrt(se2.post)
            
            ## calculate p-value
            p <- 2*pt(-abs(t), df = df.post)
            res[count, "pvalue"] <- p
            
            ## save testing results
            res[count, "log2FC"] <- FC
            res[count, "SE"] <- sqrt(se2.post)
            res[count, "DF"] <- df.post
            
            if(s2_df == 0){
              res[count, "issue"] <- "SingleMeasurePerCondition"
              
            } else{
              res[count, "issue"] <- NA
              
            }
            
          } else{
            # at least one condition is missing
            out <- .issue.checking(data = sub_data, 
                                   contrast.matrix = sub.contrast.matrix[j,])

            res[count, "log2FC"] <- out$logFC
            res[count, "pvalue"] <- NA
            res[count, "SE"] <- NA
            res[count, "DF"] <- NA
            res[count, "issue"] <- out$issue
            
          }
        } # for constrast matrix
      } else {
        # very few measurements so that the model is unfittable
        for (j in 1:nrow(sub.contrast.matrix)) {
          count <- count + 1
          res[count, "Protein"] <- proteins[i] ## protein names
          res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison
          
          out <- .issue.checking(data = sub_data, 
                                 contrast.matrix = sub.contrast.matrix[j,])
          
          res[count, "log2FC"] <- out$logFC
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          res[count, "issue"] <- out$issue
          
        } # end loop for comparison    
      } # if the linear model is fittable
      
      ## progress
      setTxtProgressBar(pb, i)
      
    } # for each protein
    
    close(pb)
    
    res <- as.data.frame(res[seq_len(count),])
    res$Protein <- as.factor(res$Protein)
    res$log2FC <- as.numeric(as.character(res$log2FC))
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$adjusted.pvalue <- NA
    comps <- unique(res$Comparison)

    ## Adjust multiple tests for each comparison
    for(i in seq_along(comps)){
        res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
    }

    res <- res[, c("Protein",
                   "Comparison",
                   "log2FC",
                   "SE",
                   "DF",
                   "pvalue",
                   "adjusted.pvalue",
                   "issue")]
    return(res)
}
