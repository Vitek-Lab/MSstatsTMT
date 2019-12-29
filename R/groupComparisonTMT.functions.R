#' @import lme4
#' @import statmod
#' @importFrom nlme fixed.effects
#' @importFrom dplyr group_by
#' @importFrom stats aggregate anova coef lm median medpolish model.matrix na.omit p.adjust pt t.test xtabs
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
      eb_fit <- limma::squeezeVar(fitted.models$s2, fitted.models$df)
        
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
    df.all <- fitted.models$df  # degree freedom
    lms <- fitted.models$model # linear models
    
    num.protein <- length(proteins)
    res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) ## store the inference results
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF", "issue")
    data$Group <- as.factor(data$Group) # make sure group is factor
    data$Run <- as.factor(data$Run)
    nrun <- length(unique(data$Run)) # check the number of MS runs in the data
    count <- 0
    for(i in seq_along(proteins)){
      message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
      
      ## get the data for protein i
      sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
      ## record the contrast matrix for each protein
      sub.contrast.matrix <- contrast.matrix
      
      sub_groups <- as.character(unique(sub_data$Group)) # groups in the sub data
      sub_groups <- sort(sub_groups) # sort the groups based on alphabetic order
      
      ## get the linear model for proteins[i]
      fit <- lms[[proteins[i]]]
      MSE <- s2.all[proteins[i]]
      df <- df.all[proteins[i]]
      
      if(!is.character(fit)){ ## check the model is fittable 
        # the protein is testable
        if(class(fit) == "lm"){# single run case 
          ## Get estimated fold change from mixed model
          coeff <- coef(fit)
          coeff[-1] <- coeff[-1] + coeff[1]
          
          ## Find the group name for baseline
          names(coeff) <- gsub("Group", "", names(coeff))
          names(coeff)[1] <- setdiff(as.character(sub_groups), names(coeff))
          
          s2.post <- (s2.prior * df.prior + MSE * df)/(df.prior + df)
          df.post <- df + df.prior
          
        } else{ # multiple run case
          ## Get estimated fold change from mixed model
          coeff <- fixed.effects(fit$mixed)
          coeff[-1] <- coeff[-1] + coeff[1]
          
          # Find the group name for baseline
          names(coeff) <- gsub("Group", "", names(coeff))
          names(coeff)[1] <- setdiff(as.character(sub_groups), names(coeff))
          
          s2.post <- (s2.prior * df.prior + MSE * df)/(df.prior + df)
          df.post <- df + df.prior
          
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
            
            # if some groups not exist in the protein data
            if(!(all(positive.groups %in% sub_groups) & 
                 all(negative.groups %in% sub_groups))){
              ## tune the coefficients of positive groups so that their summation is 1
              temp <- sub.contrast.matrix[j,sub_groups][sub.contrast.matrix[j,sub_groups] > 0]
              temp <- temp*(1/sum(temp, na.rm = TRUE))
              sub.contrast.matrix[j,sub_groups][sub.contrast.matrix[j,sub_groups] > 0] <- temp
              
              ## tune the coefficients of positive groups so that their summation is 1
              temp2 <- sub.contrast.matrix[j,sub_groups][sub.contrast.matrix[j,sub_groups] < 0]
              temp2 <- temp2*abs(1/sum(temp2, na.rm = TRUE))
              sub.contrast.matrix[j,sub_groups][sub.contrast.matrix[j,sub_groups] < 0] <- temp2
              
              ## set the coefficients of non-existing groups to zero
              sub.contrast.matrix[j,setdiff(colnames(sub.contrast.matrix), sub_groups)] <- 0
            }
            
            ## calculate the size of each group
            group_df <- sub_data %>% group_by(Group) %>% dplyr::summarise(n = sum(!is.na(Abundance)))
            group_df$Group <- as.character(group_df$Group)
            
            ## variance of diff
            variance <- s2.post * sum((1/group_df$n) * ((sub.contrast.matrix[j, group_df$Group])^2))
            FC <- sum(coeff*(sub.contrast.matrix[j, names(coeff)])) # fold change
            res[count, "log2FC"] <- FC
            
            ## Calculate the t statistic
            t <- FC/sqrt(variance)
            
            ## calculate p-value
            p <- 2*pt(-abs(t), df = df.post)
            res[count, "pvalue"] <- p
            
            ## SE
            res[count, "SE"] <- sqrt(variance)
            res[count, "DF"] <- df.post
            res[count, "issue"] <- NA
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
    } # for each protein
    
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
