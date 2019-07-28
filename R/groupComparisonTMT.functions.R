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

    ## contrast matrix can be matrix or character vector.
    contrast.pairwise <- TRUE
    if(is.matrix(contrast.matrix)){
        contrast.pairwise <- FALSE
    }

    groups <- as.character(unique(data$Group)) # groups
    if(length(groups) < 2){
        stop("Please check the Condition column in annotation file. There must be at least two conditions!")
    }

    ## set up matrix for result
    if(contrast.pairwise){
        ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison

    } else {
        # # comparison come from contrast.matrix
        if (!all(colnames(contrast.matrix) %in% groups)) {

            stop("Please check the contrast.matrix. Column names of contrast.matrix. must be matched with conditions!")
        }

        ncomp <- nrow(contrast.matrix)
    }

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

    proteins <- fitted.models$protein # proteins
    num.protein <- length(proteins)
    res <- as.data.frame(matrix(rep(NA, 6 * num.protein * ncomp), ncol = 6)) ## store the inference results
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
    data$Group <- as.factor(data$Group) # make sure group is factor
    data$Run <- as.factor(data$Run)
    count <- 0

    s2.all <- fitted.models$s2  # group variance
    df.all <- fitted.models$df  # degree freedom
    lms <- fitted.models$model # linear models
    nrun <- length(unique(data$Run)) # check the number of MS runs in the data
    for(i in 1:length(lms)){
      message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
      
      ## get the data for protein i
      sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
      sub_data <- na.omit(sub_data)
      
      if(nrow(sub_data) != 0){
        sub_groups <- as.character(unique(sub_data$Group)) # groups in the sub data
        sub_groups <- sort(sub_groups) # sort the groups based on alphabetic order
        testable <- FALSE ## Indicate whether the protein is testable
        
        ## get the linear model for protein i 
        fit <- lms[[i]]
        MSE <- s2.all[i]
        df <- df.all[i]
        
        if(!is.null(fit)){ ## check the model is fittable 
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
            testable <- TRUE # mark the protein is testable
            
          } else{ # multiple run case
            ## Get estimated fold change from mixed model
            coeff <- fixed.effects(fit$mixed)
            coeff[-1] <- coeff[-1] + coeff[1]
            
            # Find the group name for baseline
            names(coeff) <- gsub("Group", "", names(coeff))
            names(coeff)[1] <- setdiff(as.character(sub_groups), names(coeff))
            
            s2.post <- (s2.prior * df.prior + MSE * df)/(df.prior + df)
            df.post <- df + df.prior
            testable <- TRUE # mark the protein is testable
            
          }
        }
        
        if(contrast.pairwise){
          ## Pairwise comparison
          # perform testing for each pair of conditions
          for(j in 1:(length(sub_groups)-1)){
            for(k in (j+1):length(sub_groups)){
              count <- count + 1
              res[count, "Protein"] <- proteins[i] ## protein names
              res[count, "Comparison"] <- paste(sub_groups[j], sub_groups[k], sep = "-") ## comparison
              
              if(testable){
                g1_df <- nrow(sub_data %>% dplyr::filter(Group == sub_groups[j])) ## size of group 1
                g2_df <- nrow(sub_data %>% dplyr::filter(Group == sub_groups[k])) ## size of group 2
                variance <- s2.post*sum(1/g1_df + 1/g2_df) ## variance of diff
                FC <- coeff[sub_groups[j]] - coeff[sub_groups[k]] ## fold change
                res[count, "log2FC"] <- FC
                ## Calculate the t statistic
                t <- FC/sqrt(variance) ## t statistic
                p <- 2*pt(-abs(t), df = df.post) ## p value
                res[count, "pvalue"] <- p
                res[count, "SE"] <- sqrt(variance) ## se
                res[count, "DF"] <- df.post
              } else{
                res[count, "log2FC"] <- NA
                res[count, "pvalue"] <- NA
                res[count, "SE"] <- NA
                res[count, "DF"] <- NA
              }
            }
          }
        } else { ## Compare one specific contrast
          # perform testing for required contrasts
          for(j in 1:nrow(contrast.matrix)){
            count <- count + 1
            res[count, "Protein"] <- proteins[i] ## protein names
            res[count, "Comparison"] <- row.names(contrast.matrix)[j] ## comparison
            
            # make sure all the groups in the contrast exist
            if(all(colnames(contrast.matrix)[contrast.matrix[j,]!=0] %in% sub_groups) & testable){
              ## calculate the size of each group
              group_df <- sub_data %>% group_by(Group) %>% dplyr::summarise(n = sum(!is.na(Abundance)))
              group_df$Group <- as.character(group_df$Group)
              
              ## variance of diff
              variance <- s2.post * sum((1/group_df$n) * ((contrast.matrix[j, group_df$Group])^2))
              FC <- sum(coeff*(contrast.matrix[j, names(coeff)])) # fold change
              res[count, "log2FC"] <- FC
              
              ## Calculate the t statistic
              t <- FC/sqrt(variance)
              
              ## calculate p-value
              p <- 2*pt(-abs(t), df = df.post)
              res[count, "pvalue"] <- p
              
              ## SE
              res[count, "SE"] <- sqrt(variance)
              res[count, "DF"] <- df.post
              
            } else{
              res[count, "log2FC"] <- NA
              res[count, "pvalue"] <- NA
              res[count, "SE"] <- NA
              res[count, "DF"] <- NA
              
            }
          }
        }
      } # if the protein data is empty
    } # for each protein
    
    res <- as.data.frame(res[1:count,])
    res$log2FC <- as.numeric(as.character(res$log2FC))
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$adjusted.pvalue <- NA
    comps <- unique(res$Comparison)

    ## Adjust multiple tests for each comparison
    for(i in 1:length(comps)){
        res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
    }

    res <- res[, c("Protein",
                   "Comparison",
                   "log2FC",
                   "SE",
                   "DF",
                   "pvalue",
                   "adjusted.pvalue")]
    return(res)
}
