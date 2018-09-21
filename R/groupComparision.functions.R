#' @import lme4
#' @importFrom nlme fixed.effects
#' @importFrom dplyr group_by
#' @importFrom stats aggregate anova coef lm median medpolish model.matrix na.omit p.adjust pt t.test xtabs
#' @keywords internal

proposed.model <- function(data,
                           moderated = TRUE,
                           contrast.matrix = "pairwise",
                           adj.method = "BH") {

    Abundance <- Group <- Protein <- NULL

    if(moderated){ ## moderated t statistic
        ## Estimate the prior variance and degree freedom
        para <- estimate.prior.var(data)
        s2.prior <- para$s2.prior
        df.prior <- para$df.prior
    } else { ## ordinary t statistic
        s2.prior <- 0
        df.prior <- 0
    }

    data$Protein <- as.character(data$Protein) ## make sure protein names are character
    proteins <- as.character(unique(data$Protein)) ## proteins
    num.protein <- length(proteins)

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

    res <- as.data.frame(matrix(rep(0, 6 * length(proteins) * ncomp), ncol = 6)) ## store the inference results
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
    data <- as.data.table(data) ## make suree the input data is with data table format
    count <- 0
    ## do inference for each protein individually
    for(i in 1:length(proteins)) {
        message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
        sub_data <- data[Protein == proteins[i]] ## data for protein i
        sub_data <- na.omit(sub_data)
        tag <- FALSE ## Indicate whether there are enough measurements to train the linear model

        ## linear mixed model
        fit.mixed <- try(lmer(Abundance ~ 1 + (1|Mixture) + Group, data <- sub_data), TRUE)

        if(!inherits(fit.mixed, "try-error")){

            ## train linear model
            fit.fixed <- lm(Abundance ~ 1 + Mixture + Group, data <- sub_data)

            ## Get estimated fold change from mixed model
            coeff <- fixed.effects(fit.mixed)
            coeff[-1] <- coeff[-1] + coeff[1]

            # Find the group name for baseline
            names(coeff) <- gsub("Group", "", names(coeff))
            names(coeff)[1] <- setdiff(as.character(groups), names(coeff))

            # Estimate the group variance from fixed model
            av <- anova(fit.fixed)
            varcomp <- as.data.frame(VarCorr(fit.mixed))
            MSE <- varcomp[varcomp$grp == "Residual", "vcov"]
            df <- av$Df[3]
            s2.post <- (s2.prior * df.prior + MSE * df)/(df.prior + df)
            df.post <- df + df.prior

        } else {
            ## if there is only one run in the data, then train one-way anova
            fit.fixed <- try(lm(Abundance ~ Group, data <- sub_data), TRUE)

            if(!inherits(fit.fixed, "try-error")){
                ## Get estimated fold change from mixed model
                coeff <- coef(fit.fixed)
                coeff[-1] <- coeff[-1] + coeff[1]

                ## Find the group name for baseline
                names(coeff) <- gsub("Group", "", names(coeff))
                names(coeff)[1] <- setdiff(as.character(groups), names(coeff))

                ## Estimate the group variance from fixed model
                av <- anova(fit.fixed)
                MSE <- av$"Mean Sq"[2]
                df <- av$Df[2]
                s2.post <- (s2.prior * df.prior + MSE * df)/(df.prior + df)
                df.post <- df + df.prior

            } else {
                tag <- TRUE
            }
        }

        ##!! Sicheng or Ting, could you change the number for column to column name?
        if(contrast.pairwise){ ## Pairwise comparison

            for(j in 1:(length(groups)-1)){
                for(k in (j+1):length(groups)){
                    count <- count + 1
                    res[count, "Protein"] <- proteins[i] ## protein names
                    res[count, "Comparison"] <- paste(groups[j], groups[k], sep = "-") ## comparison

                    if(!tag){
                        g1_df <- nrow(sub_data[Group == groups[j]]) ## size of group 1
                        g2_df <- nrow(sub_data[Group == groups[k]]) ## size of group 2
                        variance <- s2.post*sum(1/g1_df + 1/g2_df) ## variance of diff
                        FC <- coeff[groups[j]] - coeff[groups[k]] ## fold change
                        res[count, "log2FC"] <- FC
                        ## Calculate the t statistic
                        t <- FC/sqrt(variance) ## t statistic
                        p <- 2*pt(-abs(t), df <- df.post) ## p value
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

            for(j in 1:nrow(contrast.matrix)){
                count <- count + 1
                res[count, "Protein"] <- proteins[i] ## protein names
                res[count, "Comparison"] <- row.names(contrast.matrix)[j] ## comparison

                if(!tag){
                    ## size of each group
                    group_df <- sub_data %>% group_by(Group) %>% dplyr::summarise(n = sum(!is.na(Abundance)))
                    group_df$Group <- as.character(group_df$Group)

                    ## variance of diff
                    variance <- s2.post * sum((1/group_df$n) * ((contrast.matrix[j, group_df$Group]) ^ 2))
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
    }
    res <- data.frame(res)
    res$log2FC <- as.numeric(as.character(res$log2FC))
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$adjusted.pvalue <- NA
    comps <- unique(res$Comparison)

    ## Adjust multiple tests for each comparison
    for(i in 1:ncomp){
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

#' @import limma
#' @keywords internal
estimate.prior.var <- function(data){

    Subject <- Abundance <- Protein <- NULL

    ## make sure data is in data.frame format
    data <- as.data.frame(data)
    data.mat <- data[, c("Protein", "Subject", "Abundance")]

    ## long to wide format
    data.mat <- data.mat %>% spread(Subject, Abundance)

    ## Assign the row names
    rownames(data.mat) <- data.mat$Protein
    data.mat <- data.mat %>% select(-Protein)

    ## Extract the group information
    Annotation <- unique(data[, c("Subject", "Group", "Run", "Mixture")])

    ## Assign the row names
    rownames(Annotation) <- Annotation$Subject
    Annotation <- Annotation[colnames(data.mat), ]
    group <- as.character(Annotation$Group)
    biomix <- as.character(Annotation$Mixture)

    if (length(unique(biomix)) == 1) { ## if there is only one mixture in the dataset
        design <- model.matrix(~0+group)
    } else{ ## there are multiple mixtures
        design <- model.matrix(~0+group+biomix)
    }

    data.mat <- na.omit(data.mat)
    ## Fit linear model
    fit <- lmFit(data.mat, design)
    fit2 <- eBayes(fit)

    return(list(df.prior = fit2$df.prior,
                s2.prior = fit2$s2.prior))
}
