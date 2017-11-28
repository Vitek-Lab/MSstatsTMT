
#' @importFrom nlme fixed.effects
#' @export
proposed.model <-function(data, cont.matrix = "pairwise", adj.method = "BH") {
    data$Protein <- as.character(data$Protein) # make sure protein names are character
    proteins <- unique(data$Protein) # proteins
    groups <- unique(data$Group) # groups
    ncomp <- length(groups)*(length(groups)-1)/2 # Number of comparison
    res <- matrix(rep(0, 6*length(proteins)*ncomp), ncol = 6) # store the inference results
    data <- as.data.table(data) # make suree the input data is with data table format
    count = 0
    # do inference for each protein individually
    for(i in 1:length(proteins)) {
        message("Protein: ", i)
        sub_data <- data[Protein == proteins[i]] # data for protein i
        tag <- FALSE; # Indicate whether there are enough measurements to train the linear model
        if(length(unique(na.omit(sub_data)$Run)) > 1 & length(unique(na.omit(sub_data)$Group)) > 1){
            fit.fixed <-lm(Abundance ~ 1 + Run + Group + Group:Run, data=sub_data) # train linear model
            fit.mixed <-lmer(Abundance ~ 1 + (1|Run) + Group + (1|Group:Run), data=sub_data)
            # Get estimated fold change from mixed model
            coeff <- fixed.effects(fit.mixed)
            coeff[-1] <- coeff[-1] + coeff[1]
            # Find the group name for baseline
            names(coeff) <- gsub("Group", "", names(coeff))
            names(coeff)[1] <- setdiff(as.character(groups), names(coeff))
            # Estimate the group variance from fixed model
            av <- anova(fit.fixed)
            MSE <- av$"Mean Sq"[3]
            df <- av$Df[3]

        } else {
            # if there is only one run in the data, then train one-way anova
            if(length(unique(na.omit(sub_data)$Run)) == 1 & length(unique(na.omit(sub_data)$Group)) > 1){
                fit.fixed <-lm(Abundance ~ Group, data=sub_data) # train linear model
                # Get estimated fold change from mixed model
                coeff <- coef(fit.fixed)
                coeff[-1] <- coeff[-1] + coeff[1]
                # Find the group name for baseline
                names(coeff) <- gsub("Group", "", names(coeff))
                names(coeff)[1] <- setdiff(as.character(groups), names(coeff))
                # Estimate the group variance from fixed model
                av <- anova(fit.fixed)
                MSE <- av$"Mean Sq"[2]
                df <- av$Df[2]
            } else {
                tag <- TRUE;
            }
        }
        if(cont.matrix == "pairwise"){ # Pairwise comparison
            for(j in 1:(length(groups)-1)){
                for(k in (j+1):length(groups)){
                    count = count + 1
                    res[count, 1] <- proteins[i] # protein names
                    res[count, 2] <- paste(groups[j], groups[k], sep = "-") # comparison
                    if(!tag){
                        g1_df <- nrow(sub_data[Group == groups[j]]) # size of group 1
                        g2_df <- nrow(sub_data[Group == groups[k]]) # size of group 2
                        variance <- MSE*sum(1/g1_df + 1/g2_df) # variance of diff
                        FC <- coeff[groups[j]] - coeff[groups[k]] # fold change
                        res[count, 3] <- FC
                        #Calculate the t statistic
                        t <- FC/sqrt(variance) # t statistic
                        p <- 2*pt(-abs(t), df = df) # p value
                        res[count, 4] <- p
                        res[count, 5] <- sqrt(variance) # se
                        res[count, 6] <- df
                    } else{
                        res[count, 3] <- NA
                        res[count, 4] <- NA
                        res[count, 5] <- NA
                        res[count, 6] <- NA
                    }
                }
            }
        } else{ # Compare one specific contrast
            # cont.matrix<-matrix(c(-1,0,1,0),nrow=1)
            # colnames(cont.matrix) <- c("0.125", "0.5", "0.667", "1")
            # row.names(cont.matrix) <- "0.667-0.125"
            for(j in 1:nrow(cont.matrix)){
                count = count + 1
                res[count, 1] <- proteins[i] # protein names
                res[count, 2] <- row.names(cont.matrix)[j] # comparison
                if(!tag){
                    group_df <- sub_data %>% group_by(Group) %>% dplyr::summarise(n = sum(!is.na(Abundance))) # size of each group
                    variance <- MSE*sum((1/group_df$n)*((cont.matrix[j,group_df$Group])^2)) # variance of diff
                    FC <- sum(coeff*(cont.matrix[j,names(coeff)])) # fold change
                    res[count, 3] <- FC
                    #Calculate the t statistic
                    t <- FC/sqrt(variance) # t statistic
                    p <- 2*pt(-abs(t), df = df) # p value
                    res[count, 4] <- p
                    res[count, 5] <- sqrt(variance) # se
                    res[count, 6] <- df
                } else{
                    res[count, 3] <- NA
                    res[count, 4] <- NA
                    res[count, 5] <- NA
                    res[count, 6] <- NA
                }
            }
        }
    }
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
    res <- as.data.frame(res)
    res$log2FC <- as.numeric(as.character(res$log2FC))
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
    res<-res[,c(1,2,3,5,6,4,7)]
    return(res)
}

#' @import limma
#' @import data.table
#' @export
ebayes.limma <- function(data, cont.matrix = "pairwise", adj.method = "BH"){
    data.mat <- data[,c("Protein", "Subject", "Abundance")]
    data.mat <- data.mat %>% spread(Subject, Abundance) # long to wide
    rownames(data.mat) <- data.mat$Protein # Assign the row names
    data.mat <- data.mat[, colnames(data.mat)!= "Protein"]

    # Extract the group information
    Annotation <- unique(data[,c("Subject", "Group")])
    rownames(Annotation) <- Annotation$Subject # Assign the row names
    Annotation <- Annotation[colnames(data.mat), ]
    group <- as.character(Annotation$Group)

    # Generate design matrix
    design <- model.matrix(~0+group)
    fit <- lmFit(data.mat, design) # Fit linear model

    # store inference results
    resList <- list()
    if(cont.matrix == "pairwise"){ # Do pairwise comparison
        groups <- paste("group", unique(group), sep = "")
        for(j in 1:(length(groups)-1)){
            for(k in (j+1):length(groups)){
                g1 = groups[j]
                g2 = groups[k]
                comp <- paste(g1, g2, sep="-")
                cont.matrix <- makeContrasts(contrasts=comp, levels=design)
                fit2 <- contrasts.fit(fit, cont.matrix)
                fit2 <- eBayes(fit2)

                proteins <- rownames(fit2$coefficients)
                pvalue <- as.vector(fit2$p.value)
                log2FC <- as.vector(fit2$coefficients)
                DF <- as.vector(fit2$df.total)
                SE <- as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
                resList[[paste(groups[j], groups[k],sep="-")]] <- data.frame(Protein = proteins, log2FC = log2FC, pvalue = pvalue, SE = SE, DF = DF)
            }
        }
    }else { # Compare one specific contrast
        # cont.matrix<-matrix(c(-1,0,1,0, 1, -1, 0,0), nrow=2, byrow= T)
        # colnames(cont.matrix) <- c("0.125", "0.5", "0.667", "1")
        # row.names(cont.matrix) <- c("0.667-0.125", "0.5-0.125")
        cont.matrix <- t(cont.matrix)
        rownames(cont.matrix) <- paste("group", rownames(cont.matrix), sep = "")
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)

        proteins <- rownames(fit2$coefficients)
        DF <- as.vector(fit2$df.total)
        SE <- as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
        if(ncol(cont.matrix) == 1){ # only one contrast
            resList[[colnames(cont.matrix)]] <- data.frame(Protein = proteins, log2FC = as.vector(fit2$coefficients), pvalue = as.vector(fit2$p.value), SE = SE, DF = DF)
        } else{
            for(l in 1:ncol(cont.matrix)){ # multiple contrasts
                resList[[colnames(cont.matrix)[l]]] <- data.frame(Protein = proteins, log2FC = as.vector(fit2$coefficients[,l]), pvalue = as.vector(fit2$p.value[,l]), SE = SE, DF = DF)
            }
        }
    }
    # Finalize the inference results
    res <- rbindlist(resList, use.names=TRUE, idcol = "Comparison")
    res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
    res$Comparison <- gsub("group", "", res$Comparison)
    res<-res[,c(2,1,3,5,6,4,7)]
    return(res)
}
#'@export
protein.ttest <- function (data, cont.matrix = "pairwise", adj.method = "BH"){
    data$Protein <- as.character(data$Protein) # make sure protein names are character
    proteins <- unique(data$Protein) # proteins
    groups <- unique(data$Group) # groups
    ncomp <- length(groups)*(length(groups)-1)/2 # Number of comparison
    res <- matrix(rep(0, 6*length(proteins)*ncomp), ncol = 6) # store the inference results
    data <- as.data.table(data) # make suree the input data is with data table format
    count <- 0
    # Do inference for each protein individually
    for(i in 1:length(proteins)){
        sub <- data[Protein==proteins[i]]
        sub <- na.omit(sub) # remvoe na rows
        if(nrow(sub) > 0){
            if(cont.matrix == "pairwise"){ # pairwise comparison
                for(j in 1:(length(groups)-1)){
                    for(k in (j+1):length(groups)){
                        count = count + 1
                        group1 <- sub[Group == groups[j]] # data from group 1
                        group2 <- sub[Group == groups[k]] # data from group 2
                        res[count, 1] <- proteins[i] # protein name
                        res[count, 2] <- paste(groups[j], groups[k], sep = "-") # comparison
                        if(nrow(group1) > 1 & nrow(group2) > 1){ # make sure both groups have data
                            # t test
                            t <- t.test(group1$Abundance, group2$Abundance)
                            res[count, 3] <- t$estimate[1]-t$estimate[2] # Estimate fold change
                            res[count, 4] <- t$p.value
                            res[count, 5] <- abs(t$estimate[1]-t$estimate[2])/t$statistic # se
                            res[count, 6] <- t$parameter # df
                        }
                    }
                }
            } else{ # Compare one specific contrast
                # cont.matrix<-matrix(c(-1,0,1,0),nrow=1)
                # colnames(cont.matrix) <- c("0.125", "0.5", "0.667", "1")
                # row.names(cont.matrix) <- "0.667-0.125"

                for(j in 1: nrow(cont.matrix)){
                    cont <- cont.matrix[j,]
                    cont <- cont[cont != 0]
                    if(length(cont) == 2){ # t test can only be used for two sample
                        count = count + 1
                        group1 <- sub[Group == names(cont)[1]] # data from group 1
                        group2 <- sub[Group == names(cont)[2]] # data from group 2
                        res[count, 1] <- proteins[i] # protein name
                        res[count, 2] <- paste(names(cont)[1], names(cont)[2], sep = "-") # comparison
                        if(nrow(group1) > 1 & nrow(group2) > 1){ # make sure both groups have data
                            # t test
                            t <- t.test(group1$Abundance, group2$Abundance)
                            res[count, 3] <- t$estimate[1]-t$estimate[2] # Estimate fold change
                            res[count, 4] <- t$p.value
                            res[count, 5] <- abs(t$estimate[1]-t$estimate[2])/t$statistic # se
                            res[count, 6] <- t$parameter # df
                        }
                    }
                }
            }
        }
    }
    colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
    res <- as.data.frame(res)
    res <- res[res$Comparison!=0,] # remove the tests which can't be done
    res$log2FC <- as.numeric(as.character(res$log2FC))
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
    res<-res[,c(1,2,3,5,6,4,7)]
    return(res)
}
