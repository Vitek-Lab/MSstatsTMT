

#' @export
proposed.model <-function(data, adj.method = "BH") {
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
      # Pairwise comparison
      for(j in 1:(length(groups)-1)){
        for(k in (j+1):length(groups)){
          count = count + 1
          res[count, 1] <- proteins[i] # protein names
          res[count, 2] <- paste(groups[j], groups[k], sep = "-") # comparison
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
        }
      }
    } else { # if there is only one group or one run in the data
      for(j in 1:(length(groups)-1)){
        for(k in (j+1):length(groups)){
          count = count + 1
          res[count, 1] <- proteins[i]
          res[count, 2] <- paste(groups[j], groups[k], sep = "-")
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
  return(res)
}

#' @export
ebayes.limma <- function(data, adj.method = "BH"){
  data.mat <- data[,c("Protein", "Subject", "Abundance")]
  data.mat <- data.mat %>% tidyr::spread(Subject, Abundance) # long to wide
  rownames(data.mat) <- data.mat$Protein # Assign the row names
  data.mat <- data.mat[, colnames(data.mat)!= "Protein"]

  # Extract the group information
  Annotation <- unique(data[,c("Subject", "Group")])
  rownames(Annotation) <- Annotation$Subject # Assign the row names
  Annotation <- Annotation[colnames(data.mat), ]
  Group <- Annotation$Group

  # Generate design matrix
  Group <- paste("v", Group, sep="")
  groups <-unique(Group)
  design <- model.matrix(~0+factor(Group, levels = groups))
  colnames(design) <- groups
  fit <- lmFit(data.mat, design) # Fit linear model

  # Do pairwise comparison
  resList <- list()
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
  res <- rbindlist(resList, use.names=TRUE, idcol = "Comparison")
  res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
  res$Comparison <- gsub("v", "", res$Comparison)
  return(res)
}

#' @export
protein.ttest <- function (data, adj.method = "BH"){

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
    if(nrow(sub) > 0){ # pairwise comparison
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
    }
  }
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
  res <- as.data.frame(res)
  res <- res[res$Comparison!=0,] # remove the tests which can't be done
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
  return(res)
}
