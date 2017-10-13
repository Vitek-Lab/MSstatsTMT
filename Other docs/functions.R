library(MASS)
library(WRS2)
library(affy)
library(data.table)
library(dplyr)
library(tidyr)
library(lme4)
library(nlme)
library(daewr)
library(GAD)
library(multcomp)
library(agricolae)
library(pROC)
library(ROTS)
library(plyr)
library(limma)

# Summarize PSM level data to protein level
# data: PSM level data, which has columns Protein, PSM, Subject, Run, Channel, IonIntensity
# annotation: data frame with subject label information, which has three columns: Run, Channel, Group.
# Method: summarization methods. Possible options: LogSum, Median, Biweight, MedianPolish, Huber.
protein.summarization <- function(data, annotation, method){
  # Get the protein list, subjects and runs
  proteins <- unique(data$Protein)
  subjects <- unique(data$Subject)
  runs <- unique(data$Run)
  
  # Store the estimated protein abundance
  protein.abundance <- matrix(rep(NA, length(subjects)*length(proteins)), ncol = length(subjects))
  colnames(protein.abundance) <- subjects
  # For each protein and each run, do the summarization individually
  for(i in 1:length(proteins)) {
    message("Protein: ", i)
    for(j in 1:length(runs)){
      sub_data <- data %>% filter(Protein == proteins[i] & Run == runs[j])
      if(nrow(sub_data) != 0){
        nfea <- length(unique(sub_data$PSM))
        # Change the long format to wide format
        sub_data_wide <- sub_data %>% dplyr::select(IonIntensity, PSM, Subject) %>% spread(Subject, IonIntensity)
        rownames(sub_data_wide) <- sub_data_wide[,1]
        sub_data_wide <- sub_data_wide[,-1]
        # Number of negative values
        index <- which(apply(sub_data_wide, 1, function(col) any(col < 0)))
        if(length(index) != 0){
          # MC - 20170808 : replace negative values with zero.
          sub_data_wide[!is.na(sub_data_wide) & sub_data_wide < 0 ] <- 0
          message('* replace negatives with zero')
          # end MC- 20170808
        }
        
        if(nrow(sub_data_wide) != 0){
          if(nrow(sub_data_wide) == 1){ # Only one PSM for the protein
            protein.abundance[i, colnames(sub_data_wide)] <- as.matrix(sub_data_wide)
          } else{
            if(method == "LogSum"){
              #Sum
              # MC- 20170808 : change to log2 (sum of intensity)
              protein.abundance[i, colnames(sub_data_wide)] <- log2(colSums(2^sub_data_wide, na.rm = TRUE))
            }
            if(method == "Median"){
              #Median
              protein.abundance[i, colnames(sub_data_wide)] <- colMedians(as.matrix(sub_data_wide, na.rm = TRUE))
            }
            if(method == "Biweight"){
              #Biweight
              protein.abundance[i, colnames(sub_data_wide)] <- log2(generateExprVal.method.mas(as.matrix(2^sub_data_wide))$exprs)
            }
            if(method == "MedianPolish"){
              #median polish
              meddata  <-  medpolish(as.matrix(sub_data_wide), na.rm=TRUE,trace.iter = FALSE)
              tmpresult <- meddata$overall + meddata$col
              protein.abundance[i, colnames(sub_data_wide)] <- tmpresult[colnames(sub_data_wide)]
            }
            if(method == "Huber"){
              #Huber
              protein.abundance[i, colnames(sub_data_wide)] <- unlist(apply(as.matrix(sub_data_wide), 2, function(x) huber(x, k = 1.345)$mu))
            }
          }
        }
      }
    }
  }
  rownames(protein.abundance) <- proteins
  # Get the group information for each subject
  # Make the data long format and add the group information to protein level data frame
  res <- as.data.frame(protein.abundance)
  res$Protein <- rownames(res)
  res <- res %>% gather(Subject, Abundance, -Protein) # Change to long format
  res <- res %>% separate(Subject, c("Run", "Channel"), sep= "\\.", remove = FALSE) # Get the Channel and Run information from Subject
  res <- left_join(res, annotation)
  return(res)
}

# Propsoed inference model
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, IonIntensity
# adj.method: adjusted method for multiple comparison
proposed.model <-function(data, adj.method = "BH") {
  data <- data.long
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
      } else {
        # if there is only one group in the data
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
  }
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
  res <- as.data.frame(res)
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
  return(res)
}

# Limma inference model
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, Abundance
# adj.method: adjusted method for multiple comparison
ebayes.limma <- function(data, adj.method = "BH"){
  data.mat <- data[,c("Protein", "Subject", "Abundance")]
  data.mat <- data.mat %>% spread(Subject, Abundance) # long to wide
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

# t test 
# data: protein level data, which has columns Protein, Group, Subject, Run, Channel, Abundance
# adj.method: adjusted method for multiple comparison
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