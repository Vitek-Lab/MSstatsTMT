# Load the input data
source('./code/functions.R')
load("./data/test.data.rda")
data <- test.data

# spiked-in proteins
UPS.Proteins<-as.vector(as.matrix(unique(data%>%filter(grepl("ups",Protein))%>%dplyr::select(Protein))))
# background proteins
BACK.Proteins<-as.vector(as.matrix(unique(data%>%filter(!grepl("ups",Protein))%>%dplyr::select(Protein))))

# a small subset data for testing
subset.protein <- c(UPS.Proteins[c(1,2)], BACK.Proteins[c(1,2,3)])
subset.data <- data%>%filter(Protein %in% subset.protein)
  
#######################################################
#################### Summarization #########################
#######################################################

# Load the group information
annotation <- read.csv("./data/annotation.csv")
annotation$Run <- as.character(annotation$Run)

############### LogSum #########################
LogSum.abun <- protein.summarization(subset.data, annotation,  "LogSum")

#######################################################
#################### Inference #########################
#######################################################
# For this spiked-in dataset, we are not intrested in the normalization channel, so remove it
data.long <- LogSum.abun %>% filter(Group != "Norm")

############### Inference #########################
# limma
limma.res <- ebayes.limma(data.long)
# t test
ttest.res <- protein.ttest(data.long)
# proposed
proposed.res <- proposed.model(data.long)
