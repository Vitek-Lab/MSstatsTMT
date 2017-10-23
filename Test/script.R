#testing code for MSstatsTMT

#set up
dir<-"/Users/haosicheng/Desktop/MSstatsTMT"
setwd(dir)

#write .rd fils before build package
library(devtools)
library(roxygen2)
devtools::document()
#install_github(repo = "Vitek-Lab/MSstatsTMT",auth_token = "0c4ece2b7e9aab30c8b1941e319a2c97fb2eb94e")

#######################sample testing#####################
library(MSstatsTMT)
data<-MSstatsTMT::test.data
# spiked-in proteins
UPS.Proteins<-as.vector(as.matrix(unique(data%>%filter(grepl("ups",Protein))%>%dplyr::select(Protein))))
# background proteins
BACK.Proteins<-as.vector(as.matrix(unique(data%>%filter(!grepl("ups",Protein))%>%dplyr::select(Protein))))

# a small subset data for testing
subset.protein <- c(UPS.Proteins[c(1,2)], BACK.Proteins[c(1,2,3)])
subset.data <- data%>%filter(Protein %in% subset.protein)

# summarizaiton
annotation <- MSstatsTMT::annotation.data
annotation$Run <- as.character(annotation$Run)

#LogSum 
LogSum.abun <- protein.summarization(subset.data, annotation,  "LogSum")

#######################################################
#################### Inference #########################
#######################################################
# For this spiked-in dataset, we are not intrested in the normalization channel, so remove it
data.long <- LogSum.abun %>% filter(Group != "Norm")

############### Inference #########################
# limma
limma.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "limma")

# t test
ttest.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "t")
# proposed
proposed.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "proposed")


################test the entire data set####################


##################compare#################