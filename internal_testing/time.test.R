#testing code for MSstatsTMT

#set up
dir<-"/Users/haosicheng/Desktop/MSstatsTMT/Test/data"
setwd(dir)

#write .rd fils before build package
#library(devtools)
#library(roxygen2)
#devtools::document()
#install_github(repo = "Vitek-Lab/MSstatsTMT",auth_token = "0c4ece2b7e9aab30c8b1941e319a2c97fb2eb94e")

#######################sample testing#####################
library(MSstatsTMT)
test.data<-MSstatsTMT::test.data
# spiked-in proteins
sample<-FALSE
if(sample){
    # spiked-in proteins
    UPS.Proteins<-as.vector(as.matrix(unique(test.data%>%filter(grepl("ups",Protein))%>%dplyr::select(Protein))))
    # background proteins
    BACK.Proteins<-as.vector(as.matrix(unique(test.data%>%filter(!grepl("ups",Protein))%>%dplyr::select(Protein))))

    # a small subset data for testing
    test.protein <- c(UPS.Proteins[c(1,2)], BACK.Proteins[c(1,2,3)])
    test.data <- test.data%>%filter(Protein %in% test.protein)
}

# summarizaiton
annotation <- MSstatsTMT::annotation.data
annotation$Run <- as.character(annotation$Run)

time<-rep(0,5)
start<-Sys.time()
LogSum.abun <- protein.summarization(test.data, annotation,  "LogSum")
time[1]<-Sys.time()- start

start<-Sys.time()
Median.abun <- protein.summarization(test.data, annotation,  "Median")
time[2]<-Sys.time()- start

start<-Sys.time()
Biweight.abun <- protein.summarization(test.data, annotation,  "Biweight")
time[3]<-Sys.time()- start

start<-Sys.time()
MedianPolish.abun <- protein.summarization(test.data, annotation,  "MedianPolish")
time[4]<-Sys.time()- start

start<-Sys.time()
Huber.abun <- protein.summarization(test.data, annotation,  "Huber")
time[5]<-Sys.time()- start


save(time,file = "time.rda")






