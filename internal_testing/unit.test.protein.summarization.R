
#set up
dir<-"/Users/haosicheng/Desktop/MSstatsTMT/internal_testing/data"
setwd(dir)

library(MSstatsTMT)


#TRUE to set baseline FALSE to test if the results changed
set.baseline<-FALSE

if(set.baseline){
    test.data<-MSstatsTMT::input.data
    sample<-TRUE    ##make sample = FALSE if you want to test the whole dataset, might take a while

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
    annotation <- MSstatsTMT::annotation
    annotation$Run <- as.character(annotation$Run)
    LogSum.abun <- protein.summarization(test.data, annotation,  "LogSum")
    Median.abun <- protein.summarization(test.data, annotation,  "Median")
    Biweight.abun <- protein.summarization(test.data, annotation,  "Biweight")
    MedianPolish.abun <- protein.summarization(test.data, annotation,  "MedianPolish")
    Huber.abun <- protein.summarization(test.data, annotation,  "Huber")
    save(LogSum.abun,file = "LogSum.abun.rda")
    save(Median.abun,file = "Median.abun.rda")
    save(Biweight.abun,file = "Biweight.abun.rda")
    save(MedianPolish.abun,file = "MedianPolish.abun.rda")
    save(Huber.abun,file = "Huber.abun.rda")
} else{
    test.data<-MSstatsTMT::input.data
    sample<-TRUE    ##make sample = FALSE if you want to test the whole dataset, might take a while

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
    annotation <- MSstatsTMT::annotation
    annotation$Run <- as.character(annotation$Run)
    LogSum.abun.test <- protein.summarization(test.data, annotation,  "LogSum")
    Median.abun.test <- protein.summarization(test.data, annotation,  "Median")
    Biweight.abun.test <- protein.summarization(test.data, annotation,  "Biweight")
    MedianPolish.abun.test <- protein.summarization(test.data, annotation,  "MedianPolish")
    Huber.abun.test <- protein.summarization(test.data, annotation,  "Huber")
    load("LogSum.abun.rda")
    load("Median.abun.rda")
    load("Biweight.abun.rda")
    load("MedianPolish.abun.rda")
    load("Huber.abun.rda")

    if(!all.equal(LogSum.abun,LogSum.abun.test)){
        print("LogSum has different result!!")
    }
    if(!all.equal(Median.abun,Median.abun.test)){
        print("Median has different result!!")
    }
    if(!all.equal(Biweight.abun,Biweight.abun.test)){
        print("Biweight has different result!!")
    }
    if(!all.equal(MedianPolish.abun,MedianPolish.abun.test)){
        print("MedianPolish has different result!!")
    }
    if(!all.equal(Huber.abun,Huber.abun.test)){
        print("Huber has different result!!")
    }
}
