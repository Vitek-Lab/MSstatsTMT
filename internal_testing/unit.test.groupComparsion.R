
#set up
dir<-"/Users/haosicheng/Desktop/MSstatsTMT/internal_testing/data"
setwd(dir)
library(MSstatsTMT)
library(dplyr)

load("LogSum.abun.rda")
data.long <- LogSum.abun %>% filter(Group != "Norm")

set.baseline<-FALSE
if(set.baseline){
  # limma
  limma.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "limma")
  # t test
  ttest.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "t")
  # proposed
  proposed.res <- MSstatsTMT::groupComparison.TMT(data.long,model = "proposed")

  save(limma.res,file = "limma.res.rda")
  save(ttest.res,file = "ttest.res.rda")
  save(proposed.res,file = "proposed.res.rda")
} else{
  # limma
  limma.res.test <- MSstatsTMT::groupComparison.TMT(data.long,model = "limma")
  # t test
  ttest.res.test <- MSstatsTMT::groupComparison.TMT(data.long,model = "t")
  # proposed
  proposed.res.test <- MSstatsTMT::groupComparison.TMT(data.long,model = "proposed")
  load("limma.res.rda")
  load("ttest.res.rda")
  load("proposed.res.rda")
  if(!all.equal(limma.res.test,limma.res)){
    print("limma has different result!!")
  }
  if(!all.equal(ttest.res.test,ttest.res)){
    print("ttest has different result!!")
  }
  if(!all.equal(proposed.res.test,proposed.res)){
    print("proposed model has different result!!")
  }

}
