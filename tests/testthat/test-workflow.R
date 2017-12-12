require(MSstatsTMT)
require(testthat)
require(tidyr)
require(dplyr)
#format testing for raw data



#processed data
Protein<-rep("O00299",60)
PSM<-append(rep("[K].iEEFLEAVLcPPr.[Y]_3",30),rep("[R].gFTIPEAFr.[G]_2",30))
log2Intensity<-c(16.67548,17.64749,17.60403,10.61325,12.47996,13.16848,16.76335,17.62280,17.62021,11.44671,12.76955,12.9090,16.70407,17.79888,17.56217,11.29623,12.88623,13.02754,16.77021,17.41688,17.76605,11.43990,12.37762,12.87553,16.77512,17.62933,17.49535,11.35146,12.90450,13.02290,16.72687,17.69113,17.78993,11.39507,12.88172,12.97181,16.89887,17.63559,17.79449,11.74253,12.62412,13.00693,16.78048,17.72154,17.84440,11.31412,12.82034,13.19102,16.73579,17.45225,17.75884,10.37344,13.24817,12.62399,16.58207,17.72192,17.54488,10.82730,12.84803,13.45137)
Channel<-rep(c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131"),6)
Run<-rep(append(append(rep(1,10),rep(2,10)),rep(3,10)),2)
BiologicalMixture<-rep("Mixture1",60)
Group<-rep(c("Norm","0.667","0.125","0.5","1","0.125","0.5","1","0.667","Norm"),6)
t.data<-data.frame(Protein=as.factor(Protein),PSM=as.factor(PSM),log2Intensity=log2Intensity,Channel=Channel,Run=Run,BiologicalMixture=BiologicalMixture,Group = Group)
t.data<-unite(t.data,"Subject",Run,Channel,sep = ".")
t.data$Channel<-as.factor(Channel)

t.data$Run<-Run

#annotation
Run<-append(append(rep(1,10),rep(2,10)),rep(3,10))
Channel<-rep(c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131"),3)
Group<-rep(c("Norm","0.667","0.125","0.5","1","0.125","0.5","1","0.667","Norm"),3)
annotation<-data.frame(Run=as.character(Run),Channel=Channel,Group=Group)
annotation$Channel<-as.factor(annotation$Channel)
annotation$Group<-as.factor(annotation$Group)

#expect_abun
Protein<-rep("O00299",30)
Run<-append(append(rep(1,10),rep(2,10)),rep(3,10))
Channel<-rep(c("X126","X127_N","X127_C","X128_N","X128_C","X129_N","X129_C","X130_N","X130_C","X131"),3)
BiologicalMixture<-rep("Mixture1",30)
Group<-rep(c("Norm","0.667","0.125","0.5","1","0.125","0.5","1","0.667","Norm"),3)
expect_abun<-data.frame(Protein=Protein,Subject1=Run,Subject2=Channel,Run=Run,Channel=Channel,BiologicalMixture=BiologicalMixture,Group = Group)
expect_abun<-unite(expect_abun,"Subject",Subject1,Subject2,sep = ".")
Abundance<-c(18.61628,19.58435,19.61485,12.97136,14.60965,14.98837,18.74757,19.54408,19.62485,13.51706,13.69867,13.95880,17.74278,18.76073,18.71018,12.30520,13.85366,14.11159,17.75310,18.43467,18.32140,11.56200,13.43654,13.31419,17.24078,18.23532,18.07928,11.67201,13.43549,13.81194)
expect_abun$Abundance<-Abundance
expect_abun$Protein<-as.character(expect_abun$Protein)
expect_abun<-expect_abun[,c(3,1,7,4,2,6,5)]
expect_abun<-data.table::as.data.table(expect_abun)


#expect_summary
Protein<-rep("O00299",10)
Comparison<-c("Norm-0.667","Norm-0.125","Norm-0.5","Norm-1","0.667-0.125","0.667-0.5","0.667-1","0.125-0.5","0.125-1","0.5-1")
log2FC<-c(0.080237232,0.012827161,0.112204749,0.085288311,-0.067410071,0.031967517,0.005051079,0.099377587,0.072461150,-0.026916438)
pvalue<-c(0.9506057,0.9920977,0.9309768,0.9475014,0.9584928,0.9803080,0.9968882,0.9388478,0.9553864,0.9834189)
SE<-rep(1.25536312093724,10)
DF<-rep(8,10)
adjusted.pvalue<-rep(0.9968882,10)
expect_summary<-data.frame(Protein=Protein,Comparison=Comparison,log2FC=log2FC,SE=as.factor(SE),DF=as.factor(DF),pvalue=pvalue,adjusted.pvalue=adjusted.pvalue)


# test_that("Expect raw data to be equal"),{
#
# }

test_that("Expect protein abundence to be equal",{
    abun<-MSstatsTMT::protein.summarization(t.data,method="LogSum")
    abun$Abundance<-round(abun$Abundance,digits = 3)
    expect_abun$Abundance<-round(expect_abun$Abundance,digits = 3)

    expect_equal(abun,expect_abun)
})

#Fix this after the model(groupComarison) is finalized
# test_that("Expect summarization to be equal",{
#     abun<-MSstatsTMT::protein.summarization(t.data,method="LogSum")
#     summary<-MSstatsTMT::groupComparison.TMT(abun)
#     summary$pvalue<-round(summary$pvalue,digits = 3)
#     expect_summary$pvalue<-round(expect_summary$pvalue,digits = 3)
#     summary$adjusted.pvalue<-round(summary$adjusted.pvalue,digits = 3)
#     expect_summary$adjusted.pvalue<-round(expect_summary$adjusted.pvalue,digits = 3)
#     expect_equal(summary,expect_summary)
# }
# )
