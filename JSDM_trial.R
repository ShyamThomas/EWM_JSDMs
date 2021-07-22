library(tidyverse)
library(Hmsc)
library(snow)
library(miceadds)

###Folder Path "/Users/thom7552/UMNpostdoc/ProjectEWM/RProjects/JSDMs_hmsc"
###All data in the Data folder

###First load the entire JSDM dataset which contains year, lake ids, lake covariates, survey and sampling ids, and species presence absence (1/0) 
JSDM.data=read_csv("Data/JSDM.newdata.csv")
JSDM.data
colnames(JSDM.data)[20]

### Dataset of species with lake level prevalence 10 percent and above
Taxa_10Per=read_csv("Data/Taxa_10PerAbove.csv")
Taxa_10Per
TaxaNames_10Per=Taxa_10Per$Taxa
TaxaNames_10Per
### Subset JSDM dataset with the taxa of pervelance 10 percent and higher
JSDM.subdata=JSDM.data%>%select(c(1:20),all_of(TaxaNames_10Per))
JSDM.subdata

### A seperate df showing all 8 covariates and their names; rename covariates for clarity
Xdata=as.data.frame(cbind(JSDM.subdata$DEPTH_MAX,JSDM.subdata$avgSECCHI,JSDM.subdata$Road_Density,JSDM.subdata$Littoral,
                          JSDM.subdata$GDD_WTR_10c,JSDM.subdata$Depth_mts,JSDM.subdata$DD,JSDM.subdata$CummGLD))
colnames(Xdata)=c("MaxDepth", "AvgSecchi","RoadDensity","Littoral","SurfGDD","WithinLakeDepth","WithinLakeDegreeDays",
                  "WithinLakeGLD")
tail(Xdata)
dim(Xdata)

### A quick assessment of correlation among covariates
library(corrplot)
Xdata.cor=cor(Xdata, method="pearson")
Xdata.cor
corrplot.mixed(Xdata.cor, lower.col = "black", tl.pos = "l")

### The fixed effects formula that goes into the JSDM model structure
Xformula = ~MaxDepth+AvgSecchi+RoadDensity+Littoral+SurfGDD+WithinLakeDepth+WithinLakeDegreeDays+WithinLakeGLD
Xformula

### The taxa pres/abs data matrix
Ydata=as.matrix(JSDM.subdata[,c(21:61)])
head(Ydata)
dim(Ydata)

### Defining study design df and lake ids (DOWLKNUM) as a random factor
DOWLKNUM=as.factor(JSDM.subdata$DOWLKNUM)
sample.id=as.factor(1:length(JSDM.subdata$DOWLKNUM))
studyDesign.subdata=data.frame(sample=sample.id, plot=DOWLKNUM)
head(studyDesign.subdata)
rL.DOWLKNUM = HmscRandomLevel(units = levels(studyDesign.subdata$plot))
rL.DOWLKNUM
ranlevels = list("plot" = rL.DOWLKNUM)
ranlevels

### Combine study design, covariates and species data; and save as a single df
SXY_10PerTaxa=cbind(studyDesign.subdata,Xdata,Ydata)
head(SXY_10PerTaxa)
write_csv(SXY_10PerTaxa,"Data/SXY_JSDMdata_10per.csv")

### Model 1 : A very basic model, no random effects
DOW_Top10PerTaxa_simplemodel = Hmsc(Y=Ydata, XData = Xdata,
                                    XFormula = Xformula,
                                    studyDesign = studyDesign.subdata,
                                    distr = "probit")

DOW_Top10PerTaxa_simplemodel
save(DOW_Top10PerTaxa_simplemodel,file="DOW_Top10PerTaxa_unfittedsimplemodel.Rdata")

### Define the models
### Model 2: Model with lake ids as random effects
DOW_Top10PerTaxa_model = Hmsc(Y=Ydata, XData = Xdata,
                           XFormula = Xformula,
                           studyDesign = studyDesign.subdata,
                           ranLevels  = ranlevels,
                           distr = "probit")

DOW_Top10PerTaxa_model
save(DOW_Top10PerTaxa_model,file="DOW_Top10PerTaxa_unfittedmodel.Rdata")


### Finally, run the pre-defined models
load.Rdata("DOW_Top10PerTaxa_unfittedsimplemodel.Rdata", "DOW_Top10PerTaxa_unfittedsimplemodel")
ufm=DOW_10PerTaxa_unfittedsimplemodel

samples_list = c(1000,1000,1000,1000)
thin_list = c(25,50,75,100)
nChains = 3

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  fm = sampleMcmc(ufm, samples = samples, thin=thin,
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains, nParallel = nChains) 
    
  filename = paste("simp.models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")

  save(fm, file=filename)
}

### Repeat for unfitted model with lake ids as random effects
load.Rdata("DOW_Top10PerTaxa_unfittedmodel.Rdata", "DOW_Top10PerTaxa_unfittedmodel")
ufm=DOW_10PerTaxa_unfittedmodel

samples_list = c(1000,1000,1000,1000)
thin_list = c(25,50,75,100)
nChains = 3

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  fm = sampleMcmc(ufm, samples = samples, thin=thin,
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains, nParallel = nChains) 
    
  filename = paste("lakran.models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")

  save(fm, file=filename)
}
