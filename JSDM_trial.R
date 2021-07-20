library(Hmsc)
library(snow)
library(miceadds)

setwd("./JSDM")

load.Rdata("DOW_Top10PerTaxa_unfittedmodel.Rdata", "DOW_Top10PerTaxa_unfittedmodel")
ufm=DOW_10PerTaxa_unfittedmodel

samples_list = c(1000,1000,1000)
thin_list = c(25,50, 100)
nChains = 3

for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  fm = sampleMcmc(ufm, samples = samples, thin=thin,
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains, nParallel = nChains) 
    
  filename = paste("models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")

  save(fm, file=filename)
}

