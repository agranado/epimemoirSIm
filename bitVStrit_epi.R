# Testing variations in cascade recording


rm(list=ls())
library(gplots)
source("simEpimemoir.R")
source("simulation4.R")
library(doParallel)

#request number of CPUs accordingly
os=system("cat ../os.txt",intern = T)
if(os=="mac"){
  registerDoParallel(cores=8)
}else if(os=="linux"){
  registerDoParallel(cores=6)
}else if(os=="aws"){

  registerDoParallel(cores = 72)
}

#SET PARAMETERS
barcodes = c(20,40,80,100)
#barcodes= c(2,10,50,100)
integrases = c(1)
#generations=c(8)
mus = c(0.5,0.3)
#mus=c(0.4,0.1)
generations=c(4,5,6)#,7,8,9,10)#,11,12)
#mus = c(0.99,0.6,0.5,0.4,0.1,0.01)
#mus=c(0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001)
#barcodes = c(6,7)
#generations = c(3,4,5)
if(os=="aws"){
  nRepeats=72
}else{
    nRepeats=20
}



types=c('trit')
#types=c('trit')

#Main variation of all parameters including:
#Edit rate, barcodeLength , number of integrases , number of generations
#
cascadevar=list(); cascadevar.tree = list() #this is the main object
for(ca in 1:length(integrases)){
    muVariation=list(); muVariation.tree = list()
    nIntegrases = integrases[ca]
    for(m in 1:length(mus)){
      simType=list(); simType.tree = list()
      mu=mus[m]
      for(st in 1:length(types)){
        simulationType=types[st]
        barcodeData=list(); barcodeData.tree = list()
        for(bc in 1:length(barcodes)){
           barcodeLength=barcodes[bc]
           genData=list();genData.tree = list()
           for(ng in 1:length(generations)){
              nGen=generations[ng]
              # this call is parallelized:
              xxx=compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats,nIntegrases = nIntegrases,recType="epimemoir")
              genData[[ng]]= xxx[[1]] #the first element in this list is the results.matrix object that contains all the distance measures calcualted inside the funciton
              genData.tree[[ng]]= xxx[[2]]


              print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength)," N_ints=",toString(nIntegrases),sep=""))
           }
           barcodeData[[bc]]=genData
           barcodeData.tree[[bc]] =genData.tree
        }
        simType[[simulationType]]=barcodeData
        simType.tree[[simulationType]] =barcodeData.tree
      }
      muVariation[[m]]=simType
      muVariation.tree[[m]]=simType.tree
    }
    cascadevar[[ca]] = muVariation
    cascadevar.tree[[ca]]= muVariation.tree
    save(cascadevar,file=paste("epi_intVar_int_",toString(nIntegrases),"_.rda",sep=""))
    save(cascadevar.tree,file=paste("epi_intVar_int_",toString(nIntegrases),"_trees_.rda",sep=""))
}
