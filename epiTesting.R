
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
  registerDoParallel(cores=16)
}else if(os=="aws"){

  registerDoParallel(cores = 72)
}


#just one barcode length
generations = c(4,5,6)#,6)

#generations = c(5)

barcodes = c(15) #according to Amjad estimation
integrases = c(1)

#edit rate
mus= c(0.1,0.2,0.3,0.4,0.5,0.6)#,0.8)
mus =c(0.1)
#mus=c(0.2,0.3)

nRepeats =10
#its easier to assume that open regions will have max edit rate
#also this value will be a combination of the actual edit rate and max-open accessibility


closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.9,1)
closed.vals = c(0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,1)
open.vals = rep(1,length(closed.vals)) #open is always max
#they will all map to a dynamic range: open.val / closed.val
switching.prs = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1)

#constant parameters for now
ng =1
m = 1
b=1
muVariation=list()
for(m in 1:length(mus)){
  #lets test two parameters:
  genData=list()

  for (ng in 1:length(generations)){
    dynamicData = list()
    dynamicTrees = list()
    dynamicHistories = list()
    # # # # # # # #
    mean.matrix=array(0,dim = c(length(open.vals),length(switching.prs)))
    mean.matrix2=array(0,dim = c(length(open.vals),length(switching.prs)))

    for(c in 1:length(open.vals)){
        switchingData = list()
        switchingTrees = list()
        switchingHistories = list()
        # # # # # # # # # # #
        for(s in 1:length(switching.prs)){
                xxx=compareDist(nGen =generations[ng],chr.acc=open.vals[c],chr.unacc = closed.vals[c],
                                barcodeLength =barcodes[b],nRepeat =nRepeats,recType="epimemoir",Pr_switch = switching.prs[s],mu = mus[m])
                switchingData[[s]]= xxx[[1]]
                switchingTrees[[s]]= xxx[[2]]
                switchingHistories[[s]]= xxx[[3]] #these are the epigenetic histories in text format (needs to be parse)

                mean.matrix[c,s] = apply(xxx[[1]],2,mean)[1]
                mean.matrix2[c,s] =apply(xxx[[1]],2,mean)[2]




        }

          print(paste("sim: g=",toString(generations[ng])," mu",toString(mus[m])," BC=",toString(barcodes[b])," Open vals=",toString(closed.vals[c]),sep=""))
        dynamicData[[c]] = switchingData
        dynamicTrees[[c]]=switchingTrees
        dynamicHistories[[c]] = switchingHistories
    }
        save(dynamicData,file = paste("./simdata/singleTrans_epiTest_gen_",toString(generations[ng]),"_mu", toString(mus[m]) ,"_.rdata",sep=""))
        save(dynamicTrees,file = paste("./simdata/singleTrans_epiTest_gen_",toString(generations[ng]),"_mu",toString(mus[m]),"_TREES.rdata",sep=""))
        save(dynamicHistories,file = paste("./simdata/singleTrans_epiTest_gen_",toString(generations[ng]),"_mu",toString(mus[m]),"_EPI.rdata",sep=""))

        genData[[ng]] = dynamicData
  }
        muVariation[[m]]= genData
}
#




# max.RF = 2*(2^generations[ng])-6
# mean.matrix=1-mean.matrix/max.RF
# mean.matrix2=1-mean.matrix2/max.RF
#
# colnames(mean.matrix)<-switching.prs
# row.names(mean.matrix)<-closed.vals/open.vals
#
# colnames(mean.matrix2)<-switching.prs
# row.names(mean.matrix2)<-closed.vals/open.vals
#
# heatmap.compare(mean.matrix,mean.matrix2)
