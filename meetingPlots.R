
library(gplots)
source("cascadeFunctions.R")
#just one barcode length
generations = c(4,5,6)#,6)

barcodes = c(12) #according to Amjad estimation
integrases = c(1)

#edit rate
mus= c(0.4,0.6,0.8)#,0.8)

nRepeats =100

plots.dir = "./plots/"

closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.9,1)
open.vals = rep(1,length(closed.vals)) #open is always max
#they will all map to a dynamic range: open.val / closed.val
switching.prs = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1)



      max.rf =2^(generations+1)-6
      ng=1
      mu=1

#load dynamicData object for a particular combination of G/mu
#load("simdata/epiTest_gen_5_mu0.4_.rdata")
#load dynamic trees object
#load("simdata/epiTest_gen_5_mu0.4_TREES.rdata")

plot.all.heatmaps <-function(){
  for(mu in 1:length(mus)){
      for(ng in 1:length(generations)){
          dynamicData=muVariation[[mu]][[ng]]


          mat.data = matrix(0,length(dynamicData),length(dynamicData[[1]]));
          mat.data2 = matrix(0,length(dynamicData),length(dynamicData[[1]]));

          distance.idx = 1
          for(i in 1:length(dynamicData)){
            for(j in 1:length(dynamicData[[1]])){
               mat.data[i,j] = apply(dynamicData[[i]][[j]],2,mean)[1]
               mat.data2[i,j] = apply(dynamicData[[i]][[j]],2,mean)[2]
             }
           }
           mat.data = 1-mat.data/max.rf[ng]
           mat.data2 = 1-mat.data2/max.rf[ng]

           colnames(mat.data) = switching.prs
           row.names(mat.data) = closed.vals

           colnames(mat.data2) = switching.prs
           row.names(mat.data2) = closed.vals


           #heatmap range
           this.min =min(rbind(mat.data2,mat.data))
           this.min = this.min- 0.1*this.min

           this.max = max(rbind(mat.data2,mat.data)) + 0.1*this.min


          pdf(paste( plots.dir,"heatmap_compare_G",toString(generations[ng]),"_mu",toString(mus[mu]),"_.pdf" ,sep=""))
          #heatmap.compare(mat.data,mat.data2)
            heatmap.simple(mat.data,xlab="transition rate",ylab="dynamic range",this.min=this.min,this.max =this.max,
                        title = paste("Mat1 G = ",toString(generations[ng])," mu = ",toString(mus[mu]) ,sep=""))
            heatmap.simple(mat.data2,xlab="transition rate",ylab="dynamic range",this.min=this.min,this.max = this.max,
                        title = paste("Mat 2: G = ",toString(generations[ng])," mu = ",toString(mus[mu]) ,sep=""))
          dev.off()


      }
    }
}
