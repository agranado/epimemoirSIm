
library(gplots)
source("cascadeFunctions.R")

#plots for single transition epimemoir (fist analysis)
plots.dir ="./plots/single/"

#just one barcode length
generations = c(4,5,6)#,6)
#generaitons = c(4,5)
#generations = c(5)

barcodes = c(15) #according to Amjad estimation
integrases = c(1)

#edit rate
mus= c(0.1,0.2,0.3,0.4,0.5,0.6)#,0.8)
#mus=c(0.2,0.3)

nRepeats =100
#its easier to assume that open regions will have max edit rate
#also this value will be a combination of the actual edit rate and max-open accessibility


closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.9,1)
closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.9,1)
open.vals = rep(1,length(closed.vals)) #open is always max
#they will all map to a dynamic range: open.val / closed.val
switching.prs = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1)

#constant parameters for now
ng =1
m = 1
b=1
max.rf =2^(generations+1)-6

#load dynamicData object for a particular combination of G/mu
#load("simdata/epiTest_gen_5_mu0.4_.rdata")
#load dynamic trees object
#load("simdata/epiTest_gen_5_mu0.4_TREES.rdata")

plot.all.heatmaps <-function(muVariation){
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
            heatmap.simple(mat.data2,xlab="transition rate",ylab="dynamic range",this.min=this.min,this.max =this.max,
                        title = paste("Mat1 G = ",toString(generations[ng])," mu = ",toString(mus[mu]) ,sep=""))
            heatmap.simple(mat.data,xlab="transition rate",ylab="dynamic range",this.min=this.min,this.max = this.max,
                        title = paste("Mat 2: G = ",toString(generations[ng])," mu = ",toString(mus[mu]) ,sep=""))
          dev.off()


      }
    }
}

compare.epi.memoir<-function(){

    x11()
    par(mfrow=c(1,2))
    Pr_switch =c(1/7,0,1/8)
    mus = c(0.3,0.2,0.3)
    all.data = list()
    for(j in 1:length(mus)){

       b = compareDist(nGen=5,mu=mus[j],alpha=1/2,barcodeLength=12,chr.acc=1,chr.unacc = 0.2,Pr_switch = Pr_switch[j],nRepeats=200,recType="epimemoir")
       mean.dist= 1-apply(a[[1]],2,mean)/ (2^5-6)
       #[1] 0.6884615 0.7169231 1.0000000
       ground.truth =b[[2]]
       rec=reconstruct.tree.list(ground.truth,mu=0.3,nGen=4,alpha=1/2)
       all.dist= c();

       for(i in 1:length(rec)){all.dist[i] = 1-RF.dist(rec[[i]],ground.truth[[i]],normalize = T)}
       #hist(all.dist,ylim =c(0,40),xlim=c(0,1),breaks = 10,ylab ="counts",xlab = "reconstructability")
       all.data[[j]]=all.dist
     }
     #plotting histograms
     x11();
     par(mfrow=c(1,length(mus)));
     for(m in 1:length(mus)){
     hist(all.data[[m]],breaks=8,ylim=c(1,100),xlim=c(0.2,1),
        main="transition on->off",xlab="reconstructability",col="cadetblue2",cex=2);
      }
     #hist(all.data[[2]],breaks=8,ylim=c(1,100),xlim=c(0.2,1),
      #  main="no transition",xlab ="reconstructability",col="cadetblue3",cex =2)
}


convert.epihistory<-function(){

  




}
