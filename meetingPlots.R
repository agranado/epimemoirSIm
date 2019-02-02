
library(gplots)
source("cascadeFunctions.R")

#plots for single transition epimemoir (fist analysis)
plots.dir ="./plots/single/"

#just one barcode length
#just one barcode length
# generations = c(4,5,6)#,6)
#generaitons = c(4,5)
#generations = c(5)

# barcodes = c(15) #according to Amjad estimation
# integrases = c(1)

#edit rate
mus= c(0.1,0.2,0.3,0.4,0.5,0.6)#,0.8)
#mus=c(0.2,0.3)
#mus=c(0.1)

nRepeats =72
#its easier to assume that open regions will have max edit rate
#also this value will be a combination of the actual edit rate and max-open accessibility


closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.9,1)
closed.vals = c(0,0.001,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.9,1)
open.vals = rep(1,length(closed.vals)) #open is always max
#they will all map to a dynamic range: open.val / closed.val
switching.prs = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1)

#constant parameters for now
# ng =1
# m = 1
# b=1
# max.rf =2^(generations+1)-6
#
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
    Pr_switch =c(1/2,0,1/8)
    mus = c(0.3,0.2,0.3)
    alpha  = 1/2
    nRepeats = 100
    all.data = list()
    for(j in 1:length(mus)){

       b = compareDist(nGen=5,mu=mus[j],alpha=alpha,barcodeLength=12,chr.acc=1,chr.unacc = 0.2,Pr_switch = Pr_switch[j],nRepeats=nRepeats,recType="epimemoir")
       mean.dist= 1-apply(b[[1]],2,mean)/ (2^5-6)
       #[1] 0.6884615 0.7169231 1.0000000
       ground.truth =b[[2]]
       rec=reconstruct.tree.list(ground.truth,mu=0.3,nGen=nGen,alpha=alpha)
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

#                               [[closed]][[switch.pr]]
# epihistory_4_3_=dynamicHistories[[4]][[3]]
# nGen = 4  extract.main.clades(rec)[[1]]



# execute this before
#this functino directly parses the data after loading it wiht:
# load("simdata/single...")
library(data.tree)
library(dcGOR)
# RUN THIS first
load.n.parse <-function(file,nGen, closed.vals,switching.pr){
  epihistory_4_3 = load(file)
  #data was saved as a list of data.frames:
  epihistory_4_3 = dynamicHistories[[closed.vals]][[switching.pr]]
  epi.df = history.list.ToDataFrameTree(epihistory_4_3)
  all.trees= convert.epihistory(epi.df,nGen=nGen)
  return(all.trees)
}
#includes the following functions:
history.list.ToDataFrameTree<-function(epihistory_4_3_){
  a=lapply(epihistory_4_3_,FromDataFrameTable)
  return(a)
}

convert.epihistory<-function(a, nGen,i=0){

  # up to here...

  #creates a list of data.tree Node objects
  #a=lapply(epihistory_4_3_,FromDataFrameTable)
  all.trees.epi = list()
  all.trees.barcode = list()

  if(i>0) which.trees = i else which.trees = 1:length(a)


  for(i in which.trees){
    epihistory=a[[i]]$Get("epihistory")
    barcodes = a[[i]]$Get("barcode")
    leavesID<-(length(epihistory)-2^nGen+1):length(epihistory)
    nodesID <-1:(length(epihistory)-2^nGen)

    leaves.history =epihistory[which(names(epihistory) %in% leavesID)]
    leaves.barcode =barcodes[which(names(barcodes) %in% leavesID)]

    nodes.history = epihistory[ which(names(epihistory) %in% nodesID  ) ]
    nodes.barcode = barcodes[ which(names(barcodes) %in% nodesID  ) ]

    #epigenetic tree
    b.phylo = as.phylo.Node(a[[i]])
    a.phylo = as.phylo.Node(a[[i]])

    b.phylo$tip.label = leaves.history
    b.phylo$node.label = nodes.history

    a.phylo$tip.label = leaves.barcode
    a.phylo$node.label= nodes.barcode

    all.trees.epi[[i]] = b.phylo
    all.trees.barcode[[i]] = a.phylo

  }
  #plot.phyl(b.phylo,show.nodes=T)
  return(list(all.trees.barcode,all.trees.epi))
}


#Dec 10th

#plot posterior Pr(x | mu for different mu's)
#for a given clade, we can calculate the Pr of the whole clade for a given mu value
clade.edit.rate<-function(rec,depth,genON){
  mu.scan=c(0.01,0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,0.95)
  prs=c()
  #x_mul(2,2)
  par(mfrow=c(2,2))
  for(t in 1:2){
    this.clade=extract.main.clades(rec)[[t]]
    par(cex.axis=1.1,cex.lab = 1)
    plot.phyl(this.clade,main =paste("clade",toString(t)),cex=1.3)
    for(m in 1:length(mu.scan))
      prs[m]=sum(log(clade.probability(this.clade$tip.label,nGen=depth,genON=genON,mu1=0.3,mu2=mu.scan[m],1/2)))
    plot(mu.scan,prs,ylim=c(-400,-200),type="o",lwd=2,ylab="log Pr", xlab="predited edit rate in clade")
  }
}

# Dec 13th
#assuming only two possible edit rates, we want to calculate the probabilities for all possible histories
# transition happen in G1, G2, ... or never happened, hopefully this will give us an estimate of the timing

#depth is the total number of generations
clade.history<-function(clade.barcodes, depth, mu1,mu2,alpha=1/2){
 #the minimum generation at which transition can happen is G2
  pr.history = list()
  for (gen.tr in 1:depth) {
    pr.history[[gen.tr]] = clade.probability(leaves.barcode =clade.barcodes,nGen = depth, genON = gen.tr, mu1,mu2,alpha)
  }
  return(pr.history)
}

#by default we are making alpha = 0 for bit calculations. Epimemoir is bit-based recording
clade.heatmap<-function(this.tree,nGen){

  clade.his.a.phylo = clade.history(this.tree$tip.label,depth=nGen,mu1=0.4,mu2=0.04,alpha=0)

  x11();
  pr.matrix= t(do.call(rbind, clade.his.a.phylo))

  #we use the last generation to normalize the matrix
  #The null hypothesis that NO transition happen so the branches were editing at the same time (pON)
  #in log matrix: x>1 means pON was more likely

  log.matrix = log(pr.matrix)/log(pr.matrix[,nGen])
  log.matrix[log.matrix<1] = -1/log.matrix[log.matrix<1]

  #because the transition probability is low, most cells in the matrix will have x>1 meaning not transition
  # log.matrix< -1.15
  my_palette <-  colorRampPalette(brewer.pal(11,"Spectral"))(n = 100)

  row.names(pr.matrix)<-paste(names(this.tree$tip.label),this.tree$tip.label)
  pheatmap(-log10(pr.matrix),trace="none",cluster_cols = F, cluster_rows = F,dendrogram="none",
        col = my_palette,scale ="none",treeheight_row=0,treeheight_col = 0)
  return(pr.matrix)
}

#Dec 17th
#I think I have an algortithm to call events
# "simdata/singleTrans_epiTest_gen_4_mu0.4_EPI.rdata"

# # # # # #
 # # # # #
# # # # # #
 # # # # #  LAB MEETING

# the following plots are not related to epimemoir

#Characterization of trit recording using information theory
steady.state.prs<-function(mu=0.4,alpha=1/2,gen.range=0:10){
    cols = c("cadetblue4","coral3","darkred")
    alphabet= c("u","r","x")

    prs = matrix(0,length(alphabet),length(gen.range))
    for(i in 1:length(alphabet)){
      prs[i,] =unlist(lapply(gen.range,Pr_s0,a=alphabet[i],alpha=alpha  ,mu=mu))
      if(i==1)
        plot(prs[i,],type="o",lwd=3,xlab ="generations",ylab="Expected fraction",col=cols[i],ylim=c(0,1))
      else
        lines(prs[i,],type="o",lwd=3,col = cols[i])
    }

    return(prs)
}

shannon.entropy<-function(px){
  px[px == 0]<-0.0000001
  return(sum(-px *log2(px)))
}

plot.entropy<-function(mu = 0.3,alpha=1/2,gen.range=0:10){

  recording.prs = steady.state.prs(mu=mu,alpha=alpha,gen.range=gen.range)

  plot(apply(recording.prs,2,shannon.entropy),lwd=3,type="o",xlab ="N cell divisions",ylab="Entropy")

}





##
