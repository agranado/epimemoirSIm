#This script updates the simulation to investigate epigenetic transitions:
#Recording channels will be affected by a new parameter that will account for accessibility of chromatin.
#Each region will have it's own parameter and therefore, each cell will inherit this parameter from their mother.
#Accessibility is dynamic and therefore there will be transitions (user defined) in the ground truth
#The aim of this simulation is to investigate to what extent can we reconstruct the lineage and hisotry of transitions from recorded data


library(phangorn)
library(stringdist)
library(doParallel)
library(gplots)
source("simulation4.R")
source("additionalFunctions.R")
source("cascadeFunctions.R")
#MAIN REPOSITORY FOR RECONSTRUCTION




os=system("cat ../os.txt",intern = T) #Local Mac repository (laptop)
if(os=="mac"){
  git.path = "/Users/alejandrog/MEGA/Caltech/trees/GIT/"
}else if(os =="linux"){ #local linux subsystem (maybe under Windows10)
  git.path = "/home/agranado/MEGA/Caltech/trees/lineageSoftware/"

}else if(os=="aws"){ #AMAZON cloud computing server

  git.path="/home/ubuntu/alejandrog/Caltech/macbookBranch/lineageSoftware/"
}


#Load fucntions from the main GIT repository (lineageSoftware)
source(paste(git.path,"MLfunctions.R",sep=""))

registerDoParallel(8)


rand.dist<-c(10,  26,  58, 120, 250, 506)
nRepeats = 20
# tHIS IS THE FUNCTION CURRENTLY CALLED BY bitVStrit.R
compareDist <- function(simulationType='trit',nGen=5,mu=0.3,alpha_=1/2,barcodeLength=20,nRepeats=20,methods=c(),
                        recType="integrase",nIntegrases=2,chr.acc=c(),chr.unacc=c(),Pr_switch=1/8){


  results= foreach(i=1:nRepeats) %dopar% simMemoirStrdist(nGen=nGen,mu=mu,recType=recType,alpha=alpha_,
      barcodeLength=barcodeLength,methods=methods,simulationType=simulationType,
        nIntegrases = nIntegrases,chr.acc=chr.acc,chr.unacc = chr.unacc,Pr_switch=Pr_switch) #epimemoir parameters

 #let's unlist the results from the parallel calculations:
  results_=list()
  tree.list = list()
  first.cell.list = list()

  for(i in 1:length(results)){
    results_[[i]] =results[[i]][[1]]
    tree.list[[i]] = results[[i]][[2]]
    first.cell.list[[i]] = results[[i]][[3]]
  } #this takes the first element of the results list which is the array with distance calculations

  #save the simulated trees to the hard drive
  simulation.file = paste("recType_",recType,"_mu_",toString(mu),
     "_BC_",toString(barcodeLength),
     "_nG_",toString(nGen),"_Nrep_",toString(nRepeats),"_.rda",sep="")



  results.matrix=do.call(rbind,results_)
  #Optional when only interested in the mean
  #apply(results.matrix,2,mean)
  # return(results.matrix)
return(list(results.matrix,tree.list,first.cell.list))


}
# # # # # # # # # # # #
 # # # # # # # # # #
# # # # # # # # # epiMEMOIR



 #nGen=5;mu=0.4;alpha=1/2;barcodeLength=20;methods=c();simulationType='trit';recType="integrase";nIntegrases=2



     nGen=3;mu=0.4;alpha=1/2;barcodeLength=40;methods=c();simulationType='trit';recType="epimemoir";nIntegrases=4;chr.acc = c(0,0.3,0.6,1)
#NOV 15th 2018 before thanksgiving break
#Test stringdistance measures using the stringdist R library
#use the same format as before but testing different methods included in the stringdist function

#Nov20: the array chr.acc has as many elements as "chromatin" regions
#the total lenght of the recording array will be divided evenly by the number of chromatin regions
simMemoirStrdist<-function(nGen=3,mu=0.4,alpha=1/2,barcodeLength=40,methods=c(),
                          simulationType='trit',recType="epimemoir",nIntegrases=4,chr.acc = c(0.8),chr.unacc =c(0.1),Pr_switch =1/8){
  #load necessary libraries and functions
  #detection of OS

  #update (trying to create a single branch that works in AWS and in my laptop)
  os=system("cat ../os.txt",intern = T)
  if(os=="mac"){ #local Alejandro's laptop
    pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
    pathName2="/Users/alejandrog/MEGA/Caltech/trees/simulation"
  }else if(os=="linux"){ #AWS server or any other linux machine (path is for AWS)

    pathName = "/home/agranado/MEGA/Caltech/trees/simulations/"
    pathName2= "/home/agranado/MEGA/Caltech/trees/simulations"

  }else if(os=="aws"){
    #linux desktop
    file.dir = "/home/alejandrog/MEGA/Caltech/trees/GraceData/integrase-data/10mer/"
    #setwd("/home/alejandrog/MEGA/Caltech/trees/macbookBranch/lineageSoftware")
    pathName2= "/home/alejandrog/MEGA/Caltech/lineage"
    pathName= "/home/alejandrog/MEGA/Caltech/lineage/"
    #At this point we are supposed to be in the local branch already, because we read the os.txt...
    pathName="/home/ubuntu/alejandrog/Caltech/lineage/"
    pathName2="/home/ubuntu/alejandrog/Caltech/lineage"
  }
  #clear the variable (since it behaves as global)
  if(exists("firstCell")){
    rm(firstCell)
  }

#convert numeric vectors to string for the tree to work properly
 chr = paste(chr.acc,sep="",collapse="_")
 chr.closed = paste(chr.unacc, sep="",collapse= "_")

# # # # # # # # CELL OBJECT INIT
 # # # # # # # #

 firstCell<-Node$new("1"); firstCell$barcode <-paste(rep("u",barcodeLength),collapse="");
# these number will multiply the basal edit rate
 firstCell$chr.acc = chr # accesibility values for open state
 firstCell$chr.closed = chr.closed #accessibility value for closed state

 firstCell$epihistory="1" #initial state for the fouding cell of the colony
  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

# SIMULATION
# # # # # # # # #
 # # # # # # # #
# # # # # # # # #

 #Here we implment 3 recording systems: 4 intgrases using 40 barcodes/10 each ;  2 integrases using 40 barcodes/ 20 each; and a single array using 40 barcodes
 # actIntegrase=1 # start recording using integrase 1
  #nIntegrases =4 # how many integrases comprise the cascade
 #1 is the first integrase (for nomal memoir , this is the only integrase)

#We need to assign a span of time (in generation units for each integrases to be active)
#The idea is to divide the total number of generations between the N integrases:
#So for 7 generation and 4 integrases we want a vector of activity = c(1,1,2,2,3,3,4)
# cascadeActivation function takes care of this

#this means that we want as many "channels" (integrases) as chromatin regions:
if(recType=="epimemoir") if(length(chr.acc)>0) nIntegrases=length(chr.acc)

  act_time=cascadeActivation(nGen, nIntegrases)


  for (g in 1:nGen){
    #this function simulates one generation of the tree
    actIntegrase = act_time[g]
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType,recType=recType,actIntegrase,nIntegrases,Pr_switch)
  }

  trueTree<-as.phylo.Node(firstCell)


  #get the sequences from the simulated tree + names
  barcodes<-firstCell$Get("barcode")

  #now we have the patters for ALL cells, but we need only the leaves, since that is what
  #we are going to use for reconstruction.
  #The way the tree was built, only the last 2^g cells are leaves; where g is the number of generations
  #take the number ID for the leaves.
  leavesID<-(length(barcodes)-2^nGen+1):length(barcodes)
  #grab those cells from the tree

  barcodeLeaves = array()
  namesLeaves=array()
  for (l in 1:length(leavesID)){
    barcodeLeaves[l] <-barcodes[names(barcodes)==leavesID[l]]
    namesLeaves[l] <-names(barcodes)[names(barcodes)==leavesID[l]]
  }
  names(barcodeLeaves)<-namesLeaves




  #take the barcodes from internal nodes, not useful here but for saving the complete tree simulations
  #
  nodesID = as.numeric(trueTree$node.label) #THIS WORKS
  barcodeNodes=array()
  namesNodes=array()
  for(l in 1:length(nodesID)){
    barcodeNodes[l]<-barcodes[names(barcodes)==nodesID[l]]
    namesNodes[l]<-names(barcodes)[names(barcodes)==nodesID[l]]
  }
  names(barcodeNodes)<-namesNodes
  ### SAVE TREE with NAMES

  named.tree = trueTree
  named.tree$tip.label = barcodeLeaves
  named.tree$node.label = barcodeNodes






  #now barcodeLeaves has all the leaves of the tree, with their original ID from the data.tree structure.
  #create Fasta file using only the leaves of the tree (n= 2^g)
  fastaBarcodes<-convertSimToFasta(barcodeLeaves)
  #convert name of variable to string

  #3  varName<-deparse(substitute(firstCell))

  fasID =toString(runif(1))
  #fasIN <-paste(pathName,"fasta/firstCell_",fasID,".fas",sep="")
  fasIN =tempfile("fasta/firstCell",tmpdir = pathName2)
  fasIN = paste(fasIN,fasID,".fas",sep="")

  #randomize the barcodes such that order is not a factor in the lineage reconstruction
  #the real tree (because the way is constructed, has an order 1,2,3,...N), if the barcodes are not
  #randomized, the order will "help" to the reconstruction which is not good!
  rand.barcode.order<-sample(1:length(fastaBarcodes))

  #We re-order the barcodes using a fixed (but random) order
  fastaBarcodes<-fastaBarcodes[rand.barcode.order]
  barcodeLeaves<-barcodeLeaves[rand.barcode.order]
  base::write(fastaBarcodes,file=fasIN) # dcGOR package overrides the write base method (which is stupid but that is they way it is)



  sed.return<-convertMemoirToDNA(fasIN)
  #now we can use the phyDat
  if(sed.return){
    memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
    #for distance based trees
    #from the phangorn tutorial pdf:
  }
  #Apr 8th: this is where the distance comes into place

  #Apr 8th:
  #we calculate string distances only for the leaves ( which is the data we actually get)

  allDistances = array()
  # for(m in 1:length(methods)){
  #   stringTree= upgma(stringdistmatrix(barcodeLeaves,method=methods[m]))
  #   stringTree$tip.label<-trueTree$tip.label
  #   allDistances[m]= RF.dist(stringTree,trueTree)
  # }

  #control against default dist.ml function + UPGMA, which so far is the best method
  m=0
  if(sed.return==1){
    dm<-dist.ml(memoirfas)
    dm.ham=dist.hamming(memoirfas)
    if(simulationType=='binary'){
      treeUPGMA<-upgma(dm.ham)
    }else{
      treeUPGMA<-upgma(dm)
    }
    treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
    allDistances[m+1]= RF.dist(treeUPGMA_noSeq,trueTree)
  }else{
    allDistances[m+1]= NaN
  }


  matdist_ = manualDistML(barcodeLeaves,mu*(chr.acc+chr.unacc)/2,alpha,nGen)
  manualTree_ =upgma(as.dist(t(matdist_)))
  manualTree_$tip.label= treeUPGMA$tip.label

  allDistances[m+2]= RF.dist(removeSeqLabel(manualTree_),trueTree)

  #try new distance using the built in dendrogram of heatmap2

  #h=heatmap.2(matdist_+t(matdist_),trace="none",dendrogram = 'column')
  # h=heatmap.2(matdist_+t(matdist_),Colv="Rowv")
  # heatmap.tree=as.phylo(as.hclust(h$colDendrogram))
  # heatmap.tree$tip.label = treeUPGMA$tip.label
  # allDistances[m+3]= RF.dist(removeSeqLabel(heatmap.tree),trueTree)
  #

  #alternative w/o plotting the actual heatmap, only hclust method
#  hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
#  hclust.tree$tip.label = treeUPGMA$tip.label
#  allDistances[m+3]= RF.dist(removeSeqLabel(hclust.tree),trueTree)

  # CASCADE RECONSTRUCTION::::::::: Sep27
  if(nIntegrases>1 & recType=="integrase"){
    r=cascadeReconstruction(barcodeLeaves,totalInts=nIntegrases,currentInts=1,nGen,mu,alpha)
    allDistances[m+3] = RF.dist(r,named.tree)
  }else{
    allDistances[m+3] = 0
  }




  #system(paste("rm ",firstCellFile,sep=""))
  if(sed.return){
    system(paste("rm ",paste(fasIN,".bak",sep=""),sep=""))
    system(paste("rm ",fasIN,sep=""))
  }
  #system(paste("rm ",fasIN,".bak",sep=""))


  #here we save the epihistory as data.frame
  #this has to be reconstructed into a tree later in the data analysis:

  firstCell.df=ToDataFrameTree(firstCell,"pathString","barcode","epihistory")
  #USE >n=FromDataFrameTable(firstCell.df)

  #epihistory= firstCell$Get("epihistory")
  return(list(allDistances,named.tree,firstCell.df))


}
#END of simulation function
# # # # # # # # # #
 # # # # # # # # #
# # # # # # # # # #





















runThisScript <- function (){
  #run memoirSim.R once and plot both trees.
  x11()
  par(mfrow=c(1,2))
  source("/Users/alejandrog/MEGA/Caltech/trees/simulation/memoirSim.R")
  ratioMat = manualDist(barcodeLeaves,mu,alpha,nGen );ratioTree = upgma(as.dist(t(ratioMat))); ratioTree$tip.label= treeUPGMA$tip.label;

  plot(ratioTree,main=paste("Manual tree",toString(RF.dist(trueTree,removeSeqLabel(ratioTree))) ))
  plot(treeUPGMA,main=paste("UPGMA dist tree",toString(RF.dist(trueTree,removeSeqLabel(treeUPGMA))) ))

  #execute multiple trees to get statistics
  # results= foreach(i=1:200) %dopar% simMemoirStrdist(3,0.3,alpha,barcodeLength,methods); results.matrix=do.call(rbind,results); meanDist =apply(results.matrix,2,mean)
}
