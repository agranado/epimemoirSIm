#This script updates the simulation to investigate epigenetic transitions:
#Recording channels will be affected by a new parameter that will account for accessibility of chromatin.
#Each region will have it's own parameter and therefore, each cell will inherit this parameter from their mother.
#Accessibility is dynamic and therefore there will be transitions (user defined) in the ground truth
#The aim of this simulation is to investigate to what extent can we reconstruct the lineage and hisotry of transitions from recorded data

#Update Jan30th
  # we need to include multiple regions with chromatin parameters
  # The goal is to extend the model to multiple recording regions that change edit rate

  # V2.0 should include this changes and make the overall script more efficient and less clumsy
library(phangorn)
library(stringdist)
library(doParallel)
library(gplots)
source("simulation5.R")
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
compareDist <- function(simulationType='trit',nGen=5,mu=0.3,alpha_=0,barcodeLength=c(15,15),nRepeats=20,methods=c(),
                        recType="epimemoir",nIntegrases=1,chr.acc=c(1,1),chr.unacc=c(0.9,0.2),Pr_switch=c(0,1/8)){


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



nGen=3;mu=0.4;alpha=1/2;barcodeLength=c(20,15);methods=c();simulationType='trit';
recType="epimemoir";nIntegrases=2;chr.acc = c(1,1); chr.unacc = c(1,0.2);Pr_switch =c(0,1/8)
#NOV 15th 2018 before thanksgiving break
#Test stringdistance measures using the stringdist R library
#use the same format as before but testing different methods included in the stringdist function

#Nov20: the array chr.acc has as many elements as "chromatin" regions
#the total lenght of the recording array will be divided evenly by the number of chromatin regions
#Jan 30: the array chr.acc hay multiple values, each for one region
# chr.acc and chr.unacc must have the same length
# Jan 30: Pr_switch, barcodeLength, chr.acc, chr.unacc   are ALL arrays, with N elements, where N is number of regions
simMemoirStrdist<-function(nGen=4,mu=0.4,alpha=0,barcodeLength=10,methods=c(),
                          simulationType='trit',recType="epimemoir",nIntegrases=1,
                          chr.acc = c(0.8,0.8),chr.unacc =c(0.1,0.1),Pr_switch =c(0,1/8)){
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



   # # # # # # # # CELL OBJECT INIT
   # # # # # # # #

 #convert numeric vectors to string for the tree to work properly
 chr = paste(chr.acc,sep="",collapse="_")
 chr.closed = paste(chr.unacc, sep="",collapse= "_")

 #prepare the barcodes for different regions bc.regions is an array of lengths
 if(length(barcodeLength)==1) bc.regions = rep(barcodeLength,length(chr.acc)) else bc.regions = barcodeLength

 #Jan 30th   Pr_switch will now be a new attribute of firstCell
 if(length(Pr_switch)==1) switches = paste(rep(Pr_switch,length(chr.acc)),collapse="_") else switches = paste(Pr_switch,collapse="_")

 #create cell object
 firstCell<-Node$new("1");


 #replicate barcodes with same length
  ## firstCell$barcode <-paste(  rep(paste(rep("u",barcodeLength),collapse=""),length(chr.acc)),collapse="_")   ;
 #different length barcodes for each region:
 firstCell$barcode<-do.call(paste,c(lapply(   lapply( bc.regions, rep, x="u"), paste, collapse=""),sep="_" )  )

 # these number will multiply the basal edit rate
 firstCell$chr.acc = chr # accesibility values for open state
 firstCell$chr.closed = chr.closed #accessibility value for closed state

 firstCell$epihistory=paste(rep("1",length(chr.acc)),collapse="_") #initial state for the fouding cell of the colony
  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

  firstCell$switches = switches

  # SIMULATION
  # # # # # # # # #
   # # # # # # # #
  # # # # # # # # #

  #this means that we want as many "channels" (integrases) as chromatin regions:

  #if the recording type is epimemoir, the following line does nothing: nInterases=1, and act_time = c(1,1,1,...)
  if(recType =="epimemoir") nIntegrases =1
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








  allDistances = array() # just to keep the indexes in the output consistent with further analysis

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
