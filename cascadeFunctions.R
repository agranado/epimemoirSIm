os=system("cat ../os.txt",intern = T) #Local Mac repository (laptop)
if(os=="mac"){
  git.path = "/Users/alejandrog/MEGA/Caltech/trees/GIT/"
}else if(os =="linux"){ #local linux subsystem (maybe under Windows10)
  git.path = "/home/agranado/MEGA/Caltech/trees/lineageSoftware/"

}else if(os=="aws"){ #AMAZON cloud computing server

  git.path="../lineageSoftware/"
}


#Load fucntions from the main GIT repository (lineageSoftware)
source(paste(git.path,"MLfunctions.R",sep=""))


# For a given number of generations and a given number of integrases (or recording mechanisms, gRNA etc)
# we need to indicate at which point in time (generation units), the next element will become active.
cascadeActivation<-function(nGen,nIntegrases){
    c( nGen %% nIntegrases,ceiling(nGen/nIntegrases),  nGen%/%nIntegrases)
    act_time =array()
    a =nGen %% nIntegrases
    b =ceiling(nGen/nIntegrases)
    c = nGen%/%nIntegrases

    if(a>0){
      for (i in 1:a){
          act_time =c(act_time,rep(i,b))
      }
      act_time = act_time[-1]
      act_time2=array()
      for(i in (a+1):nIntegrases){
          act_time2 =c(act_time2, rep(i,c))
      }
      act_time2 = act_time2[-1]
      act_time = c(act_time,act_time2)
    }else{

      for(i in 1:nIntegrases){
        act_time = c(act_time,  rep(i,c))
      }
      act_time = act_time[-1]
    }

    return(act_time)
  }

 partitionBarcodes <- function (fullBarcode,totalInts, currentInts){
   barcodeLength = length(fullBarcode)
   editableIndx =(barcodeLength / totalInts *(currentInts-1) +1) : (barcodeLength /totalInts * currentInts);
   a = fullBarcode[editableIndx]
   return(list(a,editableIndx))
}
library(phytools) #for pasting trees together

#this function works well for nIntegrases=2
#for more integrases, it will require recursive reconstruction
#We might not need more than 4 integrases for the paper (or in reality)
cascadeReconstruction<-function(barcodeLeaves,totalInts,currentInts,nGen,mu,alpha){

    # 1 Partition barcode into the different integrase-specific elements
      ##  barcodeLeaves= this.tree$tip.label
      #let's separate the barcode string into arrays, so we can index using the integrase indexes.
      # te=apply(as.array(barcodeLeaves),1,strsplit,"")
      # all.barcodes =as.array(unlist(te))
      # dim(all.barcodes)<-c(nchar(barcodeLeaves[1]) ,length(barcodeLeaves)) #Colums are the cells,
      # # all.barcoes[,1] == barcodesLeaves[1]

      #indexes from the function are wrong, so the first 1:n-1 elements are used then n + 1 :end
      # FIX THIS

      # for(t in 1){
      #     sub.barcodes= matrix(0,dim(all.barcodes)[2],length(partitionBarcodes(all.barcodes[,1],totalInts,1)))
      #     for(c in 1:dim(all.barcodes)[2]) { #for each cell
      #         sub.barcodes[c,]= partitionBarcodes(all.barcodes[,c],totalInts,t)
      #
      #     }
      #
      #   # 2 Reconstruct each subtree
      # }
      act_time = cascadeActivation(nGen,totalInts)

      #calcute the index then use substr
      #this function returns both the sub.barcode and the indexes
      #here we start with the first subtree, using only the first generations, therefore currentInts==1
      editableIndx=partitionBarcodes(strsplit(barcodeLeaves[1],"")[[1]],totalInts,currentInts)[[2]]
    #  editableIndx= editableIndx[-1]
      #
      sub.barcodes = substr(barcodeLeaves,range(editableIndx)[1],range(editableIndx)[2])

      #generate a small subtree based on unique profiles for the first generations
      unique.sub.barcodes = unique(sub.barcodes)

      sub.nGen = sum(act_time ==currentInts)
      #for very slow recording it might be that there are not edits during the
      #first integrase, and so all cells have same sub.barcode --> problem
      #a single if should do the job

      if(length(unique.sub.barcodes)>1){
          matdist_ = manualDistML(unique.sub.barcodes,mu,alpha,sub.nGen)
          colnames(matdist_)<-unique.sub.barcodes
      #    manualTree_1 =upgma(as.dist(t(matdist_))) #now the tree has names
          manualTree_1= as.phylo(hclust(as.dist(t(matdist_))))


        }else{ #DOES NOT WORK SO FAR
          firstCell=firstCell<-Node$new("1"); firstCell$barcode <-paste(rep("u",barcodeLength),collapse="");
        #  manualTree_1 = as.phylo(firstCell) #here we need to return a root for a tree
          #we need to add the sub tree where the parent leaves correspond, if all the sub.barcodes are the sample
          # then we should just put them randombly to any node, I thought about making only one subnode but that doesn't works
          #we can make a root node wiht N childs, where each child correspond to a parent leave wiht no identity

          #PROBLEM: as.phylo.Node does not want to convert the root only from Node->phylo

          #FIX 1: 5-Oct-2018
          #Since all the subtrees can't be related to each other (they all have same edits in 1st integrase)
          #we will just put them all together as if coming from the first node directly (which is true)
          #however, at this point we don't know how many subtrees there are, so we will just aproxximate a number
          #  for(jj in 1:sub.nGen^2) { firstCell$AddChild(jj+1) } #FIX THIS
          #  as.phylo.Node(firstCell)

          #FIX 2: make dummy node as an extension of the root (WORKS)
          dummy.root = firstCell$AddChild(unique.sub.barcodes[1]) #add dummy node such that conversion to phylo works..
          dummy.root$barcode= unique.sub.barcodes[1]
          manualTree_1 = as.phylo.Node(firstCell)
          #fix the length of the edge
          manualTree_1$edge.length  =0.01 #this is the order of the reconstructed edges
        }
      #this works
      big.tree = manualTree_1
      # for each unique.sub.barcode, get all the daughters coming from that clone
      for(i in 1:length(unique.sub.barcodes)){
        daughter.index= grep(unique.sub.barcodes[i],sub.barcodes) #get all the cells that have same barcode for first integrase (theu belong to the same clone)
        daughters.barcodes = barcodeLeaves[daughter.index] #get the barcodes

        matdist_2 = manualDistML(daughters.barcodes,mu,alpha,sub.nGen) # submatrix for distances
        colnames(matdist_2)<-daughters.barcodes # name the matrix with the actual barcodes, will be used later for distance between the real and reconstructed trees
        #manualTree_2 =upgma(as.dist(t(matdist_2))) #now the tree has names
        manualTree_2 =as.phylo(hclust(as.dist(t(matdist_2))))

        # AND REPEAT (recursively)
        # 3 Join subtrees into a big tree
        which.leaf=which(big.tree$tip.label==unique.sub.barcodes[i])
        big.tree$tip.label[which.leaf]<-"NA"
        manualTree_2$root.edge<-0
        big.tree = paste.tree(big.tree,manualTree_2)
      }



      #using phylo tools to join elements into a branch/leaf



        return(big.tree)

  }



library(RColorBrewer)


  heatmap.simple<-function(mat1,xlab="",ylab="",this.min, this.max,title=""){

        # this.min =min(mat1)
        # this.min = this.min- 0.1*this.min
        #
        # this.max = max(mat1) + 0.1*this.min
        #

        how.manycolors = 5
         my_palette <-  colorRampPalette(brewer.pal(11,"Spectral"))(n = 3*how.manycolors-1)
         this.range = this.max-this.min


         # (optional) defines the color breaks manually for a "skewed" color transition
         col_breaks = c(seq(this.min, this.min+0.33*this.range,length=how.manycolors),  # for red
                        seq( this.min+0.34*this.range, this.min +  0.66*this.range,length=how.manycolors),           # for yellow
                        seq(this.min + 0.67*this.range, this.min + 1*this.range,length=how.manycolors))
        ##### WORKS: two heatmaps side by side


        heatmap.2(mat1,trace='none',dendrogram='none',col=my_palette,breaks = col_breaks,Rowv=F,
                  Colv=F,key=T,xlab = xlab, ylab =ylab,main = title,keysize=2)


  }

  heatmap.compare<-function(mat1,mat2){

    this.min =min(rbind(mat2,mat1))
    this.min = this.min- 0.1*this.min

    this.max = max(rbind(mat2,mat1)) + 0.1*this.min


     my_palette <-  colorRampPalette(brewer.pal(11,"Spectral"))(n = 299)
     this.range = this.max-this.min
     # (optional) defines the color breaks manually for a "skewed" color transition
     col_breaks = c(seq(this.min, this.min+0.33*this.range,length=100),  # for red
                    seq( this.min+0.34*this.range, this.min +  0.66*this.range,length=100),           # for yellow
                    seq(this.min + 0.67*this.range, this.min + 1*this.range,length=100))
    ##### WORKS: two heatmaps side by side
    library(gridGraphics)
    library(grid)
    heatmap.2(mat1,trace='none',dendrogram='none',col=my_palette,breaks = col_breaks,Rowv=F,Colv=F,key=T)
    library(gridGraphics)
    grab_grob <- function(){
      grid.echo()
      grid.grab()
    }

    g <- grab_grob()
    grid.newpage()

    heatmap.2(mat2,trace='none',dendrogram='none',col=my_palette,breaks = col_breaks,Rowv=F,Colv=F,key=F)
    g2 <- grab_grob()
    #grid.newpage()
    # library(gridExtra)
    # grid.arrange(g,g, ncol=2, clip=TRUE)

    lay <- grid.layout(nrow = 1, ncol=2)
    pushViewport(viewport(layout = lay))
    grid.draw(editGrob(g, vp=viewport(layout.pos.row = 1,
                                      layout.pos.col = 1, clip=TRUE)))
    grid.draw(editGrob(g2, vp=viewport(layout.pos.row = 1,
                                      layout.pos.col = 2, clip=TRUE)))
    upViewport(1)

  }

  single.integrase.reconstruction<-function(barcodeLeaves,nGen=4,mu=0.4,alpha=1/2){

          matdist_ = manualDistML(barcodeLeaves,mu,alpha,nGen)
          colnames(matdist_)<-barcodeLeaves
          hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
          return(hclust.tree)
  }

  library(data.tree)
  library(phangorn)
  reconstruct.list<-function(cascadevar.tree,integrases,barcodes,mus,generations,alpha=1/2){


      matrices=list()

      for(ca in 1:length(integrases)){
        distBitAll = array(0,dim=c(length(barcodes),length(generations),length(mus)))
        for(mIdx in 1:length(mus)){
          print( paste("int ",toString(integrases[ca])," mu=",toString(mus[mIdx]),"\\n" ))
          for(bc in 1:length(barcodes)){
            for(ng in 1:length(generations)){

                tree.set = cascadevar.tree[[ca]][[mIdx]][[1]][[bc]][[ng]]

              list.dist=matrix(0,length(tree.set),2);
              for(i in 1:length(tree.set)){
                barcodeLeaves=tree.set[[i]]$tip.label;
                if(ca==1){

                  r = single.integrase.reconstruction(barcodeLeaves,nGen=generations[ng],mu = mus[mIdx],alpha=1/2)
                }else{
                  r=cascadeReconstruction(barcodeLeaves,totalInts=2,currentInts=1,nGen=generations[ng],mu=mus[mIdx],alpha =1/2);
                }


                list.dist[i,1] = RF.dist(r,tree.set[[i]]);
                list.dist[i,2]=RF.dist(r,tree.set[[i]],normalize =T)
              }

              distBitAll[bc,ng,mIdx] = mean(list.dist[,2])
            }
          }
        }
        matrices[[ca]] = distBitAll
      }
      return(matrices)
}


reconstruct.tree.list <-function(treelist,mu,nGen,alpha=1/2){
    reconstructed.trees = list()
    for(t in 1:length(treelist)){
      this.tree = treelist[[t]]
      barcodeLeaves=this.tree$tip.label

      matdist_ = manualDistML(barcodeLeaves,mu,alpha,nGen)
      colnames(matdist_)<-barcodeLeaves
      hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
      reconstructed.trees[[t]] = hclust.tree


    }
    return(reconstructed.trees)
}



#plot two trees side by side, usually one tree is reconstruction and the other one is ground.truth
compare.trees<-function(ground.truth,rec){

    x11();
    par(mfrow = c(1,2))

    #rename with leave number
    new.names = paste(names(ground.truth$tip.label),ground.truth$tip.label)
    ground.truth$tip.label<-new.names
    rec$tip.label <-new.names

    plot.phylo(ground.truth,show.node.label=T)
    plot(rec)



}
