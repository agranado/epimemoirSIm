

#Oct 9th,2018
#First round of cascade simulations is Done!
#Here we will generate basic plots for comparing 1 vs 2 integrases as independent channels using a new reconstruction method
#the goal:show as a proof of principle that by applying smarter memory usage we can optimize lineage recording


#this script should load the object muVariation or it should be executed afeter running bitVStrit.R
os=system("cat ../os.txt",intern = T) #Local Mac repository (laptop)
if( os=="linux"){

  #data.path = "/home/alejandrog/MEGA/Caltech/lineage/simulation_data/" #Location of simulation data from AWS
  data.path = "./"
}else if(os=="mac"){

  data.path = "./"
}
#parameters from AWS: large object cascadevar 10-Sep-2018
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



#large file with all the simulations: distances and trees.
# 1 muVariation object per integrase number, for cascade recording
file.name.cascade1="epi_intVar_int_1_.rda"

load(paste(data.path,file.name.cascade1,sep=""))

#rm("cascadevar1","cascadevar2") #DATA LOADED CORRECTLY AGC., Oct 11th

# # # # # # # # # # #
 # # # # # # # # # # # START PLOTS
# # # # # # # # # # #
n=generations
rand.dist = 2^(n+1)-6 #according to RF.dist documentation is the maximum RF distance for binary-unrooted trees
# create the data structure
a =rand.dist
dist.idx= c(2,3) #for integrase==1 there is only one distance, int==2, only distance is id=3
compare.cascade= list()
BCn=7
#for each number of integrases get the matrix:
for( ca in 1:length(integrases)){

    #matrix to fill with data, for muVariation
    distBitAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))

    muVariation = cascadevar[[ca]]
    #3D matrix with variation of barcode number, edit rate & generations
    for(ng in 1:length(generations)){
      for (mIdx in 1:length(mus)){
        simType=muVariation[[mIdx]]
        for(bc in 1:length(barcodes)){
          distBitAll[bc,ng,mIdx]= 1-apply(simType[['trit']][[bc]][[ng]],2,mean)[dist.idx[ca]]/a[ng]
        }
      }
    }

    compare.cascade[[ca]]=distBitAll
}


par(mfrow=c(2,3))
for(ng in 1:length(generations)){
  distTritAll = compare.cascade[[1]]
  plot(mus,distTritAll[BCn,ng,],ylim=c(0,1),
       type="o",col="blue",ylab="norm dist",xlab="edit rate p/site p/gen",
       main=paste("g=",toString(generations[ng]),sep=""),
       cex.axis=1.5,cex.lab=1.6,cex.main=2);
#  distBitAll = compare.cascade[[1]]
#  lines(mus,distBitAll[BCn,ng,],type="o")
}

# set names for axis:
# dimnames(distBitAll)[[1]]<-barcodes
# dimnames(distBitAll)[[2]]<-generations
# dimnames(distBitAll)[[3]]<-mus

dimnames(distTritAll)[[1]]<-barcodes
dimnames(distTritAll)[[2]]<-generations
dimnames(distTritAll)[[3]]<-mus



#optimal reconstructability
optimTritAll = apply(distTritAll,c(1,2),max)
optimBitAll = apply(distBitAll,c(1,2),max)

row.names(optimBitAll)<-as.character(barcodes)
colnames(optimBitAll)<-as.character(generations)

row.names(optimTritAll)<-as.character(barcodes)
colnames(optimTritAll)<-as.character(generations)


 #####


 my_palette <-  colorRampPalette(brewer.pal(11,"Spectral"))(n = 299)

 # (optional) defines the color breaks manually for a "skewed" color transition
 col_breaks = c(seq(0.4,0.6,length=100),  # for red
                seq(0.61,0.8,length=100),           # for yellow
                seq(0.81,1,length=100))
##### WORKS: two heatmaps side by side
library(gridGraphics)
library(grid)
heatmap.2(optimBitAll,trace='none',dendrogram='none',col=my_palette,breaks = col_breaks,Rowv=F,Colv=F,key=F)
library(gridGraphics)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

g <- grab_grob()
grid.newpage()

heatmap.2(optimTritAll,trace='none',dendrogram='none',col=my_palette,breaks = col_breaks,Rowv=F,Colv=F,key=F)
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
