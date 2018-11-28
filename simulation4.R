


library(data.tree)






#this function creates a single cell-division event asumming constant division rate
#it also duplicates the barcoded scratchpad from the parent cells
#takes thisNode as input parameter
#  thisNode <- Node$new("1")
divideCell2<-function(thisNode,mu,alpha,type='trit',recType="integrase",cInts=1,tInts=1,Pr_switch =1/8){


  x<-thisNode$levelName
  #match the characters at the end of the node's name
  #by default data.tree names nodes if format:  " ¦   ¦   ¦   ¦--16"
  mx<-regexpr("\\d+$",x,perl=T)
  cellIndex<-substr(x,mx,mx + attr(mx,"match.length")-1)

  cellIndex1<-as.numeric(cellIndex)*2
  daughter<-thisNode$AddChild(cellIndex1)
  cellIndex2<-as.numeric(cellIndex1)+1
  daughter2<-thisNode$AddChild(cellIndex2)


  # we can mutate during the lifetime of the daughters.
  #first, daughters inherit
  daughter$barcode <-thisNode$barcode
  daughter2$barcode <-thisNode$barcode
  #then mutation happens:
  #independent events for each daughter. (extract barcode, mutate, update)
  if(recType=="integrase"){
    daughter$barcode <- mutationScratchpad(daughter$barcode,mu,alpha,type,recType,cInts=cInts,tInts=tInts,chr.acc = c())
    daughter2$barcode <- mutationScratchpad(daughter2$barcode,mu,alpha,type,recType,cInts=cInts,tInts=tInts ,chr.acc=c())
    # NEW SECTION EPIMEMOIR
  }else if( recType=="epimemoir"){
    #insert mutation fuction for epimemoir (in the future should be the same)
    #open accesibility value (probability that site is accessible for recording enzyme)
    daughter$chr.acc <-thisNode$chr.acc
    daughter2$chr.acc <-thisNode$chr.acc #this is an array with the same length as the N of recording sites (channels)

    #inherit the closed "accessibility" value
    daughter$chr.closed<-thisNode$chr.closed
    daughter2$chr.closed<-thisNode$chr.closed
    #inherit the epigenetic history of the mother
    daughter$epihistory  = thisNode$epihistory
    daughter2$epihistory = thisNode$epihistory
    #ADD the sate of chromating for the new generation: THIS WILL AFFECT THE RECORDING from the current generation
    #So far works for 1 region only:
  #  Pr_switch = 1/8
    daughter$epihistory =  epigenticTransition(  daughter$epihistory,Pr_switch) #mutate and record te state in the history
    daughter2$epihistory = epigenticTransition(  daughter2$epihistory,Pr_switch)
    n.events = nchar(daughter$epihistory)
    daughter.isopen = substr(daughter$epihistory,n.events,n.events) #extract the current (new)state
    daughter2.isopen = substr(daughter2$epihistory,n.events,n.events)

    #works only for 1 region ! ! # # # # # # # #
      #based on the array of state 100100101 (Nregions) etc we will generate a new arrac of chr.acc which will correspond to the Pr values of accessibility
      #form 0.8_0.1_0.1_ (Nregions)
      #These are strings:
    if(daughter.isopen =="1") new.chr = daughter$chr.acc else new.chr =daughter$chr.closed #this is the numeric value of chromatin rate (which is assigned after the transition)

    if(daughter2.isopen =="1") new2.chr = daughter2$chr.acc else new2.chr = daughter2$chr.closed

    as.numeric(strsplit(daughter$chr.acc,"_")[[1]])
    as.numeric(strsplit(daughter2$chr.acc,"_")[[1]])
    daughter$barcode <- mutationScratchpad(daughter$barcode,mu,alpha,type,recType,cInts=0,tInts=tInts,as.numeric(strsplit(new.chr,"_")[[1]]))
    daughter2$barcode <- mutationScratchpad(daughter2$barcode,mu,alpha,type,recType,cInts=0,tInts=tInts , as.numeric(strsplit(new2.chr,"_")[[1]]))

    # Function for chromatin transitions:
    #Same way as mutationScratchpad, we will transition between open and closed states with some probability

  }

  return(thisNode)
}


epigenticTransition<-function(epihistory,Pr_switch){
    single.transition = 1
    char.history = strsplit(epihistory,"")[[1]]
    switch_event<-runif(1,0,1)
    last.event = char.history[length(char.history)]
    #FOR SINGLE TRANSITIONS
    #if all the previous states are the same, means that the transition is only allowed to happen once!
    #if unique(char.history)>1 then a transition already happened in the past and we skip the rest of the function
    if(single.transition) make.transition = length(unique(char.history)) else make.transition = 1 # if single transition is not a problem then always TRUE
    #Either way make.transition MUST eq 1 for the expression to happen

    if(make.transition==1 & switch_event<Pr_switch) if(last.event =="0") new.state ="1" else new.state="0" else new.state = last.event

    updated.history =paste(paste(char.history,collapse=""),new.state,sep="")

    return ( updated.history)
}
#
# switch_epi<-function(pre.state){
#     if(pre.state =="0") new.state ="1" else new.state="0"
#
#     return new.state
# }

#this function acts on a tree.
#finds the depth of the tree and duplicates only at the leaf level
#Effectively add one generation to the tree
divideCellRecursive2<-function(thisNode,mu,alpha,type='trit',recType="integrase",cInts=1,tInts=1,Pr_switch =1/8){
  #add Child function changes the tree permanently.
  #tree works as a global variable, all pointers represent same object

  #this function will add one generation to the three
  if(length(thisNode$children)==0){
    #the node is a leaf (or root)
    divideCell2(thisNode,mu,alpha,type,recType,cInts,tInts,Pr_switch)
    #this works
  }else{
    for (i in 1:length(thisNode$children)){
      divideCellRecursive2(thisNode$children[i][[1]],mu,alpha,type,recType,cInts,tInts,Pr_switch)
    }

  }

}


#function to mutate the scratchpad (what memoir actually does)
# Different recording schemes should be easy to implement as recType and adding
#   new functions
mutationScratchpad <- function(barcode,mu,alpha,type='trit',recType="integrase",cInts=1,tInts=1,chr.acc=c()){

  originalBarcode=strsplit(barcode,"")[[1]]
  #now a is an array of char elements representing the scratchpad
  m_event<-runif(1,0,1)
  #mutation happens (rate = mutations/generation)
  #for each element in the array, ask if mutation happens
  #with constant independent rate (mutation in [1] does not affect [2])
  if(recType=="integrase"){
    mutatedBarcode=recIntegrase(originalBarcode,mu=mu,alpha=alpha,type=type,currentInts=cInts,totalInts=tInts)
  }else if(recType=="epimemoir"){ #epimemoir
    mutatedBarcode=recEpiMemoir(originalBarcode,mu=mu,alpha=alpha,type=type,totalInts=tInts,chr.acc=chr.acc)
  }else { #no mutation
    mutatedBarcode= originalBarcode
  }
  #returns the mutated string
  return( paste(mutatedBarcode,sep="",collapse=""))
}


# currentInts --> the active integrase, in a cascade they are activated sequentially
# in normal iMEMOIR, there is only one integrase and it is activated throughout the experiment
# if no additional arguments are passed, the function will perform integrase recording 1.0
recIntegrase<-function(a,mu,alpha,type,currentInts=1,totalInts=1){

  #barcodeLength=length(a)
  fullBarcode = a
  #we want to edit only the portion that corresponds to
  #the currently active integrase
  aa=partitionBarcodes(fullBarcode, totalInts,currentInts)
  a=aa[[1]] #the sub.barcode
  editableIndx = aa[[2]] #the indexes for the sub.barcode

  for (c in 1:length(a)){
    #only mutate on unchanged elements
    if(a[c]=="u"){
      m_event<-runif(1,0,1) #random number
      if(m_event<mu){
        #mutually exclusive events
        if(type=='trit'){
          trans_pr <-runif(1,0,1) #random number

          if(trans_pr<alpha){
            a[c]="r"
          }else{
            a[c]="x"
          }
        }else if(type=='binary'){
          a[c]="x"
        }
      }
    }
  }

  fullBarcode[editableIndx]= a
  return(fullBarcode)
}

#NOV 16th
# This function needs the aruments for chromatin accessibilityo for each recording element
#this function works at the single cell leves: each cell will have it's own chr.acc vector
recEpiMemoir<-function(a,mu,alpha,type,totalInts=1,chr.acc=c()){

  #barcodeLength=length(a)
  fullBarcode = a
  #we want to edit only the portion that corresponds to
  #the currently active integrase
  #DEFAULT BEHAVIOR:
  if(length(chr.acc)==0){
    chr.acc=rep(1,totalInts)
  }
  #in epimemoir we need to cycle through all the currentInts
  #each currentInt will be a gRNA target region and will have it's own associated chr.acc value
  #in epimemoir we need to assume that even the closed regions will have chr.acc > 0
  #so we can not just skip them, but instead go through all with a FOR loop
  #we can keep the same function partitionBarcodes()
  for(t in 1:totalInts){
      mu_epi = mu * chr.acc[t]

      aa=partitionBarcodes(fullBarcode, totalInts,currentInts=t)
      a=aa[[1]] #the sub.barcode
      editableIndx = aa[[2]] #the indexes for the sub.barcode

      for (c in 1:length(a)){
        #only mutate unchanged elements
        if(a[c]=="u"){
          m_event<-runif(1,0,1) #random number
          if(m_event<mu_epi){
            trans_pr <-runif(1,0,1) #random number
            #mutually exclusive events
            if(type=='trit'){
              if(trans_pr<alpha){
                a[c]="r"
              }else{
                a[c]="x"
              }
            }else if(type=='binary'){
              a[c]="x"
            }
          }
        }
      }
      fullBarcode[editableIndx]= a #save the edited part of the barcode to the original full-lenght string and iterate to the next editable region
    }

  return(fullBarcode)
}






#OLD original backup: WORKS function to mutate the scratchpad (what memoir actually does)
#mu is the main mutation rate: pr that a mutation happens
#if a mutation happens, it will be u->x with pr a & u->r with pr b
mutationScratchpad_ <- function(barcode,mu,alpha,type='trit'){

  a=strsplit(barcode,"")[[1]]
  #now a is an array of char elements representing the scratchpad
  m_event<-runif(1,0,1)
  #mutation happens (rate = mutations/generation)
  #for each element in the array, ask if mutation happens
  #with constant independent rate (mutation in [1] does not affect [2])
  for (c in 1:length(a)){
    #only mutate on unchanged elements
    if(a[c]=="u"){
      m_event<-runif(1,0,1) #random number
      if(m_event<mu){
        trans_pr <-runif(1,0,1) #random number
        #mutually exclusive events
        if(type=='trit'){
          if(trans_pr<alpha){
            a[c]="r"
          }else{
            a[c]="x"
          }
        }else if(type=='binary'){
          a[c]="x"
        }
      }
    }
  }
  #returns the mutated string
  return( paste(a,sep="",collapse=""))
}

#convert the simulated tree to fast format.
#from the fasta format we can perfom the alignment.
#this function will add a dummy character at the end of the sequence (I think this might help with aligment but might be useless)
convertSimToFastaPlusOne <- function(barcodes){
  barcodeNames<-names(barcodes)
  fas=array()
  for (i in 1:length(barcodes)){
    #each element of the array contains a line of the fasta file
    fas[i] = paste(">",barcodeNames[i],"_",barcodes[i],"\n",barcodes[i],"R",sep="")
  }
  return(fas)
}
#Original function, only converts the barcodes to fasta format (Used everywhere)
convertSimToFasta <- function(barcodes){
  barcodeNames<-names(barcodes)
  fas=array()
  for (i in 1:length(barcodes)){
    #each element of the array contains a line of the fasta file
    fas[i] = paste(">",barcodeNames[i],"_",barcodes[i],"\n",barcodes[i],sep="")
  }
  return(fas)
}

#Apr4
convertSimToBinary <- function(barcodes){
  barcodeNames<-names(barcodes)
  fas=array()
  for (i in 1:length(barcodes)){
    #each element of the array contains a line of the fasta file
    strSplit=strsplit(barcodes[i],"")
    binaryBarcode= array()
    for (c in 1:length(strSplit[[1]])){
      if(strSplit[[1]][c]=="u"){
        binaryBarcode[c]="00"
      }else if (strSplit[[1]][c]=="r"){
        binaryBarcode[c] = "11"
      }else if (strSplit[[1]][c]=="x"){
        binaryBarcode[c]= "01"
      }

    } #end of binary barcode

    fas[i] = paste(">",barcodeNames[i],"_",barcodes[i],"\n",paste(binaryBarcode,sep="",collapse=""),sep="")
  }
  return(fas)
}




#phangorn complains about the first line (which indicates the special alphabet condition)
#so we reomve that
convertSimToPhylipBinary <- function(barcodes){
  barcodeNames<-names(barcodes)
  taxN = length(names(barcodeLeaves)) #need to speficy in phylip format
  seqLen= nchar(barcodes[1]) #length of the sequences, needed for phylip
  alphabet = "01"
  #Phylip format looks like :
  #   #SPECIALALPHABET
  #   32 6 urx
  #   cell1 urxuur
  #   cell2 uuuxxr
  #   ...

  fas=array()
  #fas[1]= "#SPECIALALPHABET" #first line of the file
  #fas[2]= paste(toString(taxN)," ",toString(seqLen*2), " ", alphabet,sep = "") # second line of the file
  fas[1]= paste(toString(taxN)," ",toString(seqLen*2),sep = "") # second line of the file

  for (i in 1:length(barcodes)){
    #each element of the array contains a line of the fasta file


    strSplit=strsplit(barcodes[i],"")
    binaryBarcode= array()
    for (c in 1:length(strSplit[[1]])){
      if(strSplit[[1]][c]=="u"){
        binaryBarcode[c]="00"
      }else if (strSplit[[1]][c]=="r"){
        binaryBarcode[c] = "11"
      }else if (strSplit[[1]][c]=="x"){
        binaryBarcode[c]= "01"
      }

    } #end of binary barcode

    fas[i+1] = paste(barcodeNames[i],"_",barcodes[i],"\t",paste(binaryBarcode,sep="",collapse=""),sep="")
  }
  return(fas)
}



convertSimToPhylip <- function(barcodes){
  barcodeNames<-names(barcodes)
  taxN = length(names(barcodeLeaves)) #need to speficy in phylip format
  seqLen= nchar(barcodes[1]) #length of the sequences, needed for phylip
  alphabet = "uxr"
  #Phylip format looks like :
  #   #SPECIALALPHABET
  #   32 6 urx
  #   cell1 urxuur
  #   cell2 uuuxxr
  #   ...

  fas=array()
  fas[1]= "#SPECIALALPHABET" #first line of the file
  fas[2]= paste(toString(taxN)," ",toString(seqLen), " ", alphabet,sep = "") # second line of the file
  for (i in 1:length(barcodes)){
    #each element of the array contains a line of the fasta file
    fas[i+2] = paste(barcodeNames[i],"_",barcodes[i],"\t",barcodes[i],sep="")
  }
  return(fas)
}


#Once the fasta file is created, we substitute letter to make the file
#look like a DNA sequence (natural MEMOIR letters u,r,x)
#mac computers we can use the sed command
# sed -i .bak 's/R/t/g' protAlign.fasta
# sed -i .bak 's/R/t/g' protAlign24.fasta
# sed -i .bak 's/x/c/g' protAlign24.fasta
# sed -i .bak 's/u/a/g' protAlign24.fasta
# sed -i .bak 's/r/g/g' protAlign24.fasta
convertMemoirToDNA <-function(fastaIN){
  cmds=array()

  #OS detection: SED command behaves differently in OS and linux platforms,
  #if Mac (local branch, Alejandro's laptop)
  #The file is NOT located in the GIT folder but in a higher folder such that it is not override or cause conflicts
  #when merging the GIT repository
  os=system("cat ../os.txt",intern = T)

  if(os=="mac"){
    cmds[1]="sed -i .bak 's/R/c/g'"
    cmds[2]="sed -i .bak 's/x/t/g'"
    cmds[3]="sed -i .bak 's/u/g/g'"
    cmds[4]="sed -i .bak 's/r/a/g'"
  }else if(os=="linux" | os== "aws"){ #AWS server or any other linux machine (This should work)

    cmds[1]="sed -i.bak 's/R/c/g'"
    cmds[2]="sed -i.bak 's/x/t/g'"
    cmds[3]="sed -i.bak 's/u/g/g'"
    cmds[4]="sed -i.bak 's/r/a/g'"

  }
  #sometimes (I think because small delays in writting/reading files) the file is not found and the script might crash (eventhough the file is in the folder)
  #We can check whether the file is there and execute the system commands only if this is true
  if(file.exists(fastaIN)){
    for(i in 1:length(cmds)){
      system(paste(cmds[i],fastaIN,sep=" "))
    }
    return(1)
  }else{
    #if the file does not exist we can read this value and return NaN (or something else)
    return(0)
  }


}


# Rename the tip labels: remove the sequence from the label
removeSeqLabel <-function(treeUPGMA_noSeq){
  for (i in 1:length(treeUPGMA_noSeq$tip.label)){
    label1=treeUPGMA_noSeq$tip.label[i]
    #take the number ID from the format   23_ACGACGAGGCAAAG
    #it could be one or more digits, so we use a regular expression
    mx<-regexpr("^\\d+",label1,perl=T)
    treeUPGMA_noSeq$tip.label[i]=substr(label1,mx,mx + attr(mx,"match.length")-1)
    #treeUPGMA_noSeq$tip.label[i]=substr(label1,1,2)
  }
  return(treeUPGMA_noSeq)
}


randomSeqLabel <-function(treeUPGMA_noSeq){
  #duplicate the tree
  randomTree= treeUPGMA_noSeq
  newLabels=array()
  for (i in 1:length(treeUPGMA_noSeq$tip.label)){
    label1=treeUPGMA_noSeq$tip.label[i]
    #take the number ID from the format   23_ACGACGAGGCAAAG
    #it could be one or more digits, so we use a regular expression
    mx<-regexpr("^\\d+",label1,perl=T)
    newLabels[i] =substr(label1,mx,mx + attr(mx,"match.length")-1)
    #treeUPGMA_noSeq$tip.label[i]=substr(label1,1,2)
  }
  randLabels=sample(newLabels)
  for (i in 1:length(randomTree$tip.label)){
    randomTree$tip.label[i]=randLabels[i]
  }
  return(randomTree)
}

reverseLabels = function(treeUPGMA){

  newTips =gsub("g","u",treeUPGMA$tip.label)
  newTips =gsub("t","x",newTips)
  newTips =gsub("a","r",newTips)
  treeUPGMA$tip.label = newTips

  return(treeUPGMA)
}

randomTrueTree <-function(treeUPGMA_noSeq){
  #duplicate the tree
  randomTree= treeUPGMA_noSeq
  newLabels=array()
  for (i in 1:length(treeUPGMA_noSeq$tip.label)){
    newLabels[i]=treeUPGMA_noSeq$tip.label[i]
  }
  randLabels=sample(newLabels)
  for (i in 1:length(randomTree$tip.label)){
    randomTree$tip.label[i]=randLabels[i]
  }
  return(randomTree)
}


#Based on profiles of unique barcodes in a tree, we estimate the proportion of each letter
#We can use this for two things 1)estimate the number of "u" hence edit rate and 2) estimate bias between r/x edits
baseFreq = function(unique.BC.matrix){
  #read barcodes
  alphabet= c("u","r","x")
  letterFreq=array(0,dim=c(length(alphabet),dim(unique.BC.matrix)[1]));
  for (a in 1:length(alphabet)){
    for (d in 1:dim(unique.BC.matrix)[1]){ letterFreq[a,d]=length(which(unique.BC.matrix[d,]==alphabet[a]));}
    #average number of unedited barcodes across all cells (we consider only unique barcode profiles)

  }
  meanFreq=apply(letterFreq,1,mean)
  return(meanFreq)
}
