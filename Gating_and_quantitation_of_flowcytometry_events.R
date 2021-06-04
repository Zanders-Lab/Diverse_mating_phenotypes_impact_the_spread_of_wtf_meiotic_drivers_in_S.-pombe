#### This script contains the functions used to collect, filter and quantify the data coming from the ZE5 cytometer.

### Dependencies:
# library(Rmisc)
# library(flowCore)
# library(flowClust)
# library(flowTrans)
# library(tidyverse)
# library(reshape2)
# library(stringr)
# library(cowplot)
# library(gridExtra)
# library(ggplot2)


####This function uses the timming.data to extract the data from standards where single strains where plated.
#### per experiment. Collect only the polygon gate.
#### Inside contains functions:
####<<<<<<ReadFlowSetToDataFrame>>>>>>> Reads flow cytometer file where standards are located.
#### Timing data file has the following structure:
#### no column names but each column has the following composition
# timing.data.1 timing.data2  well.position genotype  cross(biol.rep)  tech.rep initial.freq.for.a.fluor generation
###timing.data.1  = Data produced by cytometer.
###timing.data.2  = Data produced by cytometer.
###well.position  = Well position in axygen plate.
###genotype       = mC (mCherry), G (GFP), *(wtf4+), M(Mendelian), B(Bias)
#Example
#    909984.11  8392902.93    A1          	mC-G-M	    1               	1	      0.5	                      2
#    932984.31  8392232.01    A2	          mC-G-M	    1	                2 	    0.5	                      2
#.
GetPolygonFromGenotype<-function(timing.data,genotype){
  #####Read Flow.set
  ####################################################################3
  ## Read FlowSet
  #ReadFlowSetWithDesc
  ##fcs.file="/n/core/cyto/_Data/Zanders/JLH/17cy1827/A1.fcs"
  ReadFlowSetToDataFrame<-function(fcs.file){
    #Reading fcs.file
    myframe<-read.FCS(file = fcs.file, transformation = TRUE,emptyValue = FALSE)
    #Reading and changing detector names for description. This makes easier to recognize channels with labels of flurophores
    detectors <- as.character(colnames(exprs(myframe)))
    desc <- myframe@description
    desc  <- data.frame(desc)
    desc <- desc[1,]
    desc <- data.frame(t(desc))
    dyesPnS <- data.frame(subset(desc, grepl("P\\d{1,}S", rownames(desc))))
    dyesPnS$IndexCol <- lapply(rownames(dyesPnS), str_extract, pattern = "P\\d{1,}") #get parameter numbers for indexing
    rownames(dyesPnS) <- dyesPnS[,2]
    dyesPnS[,2] <- NULL
    colnames(dyesPnS) <- c("Dye-PnS")
    res1<- as.character(dyesPnS$`Dye-PnS`)
    dyesPnS <- dyesPnS[,1]
    colnames(myframe)<-gsub("-","_",dyesPnS)
    return(myframe)
  }
  
  renamed.genotped=gsub("\\*","",genotype)
  selected.rows=timing.data[gsub("\\*","",timing.data$V4)==renamed.genotped,]
  files=rownames(selected.rows)
  files<-gsub("\\.\\d+","",files)
  files=paste(files,"/",as.vector(selected.rows$V3),".fcs",sep="")
  FlowDataSet<-NULL
  for(i in files){
    print(i)
    FlowDataSet<-rbind(FlowDataSet,
                       as.data.frame(exprs(ReadFlowSetToDataFrame(fcs.file = i)))[,c("GFP_A","mCherry_A")])
  }
  FlowDataSet<-log(FlowDataSet+1,10)
  #### This limits were stablished after running pilot experiments and seeing a
  #### constant channel intensity for each fluorophore.
  ####Defining the limits noColor= no fluorophore
  if(renamed.genotped=="NC"){
    FlowDataSet=FlowDataSet[ (FlowDataSet$GFP_A<6.5 & FlowDataSet$mCherry_A<7)
                             ,]
    hull.points=FlowDataSet[chull(x=FlowDataSet$mCherry_A, y = FlowDataSet$GFP_A),]
  }
  ###Defining the limits GFP= GFP
  if(renamed.genotped=="GFP"){
    FlowDataSet=FlowDataSet[FlowDataSet$GFP_A>6.5 & FlowDataSet$mCherry_A<7,]
    hull.points=FlowDataSet[chull(x=FlowDataSet$mCherry_A, y = FlowDataSet$GFP_A),]
  }
  ###Defining the limits GFP= GFP
  if(renamed.genotped=="mCherry"){
    FlowDataSet=FlowDataSet[FlowDataSet$GFP_A<7 & FlowDataSet$mCherry_A>6.5,]
    hull.points=FlowDataSet[chull(x=FlowDataSet$mCherry_A, y = FlowDataSet$GFP_A),]
  }
  as.list(hull.points)
}
###########################################################################
## ***************Usage example****************##
# GFP<-GetPolygonFromGenotype(timing.data = timing.data,genotype = "GFP")
# mCherry<-GetPolygonFromGenotype(timing.data = timing.data,genotype = "mCherry")
###########################################################################


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

### This function takes everythin into accout.
### This function is the master function that does
### 1.- Take FCS file.
### 2.- Separates events that are unicellular. Defined by limits in 
###     the file for each generation (gates).
### 3.- Separates living and dead cells by DAPI (dead cells).
### 4.- Living cells are then separated by polygons that have 
### GFP,mCherry and Nofluorophores.
### Inside contains functions:
### ####<<<<<<ReadFlowSetToDataFrame>>>>>>> Read flow cytomter file.
### ####<<<<<<ExtractPointsFromPolygon>>>>>>> Collects events insidethe polygon where round cells are located.
### ####<<<<<<ExtractNegativePointsFromOuterPolygon>>>>>>> Collects events outside the polygon where round cells are located.
### ####<<<<<<TransformChannelsToLog>>>>>>> Log transformation for better visual representation of channel intensity.
ProcessData<-function(x,gfp.polygon,mCherry.polygon,negative.polygon,timing.data){
  path<-gsub("\\/\\w\\d+.fcs","",x)
  FilePos<-grep(path,paths)
  Plot<-as.list(NULL)
  #####Read Flow.set
  ####################################################################3
  ## Read FlowSet
  #ReadFlowSetWithDesc
  ##fcs.file="/n/core/cyto/_Data/Zanders/JLH/17cy1827/A1.fcs"
  ReadFlowSetToDataFrame<-function(fcs.file){
    #Reading fcs.file
    myframe<-read.FCS(file = fcs.file, transformation = TRUE,emptyValue = FALSE)
    #Reading and changing detector names for description. This makes easier to recognize channels with labels of flurophores
    detectors <- as.character(colnames(exprs(myframe)))
    desc <- myframe@description
    desc  <- data.frame(desc)
    desc <- desc[1,]
    desc <- data.frame(t(desc))
    dyesPnS <- data.frame(subset(desc, grepl("P\\d{1,}S", rownames(desc))))
    dyesPnS$IndexCol <- lapply(rownames(dyesPnS), str_extract, pattern = "P\\d{1,}") #get parameter numbers for indexing
    rownames(dyesPnS) <- dyesPnS[,2]
    dyesPnS[,2] <- NULL
    colnames(dyesPnS) <- c("Dye-PnS")
    res1<- as.character(dyesPnS$`Dye-PnS`)
    dyesPnS <- dyesPnS[,1]
    colnames(myframe)<-gsub("-","_",dyesPnS)
    return(myframe)
  }
  ####After collection from polygon. I uses a function to take level% of the sample.
  ####The algorith takes level% based on a Box-Cox transformation.
  ####Basically from the whole, takes level% more in the  "center", leaves
  ### outliers.
  ####The function uses the polygon data processed by GetPolygonFromGenotype function
  ####Uses the FSC files collected and the level(% of cells in that polygon)   
  ExtractPointsFromPolygon<-function(polygon.gate,log.transformed.FSC,level){
    ###Polygon Gate
    polygate <- polygonGate(filterId = "Polygon", polygon.gate)
    #FilterFluorPoly <- filter(log.transformed.FSC, polygate)
    SubSetFluor=Subset(x = log.transformed.FSC,subset =  polygate)
    if(nrow(SubSetFluor)>0){
      ### Taking 99% of sample
      xfilter <- tmixFilter("xfilter", c("GFP_A","mCherry_A"), K=1, B=50,level=level) #create filter for 'gating'
      #xf <- filter(SubSetFluor, xfilter)
      ### Subseting to get core color
      Core.Size.Fluor<-Subset(SubSetFluor,xfilter)
      if(nrow(Core.Size.Fluor)>0){
        return(Core.Size.Fluor)
      }else(return(NULL))
    }else(return(NULL))
    
  }
  
  ####After collection from polygon. I uses a function to take level% of the sample.
  ####The algorith takes level% based on a Box-Cox transformation.
  ####Basically from the whole, takes level% more in the  "center", leaves
  ### outliers.
  ####The function uses the polygon data processed by GetPolygonFromGenotype function
  ####Uses the FSC files collected and the level(% of cells in that polygon)   
  ExtractNegativePointsFromOuterPolygon<-function(polygon.gate,log.transformed.FSC,level,outer){
    ###Polygon Gate
    polygate <- polygonGate(filterId = "Polygon", polygon.gate)
    #FilterFluorPoly <- filter(log.transformed.FSC, polygate)
    if(outer=="TRUE"){
      SubSetFluor=Subset(x = log.transformed.FSC,subset =  !polygate)
    }else{SubSetFluor=Subset(x = log.transformed.FSC,subset =  !polygate)}
    if(nrow(SubSetFluor)>0){
      ### Taking 99% of sample
      xfilter <- tmixFilter("xfilter", c("GFP_A","mCherry_A"), K=1, B=50,level=level) #create filter for 'gating'
      #xf <- filter(SubSetFluor, xfilter)
      ### Subseting to get core color
      Core.Size.Fluor<-Subset(SubSetFluor,xfilter)
      if(nrow(Core.Size.Fluor)>0){
        return(Core.Size.Fluor)
      }else(return(NULL))
    }else(return(NULL))
    
  } 
  
  
  FlowSetData<-ReadFlowSetToDataFrame(fcs.file = x)
  list.channels.and.limits<-list("FSC 488/10_A"=c(fsc.488.min,fsc.488.max),
                                 "SSC 488/10_A"=c(ssc.488.min,ssc.488.max))
  eg <- rectangleGate(filterId= "SizeRectGate",list.channels.and.limits)
  SizeGatedFlowSet=Subset(FlowSetData, eg)
  Plot[[1]]<-xyplot(`SSC 488/10_A`~`FSC 488/10_A`, FlowSetData,xbin = 128, smooth = F,
                    filter=eg,stat=TRUE,outline=T,
                    par.settings=list(gate=list(fill="black", alpha=0.2)),
                    gate.text=list(col="Black", alpha=0.7, cex=1),
                    flow.symbol=list(alpha=0.04, pch=20, cex=0.7),
                    xlab=c("FSC 488/10_A"),ylab=c("SSC 488/10_A"))
  xfilter <- tmixFilter("xfilter", c("FSC 488/10_A","SSC 488/10_A"), K=1, B=50,level=.5) #create filter for 'gating'
  Core.Size.Gate<-Subset(SizeGatedFlowSet,xfilter)
  Plot[[2]]<-  xyplot(`SSC 488/10_A`~`FSC 488/10_A`, Core.Size.Gate,xbin = 128, smooth = F,
                      stat=TRUE,outline=T,
                      par.settings=list(gate=list(fill="black", alpha=0.2)),
                      gate.text=list(col="Black", alpha=0.7, cex=1),
                      flow.symbol=list(alpha=0.04, pch=20, cex=0.7),
                      xlab=c("FSC 488/10_A"),ylab=c("SSC 488/10_A"))
  ####################################################################
  ##Transform data to logaritmic
  #Channels<-c("FSC 488/10_A","SSC 488/10_A","SYTOX Red_A",mCherry_A","GFP_A")
  #log.base<-log10
  #flow.set<-SizeGatedSet
  #cut<-1
  TransformChannelsToLog<-function(channels,log.base,flow.set,cut){
    truncateTrans <- truncateTransform(transformationId="Truncate-transformation", a=1)
    dataTransform <- transform(flow.set,
                               transformList(channels,truncateTrans))
    translist <- transformList(channels, log.base)
    return(transform(dataTransform, translist))
  }
  ####################################################################
  Transfomed_SizeGatedFSD<-TransformChannelsToLog(channels = channels.to.transform,
                                                  log.base = log.base,
                                                  flow.set = Core.Size.Gate, 
                                                  cut = Cut)
  
  rg <- rectangleGate(filterId  ="Rectangle",list.channels.and.limits.DAPI)
  
  DAPIFree_SizeGatedFSC=Subset(Transfomed_SizeGatedFSD, rg)
  Plot[[3]]<-xyplot(`FSC 488/10_A`~`DAPropidium Iodide_A`,Transfomed_SizeGatedFSD,xbin = 128, smooth = F,
                    xlim=c(5,10),ylim=c(log(fsc.488.min,10),log(fsc.488.max,10)),filter=rg,
                    stat=TRUE,outline=T,
                    par.settings=list(gate=list(fill="black", alpha=0.2)),
                    gate.text=list(col="Black", alpha=0.7, cex=1),
                    flow.symbol=list(alpha=0.04, pch=20, cex=0.7),
                    xlab=c("DAPI_A"),ylab=c("FSC 488/10_A"))
  ####Plot filtering gates
  ggsave(path = "../images/CellGate/",
         filename = paste("Size_Cluster_DAPI.",DAPIFree_SizeGatedFSC@description$GUID,
                          Generation,".",FilePos,".png",sep=""),
         width = 13,height =30,units = "cm",scale=1,dpi = 200,
         plot=grid.arrange(Plot[[1]],Plot[[2]],Plot[[3]],nrow=3,ncol=1))
  ###Assigning populations
  scatter.type<-as.vector(timing.data[as.vector(timing.data$V3)%in%
                                        sub(".fcs","",
                                            DAPIFree_SizeGatedFSC@description$GUID),4])
  
  #####Uses the polygons defines before of GFP, mCherry and no fluorophore.
  GFPPositive<-ExtractPointsFromPolygon(polygon.gate =gfp.polygon,
                                        log.transformed.FSC = DAPIFree_SizeGatedFSC,
                                        level=0.99)
  
  mCherryPositive<-ExtractPointsFromPolygon(polygon.gate =mCherry.polygon,log.transformed.FSC = DAPIFree_SizeGatedFSC,
                                            level=0.99)
  #NegativeCells<-ExtractPointsFromPolygon(polygon.gate =negative.polygon,log.transformed.FSC = DAPIFree_SizeGatedFSC,
  #                                       level=0.99)
  Negative.mCherry<-ExtractNegativePointsFromOuterPolygon(polygon.gate =mCherry.polygon,
                                                          log.transformed.FSC = DAPIFree_SizeGatedFSC,
                                                          level=0.99,outer = TRUE)
  NegativeCells<-ExtractNegativePointsFromOuterPolygon(polygon.gate =gfp.polygon,
                                                       log.transformed.FSC = Negative.mCherry,
                                                       level=0.99,outer=TRUE)
  
  ####Those cells that are GFP positive.
  DF.GFPpos=NULL
  if(!is.null(GFPPositive)){
    DF.GFPpos<-as.data.frame(exprs(GFPPositive)[,channels])
    DF.GFPpos$Classification<-"GFP"
  }
  ####Those cells that are mCherry positive.
  DF.mCherrypos=NULL
  if(!is.null(mCherryPositive)){
    DF.mCherrypos<-as.data.frame(exprs(mCherryPositive)[,channels])
    DF.mCherrypos$Classification<-"mCherry"
  }
  ####Negative cells as well.
  Negative.Cells=NULL
  if(!is.null(NegativeCells)){
    Negative.Cells<-as.data.frame(exprs(NegativeCells)[,channels])
    Negative.Cells$Classification<-"NotAssigned"
  }
  ### Concatenate data.
  TotalTableClassified<-do.call("rbind",list(DF.GFPpos,DF.mCherrypos,Negative.Cells))
  ### Plot the entire number of cells separated in gates.
  ClusteredPlot<-ggplot(TotalTableClassified)+
    geom_point(aes_string(x="GFP_A",y="mCherry_A",color="Classification"),
               alpha=.2,size=.01)+coord_cartesian(xlim=c(4,9),ylim=c(4,9)) +
    scale_color_manual(values = c("mCherry"="#BA5798", "GFP"="#71BB52","NotAssigned"="black"))+
    theme_minimal()
  ggsave(path = "../images/CellGate",
         filename = paste("ColorGatesClustered",
                          DAPIFree_SizeGatedFSC@description$GUID,
                          "Generation",Generation,FilePos,".png",sep="_"),
         width = 13,height = 13,units = "cm",scale=1,dpi = 100,
         plot=ClusteredPlot)
  
  ######Adding details to table 
  if(sum(timing.data$V3%in%gsub(paste(path,"",sep="/"),"",
                                gsub(".fcs","",x)))>1){
    path<-gsub("\\/\\w\\d+.fcs","",x)
    FilePos.1<-grep(path,rownames(timing.data[timing.data$V3 %in%gsub(paste(path,"",sep="/"),"",
                                                                      gsub(".fcs","",x)),]))
    SampleInfo<- timing.data[timing.data$V3 %in%gsub(paste(path,"",sep="/"),"",
                                                     gsub(".fcs","",x)),][FilePos.1,]
    
    
  }else{SampleInfo<- timing.data[timing.data$V3 %in%gsub(paste(path,"",sep="/"),"",
                                                         gsub(".fcs","",x)),]}
  
  SampleInfo$V4<-as.vector(SampleInfo$V4)
  ##Kind of segregation
  TotalTableClassified$Segregation<-SampleInfo$V4
  TotalTableClassified$Cross<-SampleInfo$V5
  TotalTableClassified$Replicate<-SampleInfo$V6
  TotalTableClassified$InitFreqClass<-SampleInfo$V7
  TotalTableClassified$Generation<-SampleInfo$V8
  TotalTableClassified$MatingTypeCross<-SampleInfo$V9
  if(nrow(TotalTableClassified)>SampleSize){
    return(TotalTableClassified[sample(x = nrow(TotalTableClassified),size = SampleSize,replace = F),])
  }else(return((TotalTableClassified)))
}
### The output of this file, includes a report of cells gated using forward and side scatter plots.
### I also produces figures where quantitation of fluorophores for each gate is produced.
### Information of each analyzed well is saved a data.frame. 
### This fuctions may use an exceeding ammount of memory. Use with caution.

###########################################################################
## ***************Usage example****************##
#ProcessData(x=fcs.file,gfp.polygon=GFP,mCherry.polygon=mCherry,timing.data=timing.data)
###########################################################################


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################





