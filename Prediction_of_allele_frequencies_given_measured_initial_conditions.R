#### This script contains the functions used to calculate the cost of carrying a gfp performs a sumary of each experiment perfromed, 
#### deleterious allele cost, and the predicted values of alleles. 
####The prediction considers given  genotype, inbreeding coefficient and linked deleterious alleles.
### 

### Dependencies:
# library(reshape2)
# library(dplyr)
# library(stats4)
# library(Rmisc)
# library(ggplot2)
# library(ggpubr)
# library(tidyverse)
# library(rstatix)
# library(ggpubr)



######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
### StatisticsPerGeneration. This function gathers the information from the gated and classified events produced by Gating_and_quantitation_of_flowcytometry_events.R
### This functions uses as an input the list with all the measured samples.
### The output of this function is a table with the flurophore frequencies, and the descripcition of each sample.
StatisticsPerGeneration<-function(generation){
  ## First we separate by the type of cross. 
  ## The segreations are either the crosses we generation use:
  ## Mendelian as mC-G-M and mC*-G*-M or biased with drive:
  ## mC*-G-B and mC-G*-B.
  StatisticsPerGeneration
  Segregations<-unique(sub("_.+","",sub("\\w\\d+_","",names(generation))))
  Names<-sub("_.+","",sub("\\w\\d+_","",names(generation)))
  
  GenerationData<-sapply(Segregations,function(x){
    ## For each genotype present in this generation.
    Genotype<-x
    seg.type.set<- generation[Names%in%x]
    seg.type.sum<-lapply(
      seg.type.set,
      function(x){
        if(sum(grep("mCherry",unique(x$Segregation)),
               grep("GFP",unique(x$Segregation)),
               grep("NC",unique(x$Segregation)))>0)
        {
          return(NULL)
        }
        ## The segregation names determine what fluorophore is the one that drives.
        if(length(grep("mC\\*",unique(x$Segregation)))==1){driver.pop<-"mCherry";Segregation="Biased"}
        if(length(grep("G\\*",unique(x$Segregation)))==1){driver.pop<-"GFP";Segregation="Biased"}
        if(length(grep("G\\*",unique(x$Segregation)))==1 &
           length(grep("mC\\*",unique(x$Segregation)))==1){driver.pop<-"mCherry";Segregation="Mendelian"}
        if(length(grep("mC-",unique(x$Segregation)))==1 &
           length(grep("G-",unique(x$Segregation)))==1){driver.pop<-"mCherry";Segregation="Mendelian"}
        if(sum(c("GFP","mCherry")%in%names(table(x$Classification)))<2){
          NotDetected<-(c("GFP","mCherry")[!c("GFP","mCherry")%in%
                                             names(table(x$Classification))])
          tot.pop<-table(x$Classification)[c("GFP","mCherry")]
          names(tot.pop)[is.na(names(tot.pop))]<-NotDetected
          tot.pop[NotDetected]<-0
        }else{
          tot.pop<-table(x$Classification)[c("GFP","mCherry")]
        }
        #### Here I gather all the information necesary.
        Populations<-tot.pop/sum(tot.pop) ## fluorophores
        Cross<-unique(x$Cross)            ## The experiment or cross performed.
        Replicate<-unique(x$Replicate)    ## Technical replicate during mitotic growth. 
        InitialFreq<-Populations[driver.pop] ## This colums is deprecate, but remains in the script.
        Generation<-unique(x$Generation)  ## The generation at the experiment.
        InitFreqClass<-unique(x$InitFreqClass) ## This is the intended initial frequency at which crosses started.
        MatingTypeCross<-unique(x$MatingTypeCross) ## Importan while distinguishing different mating type of different strains.
        Data<-data.frame(t(as.data.frame(as.numeric(as.vector(as.matrix(Populations))))),
                         Cross,Replicate,InitialFreq,Segregation,Genotype,
                         Generation,
                         InitFreqClass,MatingTypeCross)
        colnames(Data)[1:2]<-c("GFP","mCherry")
        rownames(Data)<-Genotype
        return(Data)
      })
    return(list(do.call("rbind",seg.type.sum)))
  })
  GenerationData<-GenerationData[names(unlist(lapply(X = GenerationData,nrow)))]
  GenerationData<-do.call("rbind",GenerationData)
  print("Generation Done")
  flush.console()
  GenerationData
}
###########################################################################
## ***************Usage example****************##
#StatisticsPerGeneration(generation = Table_from_ProcessData_function)
###########################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
####Calculate.FlurescentvsNonFlurescentEvents. This function checks for the frequency of fluorescent
#### and non fluorescent cells.
Calculate.FlurescentvsNonFlurescentEvents<-function(generation){
  list.table<-lapply(generation,function(X){
    table(X$Classification)
  })
  sample.info<-lapply(generation,function(X){
    sample<-c(Experiment=unique(X$Cross),Replicate=unique(X$Replicate))
    sample
  })
  list.clean.names<-lapply(list.table,function(X){
    matched<-match(names(X),c("GFP" ,"mCherry" ,"NotAssigned"))
    new.table<-c(0,0,0)
    if(length(matched)>0){
      new.table[matched]<-X
    }
    names(new.table)<-c("GFP" ,"mCherry" ,"NotAssigned")
    new.table
  })
  
  abs.freq<-do.call("rbind",list.clean.names)
  samp.info<-do.call("rbind",sample.info)
  rel.freq<-data.frame((t(apply(X = abs.freq,1,function(X){
    (X/sum(X))
  }))))
  freqs<-cbind(abs.freq,rel.freq)
  colnames(freqs)<-paste(colnames(freqs),c(rep(".Abs",3),rep(".Rel",3)),sep="")
  freqs<-cbind(freqs,samp.info)
  freqs
}
###########################################################################
## ***************Usage example****************##
#Calculate.FlurescentvsNonFlurescentEvents(generation = Table_from_ProcessData_function)
###########################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
### calculateGenerationsDriveAndInbreeding: This function calculates drive based on JCrow and we adapted ti to include inbreeding.
### This funtion requires the Wright parameter for drive, also associated with the JCrow notation. These values can be
## Change to include the cost of bearing a meiotic drivers.
calculateGenerationsDriveAndInbreeding<-function(w1,w2,w3,p,q,k,Generations,f){
  W.inb=NULL
  Total=NULL
  for(i in 1:(Generations)){
    W.inb      =c(W.inb    , p[i]^2*w1 + 2*(p[i]*q[i])*w2 + q[i]^2*w3 + f*p[i]*q[i]*(w1+w3-2*w2))       # W calculation
    p      =c(p    ,(p[i]^2*w1 + f*p[i]*q[i]*w1 + 2*p[i]*q[i]*w2*k*(1-f))       /W.inb[i])   # p calculation
    q      =c(q    ,(q[i]^2*w3 + f*p[i]*q[i]*w3 + 2*p[i]*q[i]*w2*(1-k)*(1-f))   /W.inb[i])      # q calculation
  }
  W.inb      =c(W.inb    , p[Generations]^2*w1 +
                  2*(p[Generations]*q[Generations])*w2 +
                  q[Generations]^2*w3 + 
                  f*p[Generations]*q[Generations]*(w1+w3-2*w2)) #Adding the last W to fulfill table
  MeioticDrive<-data.frame(DriverInbreeding=p,FitnessInbreeding=W.inb)
  return(MeioticDrive)
}
###########################################################################
## ***************Usage example****************##
# calculateGenerationsDriveAndInbreeding(w1 = W1 ,w2 = W2,w3 = W3,p = P,
#                                        q = 1-P,k = K,f = inbreeding.coeff,
#                                        Generations =generations)
###########################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################



######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
##GetRelativeFitnessFromMLE: This function estimates the expected relative fitness of GFP given fixed parameters using LGBS
## The algorithm calculates the expected cost of being a GFP homozygous diploid.
Get.cost.carrying.a.gfp.allele<-function(frequencies.per.generation,genotype,
                                         inbreeding.coeff,mating.type,generations){
  
  costs.gfp<-sapply(sort(unique(frequencies.per.generation$InitFreqClass)),function(X){
    dplyr::filter(frequencies.per.generation,InitFreqClass==X)->Init.Freq.Class.Selected
    frequencies.per.generation.selected<-as_tibble(summarySE(dplyr::filter(Init.Freq.Class.Selected,MatingTypeCross==mating.type,
                                                                           Genotype==genotype),
                                                             groupvars = c("Generation","Cross","Segregation","Genotype",
                                                                           "InitFreqClass","MatingTypeCross"),measurevar = "mCherry"))
    
    gfp.allele.cost<-sapply(unique(frequencies.per.generation.selected$Cross),function(X){
      ### Here we are using the observed values measured. I used mCherry as the objective
      ### frequency to measure. We do this since GFP is the one that showed different rel.fit.
      y=dplyr::filter(frequencies.per.generation.selected,Cross==X, Generation<=generations)%>%
        dplyr::select(mCherry)%>%unlist()
      ### This is the initial frequency then names the p initial frequency.
      p<-unlist(dplyr::filter(frequencies.per.generation.selected,Generation==0,Cross==X)%>%
                  dplyr::select(mCherry))
      LL <- function(c, mu, sigma,f) {
        # Find residuals
        R = y - calculateGenerationsDriveAndInbreeding(w1 = 1-c,w2 = (1-c*0.5),w3 =1-2*c,
                                                       p = p,q = 1-p,k = .5,Generations = generations,f=f)$Driver
        # Calculate the likelihood for the residuals (with mu and sigma as parameters
        R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
        # Sum the log likelihoods for all of the data points
        -sum(R)
      }
      fit<-stats4::mle(LL, start = list(c = 0), fixed = list(mu = 0,sigma=1,f=inbreeding.coeff),
                       nobs = length(y),method = "L-BFGS-B", lower = c(0),
                       upper = c(Inf)) 
      relatve.fitness.MLE<-stats4::coef(fit)["c"]
      return(relatve.fitness.MLE)
    })
    return(list(data.frame(GFPCostInDiploid=gfp.allele.cost,Genotype=genotype,
                           MatingType=mating.type,I.C=inbreeding.coeff,Init.Freq.Class=X,
                           row.names = c(unique(frequencies.per.generation.selected$Cross)))
    ))
  })
  do.call("rbind",costs.gfp)
}
###########################################################################
## ***************Usage example****************##
####### Using Maximum Likelihood Estimation to get relative fitness
#init.freq.class<-0.5
#mating.type<-"h90"
#inbreeding.coeff<-0.5
#Genotype<-"mC-G-M"
#Generations.to.fit<-6
#Get.cost.carrying.a.gfp.allele(frequencies.per.generation = frequencies.per.generation.E1,
#                               genotype = Genotype,inbreeding.coeff=inbreeding.coeff,
#                               mating.type = mating.type,generations = Generations.to.fit)
###########################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################


######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
##CalculatePredictedFrequencies: This function calculates the predicted frequencies given an inbreeding coefficient and
##the cost of carrying a deleterious allele, in this case gfp.
CalculatePredictedFrequencies<-function(frequencies.per.generation,
                                        generations,genotype,inbreeding.coeff,
                                        mating.type,gfp.cost.diploid,measure.allele,driver.allele){
  
  ##  This cross uses the GFP no wtf homozygote relative fitness
  RandMat.Sim<-SimulateFrequencies(frequencies.per.generation = frequencies.per.generation,generations = generations,
                                   genotype = genotype,inbreeding.coeff=0,mating.type = mating.type,
                                   gfp.cost.diploid=gfp.cost.diploid,
                                   measure.allele =measure.allele ,driver.allele = driver.allele)
  RandMat.Sim$Class<-"Exp.RandomMating"
  inbreeding.coeff<-inbreeding.coeff
  Inbreeding.Sim<-SimulateFrequencies(frequencies.per.generation = frequencies.per.generation,generations = generations,
                                      genotype = genotype,inbreeding.coeff=inbreeding.coeff,mating.type = mating.type,
                                      gfp.cost.diploid=gfp.cost.diploid,
                                      measure.allele =measure.allele,driver.allele = driver.allele)
  Inbreeding.Sim$Class<-"Exp.Inbreeding"
  dplyr::filter(frequencies.per.generation,MatingTypeCross==mating.type,Genotype==genotype)%>%
    select_(measure.allele,"Generation","Cross","InitFreqClass","MatingTypeCross") %>%
    mutate(InbreedingCoefficient=0) %>%
    mutate(Class="Measured") %>%
    mutate(Cross=as.character(Cross),
           InitFreqClass=as.character(InitFreqClass),
           MatingTypeCross=as.character(MatingTypeCross))->
    Measured
  SimulationsAndMesured<-bind_rows(as_tibble(RandMat.Sim),as_tibble(Inbreeding.Sim),Measured)
  mutate(SimulationsAndMesured,Genotype=as.character(genotype))
}
###########################################################################
## ***************Usage example****************##
####### Using Maximum Likelihood Estimation to get relative fitness
#Generations<-10
#Genotype<-"mC-G-M"
#mating.type<-"h90"
#inbreeding.coeff<-0.5
#measure.allele<-"mCherry"
#driver.allele = "mCherry"
# CalculatePredictedFrequencies(frequencies.per.generation = Table_produced_in_StatisticsPerGeneration,
#                               generations = Generations,genotype = Genotype,inbreeding.coeff=inbreeding.coeff,
#                               mating.type = mating.type,gfp.cost.diploid =AvergareGFPCost,
#                               measure.allele = measure.allele,driver.allele = driver.allele)
###########################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################







