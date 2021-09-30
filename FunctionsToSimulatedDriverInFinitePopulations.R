#Function that rounds real number such that all of the number have an approximate 
# total value.
# Useful when populations are small.
smart.round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}


##### Function:
##### SimulateDrift_and_MeioticDriveWithCloneSelfing
##### This function is the core of the stochastic simulation. This function is called
##### upon other functions to simulate the evolution of a driver under the set conditions.
##### The output is a table
## The function requires the following arguments
## N.individuals, total individuals selected each generation.
## driving.individuals, Number of individuals that carry a driver.
## meiotic.generations, are the total number of sexual cycles.
## selfing.portion, Inbreeding coefficient F.
## Iterations, times the simulation is run.
## Replacing, TRUE=Infinite gametes pool or finite=FALSE.
## associated.cot, Value of cost linked to the driver.
## cost.dominance, Value of dominance coefficient of associated cost to th elinked driver.
## driver.birth, The driver can be born durin meiosis or after haploids have been selected
# Usage:
#   SimulateDrift_and_MeioticDriveWithCloneSelfing(
#     associated.cost = 0 ,
#     iterations=Iterations,
#     initial.driving.individuals = InitialDrivingIndividuals,
#     generation.mitotic.growth = MitoticGenerations,
#     meiotic.generations = MeioticGenerations,
#     replacing = TRUE,pop.sizes = PopulationSizes,
#     clonal.selfing = Clonal.selfing.F,
#     driver.birth = MomentOfDriverBirth,cost.dominance = DominanceOfLinkedCost)
SimulateDrift_and_MeioticDriveWithCloneSelfing<-function(N.individuals,driving.individuals,
                                                         generation.mitotic.growth,
                                                         meiotic.generations,
                                                         selfing.portion,
                                                         replacing,associated.cost,
                                                         driver.birth,
                                                         cost.dominance){
  Cost<-associated.cost
  N <- N.individuals #N is number of individuals in population
  mi.g <- generation.mitotic.growth # Mitotic generations in between sexual 
  me.g <- meiotic.generations #number of meiotic generations 
  selfing.portion<-selfing.portion
  Generations<-data.frame(ND=N-driving.individuals,D=driving.individuals)
  for(g.me in 1:me.g){
    N.D.A <- rep(0, Generations[g.me,"ND"]) # number of copies of allele 0
    D.A <- rep(1, Generations[g.me,"D"]) 
    
    pop<-c(N.D.A,D.A)
    #### Mitotic growth assuming equal chances of cells growing
    pop.grown<-rep(x = pop,2^mi.g)
    
    #### Sampling for meiosis
    if(length(pop.grown)!=0){
      
      if(driver.birth=="Meiosis"){
        SampledPop<-sample(pop.grown,N,replace = replacing)
      }
      if(driver.birth=="SelectedHaploids"){
        if(g.me==1){
          SampledPop<-pop
        }else SampledPop<-sample(pop.grown,N,replace = replacing)
      }
      if(driver.birth!="Meiosis" &driver.birth!="SelectedHaploids" ){
        print("No driver birth")
        flush.console()
        return(NULL)
      }
      
      
      #### Selecting fraction that inbreeds
      Selfing<-sample(SampledPop,length(SampledPop)*selfing.portion,
                      replace = FALSE)
      ##### Defining the total values of driving and non driving alleles 
      ##### in selfing portion.
      Selfing.v<-c(ND=sum(Selfing==0),D=sum(Selfing==1))
      Selfing.v.Progeny<-Selfing.v
      ##Selfing always doubles for non drivers
      Selfing.v.Progeny["ND"]<-Selfing.v["ND"]*2
      ##Selfing always doubles but loses the portion by the cost. It is rounded
      Selfing.v.Progeny["D"]<-round(Selfing.v["D"]*2*(1-Cost))
      
      ##### Defining the total values of driving and non driving alleles 
      ##### in non-selfing portion.
      Non.Selfing.v<-c(ND=unname(sum(SampledPop==0)-Selfing.v["ND"]),
                       D=unname(sum(SampledPop==1)-Selfing.v["D"]))
      
      if(sum(Non.Selfing.v)!=0){
        #### Non selfing can form homozygotes or heterozygotes.
        Non.Selfing.freq<-(Non.Selfing.v)/sum(Non.Selfing.v)
        Non.Selfing.freq.HW<-c(Hom.ND=unname(Non.Selfing.freq["ND"]^2),
                               Het=unname(2*Non.Selfing.freq["ND"]*Non.Selfing.freq["D"])*(1-Cost*cost.dominance),
                               Hom.D=unname(Non.Selfing.freq["D"]^2)*(1-Cost))
        if(sum(Non.Selfing.freq.HW)>0){ 
          Non.Selfing.freq.Norm<-Non.Selfing.freq.HW/sum(Non.Selfing.freq.HW)}
        
        if(sum(Non.Selfing.freq.HW)==0){ 
          Non.Selfing.freq.Norm[1:3]<-0}
        ##### Since they portion that mates randomly is in realnumbers, I use
        #### round smart to select to apporximate round numbers to sum a total of the
        #### n individuals.
        Non.Selfing.v.Progeny<-smart.round(sum(Non.Selfing.v)*Non.Selfing.freq.Norm)
        #### Zygote numbers
        Non.Selfing.v.Progeny.Zygotes<-c(Hom.ND=unname(Non.Selfing.v.Progeny["Hom.ND"]*2),
                                         Het=unname(Non.Selfing.v.Progeny["Het"]),
                                         Hom.D=unname(Non.Selfing.v.Progeny["Hom.D"]*2))
        #### Progeny getting final frequencies of driving and non driving.
        Non.Selfing.v.Progeny<-c(ND=unname(Non.Selfing.v.Progeny.Zygotes["Hom.ND"]),
                                 D=unname(sum(Non.Selfing.v.Progeny.Zygotes[c("Het","Hom.D")])))
      }
      
      if(sum(Non.Selfing.v)==0){
        Non.Selfing.v.Progeny<-c(ND=0,D=0)
      }
      Generations<-rbind(Generations,Selfing.v.Progeny+Non.Selfing.v.Progeny)
      
    }else{Generations<-rbind(Generations,data.frame(ND=0,D=0))}
  }
  Generations$MeioticGeneration<-0:me.g
  Generations$ClonalSelfing<-selfing.portion
  Generations$AssociatedCost<-associated.cost
  return(Generations)
}


#####Function:
##### MeioticDriveWithDrift_Iterations
##### This function uses the that simulates the evolution of a driver under multiple inbreeding coefficients
##### It calls the previously declared function and  parameters are parsed:
#####
#####  SimulateDrift_and_MeioticDriveWithCloneSelfing
#####
##### The added argument for this fuction parese to SimulateDrift_and_MeioticDriveWithCloneSelfing is:
#####  Iterations
#####
##### Usage:
# Iterations<-1000
#   MeioticDriveWithDrift_Iterations(
#     N.individuals=PopulationSizes
#     associated.cost = 0 ,
#     iterations=Iterations,
#     initial.driving.individuals = InitialDrivingIndividuals,
#     generation.mitotic.growth = MitoticGenerations,
#     meiotic.generations = MeioticGenerations,
#     replacing = TRUE,pop.sizes = PopulationSizes,
#     clonal.selfing = Clonal.selfing.F,
#     driver.birth = MomentOfDriverBirth,cost.dominance = DominanceOfLinkedCost)
MeioticDriveWithDrift_Iterations<-function(N.individuals,
                                           driving.individuals,
                                           generation.mitotic.growth,
                                           meiotic.generations,
                                           selfing.portion,Iterations,replacing,
                                           associated.cost,cost.dominance,
                                           driver.birth){
  
  All.Iterations<-data.frame()
  for(I in 1:Iterations){
    if((I%%500)==0){print(I)}
    DrivingSimulation<-
      SimulateDrift_and_MeioticDriveWithCloneSelfing(N.individuals = N.individuals,
                                                     driving.individuals = driving.individuals,
                                                     generation.mitotic.growth = generation.mitotic.growth,
                                                     meiotic.generations = meiotic.generations,
                                                     selfing.portion = selfing.portion,replacing = replacing,
                                                     associated.cost=associated.cost,
                                                     driver.birth=driver.birth,
                                                     cost.dominance=cost.dominance)
    
    DrivingSimulation$Iteration<-I
    if(tail(DrivingSimulation,1)$ND==0 & tail(DrivingSimulation,1)$D!=0){
      DrivingSimulation$FixedDriver<-"Fixed"
    }
    if(tail(DrivingSimulation,1)$D==0 & tail(DrivingSimulation,1)$ND!=0){
      DrivingSimulation$FixedDriver<-"Lost"
    }
    if(tail(DrivingSimulation,1)$ND!=0 & tail(DrivingSimulation,1)$D!=0){
      DrivingSimulation$FixedDriver<-"Drifting"
    }
    if(tail(DrivingSimulation,1)$ND==0 & tail(DrivingSimulation,1)$D==0){
      DrivingSimulation$FixedDriver<-"Extinct"
    }
    All.Iterations<-rbind(All.Iterations,DrivingSimulation)
  }
  return(All.Iterations)
}


#####


#####Function:
##### MeioticDriveWithDrift_Iterations_Selfing_PopulationSize.Cost_Dominance
##### This function calls the function MeioticDriveWithDrift_Iterations
##### It asssigns the correct values that are parsed down to the core function:
#####
##### SimulateDrift_and_MeioticDriveWithCloneSelfing
#####
##### Function uses same arguments as: StochasticSimulationOfaDriver
# Usage:
# StochasticSimulationOfaDriver<-
#   MeioticDriveWithDrift_Iterations_Selfing_PopulationSize.Cost_Dominance(
#     N.individuals=PopulationSizes
#     associated.cost = 0 ,
#     iterations=Iterations,
#     initial.driving.individuals = InitialDrivingIndividuals,
#     generation.mitotic.growth = MitoticGenerations,
#     meiotic.generations = MeioticGenerations,
#     replacing = TRUE,pop.sizes = PopulationSizes,
#     clonal.selfing = Clonal.selfing.F,
#     driver.birth = MomentOfDriverBirth,cost.dominance = DominanceOfLinkedCost)
MeioticDriveWithDrift_Iterations_Selfing_PopulationSize.Cost_Dominance<-function(associated.cost,
                                                                                 iterations,
                                                                                 initial.driving.individuals,
                                                                                 generation.mitotic.growth,
                                                                                 meiotic.generations,replacing,
                                                                                 pop.sizes,clonal.selfing,
                                                                                 driver.birth,
                                                                                 cost.dominance){
  Iterations=iterations
  #######################################
  FixProbs.Selfing.PopSize<-data.frame()
  for(pop.size in pop.sizes){
    print(paste("***********pop.size =",pop.size,"*******"))
    FixProbs.Selfing<-data.frame()
    for(s in clonal.selfing){
      print(paste("***********selfing =",s,"*******"))
      FixProbs.Cost<-data.frame()
      for(c in associated.cost){
        ValuesSelfing<- MeioticDriveWithDrift_Iterations(N.individuals = pop.size,
                                                         driving.individuals = initial.driving.individuals,
                                                         generation.mitotic.growth = generation.mitotic.growth,
                                                         meiotic.generations = meiotic.generations,
                                                         selfing.portion = s,
                                                         Iterations=iterations,replacing = replacing,
                                                         associated.cost=c,
                                                         driver.birth=driver.birth,
                                                         cost.dominance=cost.dominance)
        
        data.frame(Fixed=sum(ValuesSelfing$FixedDriver=="Fixed")/(meiotic.generations+1),
                   Lost=sum(ValuesSelfing$FixedDriver=="Lost")/(meiotic.generations+1),
                   Drifting=sum(ValuesSelfing$FixedDriver=="Drifting")/(meiotic.generations+1),
                   Extinct=sum(ValuesSelfing$FixedDriver=="Extinct")/(meiotic.generations+1)
        )->
          Values.Table
        Values.Table$Selfing<-s
        Values.Table$AssociatedCost<-c
        FixProbs.Cost<-rbind(FixProbs.Cost,Values.Table)
      }
      FixProbs.Selfing<-rbind(FixProbs.Selfing,FixProbs.Cost)
    }
    FixProbs.Selfing$Pop.Size<-pop.size
    FixProbs.Selfing.PopSize<-rbind(FixProbs.Selfing.PopSize,FixProbs.Selfing)
  }
  return(FixProbs.Selfing.PopSize)
}




##library(ggplot2)
##library(dplyr)
##library(Rmisc)
##library(scales)
##library(scico)
### ### ### ### ### ### 
### Example of input.### 
### ### ### ### ### ### 
### ### ### ### ### ### 
# Iterations<-100
# InitialDrivingIndividuals<-1
# MitoticGenerations<-0
# MeioticGenerations<-100
# PopulationSizes<-10^(1:5)
# Clonal.selfing.F<-seq(0,1,0.05)
# MomentOfDriverBirth<-"Meiosis"
# DominanceOfLinkedCost<-0.5
# CostLinked<-0

# StochasticSimulationOfaDriver<-
#   MeioticDriveWithDrift_Iterations_Selfing_PopulationSize.Cost_Dominance(
#     N.individuals=PopulationSizes
#     associated.cost = CostLinked ,
#     iterations=Iterations,
#     initial.driving.individuals = InitialDrivingIndividuals,
#     generation.mitotic.growth = MitoticGenerations,
#     meiotic.generations = MeioticGenerations,
#     replacing = TRUE,pop.sizes = PopulationSizes,
#     clonal.selfing = Clonal.selfing.F,
#     driver.birth = MomentOfDriverBirth,cost.dominance = DominanceOfLinkedCost)



### ### ### ### ### ### 
### Example of output### 
### ### ### ### ### ### 
### ### ### ### ### ### 
#StochasticSimulationOfaDriver

#     Fixed Lost Drifting Extinct Selfing AssociatedCost  Pop.Size
# 1   3     6    1        0       0.00    0               10
# 2   2     8    0        0       0.05    0               10
# 3   2     8    0        0       0.10    0               10
# 4   2     8    0        0       0.15    0               10
# 5   6     4    0        0       0.20    0               10
# .....................


