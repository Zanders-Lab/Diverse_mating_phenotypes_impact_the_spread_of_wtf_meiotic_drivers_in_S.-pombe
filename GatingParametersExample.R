
##### Example of parameters used to gate round cells
#
#Size gate
#X
fsc.488.min<<-6e8 
fsc.488.max<<-1.2e9
#Y
ssc.488.min<-2e8
ssc.488.max<-7e8

#### Example of gates to optain living cells
list.channels.and.limits.DAPI<<-list("DAPropidium Iodide_A"=c(5.5, 7),"FSC 488/10_A"=c(log(fsc.488.min,10),log(fsc.488.max,10)))
channels.to.transform<<-c("FSC 488/10_A","SSC 488/10_A","DAPropidium Iodide_A","GFP_A","mCherry_A")
log.base<<-log10
Cut<<-1
##Channels to classify
channels<<-c("GFP_A","mCherry_A")
#Number of cells to plot
SampleSize<<-100000