## In this file we generate trajectories of length W*60 i.e. W windows of 1 hour.
## We match seeds to the waFm CUSUM experiment so we're considering the same data.
## We identify departure events and use them within the CUSUM-C statistic.

## We also optimise the control limit parameter used within the CUSUM-C chart so that the in-control
## ARL is approximately 371 (we used this value in the initial MtG1 experiment)
prob_i<-function(t,i,lam){
  return((exp(-lam*t)*(lam*t)^i)/factorial(i))
}
PhaseIseeds<-read.csv(file = "PhaseIseeds_ARL0.csv")
PhaseIIseeds<-read.csv(file = "PhaseIIseeds_ARL0.csv")

W<-10000  #windows

source("MGt1_withBurstOption_lam.R")  # generates trajectories

#system parameters
window_len<-60
lam_orig<-0.7             

n<-1  # lets just use a single subgroup  

Tend<-W*60 # as we want a trajectory W windows long
c=1

## ARL0 calculation
# 1. calculate the control limits using W windows
# 2. record length of time to detect an out of control signal

G<-100      #macro replications
max_len<-1200 # max number of windows in ARL0 experiment 
lam_0<-0.7


for (kk in 1:G){
    print("macro replication")
    print(kk)
    
    #Phase I
    
    h<-c() # this is our store for the CUSUM stat
    hh<-c()
    hhh<-c()
    hhhh<-c()

    h_last<-0
    hh_last<-0
    hhh_last<-0
    hhhh_last<-0

    # we run the simulation and collect event times and NIS
    # note that we run this simulation at 0.7 (lam_orig)
    run<-MGt1(PhaseIIseeds[kk,1],c,Tend,burst=FALSE,when=0,len=0,lam_orig,lam_orig)
    NIS_times<-run[[5]]
    NIS_queue<-run[[4]]

    #identify departure events
    dept_ind<-c()
      for(p in 1:(length(NIS_queue)-1)){
        if( NIS_queue[p+1] == (NIS_queue[p]-1) ){dept_ind<-c(dept_ind,p+1)}
      }
    T_k<-NIS_times[dept_ind]   #times of departures
    Q_k<-NIS_queue[dept_ind]   #NIS at departures
    Z=0 

      #we're interested in shifts from 0.7 up and down by 0.1 and 0.2 so we collect 4 sets of CUSUM-C stats
      # now for the cusum stat we're tracking for each lam_1 parameter
      for(o in 1:(length(T_k)-1)){
        
        if(Q_k[o]==0){Z=1}else{Z=0}
        qu<- Q_k[o+1] - Q_k[o] + 1 - Z
        lam_0_int <- integrate(prob_i,lower=0.5,upper=1.5,i=qu,lam=lam_0)$value
        
        lam_1_int <- integrate(prob_i,lower=0.5,upper=1.5,i=qu,lam=0.8)$value
        h<-c(h, max(0,h_last + log(lam_1_int) - log(lam_0_int)))
        h_last <- max(0,h_last + log(lam_1_int) - log(lam_0_int))
        
        lam_1_int <- integrate(prob_i,lower=0.5,upper=1.5,i=qu,lam=0.9)$value
        hh<-c(hh, max(0, hh_last + log(lam_1_int) - log(lam_0_int)))
        hh_last <- max(0,hh_last + log(lam_1_int) - log(lam_0_int))
        
        lam_1_int <- integrate(prob_i,lower=0.5,upper=1.5,i=qu,lam=0.6)$value
        hhh<-c(hhh, max(0, hhh_last + log(lam_1_int) - log(lam_0_int)))
        hhh_last <- max(0,hhh_last + log(lam_1_int) - log(lam_0_int))
        
        lam_1_int <- integrate(prob_i,lower=0.5,upper=1.5,i=qu,lam=0.5)$value
        hhhh<-c(hhhh, max(0, hhhh_last + log(lam_1_int) - log(lam_0_int)))
        hhhh_last <- max(0,hhhh_last + log(lam_1_int) - log(lam_0_int))
      }
    
    # all stats have now been collected in this macro rep so let's print to file

  write.table(t(h),"CUSUM-P_pt7_pt8_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(hh),"CUSUM-P_pt7_pt9_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(hhh),"CUSUM-P_pt7_pt6_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(hhhh),"CUSUM-P_pt7_pt5_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)

  # storing the departure times so we can find the detection time 
  write.table(t(T_k),"CUSUM-P_pt7_Tk.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
}



# The following function is for optimizing the control limit parameter to some value (here 371) when the system is in-control 
# i.e. we find CL_para that gives us ARL_0 (in control average run length) of 371

find_param<-function(CL_para,data,Tk_vec){
  G<-100
  ooc<-c()
  max_len<-1200
  print(CL_para)
  
  for (i in 1:G){ # in the G macro replications
    phaseIIdata<-as.numeric(data[i,])  # cusum-c stats given lam_1 of 0.5,0.6,0.8 or 0.9
    T_k<-as.numeric(Tk_vec[i,]) # departure times rep i

    h_ind <- as.numeric(min(which( phaseIIdata > CL_para ))) # which index surpasses CL_pata
    if(h_ind != Inf){
      dept_alert<-T_k[h_ind]     # this is the departure the alert occured at
      ooc[i] <- dept_alert/60    # this tells us the run length in units of hour long windows (so we can compare ARL_0 across methods and ensure they're similar)
    }else{ooc[i]=Inf}       # if no alert put Inf 

  }
  ooc[which(ooc==Inf)]=max_len    # we now treat any Infs as the maximum of the windows so 1200 i.e. max_len
  #print(ooc)
  print(mean(ooc)-371)       # uniroot is looking for a 0 so we consider ARL_0 - target ARL_0
  return(mean(ooc) - 371)
}

indFile="CUSUM-P_pt7_Tk.csv"
Tk_vec<-read.csv(file = indFile,header=FALSE)  #only needs reading once

filename="CUSUM-P_pt7_pt5_W1200.csv"
data<-read.csv(file = filename,header = FALSE)
uniroot(find_param,data=data,Tk_vec=Tk_vec, lower=2,upper=6)

filename="CUSUM-P_pt7_pt6_W1200.csv"
data<-read.csv(file = filename,header = FALSE)
uniroot(find_param,data=data,Tk_vec=Tk_vec, lower=2,upper=6)

filename="CUSUM-P_pt7_pt8_W1200.csv"
data<-read.csv(file = filename,header = FALSE)
uniroot(find_param,data=data,Tk_vec=Tk_vec, lower=2,upper=6)

filename="CUSUM-P_pt7_pt9_W1200.csv"
data<-read.csv(file = filename,header = FALSE)
uniroot(find_param,data=data,Tk_vec=Tk_vec, lower=2,upper=6)