## In this file we generate trajectories of length W*60 i.e. W windows of 1 hour.
## We match seeds to the waFm CUSUM experiment so we're considering the same data.
## We identify departure events and use them within the CUSUM-C statistic.

## We also optimise the control limit parameter used within the CUSUM-C chart so that the in-control
## ARL is approximately 371 (we used this value in the initial MtG1 experiment)

PhaseIseeds<-read.csv(file = "PhaseIseeds_ARL0.csv")
PhaseIIseeds<-read.csv(file = "PhaseIIseeds_ARL0.csv")



source("MGt1_withBurstOption_lamMod2.R")  # generates trajectories

#system parameters
window_len<-60
lam_orig<-0.7             
n<-1  # lets just use a single subgroup  
c=1 # single server

## ARL0 calculation
# 1. calculate the control limits using W windows
# 2. record length of time to detect an out of control signal



G<-100      #macro replications
W<-1200  #1200 windows for phase I
Tend<-W*window_len # as we want a trajectory W windows long 
warmup = 3*window_len # delete first three hours from M/G/ 1 simulation

# since each macroreplication has different number of events, set
#   oversize array to hold each (h8, h9, h5, h6) and the actual 
#   number of values in each (n_kk)

h8 = matrix(0,G,as.integer((lam_orig+.2)*Tend)) #larger than max departures by Tend
h9 = matrix(0,G,as.integer((lam_orig+.2)*Tend)) #larger than max departures by Tend
h5 = matrix(0,G,as.integer((lam_orig+.2)*Tend)) #larger than max departures by Tend
h6 = matrix(0,G,as.integer((lam_orig+.2)*Tend)) #larger than max departures by Tend
T_k = matrix(0,G,as.integer((lam_orig+.2)*Tend)) #larger than max departures by Tend
l_T_k_m1 = rep(0,G) # number of departures - 1 for each macroreplication

for (kk in 1:G){
    print("macro replication")
    print(kk)
    
    #Phase I
    
    #h<-c() # this is our store for the CUSUM stat
    #hh<-c()
    #hhh<-c()
    #hhhh<-c()

    h_last<-0
    hh_last<-0
    hhh_last<-0
    hhhh_last<-0

    # we run the simulation and collect event times and NIS
    # note that we run this simulation at 0.7 (lam_orig)
    run<-MGt1(PhaseIseeds[kk,1],c,Tend,burst=FALSE,when=warmup,len=0,lam_orig,lam_orig)
    NIS_times<-run[[5]]
    NIS_queue<-run[[4]]

    #identify departure events
    dept_ind<-c()
      for(p in 1:(length(NIS_queue)-1)){
        if( NIS_queue[p+1] == (NIS_queue[p]-1) ){dept_ind<-c(dept_ind,p+1)}
      }
    l_NIS_times = length(NIS_times[dept_ind])
    T_k[kk,1:l_NIS_times]<-NIS_times[dept_ind]   #times of departures
    l_T_k_m1[kk] = l_NIS_times - 1
    Q_k<-NIS_queue[dept_ind]   #NIS at departures
      

      #we're interested in shifts from 0.7 up and down by 0.1 and 0.2 so we collect 4 sets of CUSUM-C stats
      # now for the cusum stat we're tracking for each lam_1 parameter
    # NOTE: to shorten formula, use lam_0 rather than lam_orig
    lam_0 = lam_orig
      for(o in 1:l_T_k_m1[kk]){
        cusum_stat <- h_last + (Q_k[o+1] - Q_k[o] + 1)*log(0.8/lam_0) - (0.8 - lam_0)*(T_k[kk,o+1] - T_k[kk,o])
        #print(cusum_stat)
        h8[kk,o] = max(0,cusum_stat)
        h_last <- max(0,cusum_stat)
      }
      
      for(o in 1:l_T_k_m1[kk]){
        cusum_stat <- hh_last + (Q_k[o+1] - Q_k[o] + 1)*log(0.9/lam_0) - (0.9 - lam_0)*(T_k[kk,o+1] - T_k[kk,o])
        #print(cusum_stat)
        h9[kk,o] = max(0,cusum_stat)
        hh_last <- max(0,cusum_stat)
      }

      for(o in 1:l_T_k_m1[kk]){
        cusum_stat <- hhh_last + (Q_k[o+1] - Q_k[o] + 1)*log(0.6/lam_0) - (0.6 - lam_0)*(T_k[kk,o+1] - T_k[kk,o])
        h6[kk,o] = max(0,cusum_stat)
        hhh_last <- max(0,cusum_stat)
      }

      for(o in 1:l_T_k_m1[kk]){
        cusum_stat <- hhhh_last + (Q_k[o+1] - Q_k[o] + 1)*log(0.5/lam_0) - (0.5 - lam_0)*(T_k[kk,o+1] - T_k[kk,o])
        h5[kk,o] = max(0,cusum_stat)
        hhhh_last <- max(0,cusum_stat)
      }

    # all stats have now been collected in this macro rep so let's print to file
    # DO NOT DO THIS - h, hh, hhh, hhhh different lengths in each macro rep
    # PLUS writing and reading such large files takes very long time

  #write.table(t(h),"CUSUM-C_pt7_pt8_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  #write.table(t(hh),"CUSUM-C_pt7_pt9_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  #write.table(t(hhh),"CUSUM-C_pt7_pt6_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  #write.table(t(hhhh),"CUSUM-C_pt7_pt5_W1200.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
  
  # storing the departure times so we can find the detection time 
    # DO NOT DO THIS - T_k different lengths in each macro rep
    # PLUS writing and reading such large files takes very long time
  #write.table(t(T_k),"CUSUM-C_pt7_Tk.csv",sep=",",row.names = FALSE,col.names = FALSE,append = TRUE)
}



# The following function is for optimizing the control limit parameter to some value (here 371) when the system is in-control 
# i.e. we find CL_para that gives us ARL_0 (in control average run length) of 371

find_param<-function(CL_para,data,Tk_vec,len_Tk){
  G<-100
  ooc<-c()
  max_len<-1200
  print(CL_para)
  
  for (i in 1:G){ # in the G macro replications
    phaseIdata<-as.numeric(data[i,1:len_Tk[i]])  # cusum-c stats given lam_1 of 0.5,0.6,0.8 or 0.9
    T_i<-as.numeric(Tk_vec[i,1:(len_Tk[i]+1)]) # departure times rep i

    if(length(which( phaseIdata > CL_para ))>0) {
      h_ind <- as.numeric(min(which( phaseIdata > CL_para ))) # which index surpasses CL_para
      dept_alert<-T_i[h_ind+1]     # this is the departure the alert occurred at
      ooc[i] <- dept_alert/60    # this tells us the run length in units of hour long windows
                                 # (so we can compare ARL_0 across methods and ensure they're similar)
    }else{
      ooc[i]=max_len
      }
}
  #print(ooc)
  print(mean(ooc)-371)       # uniroot is looking for a 0 so we consider ARL_0 - target ARL_0
  return(mean(ooc) - 371)
}

# no need to read data - held in h8, h9, h5, h6 and T_k
#filename="CUSUM-C_pt7_pt8_W1200.csv"
#indFile="CUSUM-C_pt7_Tk.csv"
Tk_vec = T_k

#data<-read.csv(file = filename,header = FALSE)
#Tk_vec<-read.csv(file = indFile,header=FALSE)
uniroot(find_param,data=h8,Tk_vec=T_k,len_Tk=l_T_k_m1, lower=1,upper=7)

#filename="CUSUM-C_pt7_pt9_W1200.csv"
#indFile="CUSUM-C_pt7_Tk.csv"
#data<-read.csv(file = filename,header = FALSE)
#Tk_vec<-read.csv(file = indFile,header=FALSE)
uniroot(find_param,data=h9,Tk_vec=T_k,len_Tk=l_T_k_m1, lower=1,upper=7)

#filename="CUSUM-C_pt7_pt6_W1200.csv"
#indFile="CUSUM-C_pt7_Tk.csv"
#data<-read.csv(file = filename,header = FALSE)
#Tk_vec<-read.csv(file = indFile,header=FALSE)
uniroot(find_param,data=h6,Tk_vec=T_k,len_Tk=l_T_k_m1, lower=1,upper=7)

#filename="CUSUM-C_pt7_pt5_W1200.csv"
#indFile="CUSUM-C_pt7_Tk.csv"
#data<-read.csv(file = filename,header = FALSE)
#Tk_vec<-read.csv(file = indFile,header=FALSE)
uniroot(find_param,data=h5,Tk_vec=T_k,len_Tk=l_T_k_m1, lower=1,upper=8)

# At the end of running this file we should have 4 control limit parameters to work with that give ARL_0 close to 371. One for each shift.



