#### This is the CUSUM-C Chen and Zhou Experiment ARL0 ####

## to be edited but we'll do the CL param experiment first
## run for 1 to 8 ## COMPLETED 1, 2, 3, 4, 5, 6

data = read.csv("Experiment_Settings_CUSUM-P.csv",header=TRUE)
#id = as.numeric(Sys.getenv("SGE_TASK_ID"))
for(id in 1:16){
  print(id)
  read = data[id,]
  
  seed<-123456
  set.seed(seed)
  
  PhaseIseeds<-read.csv(file = "PhaseIseeds_ARL1.csv")
  PhaseIIseeds<-read.csv(file = "PhaseIIseeds_ARL1.csv")
  
  lam_1<-as.numeric(read[1])            # what is the target shift
  lam_0<-0.7
  lam_true<-as.numeric(read[2])           # what is the true shift
  CL_parameter<-as.numeric(read[3])
  
  source("MGt1_withBurstOption_lam.R")
  
  
  #system parameters
  window_len<-60
  lam_orig<-0.7
  c=1
  
  G=100 #macro reps
      
  max_len<-1200
  Tend<-max_len*60 #trajectories of length 1200 windows
  
  store_ARL1<-c()
  store_ARL1_specific<-c()
  
  for(kk in 1:G){  
    ww<-0
    h_last<-0
    h<-c()
    dept_times<-c()
    detect_time<-0
    
    run<-MGt1(PhaseIIseeds[kk,3],c,Tend,burst=TRUE,when=0,len=0,lam_true,lam_true)  # running at the true shift parameter
    NIS_times<-run[[5]]
    NIS_queue<-run[[4]]
  
    dept_ind<-c()
    for(p in 1:(length(NIS_queue)-1)){
        if( NIS_queue[p+1] == (NIS_queue[p]-1) ){dept_ind<-c(dept_ind,p+1)}
    }
    T_k<-NIS_times[dept_ind]
    Q_k<-NIS_queue[(dept_ind)]
    Z=0 
    
    #we're interested in shifts from 0.7 up and down by 0.1 and 0.2 so we collect 4 sets of CUSUM-C stats
    # now for the cusum stat we're tracking for each lam_1 parameter
    for(o in 1:(length(T_k)-1)){
      if(Q_k[o]==0){Z=1}else{Z=0}
      cusum_stat <- h_last + log((lam_1*(1+lam_0))/(lam_0*(1+lam_1)))*(Q_k[o+1] - Q_k[o] + 1 - Z) - log((1+lam_1)/(1+lam_0)) #CUSUM-P
      #CUSUM-C  #h_last + (Q_k[o+1] - Q_k[o] + 1)*log(0.8/lam_0) - (0.8 - lam_0)*(T_k[o+1] - T_k[o])
      h<-c(h, max(0,cusum_stat))
      h_last <- max(0,cusum_stat)
    }    
    #for(o in 1:length(T_k)-1){
    #    cusum_stat <- h_last + (Q_k[o+1] - Q_k[o] + 1)*log(lam_1/lam_0) - (lam_1 - lam_0)*(T_k[o+1] - T_k[o])    # collecting CUSUM-C at target shift parameter lam_1
    #    h<-c(h, max(0, cusum_stat ))
    #    h_last <- max(0,cusum_stat)
    #}
        
    dept_times<-T_k
    detect_ind<-min(which(h>CL_parameter))
    detect_time<-dept_times[detect_ind]/60  # in units of windows
  
    if(is.na(detect_time)==TRUE){detect_time=max_len}
               # we'll get a window of detection this way
    ARL_1_specific<-detect_time
    store_ARL1_specific[kk]<-ARL_1_specific
  }  
  
  write.table(t(store_ARL1_specific),paste("BurstyArrivals/CUSUM-P_bursty_ARL1_lam1_",lam_1,"_lam_true",lam_true,".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
  
}