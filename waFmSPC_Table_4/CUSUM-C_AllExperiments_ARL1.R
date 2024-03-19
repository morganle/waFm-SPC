#### This is the CUSUM-C Chen and Zhou Experiment ARL0 ####

## to be edited but we'll do the CL param experiment first

# now all done on one computer so id aspects deleted
#id = as.numeric(Sys.getenv("SGE_TASK_ID"))
#print(id)
expt_data = read.csv("Experiment_Settings_CUSUM-C.csv",header=TRUE)




seed<-123456
set.seed(seed)

PhaseIseeds<-read.csv(file = "PhaseIseeds_ARL1.csv")
PhaseIIseeds<-read.csv(file = "PhaseIIseeds_ARL1.csv")

source("MGt1_withBurstOption_lamMod2.R")

#system parameters
window_len<-60
lam_0<-0.7
lam_orig<-lam_0
c=1
G=100 #macro reps
max_len<-1200
Tend<-max_len*window_len # trajectories of length 1200 windows
warmup = 3*window_len    # initial portion deleted, optional change in lambda

CUSUMC_results = as.data.frame(matrix(0,nrow(expt_data),(3+G)))

for(i_expt in 1:nrow(expt_data)){
  expt = expt_data[i_expt,]

lam_1<-as.numeric(expt[1])            # what is the target shift

lam_true<-as.numeric(expt[2])           # what is the true shift
CL_parameter<-as.numeric(expt[3])


store_ARL1<-c()
store_ARL1_specific<-c()

# since each macroreplication has different number of events, set
#   oversize array to hold each (h8, h9, h5, h6) and the actual 
#   number of values in each (n_kk)

h_results = matrix(0,G,as.integer((lam_orig+.5)*Tend)) # larger than max departures by Tend
                                                       # since changing lam, need +.5

T_k = matrix(0,G,as.integer((lam_orig+.5)*Tend)) #larger than max departures by Tend
l_T_k_m1 = rep(0,G) # number of departures - 1 for each macroreplication

for(kk in 1:G){ 
  print(kk)
  ww<-0
  h_last<-0
  h<-c()
  dept_times<-c()
  detect_time<-0
  
  #NOTE: len=0 not used, 5 window warmup
  #run<-MGt1(PhaseIIseeds[kk,1],c,Tend,burst=TRUE,when=warmup,
  #          len=0,lam_orig,lam_true)  # running at the true shift parameter
  run<-MGt1(PhaseIIseeds[kk,2],c,Tend,burst=TRUE,when=warmup,
            len=0,lam_orig,lam_true)  # running at the true shift parameter
  
  NIS_times<-run[[5]]
  NIS_queue<-run[[4]]

  dept_ind<-c()
  for(p in 1:(length(NIS_queue)-1)){
      if( NIS_queue[p+1] == (NIS_queue[p]-1) ){dept_ind<-c(dept_ind,p+1)}
  }
  l_NIS_times = length(NIS_times[dept_ind])
  T_k[kk,1:l_NIS_times]<-NIS_times[dept_ind]   #times of departures
  l_T_k_m1[kk] = l_NIS_times - 1
  Q_k<-NIS_queue[(dept_ind)]
      
  for(o in 1:l_T_k_m1[kk]){
      cusum_stat <- h_last + (Q_k[o+1] - Q_k[o] + 1)*log(lam_1/lam_0) - (lam_1 - lam_0)*(T_k[kk,o+1] - T_k[kk,o])    # collecting CUSUM-C at target shift parameter lam_1
      h_results[kk,o] = max(0,cusum_stat)
      h_last <- max(0,cusum_stat)
  }
      
  dept_times<-T_k[kk,1:l_NIS_times]
  if(length(which(h_results[kk,1:l_T_k_m1[kk]]>CL_parameter))>0) {
  detect_ind<-min(which(h_results[kk,1:l_T_k_m1[kk]]>CL_parameter))+1 # this is the departure the alert occurred at
  detect_time<-dept_times[detect_ind]/60  # in units of windows
  # we'll get a window of detection this way
  }else{
    detect_time=max_len
  }

  store_ARL1_specific[kk]<-detect_time
  #print(kk)
  #print(detect_time)
}  
CUSUMC_results[i_expt,1]=lam_0     # value from which departure happens
CUSUMC_results[i_expt,2]=lam_true  # true departure
CUSUMC_results[i_expt,3]=lam_1     # CUSUM-C planned departure
CUSUMC_results[i_expt,4:(G+3)]=store_ARL1_specific[1:G]
# all experiments done
}

write.table(CUSUMC_results,paste("BurstyArrivals/CUSUM-C_MG1_ARL1",".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=FALSE)

## we get 16 output rows looking at the 4 different shifts lam_true = 0.5,0.6,0.8 and 0.9 
## and we consider for each shift lam_1 of any value lam_1=0.5, 0.6, 0.8, 0.9
## 16 can be seen as nrow(expt_data) from the experiment file we read in at the top

