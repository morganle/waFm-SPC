## In this script we calculate ARL_1 (out-of-control run length)
## we use the optimized CL parameter within the CUSUM method
## we consider all possible shifts of the arrival rate.

library("qcc")
library(tictoc)

#id = as.numeric(Sys.getenv("SGE_TASK_ID"))
#print(id)
data = read.csv("Experiment_Settings_ARL1_CUSUM.csv",header=TRUE)
#read = data[id,]
read = data[1,]

print(read)
seed<-123456
set.seed(seed)

PhaseIseeds<-read.csv("PhaseIseeds_ARL1.csv")
PhaseIIseeds<-read.csv("PhaseIIseeds_ARL1.csv")

W<-as.numeric(read[1])
mag_L<-as.numeric(read[2])
mag_U<-as.numeric(read[3])
CL_para<-as.numeric(read[4])

## CUSUM Chart
source("MGt1_withBurstOption_lamMod2.R")
source("WindowMod.R")

window_len<-60


n<-1  # lets just use a single subgroup    # batch num
batch_len <- window_len/n
Tend<-W*window_len 
warmup = 3*window_len
c=1

## break down into windows and batches
waFm_store<-numeric(W*n)
window_mean<-numeric(W)


## PHASE I a single really long run to try to remove error in the control limits W=10,000
print("Start of Phase I")
lam_orig<-0.7
run<-MGt1(PhaseIseeds[1,1],c,Tend,burst=FALSE,when=warmup,len=0,lam_orig,lam_orig)

for (ii in 1:W){  # window 

  ind_min<-min(which(run[[5]] > (ii-1)*window_len )) 
  ind_max<-max(which(run[[5]] < ii*window_len ))  
  NIS_times<-run[[5]][ind_min:ind_max] - (ii-1)*window_len
  NIS_queue<-run[[4]][ind_min:ind_max] 
  
  batch_means<-numeric(n)
  for (jj in 1:n){ # batch
    
    min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
    max_i<-max(which(NIS_times <  (jj*batch_len)))
    
    times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len)
    times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
    queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
    
    N_batch<-2^(16)
    N<- (N_batch+2)/2  # appropriate for mirroring 
    min_int<-batch_len/N_batch
    # special case: interval data - sampling interval 30 seconds (.5 minute)
    samp_int = .5
    N_int = floor(batch_len/samp_int)
    frac_samp = N_int/N_batch
    
    c_brkpt <- numeric(length(times))
    
    j = 1
    for (i in 1:length(times))
    {
      c_brkpt[j] <- floor(times[i]/min_int)  #Calculate max order of; pts w/ each int.
      j = j+1
    }
    
    Fqueue <- numeric(c_brkpt[length(c_brkpt)])
    Ftime <- numeric(c_brkpt[length(c_brkpt)])
    
    for (i in 1:(length(times)-1))
    {
      for (j in (c_brkpt[i]+1):c_brkpt[i+1])
      {
        Fqueue[j] <- queue[i]
        Ftime[j] <- min_int * j
      }
    }
    
    # interval sampling here
    ic_index_part = 1:N_int
    ic_index = floor((1/(2*frac_samp)) + ((ic_index_part - 1) / frac_samp))
    # reduce ic_index to only points of change in number in system
    ic_index = ic_index[(which(diff(Fqueue[ic_index])!=0))]
    ic_queue = Fqueue[ic_index]
    ic_times = Ftime[ic_index]
    
    # interval-censored batch high-frequency samples
    
    ic_times<- c(0,ic_times,batch_len)
    queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
    c_brkpt <- numeric(length(ic_times))
    
    j = 1
    for (i in 1:length(ic_times))
    {
      c_brkpt[j] <- floor(ic_times[i]/min_int)  #Calculate max order of; pts w/ each int.
      j = j+1
    }
    
    Fqueue <- numeric(c_brkpt[length(c_brkpt)])
    Ftime <- numeric(c_brkpt[length(c_brkpt)])
    
    for (i in 1:(length(ic_times)-1))
    {
      for (j in (c_brkpt[i]+1):c_brkpt[i+1])
      {
        Fqueue[j] <- queue[i]
        Ftime[j] <- min_int * j
      }
    }
    
    ## perform mirroring
    Fqueue = WindowMod(w_method ="baseline",w_vec=Fqueue)
    Y<-fft(Fqueue)/N_batch
    
    ##
    mag <- sqrt(Re(Y)^2+Im(Y)^2)
    if((mag_L == 1) & (mag_U==1)){
      waFm<-mag[1]
    }else{
      waFm<-mean(mag[mag_L:mag_U]*seq((mag_L-1),(mag_U-1),by=1)*2*pi)     
    }

    waFm_store[(ii-1)*n+jj]<-waFm  # store for the waFm stats before we transform them
  }
}


for (ii in 1:W){  # window
  window_ii_waFm<-waFm_store[((ii-1)*n+1):((ii-1)*n+n)]
  window_mean[ii]<-mean(window_ii_waFm)
}

PhaseIdata<-window_mean

#===============================================================
# Phase II runs
print("Start Phase II")
pre_runs<-0 
G<-100
max_len<-1200
Tend = max_len*window_len

# ARL Down pt 1 -----------------------------------------------------------
tic()
# ## Step 2 - ARL_1 experiment
lam_orig = 0.7
lam_burst<-0.6
store_ARL<-numeric(G)

for (kk in 1:G){
  print(kk)
  window_mean = NULL
  window_mean<-numeric(max_len)
  
  run<-MGt1(PhaseIIseeds[kk,1],c,Tend,burst=TRUE,when=warmup,len=0,lam_orig,lam_burst) 
  
  for(ww in 1:max_len){

    ind_min<-min(which(run[[5]] > (ww-1)*window_len ))
    ind_max<-max(which(run[[5]] < ww*window_len ))
    NIS_times<-run[[5]][ind_min:ind_max] - (ww-1)*window_len
    NIS_queue<-run[[4]][ind_min:ind_max] 
    
    batch_means<-numeric(n)
    for (jj in 1:n){ # batch
      
      min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
      max_i<-max(which(NIS_times <  (jj*batch_len)))
      
      times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len)
      times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
      queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
      
      N_batch<-2^(16)
      N<- (N_batch+2)/2   # appropriate for mirroring
      min_int<-batch_len/N_batch
      # special case: interval data - sampling interval 30 seconds (.5 minute)
      samp_int = .5
      N_int = floor(batch_len/samp_int)
      frac_samp = N_int/N_batch
      
      c_brkpt <- numeric(length(times))
      
      j = 1
      for (i in 1:length(times))
      {
        c_brkpt[j] <- floor(times[i]/min_int)  #Calculate max order of; pts w/ each int.
        j = j+1
      }
      
      Fqueue <- numeric(c_brkpt[length(c_brkpt)])
      Ftime <- numeric(c_brkpt[length(c_brkpt)])
      
      for (i in 1:(length(times)-1))
      {
        for (j in (c_brkpt[i]+1):c_brkpt[i+1])
        {
          Fqueue[j] <- queue[i]
          Ftime[j] <- min_int * j
        }
      }
      
      # interval sampling here
      ic_index_part = 1:N_int
      ic_index = floor((1/(2*frac_samp)) + ((ic_index_part - 1) / frac_samp))
      # reduce ic_index to only points of change in number in system
      ic_index = ic_index[(which(diff(Fqueue[ic_index])!=0))]
      ic_queue = Fqueue[ic_index]
      ic_times = Ftime[ic_index]
      
      # interval-censored batch high-frequency samples
      
      ic_times<- c(0,ic_times,batch_len)
      queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
      c_brkpt <- numeric(length(ic_times))
      
      j = 1
      for (i in 1:length(ic_times))
      {
        c_brkpt[j] <- floor(ic_times[i]/min_int)  #Calculate max order of; pts w/ each int.
        j = j+1
      }
      
      Fqueue <- numeric(c_brkpt[length(c_brkpt)])
      Ftime <- numeric(c_brkpt[length(c_brkpt)])
      
      for (i in 1:(length(ic_times)-1))
      {
        for (j in (c_brkpt[i]+1):c_brkpt[i+1])
        {
          Fqueue[j] <- queue[i]
          Ftime[j] <- min_int * j
        }
      }
      
      
      ## perform mirroring
      Fqueue = WindowMod(w_method ="baseline",w_vec=Fqueue)
      Y<-fft(Fqueue)/N_batch
      
      ##
      mag <- sqrt(Re(Y)^2+Im(Y)^2)
      if((mag_L == 1) & (mag_U==1)){
        waFm<-mag[1]
      }else{
        waFm<-mean(mag[mag_L:mag_U]*seq((mag_L-1),(mag_U-1),by=1)*2*pi)     
      }

      batch_means[jj]<-waFm  # store for the waFm stats before we transform them
    }
    ##############
    window_mean[ww] <- mean(batch_means)
  }  
  PhaseIIdata<-window_mean
  
  q<-cusum(PhaseIdata,newdata = PhaseIIdata,decision.interval = CL_para, plot = FALSE)   # no we perform waFm CUSUM
  violations<-sort(c(q$violations$lower,q$violations$upper))
  
  if(length(which(violations>W))>0) {
    RL <- violations[min(which(violations>W))] - W
  }else{
    RL <- max_len
  }
  store_ARL[kk] <- RL
}

 print(store_ARL)
 mean(store_ARL)
 toc()
 OO<-c(W,mag_L,mag_U)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerDown_ARL1_pt1_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 OO<-c(store_ARL)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerDown_ARL1_pt1_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 

# ARL Down pt 2 -----------------------------------------------------------
 
 # ## Step 2 - ARL_1 experiment
 lam_orig<-0.7
 lam_burst<-0.5
 store_ARL<-numeric(G)

 for (kk in 1:G){
   print(kk)
   ww<-0
   window_mean<-c()

   run<-MGt1(PhaseIIseeds[kk,2],c,Tend,burst=TRUE,when=warmup,len=0,lam_orig,lam_burst)
   
   for(ww in 1:max_len){
     
     ind_min<-min(which(run[[5]] > (ww-1)*window_len ))  # 1PM
     ind_max<-max(which(run[[5]] < ww*window_len))  # 2PM
     NIS_times<-run[[5]][ind_min:ind_max] - (ww-1)*window_len
     NIS_queue<-run[[4]][ind_min:ind_max] 
     
     batch_means<-numeric(n)
     for (jj in 1:n){ # batch
       
       min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
       max_i<-max(which(NIS_times <  (jj*batch_len)))
       
       times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len)
       times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
       queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
       
       N_batch<-2^(16)
       N<- (N_batch+2)/2   # useful for mirroring 
       min_int<-batch_len/N_batch
       # special case: interval data - sampling interval 30 seconds (.5 minute)
       samp_int = .5
       N_int = floor(batch_len/samp_int)
       frac_samp = N_int/N_batch
       
       c_brkpt <- numeric(length(times))
       
       j = 1
       for (i in 1:length(times))
       {
         c_brkpt[j] <- floor(times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       
       # interval sampling here
       ic_index_part = 1:N_int
       ic_index = floor((1/(2*frac_samp)) + ((ic_index_part - 1) / frac_samp))
       # reduce ic_index to only points of change in number in system
       ic_index = ic_index[(which(diff(Fqueue[ic_index])!=0))]
       ic_queue = Fqueue[ic_index]
       ic_times = Ftime[ic_index]
       
       # interval-censored batch high-frequency samples
       
       ic_times<- c(0,ic_times,batch_len)
       queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
       c_brkpt <- numeric(length(ic_times))
       
       j = 1
       for (i in 1:length(ic_times))
       {
         c_brkpt[j] <- floor(ic_times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(ic_times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       
       ## perform mirroring
       Fqueue = WindowMod(w_method ="baseline",w_vec=Fqueue)
       Y<-fft(Fqueue)/N_batch
       
       ##
       mag <- sqrt(Re(Y)^2+Im(Y)^2)
       if((mag_L == 1) & (mag_U==1)){
         waFm<-mag[1]
       }else{
         waFm<-mean(mag[mag_L:mag_U]*seq((mag_L-1),(mag_U-1),by=1)*2*pi)     
       }
       batch_means[jj]<-waFm  # store for the waFm stats before we transform them
     }
     ##############
     window_mean[ww] <-mean(batch_means)
   }  
   PhaseIIdata<-window_mean
   q<-cusum(PhaseIdata,newdata = PhaseIIdata,decision.interval = CL_para, plot = FALSE)
   violations<-sort(c(q$violations$lower,q$violations$upper))
   
   if(length(which(violations>W))>0) {
     RL <- violations[min(which(violations>W))] - W
   }else{
     RL <- max_len
   }

   store_ARL[kk] <- RL
 }
 
 print(store_ARL)
 mean(store_ARL) 

 OO<-c(W,mag_L,mag_U)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerDown_ARL1_pt2_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 OO<-c(store_ARL)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerDown_ARL1_pt2_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 

# ARL Up pt 1 -----------------------------------------------------------
 
 # ## Step 2 - ARL_1 experiment
 lam_orig<-0.7
 lam_burst<-0.8
 store_ARL<-numeric(G)

 for (kk in 1:G){
   print(kk)
   window_mean<-c()
   
   run<-MGt1(PhaseIIseeds[kk,3],c,Tend,burst=TRUE,when=warmup,len=0,lam_orig,lam_burst) 
   
   for(ww in 1:max_len){
     
     ind_min<-min(which(run[[5]] > (ww-1)*window_len ))  
     ind_max<-max(which(run[[5]] < ww*window_len ))  
     NIS_times<-run[[5]][ind_min:ind_max] - (ww-1)*window_len
     NIS_queue<-run[[4]][ind_min:ind_max] 
     
     batch_means<-numeric(n)
     for (jj in 1:n){ # batch
       
       min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
       max_i<-max(which(NIS_times <  (jj*batch_len)))
       
       times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len)
       times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
       queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
       
       N_batch<-2^(16)
       N<- (N_batch+2)/2   # useful for mirroring
       min_int<-batch_len/N_batch
       # special case: interval data - sampling interval 30 seconds (.5 minute)
       samp_int = .5
       N_int = floor(batch_len/samp_int)
       frac_samp = N_int/N_batch
       
       c_brkpt <- numeric(length(times))
       
       j = 1
       for (i in 1:length(times))
       {
         c_brkpt[j] <- floor(times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       # interval sampling here
       ic_index_part = 1:N_int
       ic_index = floor((1/(2*frac_samp)) + ((ic_index_part - 1) / frac_samp))
       # reduce ic_index to only points of change in number in system
       ic_index = ic_index[(which(diff(Fqueue[ic_index])!=0))]
       ic_queue = Fqueue[ic_index]
       ic_times = Ftime[ic_index]
       
       # interval-censored batch high-frequency samples
       
       ic_times<- c(0,ic_times,batch_len)
       queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
       c_brkpt <- numeric(length(ic_times))
       
       j = 1
       for (i in 1:length(ic_times))
       {
         c_brkpt[j] <- floor(ic_times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(ic_times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       ## perform mirroring
       Fqueue = WindowMod(w_method ="baseline",w_vec=Fqueue)
       Y<-fft(Fqueue)/N_batch
       
       ##
       mag <- sqrt(Re(Y)^2+Im(Y)^2)
       if((mag_L == 1) & (mag_U==1)){
         waFm<-mag[1]
       }else{
         waFm<-mean(mag[mag_L:mag_U]*seq((mag_L-1),(mag_U-1),by=1)*2*pi)     
       }
       batch_means[jj]<-waFm  # store for the waFm stats before we transform them
     }
     ##############
     window_mean[ww] <-mean(batch_means)
     
   }  
   PhaseIIdata<-window_mean
   q<-cusum(PhaseIdata,newdata = PhaseIIdata,decision.interval = CL_para, plot = FALSE)
   violations<-sort(c(q$violations$lower,q$violations$upper))
   
   if(length(which(violations>W))>0) {
     RL <- violations[min(which(violations>W))] - W
   }else{
     RL <- max_len
   }
   
   store_ARL[kk] <- RL
 }
 
 print(store_ARL)
 mean(store_ARL) 
 OO<-c(W,mag_L,mag_U)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerUp_ARL1_pt1_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 OO<-c(store_ARL)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerUp_ARL1_pt1_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 
# ARL Up pt 2 -----------------------------------------------------------
 
 # ## Step 2 - ARL_1 experiment
 lam_orig<-0.7
 lam_burst<-0.9
 store_ARL<-numeric(G)

 for (kk in 1:G){
   print(kk)
   window_mean<-c()

   run<-MGt1(PhaseIIseeds[kk,4],c,Tend,burst=TRUE,when=warmup,len=0,lam_orig,lam_burst) 
   
   for(ww in 1:max_len){
     
     ind_min<-min(which(run[[5]] > (ww-1)*window_len ))  
     ind_max<-max(which(run[[5]] < ww*window_len ))  
     NIS_times<-run[[5]][ind_min:ind_max] - (ww-1)*window_len
     NIS_queue<-run[[4]][ind_min:ind_max] 
     
     batch_means<-numeric(n)
     for (jj in 1:n){ # batch
       
       min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
       max_i<-max(which(NIS_times <  (jj*batch_len)))
       
       times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len)
       times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
       queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
       
       N_batch<-2^(16)
       N<- (N_batch+2)/2   # useful when mirroring
       min_int<-batch_len/N_batch
       # special case: interval data - sampling interval 30 seconds (.5 minute)
       samp_int = .5
       N_int = floor(batch_len/samp_int)
       frac_samp = N_int/N_batch
       
       c_brkpt <- numeric(length(times))
       
       j = 1
       for (i in 1:length(times))
       {
         c_brkpt[j] <- floor(times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       # interval sampling here
       ic_index_part = 1:N_int
       ic_index = floor((1/(2*frac_samp)) + ((ic_index_part - 1) / frac_samp))
       # reduce ic_index to only points of change in number in system
       ic_index = ic_index[(which(diff(Fqueue[ic_index])!=0))]
       ic_queue = Fqueue[ic_index]
       ic_times = Ftime[ic_index]
       
       # interval-censored batch high-frequency samples
       
       ic_times<- c(0,ic_times,batch_len)
       queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
       c_brkpt <- numeric(length(ic_times))
       
       j = 1
       for (i in 1:length(ic_times))
       {
         c_brkpt[j] <- floor(ic_times[i]/min_int)  #Calculate max order of; pts w/ each int.
         j = j+1
       }
       
       Fqueue <- numeric(c_brkpt[length(c_brkpt)])
       Ftime <- numeric(c_brkpt[length(c_brkpt)])
       
       for (i in 1:(length(ic_times)-1))
       {
         for (j in (c_brkpt[i]+1):c_brkpt[i+1])
         {
           Fqueue[j] <- queue[i]
           Ftime[j] <- min_int * j
         }
       }
       
       ## perform mirroring
       Fqueue = WindowMod(w_method ="baseline",w_vec=Fqueue)
       Y<-fft(Fqueue)/N_batch
       
       ##
       mag <- sqrt(Re(Y)^2+Im(Y)^2)
       if((mag_L == 1) & (mag_U==1)){
         waFm<-mag[1]
       }else{
         waFm<-mean(mag[mag_L:mag_U]*seq((mag_L-1),(mag_U-1),by=1)*2*pi)     
       }
       batch_means[jj]<-waFm  # store for the waFm stats before we transform them
     }
     ##############
     window_mean[ww] <-mean(batch_means)

   }  
   PhaseIIdata<-window_mean
   q<-cusum(PhaseIdata,newdata = PhaseIIdata,decision.interval = CL_para, plot = FALSE)
   violations<-sort(c(q$violations$lower,q$violations$upper))
   
   if(length(which(violations>W))>0) {
     RL <- violations[min(which(violations>W))] - W
   }else{
     RL <- max_len
   }
   
   store_ARL[kk] <- RL
 }
 
 print(store_ARL)
 mean(store_ARL) 
 OO<-c(W,mag_L,mag_U)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerUp_ARL1_pt2_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
 OO<-c(store_ARL)
 write.table(t(OO),paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_SerUp_ARL1_pt2_interval.csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

 