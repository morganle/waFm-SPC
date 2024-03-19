#### This is data generation to speed up analysis ####
## Here we generate data for the ARL0 experiment - so Phase I and Phase II parameters are the same. 
## We then feed this data into "Optimising CL parameters.R" to aim to force ARL_0 to equal 371
library(tictoc)

PhaseIseeds<-read.csv(file = "PhaseIseeds_ARL0.csv")
PhaseIIseeds<-read.csv(file = "PhaseIIseeds_ARL0.csv")

WI<-1000 # number of windows of window_len minutes we're considering in Phase I

source("MGt1_withBurstOption_lamMod2.R")
source("WindowMod.R")
#system parameters
window_len<-60
lam_orig<-0.7 

# set matrices to store waFm data
mag_store_all_matrix<-NULL
mag_store_all_p2_matrix<-NULL

n<-1  # lets just use a single subgroup   #as.numeric(read[2])  # batch num
batch_len <- window_len/n
Tend<- WI*window_len # WI contiguous windows # WI*window_len minutes total
warmup = 3*window_len

c=1

## ARL0 calculation
# 1. calculate the control limits using WI windows
# 2. record length of time to detect an out of control signal

## break down into windows and batches
waFm_store<-numeric(WI*n)
window_mean<-numeric(WI)

G<-100       #macro replications


for (kk in 1:G){
  tic()
    print("macro replication")
    print(kk)
    
    magnitude_store<-NULL
    magnitude_store_p2<-NULL
    
    mag_store_all<-NULL
    mag_store_all_p2<-NULL
    
    # in each rep we run a new phase I to calculate the control limits
    # for contiguous windows we run one long WI*window_len window and then split it
    run<-MGt1(PhaseIseeds[kk,3],c,Tend,burst=FALSE,when=warmup,len=0,lam_orig,lam_orig)
    
    for (ii in 1:WI) {  # window 

      ind_min<-min(which(run[[5]] > (ii-1)*window_len ))  
      ind_max<-max(which(run[[5]] < ii*window_len ))  
      NIS_times<-run[[5]][ind_min:ind_max] - (ii-1)*window_len
      NIS_queue<-run[[4]][ind_min:ind_max] 
      
      batch_means<-numeric(n)
      for (jj in 1:n){ # batch
        
        min_i<-min(which(NIS_times > ((jj-1)*batch_len)))
        max_i<-max(which(NIS_times <  (jj*batch_len)))
        
        times<-c((jj-1)*batch_len,NIS_times[min_i:max_i],jj*batch_len) ## adding start and end point
        times<-times-rep((jj-1)*batch_len,length(times))   #translating the batch back to time 0
        queue<-c(NIS_queue[(min_i-1)],NIS_queue[min_i:max_i],NIS_queue[(max_i)])
        
        N_batch<-2^(16)      ## choice of sampling frequency > 20*2500 
        #N<- (N_batch+2)/2   ## appropriate for mirroring
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
        ic_queue = Fqueue[ic_index]
        # reduce ic_index to only points of change in number in system
        ic_index = c(1,which(diff(Fqueue)==1)+1)
        ic_queue = Fqueue[ic_index]
        ic_times = Ftime[ic_index]
        
        # interval-censored batch high-frequency samples
        
        times<- c(0,ic_times,batch_len)
        queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
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
        
        
        # perform the baseline tilt transformation
        Fqueue = WindowMod(w_method = "baseline",w_vec = Fqueue)
        
        Y <-fft(Fqueue)/N_batch 
        ##
        mag <- sqrt(Re(Y)^2+Im(Y)^2)

      }
      magnitudes<-mag[1]  # just looking at the 0th freq  
      magnitude_store<-rbind(magnitude_store,magnitudes)
      
      mags<-mean(mag[2:2501]*seq(1,2500,by=1)*2*pi)    # waFm for this window in Phase I
      mag_store_all<-rbind(mag_store_all,mags)
    }
    mag_store_all_matrix = rbind(mag_store_all_matrix,t(mag_store_all))
    
   ## Phase 2 - ARL_0 experiment
    print("Start of Phase 2")
    WII<-1200     #  max number of windows in Phase II of ARL0 experiment
    Tend <- WII*window_len
    
    run<-MGt1(PhaseIIseeds[kk,3],c,Tend,burst=FALSE,when=warmup,len=0,lam_orig,lam_orig)   # running an in-control experiment at original lambda value
    for (ww in 1:WII){
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
        ic_queue = Fqueue[ic_index]
        # reduce ic_index to only points of change in number in system
        ic_index = c(1,which(diff(Fqueue)==1)+1)
        ic_queue = Fqueue[ic_index]
        ic_times = Ftime[ic_index]
        
        # interval-censored batch high-frequency samples
        
        times<- c(0,ic_times,batch_len)
        queue<- c(NIS_queue[(min_i-1)],ic_queue,NIS_queue[(max_i)])
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
        
        #perform baseline transform
        Fqueue = WindowMod(w_method = "baseline",w_vec= Fqueue)
        
        Y <-fft(Fqueue)/N_batch 
        ##
        mag <- sqrt(Re(Y)^2+Im(Y)^2)
        #waFm<-mean(mag[mag_L:mag_U]*seq(mag_L,mag_U,by=1)*2*pi)    
      }
      magnitudes<-mag[1] 
      magnitude_store_p2<-rbind(magnitude_store_p2,magnitudes)
      
      mags<-mean(mag[2:2501]*seq(1,2500,by=1)*2*pi )   # all magnitudes
      mag_store_all_p2 <- rbind(mag_store_all_p2,mags)
    } 
    mag_store_all_p2_matrix = rbind(mag_store_all_p2_matrix,t(mag_store_all_p2))
 
    toc()
}


write.table(mag_store_all_matrix,"PhaseI_magnitudes_lambda_pt7interval.csv",
            sep=",",row.names = FALSE,col.names = FALSE,append = FALSE)
write.table(mag_store_all_p2_matrix,"PhaseII_magnitudes_lambda_pt7interval.csv",
            sep=",",row.names = FALSE,col.names = FALSE,append = FALSE)

## we end up with two files one including the waFm stat for each window in Phase I for G macro-reps and the other containing the same data but for Phase II

