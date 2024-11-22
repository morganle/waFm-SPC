## analysis of experiments from the SPC bursty experiment

# CUSUM -------------------------------------------------------------------
## ARL0
G<-100 # number macro reps
par(mfrow=c(3,2))
W<-c(1000)
magL<-c(1)
magU<-c(2501)
All<-c()
max_len<-1200    # if we didn't detect a change then we cap this at the threshold of the length of Phase II(line 18)
for(i in 1:length(W)){
  for(j in 1:length(magU)){
    ##e.g. CUSUM_mirrWaFm_bursty_W3000_magL1_magU1_ARL0_lambda_pt7
    data<-read.csv(paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_W",W[i],"_magL",magL[j],"_magU",magU[j],"_ARL0_lambda_pt7_adjustedCLs.csv",sep=""),sep=",",header = FALSE)
    RL<-as.numeric(data[1:100,1])
    inds<-which(is.na(RL)==TRUE)
    if(length(inds)>0){ RL[inds]<-max_len }
    
    mean_ARL0<-mean(RL)
    se_ARL0<-sqrt(var(RL)/length(RL))
    
    OO<-t(c(W[i],magL[j],magU[j],mean_ARL0,se_ARL0))
    All<-rbind(All,OO)
  }
}
print(All)

write.table(All,"BurstyArrivals/CUSUM_ARL0_experiment_pt7.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)

## ARL1
##repeat for pt1 and pt2

par(mfrow=c(3,2))
ud<-c("SerDown","SerUp")
magL<-c(1,2,2,2,2,4)
magU<-c(1,2,3,4,2501,2501)
preruns<-0
max_runlen<-1000

W<-10000
All<-c()
for(i in 1:length(ud)){
  for(j in 1:length(magU)){
    #e.g. CUSUM_mirrWaFm_bursty_SerUp_ARL1_W10000_magL2_magU2_pt2
    data<-read.csv(paste("BurstyArrivals/CUSUM_mirrWaFm_bursty_",ud[i],"_ARL1_W",W,"_magL",magL[j],"_magU",magU[j],"_pt2_adjustedCLs.csv",sep=""),sep=",",header = FALSE)
    
    RL1<-as.numeric(data[2,])
    inds<-which(is.na(RL1==TRUE))
    #RL1[inds]<-max_runlen
    if(length(inds)>0){ RL1 <- RL1[-inds] }
    #hist(as.numeric(data[2,]),main=paste0("Magnitudes - ",magL[j],":",magU[j]),xlim=c(0,50))
    #abline(v=mean(RL1),col="red")
    
    RL1<-sort(RL1)
    trimmed_mean<-mean(RL1[which(RL1>=quantile(RL1,0.9))])
    #qqnorm(RL1,main=paste(ud[i],"_ARL1_magL",magL[j],"_magU",magU[j],"_pt2",sep=""))
    
    OO<-t(c(ud[i],W,magL[j],magU[j],mean(RL1),median(RL1),trimmed_mean))
    All<-rbind(All,OO)
  }
  write.table(All,paste("CUSUM_ARL1_pt2_",ud[i],".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
}
print(All)


### Experiment analysis CUSUM-C

data<-read.csv("BurstyArrivals/CUSUM-C_bursty_ARL1_lam1_0.9_lam_true0.9.csv",header = FALSE)
ARL1<-mean(as.numeric(data[1,]))
MRL1<-median(as.numeric(data[1,]))
TMRL1<-mean(sort(as.numeric(data[1,]),decreasing=TRUE)[1:10])

print(ARL1)
print(MRL1)

print(TMRL1)

### Experiment analysis CUSUM-P
##
data<-read.csv("BurstyArrivals/CUSUM-P_bursty_ARL1_lam1_0.6_lam_true0.9.csv",header = FALSE)
dim(data)
ARL1<-mean(as.numeric(data[1,]))

print(ARL1)


MRL1<-median(as.numeric(data[1,]))
TMRL1<-mean(sort(as.numeric(data[1,]),decreasing=TRUE)[1:10])
print(ARL1)
print(MRL1)
print(TMRL1)

