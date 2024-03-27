## we read in the Phase I and Phase II data and then use the uniroot function to optimise our ARL_0 to the value we want it to be so waFm CUSUM and CUSUM-C have 
## comparable ARL_0 (in-control run-length)

data =read.csv("PhaseII_magnitudes_lambda_pt7NoGibbs.csv",header = FALSE)  # G=100 macro reps W=1000 contiguous windows waFm stats taking in 0th to 2500th mags
data2=read.csv("PhaseII_magnitudes_lambda_pt7NoGibbs.csv",header=FALSE)    # G=100 macro reps W=1200 contiguous windows
acf(as.vector(as.numeric(data[1,])),1:20,pl=TRUE) # check autocorrelation

#CUSUM
# needs library qcc
library(qcc)
library(tictoc)

tic()
find_param<-function(CL_para,W){
  G<-100
  ooc<-c()

  max_len<-1200 # phase II run length
  for (i in 1:G){
      phaseIdata<-as.numeric(data[i,]) 
      phaseIIdata<-as.numeric(data2[i,])

      q<-cusum(phaseIdata,newdata = phaseIIdata,decision.interval = CL_para,plot=FALSE)   # cusum function is inbuilt within qcc package
      violations<-sort(c(q$violations$lower,q$violations$upper))              # when did violations occur
      if(length(violations)>0 && violations[length(violations)] > W){
        ooc[i] <- violations[min(which(violations>W))] - W
      }else{
        ooc[i]<-max_len    # if no violations then return 1200 max len of Phase II
      }

  }

  ooc[which(ooc==Inf)]=max_len
#  if(length(which(is.na(ooc)==TRUE))>0){
#  ooc<-ooc[-which(is.na(ooc)==TRUE)]
#  }

  print(CL_para)
  print(mean(ooc))
  return(mean(ooc) - 371)   ## we want ARL_0 of 371
}

W<-1000  # length of phase I data

uniroot(find_param,W=W,lower=4,upper=12)  # searches for CL parameter 

toc()




