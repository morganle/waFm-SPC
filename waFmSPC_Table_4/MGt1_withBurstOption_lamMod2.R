### note: this file is now simulating an M/U/1 queue

##Slength and Stimes hold the system state changes and times they changed.

#original function starts from empty
MGt1 <- function(Starting.seed,c,Tend,burst,when,len,lam_orig,lam_burst){
  ## if burst = TRUE then we have a burst in service times at time point "when" for length "len"
  # 'when' is used to delete warmup and set lambda change if any
  # Main simulation loop; execution starts here
  set.seed(Starting.seed)
  Queue <<- NULL
  #S <<- 0
  Busy <<- 0      ## is the system empty?
  Clock <<- 0
  Ndept<<-0
  Waits <<-  0 #### not used
  Slength <<- numeric(length = as.integer(2*Tend))
  Stimes <<- numeric(length = as.integer(2*Tend))
  sy_idx <<- 1 # index to keep track of number in system (different from s_idx)
  Qlength <<- numeric(length = as.integer(2*Tend))
  Qtimes <<- numeric(length = as.integer(2*Tend))
  q_idx <<- 1
  IAmean <<- as.double(0)  ### not used
  IAmean_last <<- as.double(0) ### not used
  IAvar <<- as.double(0)   ### not used
  IAsqsum <<- as.double(0) ### not used
  IAsqsum_last <<- as.double(0) ### not used
  IAn <<- as.double(0)     ### not used
  
  # make Tend global
  Tend <<- Tend
  
  ## arrival distribution parameters, with values generated in advance
  lam_orig <<- lam_orig
  lam_burst <<- lam_burst 
  arr_vec <<- rexp(as.integer(1.5*Tend),rate=lam_orig)
  arr_vec_b <<- rexp(as.integer(1.5*Tend),rate=lam_burst)
  a_idx <<- 1 #index for using interarrival times
  
  ## service distribution parameters, with values generated in advance
  s_min_orig <<- 0.5
  s_max_orig <<- 1.5
  ser_vec <<- runif(as.integer(1.5*Tend),min=s_min_orig,max=s_max_orig)
  s_idx <<- 1
  
  SSmean<<-as.double(0) ### not used
  SSmean_last<<-as.double(0) ### not used
  SSn<<-as.double(0) ### not used
  
  bursts <<-when   # when the burst starts ### not used - just use 'when'
  period <<-when+len  # when the burst ends ### not used
  
  # create event calendar and schedule initial events need c+3 for clear event
  EventTime <<- rep(Inf, c+3)
  EventType <<- rep("", c+3)
  Schedule("Arrival", ArrivalDistribution())
  Schedule("ClearAndChangeLambda", when)
  Schedule("EndSimulation", (Tend + when)) # add to length to discard warmup
  NextEvent <- ""
  
  #period_length<<-len
  while(NextEvent != "EndSimulation"){
    #print(NextEvent)
    Qlength[q_idx] <<- length(Queue)
    Qtimes[q_idx] <<- Clock
    q_idx <- q_idx + 1
    
    if(Busy==0){Slength[sy_idx]<<-0}else{Slength[sy_idx]<<-length(Queue)+Busy}
    Stimes[sy_idx] <<- Clock
    sy_idx <<- sy_idx + 1
    
    NextEvent <- TimerQ()
    switch(NextEvent,
           Arrival = Arrival(),
           Departure = Departure(),
           ClearAndChangeLambda = ClearAndChangeLambda()
    )
  } 
  OO<-list(Waits,Qlength[1:(q_idx-1)],Qtimes[1:(q_idx-1)]-when,
           Slength[1:(sy_idx-1)],Stimes[1:(sy_idx-1)]-when)
  return(OO)
}


# Define arrival and service distributions here  
ArrivalDistribution <- function(){
    #arr<-rexp(1,rate=lam_orig)
    arr = arr_vec[a_idx]
    a_idx <<- a_idx + 1
    return(arr)
}

ServiceDistribution <- function(){
    #ser<-runif(1,min=s_min_orig,max=s_max_orig)
    ser = ser_vec[s_idx]
    s_idx <<- s_idx + 1
    return(ser)
  
}

# Event calendar management functions
Schedule <- function(Event, Time){
  # inserts events into the event calendar
  CalIndex <- which.max(EventTime)  
  EventTime[CalIndex] <<- Clock + Time
  EventType[CalIndex] <<- Event
}

TimerQ <- function(){
  # find next event, update Clock, and return event type
  CalIndex <- which.min(EventTime)
  Clock <<- EventTime[CalIndex]
  #if(Clock > (period)){
   # burst<<-burst+period_length
   # period<<-period+period_length
  #  }
  EventTime[CalIndex] <<- Inf
  return(EventType[CalIndex])
}


# Events

Arrival <- function(){
  # models the arrival of a customer event
  Schedule("Arrival", ArrivalDistribution())
  #S<<-S+1
  if (Busy < c){
    Busy <<- Busy + 1
    Waits <<- c(Waits, 0)
    Schedule("Departure", ServiceDistribution())
  }
  else{
    Queue <<- c(Queue, Clock)
  }
}

ClearAndChangeLambda <- function(){
  #print("Clear change lambda")
  Waits <<-  0 #### not used
  Slength <<- NULL
  Slength <<- numeric(length = as.integer(2*Tend))
  Stimes <<- NULL
  Stimes <<- numeric(length = as.integer(2*Tend))
  sy_idx <<- 1
  Qlength <<- NULL
  Qlength <<- numeric(length = as.integer(2*Tend))
  Qtimes <<- NULL
  Qtimes <<- numeric(length = as.integer(2*Tend))
  q_idx <<- 1
  arr_vec <<- arr_vec_b #update to pull interarrival times from burst data
  
}

Departure <- function(){
  # models departure of a customer event
  NumQ <- length(Queue)
  if (NumQ > 0){
    Schedule("Departure", ServiceDistribution())
    if (NumQ == 1){
      Queue <<- NULL
    }
    else {
      Queue <<- Queue[2:NumQ]
    }
  }
  else{
    Busy <<- Busy - 1
  }   

  #S<<-S-1
} 


