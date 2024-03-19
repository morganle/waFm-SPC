# WindowMod.R
# This function modifies windowed data for Fourier transform
# Requires as input the window data (w_vec), the method (w_method)
# Output is the Fourier coefficients for frequencies up to
#   1/20 the length of w_vec
WindowMod = function (w_method,w_vec){
  
  # length of w_vec
  N = length(w_vec)
  a0 = .54 # needed for Hamming
  
  # case for windowing method implemented via 'switch'
  result = switch(w_method,
  # no change
  rect = w_vec,
  
  # baseline tilt
  baseline = w_vec - ((w_vec[N]-w_vec[1])*(0:(N-1))/(N-1)),
  
  # mirror/flip
  flip = c(w_vec,w_vec[N:1]),
  
  # Hamming
  Hamming = w_vec*(a0-(1-a0)*cos(2*pi*0:(N-1)/(N-1))))
  return(result)
}

w_vec = sin(.14*(1:128)) + .05*(1:128)+ runif(128,-.2,.2)
plot(w_vec)
w_type = list("rect","baseline","flip","Hamming")
window_mods = sapply(w_type,WindowMod,w_vec=w_vec)
for (i in 1:4) {plot(window_mods[[i]])}


