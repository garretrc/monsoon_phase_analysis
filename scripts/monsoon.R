#### Phase alignment and monsoon timing ####

#ndims and trapz from fdasrvf 2.2.0 (J. Derek Tucker; GPL-3)
ndims <- function(x){
  return(length(dim(x)))
}

trapz <- function(x,y,dims=1){
  if ((dims-1)>0){
    perm = c(dims:max(ndims(y),dims), 1:(dims-1))
  } else {
    perm = c(dims:max(ndims(y),dims))
  }
  
  if (ndims(y) == 0){
    m = 1
  } else {
    if (length(x) != dim(y)[dims])
      stop('Dimension Mismatch')
    y = aperm(y, perm)
    m = nrow(y)
  }
  
  if (m==1){
    M = length(y)
    out = sum(diff(x)*(y[-M]+y[-1])/2)
  } else {
    slice1 = y[as.vector(outer(1:(m-1), dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+')) ]
    dim(slice1) = c(m-1, length(slice1)/(m-1))
    slice2 = y[as.vector(outer(2:m, dim(y)[1]*( 1:prod(dim(y)[-1])-1 ), '+'))]
    dim(slice2) = c(m-1, length(slice2)/(m-1))
    out = t(diff(x)) %*% (slice1+slice2)/(2.0)
    siz = dim(y)
    siz[1] = 1
    out = array(out, siz)
    perm2 = rep(0, length(perm))
    perm2[perm] = 1:length(perm)
    out = aperm(out, perm2)
    ind = which(dim(out) != 1)
    out = array(out, dim(out)[ind])
  }
  
  out
}

efda_warp = function(f1, f2, time, lambda = 0, pen = "roughness", method = 'DP',nbhd_dim=7,t_arrive,t_peak,t_retreat){
  q1 <- f_to_srvf(f1, time)
  q2 <- f_to_srvf(f2, time)
  gam <- optimum.reparam(q1, time, q2, time, lambda, pen, method = method,nbhd_dim = nbhd_dim)
  fw <- warp_f_gamma(f2, time, gam)
  qw <- warp_q_gamma(q2, time, gam)
  Dy <- sqrt(trapz(time, (q1 - qw)^2))
  time1 <- seq(0, 1, length.out = length(time))
  binsize <- mean(diff(time1))
  psi <- sqrt(gradient(gam, binsize))
  q1dotq2 = trapz(time1, psi)
  if (q1dotq2 > 1) {
    q1dotq2 = 1
  }
  else if (q1dotq2 < -1) {
    q1dotq2 = -1
  }
  Dx <- acos(q1dotq2)
  
  list(c(Dy, Dx, f1[1]-f2[1], mean(f2)-mean(f1),
         gam[t_arrive], gam[t_peak], gam[t_retreat]),
       fw,
       gam)   
}

efda_timing_bias = function(c1,c2,mt_a=0.5,mt_r=0.5){
  t_onset = min(which(c1>=mt_a*max(c1)))
  t_peak = which.max(c1)
  t_retreat = max(which(c1>=mt_r*max(c1)))
  time=1:365
  res = efda_warp(c1,c2,time, 
                  lambda = 0, pen = "roughness", method = 'DP',nbhd_dim=7,
                  t_onset,t_peak,t_retreat)
  ob = 364*res[[1]][5]+1
  pb = 364*res[[1]][6]+1
  rb = 364*res[[1]][7]+1
  c(ob-t_onset,pb-t_peak,rb-t_retreat)
}

calc_tb = function(x){
  tb = array(NA,c(3,dim(x)[2:3]))
  for(i in 1:dim(x)[2]){
    for(j in 1:dim(x)[3]){
      tb[,i,j] = efda_timing_bias(gpcp_ismr[,i,j],x[,i,j])
    }
  }
  tb
}


gpcp_ismr = climatologies$slices$GPCP
onset = apply(gpcp_ismr,2:3,function(x){min(which(x>=0.5*max(x)))})
retreat = apply(gpcp_ismr,2:3,function(x){max(which(x>=0.5*max(x)))})
mods_tb = lapply(mods,function(model){
  message('Monsoon timing: ',model)
  calc_tb(climatologies$slices[[model]])
})
names(mods_tb) = mods
