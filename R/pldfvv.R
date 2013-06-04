pldfvv <- function(){
# 
# plasma dispersion function interpolation table for spectrum calculations
# 
# I. Virtanen 2011
# 

  Q <- matrix( c( -1.0/6.0 , 0.5 , -1.0/3.0 , 0.0 , 0.5 , -1.0 , -0.5 , 1.0 , -0.5 , 0.5 , 1.0 , 0.0 , 1.0/6.0 , 0.0 , -1.0/6.0 , 0.0 ),
       	       ncol=4 , nrow=4)

  f <- matrix(0+0i,nrow=67,ncol=65)

  for(m in seq(-1,(16*4+1))){
    for(n in seq(0,(16*4))){
      z <- (m-1i*n)/16.0
      f[(m+2),(n+1)] <- pldf(matrix(z))
    }
  }

  pldfvv <- rep(0+0i,16640)

  for(n in seq(0,16*4)){
    for(m in seq(0,(16*4-1))){
      apu <- Q%*%f[(m+1):(m+4),(n+1)]
      pldfvv[(64*4*n + 4*m)+seq(4)] <- apu
    }
  }

  return(list(Re=Re(pldfvv),Im=Im(pldfvv)))

}
