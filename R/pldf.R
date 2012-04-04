# pldf.m: complex plasma dispersion function  
# GUISDAP v.1.60 96-05-27 Copyright Asko Huuskonen and Markku Lehtinen
#
# Conversion to R by Ilkka Virtanen 2011
#
#
# accuracy:  8 numbers in the whole complex plane
#            (main algorithm accurate to 1e-12, accuracy is
#            limited by subroutine cerfexp and numerical accuracy)
#
# method:    approximation described by Salzer in
#            Math. Tables Aids Comput., Vol. 5 (1951) pp. 67-70.
#            The same formula appears in
#            Abramovitch-Stegun p.299 formula 7.1.29.
#            (the exponentials in the formulas have been written
#            in a slightly modified way to get rid of unnecessary
#            overflows)
#                                           Markku Lehtinen 10.3.1979
#
# res=pldf(z)
#
#  function res=pldf(z)
##
##  if max(size(z))==1
##    res=0;
##  else
##    res=zeros(z); 
##  end
#  [nrow ncol]=size(z);
#  res=zeros(nrow, ncol);
#  j=sqrt(-1);
#  for irow=1:nrow
#    for icol=1:ncol
#      z1=-j*z(irow,icol);
#      x=real(z1);
#      y=imag(z1);
#      if(abs(2*x*y)<30000)
#        cs=cos(2*x*y);
#        sn=sin(2*x*y);
#      else
#        cs=1; sn=0;
#      end
#      pisqr=sqrt(pi);
#      fn=max(1,floor(abs(2*y))-11):floor(abs(2*y))+11;
#      term1=exp(-fn.*fn/4-y*y);
#      if max(abs(term1))<1e-100, term1=zeros(fn); end
#      term2=exp(-(fn/2-y).*(fn/2-y))/2;
#      term3=exp(-(fn/2+y).*(fn/2+y))/2;
#      factor=ones(size(fn))./(fn.*fn+4*x*x);
#      sume=sum(factor.*(term1*cs-term2-term3));
#      sumf=sum(factor.*(term1*2*x*sn+fn.*(term2-term3)));
#      if (abs(x*y)<1e-4)
#        e1=(-exp(-y*y)*x*y*y/2+2*x*sume)/pisqr;
#        f1=(exp(-y*y)*y/2+sumf)/pisqr;
#      else
#         e1=(exp(-y*y)*(cs-1)/4/x+2*x*sume)/pisqr;
#         f1=(exp(-y*y)*sn/4/x+sumf)/pisqr;
#      end
#      if abs(e1)<1e-100, e1=0; end
#      z1=2*j*(e1+j*f1)+j*pisqr*exp(-y*y)*(cs+j*sn)*cerfexp(-x);
#      res(irow,icol)=z1;
#    end 
#  end

pldf <- function(z){
  dims <- dim(z)
  nrow <- dims[1]
  ncol <- dims[2]
  res  <- matrix(0,nrow=nrow,ncol=ncol)
  for(irow in seq(nrow)){
    for(icol in seq(ncol)){
      z1 <- -1i*z[irow,icol]
      x  <- Re(z1)
      y  <- Im(z1)
      if(abs(2*x*y)<30000){
        cs <- cos(2*x*y)
        sn <- sin(2*x*y)
      }else{
        cs <- 1
        sn <- 0
      }
      pisqr <- sqrt(pi)
      fn <- seq( max(1,floor(abs(2*y))-11) , (floor(abs(2*y))+11) )
      term1 <- exp(-fn*fn/4-y*y)
      if( max(abs(term1)) < 1e-100 ) term1 <- rep(0,length(fn))
      term2 <- exp( -(fn/2-y) * (fn/2-y) )/2
      term3 <- exp( -(fn/2+y) * (fn/2+y) )/2
      factor <- rep(1,length(fn)) / (fn*fn+4*x*x)
      sume <- sum(factor*(term1*cs-term2-term3))
      sumf <- sum(factor*(term1*2*x*sn+fn*(term2-term3)))
      if(abs(x*y)<1e-4){
        e1 <- (-exp(-y*y)*x*y*y/2+2*x*sume ) / pisqr
        f1 <- ( exp(-y*y)*y/2+sumf ) / pisqr
      }else{
        e1 <- ( exp(-y*y)*(cs-1)/4/x + 2*x*sume ) / pisqr
        f1 <- ( exp(-y*y)*sn/4/x+sumf ) / pisqr
      }
      if(abs(e1)<1e-100) e1 <- 0
      z1 <- 2i*(e1+1i*f1) + 1i*pisqr*exp(-y*y)*(cs+1i*sn)*cerfexp(-x) 
      res[irow,icol] <- z1
    }
  }

  return(res)

} # pldf