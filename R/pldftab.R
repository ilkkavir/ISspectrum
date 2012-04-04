# pldftab.m : function to create the plasma dispersion function table
# GUISDAP v.1.60 96-05-27 Copyright Asko Huuskonen and Markku Lehtinen
#
# Conversion to R by Ilkka Virtanen 2011
# There is no R version of pldfvinit.m, but pldftab should be called with default arguments instead
#
# res=pldftab(dx,dy,nx,ny,nx1,ny1) 
#
#  function res=pldftab(dx,dy,nx,ny,nx1,ny1) 
#
#  res=zeros(nx*ny,1);
#  for i=1:nx  
#    disp(i);
#    for j=1:ny 
#      res(i+(j-1)*nx,1)=pldf( (nx1+i-1)*dx +sqrt(-1)*(ny1+j-1)*(-dy) );
#    end  
#  end

pldftab <- function(dx=.075,dy=.075,nx=61,ny=48,nx1=-3,ny1=0){
  res <- matrix(0,nrow=nx*ny,ncol=1)
  for(i in seq(nx)){
#    print(i)
    for(j in seq(ny)){
      res[i+(j-1)*nx,1] <- pldf( matrix((nx1+i-1)*dx + 1i*(ny1+j-1)*(-dy)))
    }
  }
  return(res)
}