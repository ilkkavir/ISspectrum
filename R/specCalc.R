specCalc <- function(pldfv,nin0,tit0,mim0,psi,vi,kd2,om){

  nom  <- length(om)
  nion <- length(mim0)-1

  tmp <- .C("specCalc",
             pldfvPr = as.double(Re(pldfv)),
	     pldfvPi = as.double(Im(pldfv)),
	     nin0Pr  = as.double(nin0),
	     tit0Pr  = as.double(tit0),
	     nion    = as.integer(nion),
	     mim0Pr  = as.double(mim0),
	     psiPr   = as.double(psi),
	     viPr    = as.double(vi),
	     kd2Pr   = as.double(kd2),
	     scr     = as.double(rep(0,(3*(nion+1)*(1+nom)))),
	     nom     = as.integer(nom),
	     omPr    = as.double(om),
	     resPr   = as.double(rep(0,nom)),
	     ifref   = as.integer(0)
	   )

  return(tmp$resPr)

} # specCalc