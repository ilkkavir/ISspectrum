adminvec <- function(pldfv,om,vi,psi,gam){
  
  nom <- length(om)
  tmp <- .C("adminvec"   ,
              nom     = as.integer(nom),
  	      pldfvPr = as.double(Re(c(pldfv))),
      	      pldfvPr = as.double(Im(c(pldfv))),
      	      apuprv  = as.double(rep(0,nom)),
      	      apupiv  = as.double(rep(0,nom)),
      	      omv     = as.double(om),
      	      vi      = as.double(vi),
      	      psi     = as.double(psi),
      	      gam     = as.double(gam)
    	    )
  
  return(tmp$apuprv + 1i*tmp$apupiv)

} # adminvec