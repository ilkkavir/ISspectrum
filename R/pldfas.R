pldfas <- function(z){

  tmp <- .C("pldfas"                 ,
              retR = as.double(0.0)  ,
     	      retI = as.double(0.0)  ,
	      zR   = as.double(Re(z)),
	      zI   = as.double(Im(z))
	    )

  return(tmp$retR+1i*tmp$retI)

}