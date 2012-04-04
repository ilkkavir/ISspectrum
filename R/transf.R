transf <- function(p,pm0){

  nion  <- length(pm0)

  nin0 <- tit0 <- mim0 <- psi <- vi <- rep(0,(nion+1))

  return(.C("Transf",
	     pPr    = as.double(p),
	     nin0Pr = as.double(nin0),
	     tit0Pr = as.double(tit0),
	     mim0Pr = as.double(mim0),
	     psiPr  = as.double(psi),
	     viPr   = as.double(vi),
	     p_m0   = as.double(pm0),
	     nion   = as.integer(nion)
	   )
	 )

} # transf 