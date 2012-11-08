ISspectrum.general <- function(ele=c(1e11,300,0,0),ion=list(c(30.5,.7e11,300,0,0),c(16,.3e11,300,0,0),c(1,0,300,0,0)),fradar=933e6,scattAngle=180,freq=seq(-1000,1000)*10){
# 
# a general incoherent scatter spectrum calculation, with max 7  singly-charged ions. 
# 
# INPUT:
#  ele         c(Ne,Te,nu_en,ve)
#  ion         list( c(m1,N1,T1,nu_1n,v1) , c(m2,N2,T2,nu_2n,v2) , ... )
#  fradar      radar frequency
#  scattAngle  scattering angle (the angle between incident and scattered wave vectors, 180 for backscattering)
#  freq        frequency axis
# 
#  ion masses in amu
#  ion densities can be given as absolute values, but they are treated as relative abundances, so that sum(Ni) = Ne
#  only singly-charged ions at the moment
#  radar frequency in Hz
#  scattering angle in degrees
#  frequency axis points in Hz
#
#  OUTPUT:
#   a vector of power spectral densities at the given frequencies
# 
# I. Virtanen 2011
# 


 # transform the parameters
 if(is.list(ion)){
    nIon <- length(ion)
  }else{
    ion <- list(ion)
    nIon <- 1
  }

  if(nIon>7) stop('Max 7 ion species.')

  # allocate vectors for specCalc input parameters
  nin0 <- tit0 <- mim0 <- psi <- vi <- rep(0,nIon+1)

  # Boltzmann constant
  kb <- 1.3806503e-23

  # atomic mass unit
  amu <- 1.66053886e-27

  # electron mass in amy
  me <- 0.000548579867

  # electron charge
  q <- 1.60217646e-19

  # permittivity of free space
  eps0 <- 8.85418782e-12
 
  # debye length
  D <- sqrt(eps0*kb/q**2)
 
  # speed of light
  c <- 299792458.

  # scattering wave number for backscattering
  ks0 <- 4*pi*fradar/c

  # actual scattering wave number
  ks <- ks0 * sin(scattAngle/2*pi/180)

  # (k*D)^2
  kd2 <- ks**2 * D**2

  # angular frequency axis in normalized units
  om0 <- (ks * sqrt(2*kb/amu))
  om <- freq*2*pi / om0

  # normalized parameters
  nin0[nIon+1] <- ele[1]		# electron density
  tit0[nIon+1] <- ele[2]     		# electron temperature
  mim0[nIon+1] <- me         		# electron mass 
  psi[nIon+1]  <- ele[3]     		# electron-neutral collision frequency
  vi[nIon+1]   <- 4*pi*fradar*ele[4]/c  # Doppler shift corresponding the electron velocity

  nsum <- 0
  for (k in seq(nIon)){nsum <- nsum + ion[[k]][2]}
  if(nsum==0){
    nscale <- 1
  }else{
    nscale <- nsum/ele[1]
  }
  for(k in seq(nIon)){
    nin0[k] <- ion[[k]][2]/nscale
    tit0[k] <- ion[[k]][3]
    mim0[k] <- ion[[k]][1]
    psi[k]  <- ion[[k]][4]
    vi[k]   <- 4*pi*fradar*ion[[k]][5]/c
  }
  psi <- psi / om0
  vi  <- vi  / om0

 # plasma dispersion function interpolation table
 # the table is stored in the global workspace during first call in each R session
 if(!exists('pldfIntTab')){
   pldfIntTab <<- pldfvv()
  }

 # the spectrum  (4*pi*r_e**2 = 9.978688e-29 m**2, Vallinkoski 1989)
 s <- specCalc(pldfIntTab,nin0,tit0,mim0,psi,vi,kd2,om) * 9.978688e-29 / om0

 return(s)

} # ISspectrum.general
