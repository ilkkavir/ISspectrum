ISspectrum.guisdap <- function(p=c(1e11,300,1,0,0,.3),pm0=c(30.5,16),fradar=933e6,scattAngle=180,freq=seq(-1000,1000)*10){
# 
# An incohrent scatter spectrum with similar input that is used in GUISDAP output files
#
# INPUT:
#   p          c(Ne,T1,Te/Ti,coll,Vi,pm2)
#   pm0        c(m1,m2), in amu
#   fradar     radar carrier frequency in Hz
#   scattAngle scattering angle in degrees, 180 for backscattering
#   freq       frequency axis in Hz
#
# OUTPUT:
#   a vector of power spectral densities at the given frequencies
#
# I. Virtanen 2011
# 

 # Boltzmann constant
 kb <- 1.3806503e-23

 # atomic mass unit
 amu <- 1.66053886e-27

 # electron charge
 q <- 1.60217646e-19

 # permittivity of free space
 eps0 <- 8.85418782e-12
 
 # debye length for Te=Ti, specCalc will add a correction according to the temperature ratio
 D <- sqrt(eps0*kb/q**2)
 
 # speed of light
 c <- 299792458.

 # scattering wave number
 ks <- 4*pi*fradar/c * sin(scattAngle/2*pi/180)
 
 # (k*D)^2
 kd2 <- ks**2 * D**2

 # angular frequency axis in normalized units
 om0 <- (ks * sqrt(2*kb/pm0[1]/amu))
 om <- freq*2*pi / om0

 # transform the parameters, collision frequency is normalised with om0, and ion velocity also converted to Doppler shift
 tpar <- transf( (p * c( 1 , 1 , 1 , 1/om0 , ks/om0 , rep(1,length(p)-5) ) ) , pm0 )

 # plasma dispersion function interpolation table
 # the table is stored in the global workspace during first call in each R session
 if(!exists('pldfIntTab')){
   pldfIntTab <<- pldfvv()
  }

 # the spectrum, no multiplication with p[1], because specCalc automatically multiplies with it
 s <- specCalc(pldfIntTab,tpar$nin0Pr,tpar$tit0Pr,tpar$mim0Pr,tpar$psiPr,tpar$viPr,kd2,om) * 9.978688e-29 / om0

 return(s)

}

