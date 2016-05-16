ISspectrum.3D <- function(ele=c(1e11,300,300,0,0,0,0),ion=list(c(30.5,.7e11,300,300,0,0,0,0),c(16,.3e11,300,300,0,0,0,0),c(1,0,300,300,0,0,0,0)),Bdir=c(0,0,1),kdir=c(0,0,1),fradar=993e6,scattAngle=180,freq=seq(-1000,1000)*100){
#
# a general three-dimensional incoherent scatter model, with max 7 singly-charged ions
#
# INPUT:
#  ele        c(Ne,Te_par,Te_perp,nu_en,ve_x,ve_y,ve_z)
#  ion        list( c(m1,N1,T1_par,T1_perp,nu_1n,v1_x,v1_y,v1_z) , c(m2,N2,T2_par,T2_perp,nu_2n,v2_x,v2_y,v2_z) , ... )
#  Bdir       c(x,y,z) Magnetic field direction
#  kdir       c(x,y,z) Scattering wave vector direction
#  fradar     radar frequency in Hz
#  scattAngle scattering angle (the angle between incident and scattered wave vectors, 180 for backscattering)
#  freq       frequency axis
#
#  ion masses in amu
#  ion densities can be given as absolute values, but they are treated as relative abundances, so that sum(Ni) = Ne
#  only singly-charged ions
#  radar frequency in Hz
#  scattering angle in degrees
#  frequency axis in Hz
#
# OUTPUT:
#  a vector of power spectral densities at the given frequencies
#
# I. Virtanen 2012
#

  # conversion to unit vectors
  bunit <- Bdir/sqrt(sum(Bdir**2))
  kunit <- kdir/sqrt(sum(kdir**2))

  # angle between magnetic field direction and scattering wave vector
  bangle <- acos(sum(bunit*kunit))

  # projected parameters
  # velocity seen at a site is the projection of the true velocity to the scattering wave vector direction,
  ele2   <- c( ele[1], (ele[2]*cos(bangle)**2 + ele[3]*sin(bangle)**2) , ele[4] , sum(ele[5:7]*kunit) )
  ion2   <- list()
  for(k in seq(length(ion))){
    ion2[[k]] <- c( ion[[k]][1:2] , (ion[[k]][3]*cos(bangle)**2 + ion[[k]][4]*sin(bangle)**2) , ion[[k]][5] , sum(ion[[k]][6:8]*kunit) )
  }

  return(ISspectrum.general(ele=ele2,ion=ion2,fradar=fradar,scattAngle=scattAngle,freq=freq))

} # ISspectrum.3D
