\name{ISspectrum.general}
\title{ISspectrum.general}
\alias{ISspectrum.general}
\description{General incoherent scatter spectrum}
\usage{ISspectrum.general(ele=c(1e11,300,0,0),ion=list(c(30.5,.7e11,300,0,0),c(16,.3e11,300,0,0),c(1,0,300,0,0)),
                   fradar=933e6,scattAngle=180,freq=seq(-1000,1000)*10)}

\arguments{
\item{ele}{c( Ne [m^-3] , Te [K] , nu_en [s^-1]  , ve [ms^-1] )}
\item{ion}{list( c( m1 [amu] , N1 [m^-3] , T1 [k] , nu_1n [s^-1] , v1 [ms^-1] ) , c(m2,N2,T2,nu_2n,v2) , ... ) }
\item{fradar}{Radar carrier frequency in Hz}
\item{scattAngle}{Scattering angle in degrees, 180 for backscattering.}
\item{freq}{Frequency axis in Hz}
}

\value{
\item{s}{power spectral densities at the given frequencies}
}

\author{Ilkka Virtanen (University of Oulu, Finland) \cr
  \email{ilkka.i.virtanen@oulu.fi}}

\references{The function uses C routines of GUISDAP (Lehtinen and
  Huuskonen, 1996) for calculating the incoherent scatter spectrum}