\name{ISspectrum.guisdap}
\title{ISspectrum.guisdap}
\alias{ISspectrum.guisdap}
\description{GUISDAP-style incoherent scatter spectrum}
\usage{ISspectrum.guisdap(p=c(1e11,300,1,0,0,.3),pm0=c(30.5,16),fradar=933e6,scattAngle=180,freq=seq(-1000,1000)*10)}


\arguments{
\item{p}{c( Ne [m^-3] , Ti [K] , Te/Ti , nu_in [s^-1]  , vi [ms^-1] ,
  pm2 , pm3 , ...)}
\item{pm0}{list( c( m1 [amu] , m2 , m3 , ... ) }
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