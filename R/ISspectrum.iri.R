ISspectrum.iri <- function( time = c(2000,1,1,11,0,0) , latitude=69.58864 , longitude=19.2272 , heights=seq(1,1000) , fradar=233e6,scattAngle=180,freq=seq(-1000,1000)*4,savePlasmaParams=FALSE)
  {

    
    iripar <- IRI2016::iriParams( time=time , latitude=latitude , longitude=longitude , heights=heights)


    # the iri model returns -1 at heights where densities are not calculated
    ions <- c('O+','H+','He+','O2+','N+','NO+')
    neutrals <- c('O','H','He','O2','N','N2')

    for(n in c(ions,neutrals)) iripar[n,iripar[n,]<0] <- 0
      
    if(savePlasmaParams) save(iripar,file='ISspectrum.iri.PlasmaParam.Rdata')

    nh <- length(heights)
    nf <- length(freq)
    spectra <- matrix(nrow=nh,ncol=nf)
    for(k in seq(nh)){

      cfreq <- IRI2016::ionNeutralCollisionFrequency(iripar[,k])
      ele <- c(iripar[c('e-','Te'),k],0,0)
      ion <- list(
                  c(16,iripar[c('O+','Ti') , k] ,sum(cfreq['O+',]) ,0),
                  c(1 ,iripar[c('H+','Ti') , k] ,sum(cfreq['H+',]) ,0),
                  c(4 ,iripar[c('He+','Ti'), k] ,sum(cfreq['He+',]),0),
                  c(32,iripar[c('O2+','Ti'), k] ,sum(cfreq['O2+',]),0),
                  c(14,iripar[c('N+','Ti') , k] ,sum(cfreq['N+',]) ,0),
                  c(30,iripar[c('NO+','Ti'), k] ,sum(cfreq['NO+',]),0)
                  )

      spectra[k,] <- ISspectrum.general( ele=ele , ion=ion , fradar=fradar , scattAngle=scattAngle , freq=freq)

    }

    return(spectra)
  }

