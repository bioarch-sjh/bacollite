
classify <- function(froot,spots,pepnames,pepdata,name,manualID="unknown",gauss = 0.3,verbose=F,doplot=T){

  sample <- load.sample(froot = froot,spots = spots,name=name)

  cd <- list()
  for(spp in 1:length(pepnames))
    cd[[spp]] <- ms_fit(peptides = pepdata[[spp]], sample = sample, doplot = F, force= T, gauss = gauss)

  labs <- pepnames

  #Function
  corlim = seq(0,1,0.05)
  scores <- vector(length = length(cd))
  scores[] <- 0
  laglim <- 0.6

  #massage the raw cordata into a form we can work with:
  cld <- list()
  for(cc in 1:length(cd)){
    cld[[cc]] <- corlim_data(cd[[cc]],laglim)
    #initialise the scores for each sample
    cld[[cc]]$cumscore <- 0
  }

  for(cl in 1:length(corlim)){
    for(ss in 1:length(cld)){
      nh <- cld[[ss]]$nh[cl]
      maxonh<-0
      for(tt in 1:length(cld)){
        if(ss != tt){
          maxonh <- max(maxonh,cld[[tt]]$nh[cl])
        }
      }

      if(nh > maxonh){
        cld[[ss]]$cumscore[cl] <- (nh-maxonh)*corlim[cl]
      }
    }
  }

  result <- data.frame("id" = labs, "score" = 0, stringsAsFactors = F)

  for(ss in 1:length(cld)){
    if(verbose)message(sprintf("Score for species %d (%s) = %f" ,ss,labs[ss],sum(cld[[ss]]$cumscore)))
    result$score[ss] <- sum(cld[[ss]]$cumscore)
  }

  if(doplot){
    plot.classification(result,sample$name,manualID,corlim,cld,labs)
  }

  return(result)

}
