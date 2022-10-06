
#' @export
plot.classification <- function(result,samplename,manualID,corlim,cld,labs){

  #if(max(result$score)>0.001){
  #  title = sprintf("Sample '%s':\n manual ID: '%s'; Calc ID: '%s'\nscores ",
  #                  samplename,manualID,result$id[result$score == max(result$score)])
  #}
  #else{
  #  title = sprintf("Sample '%s':\n manual ID: '%s'; Calc ID: 'unknown'\nscores ",samplename,manualID)
  #}

  title = ""
  for(ss in 1:nrow(result)){
    title = sprintf("%s %s = %0.3f",title,result$id[ss],result$score[ss])
  }

  par(mar = c(4,4,5,4))
  plot(NA,xlab="Correlation Threshold",ylab = "Number of Hits",ylim=c(0,3*nrow(pepdata[[1]])),xlim = c(0,1), main = title)

  pcols = c("orange","blue","red","black")

  #todo: this is a hacky way of handling up to four candidate species
  plabs = c("","","","")

  for(cc in 1:length(cld)){
    plabs[cc] <- labs[cc]
    points(x=corlim,y=cld[[cc]]$nh,col=pcols[cc])
    lines (x=corlim,y=cld[[cc]]$nh,col=pcols[cc])

  }


  legend("topright",legend = labs,col=pcols[1:length(cld)],lty = 1,pch=1)
}



#' @export
ppclassify <- function(froot,pdrow,pepdata,pepnames,gauss = 0.3,verbose=F,doplot=T,cormin=0){

  spots <- c(pdrow$spot1,pdrow$spot2,pdrow$spot3)
  #reps<- load.sample(fr,pd$sampleID[rr],spots)
  sample <- load.sample(froot = sprintf("%s%s",froot,pdrow$froot),spots = spots,name=pdrow$sampleID)

  cd <- list()
  for(spp in 1:length(pepdata))
    cd[[spp]] <- ms_fit(peptides = pepdata[[spp]], sample = sample, doplot = F, force= T, gauss = gauss)

  labs <- pepnames

  #Function
  corlim = seq(cormin,1,0.05)
  scores <- vector(length = length(cd))
  scores[] <- 0
  laglim <- 0.6

  #massage the raw cordata into a form we can work with:
  cld <- list()
  for(cc in 1:length(cd)){
    cld[[cc]] <- corlim_data(cd[[cc]],laglim,corlim)
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
    plot.classification(result,sample$name,pdrow$manualID,corlim,cld,labs)
  }

  return(result)

}
