
#' get the hits for a range of correlation limit values
#'
#' @param cd list of outputs from the ms_fit function for a set of peptides
#' @param corlim a sequence of correlation limit values for thresholding the data
#' @param laglim maximum acceptable value of lag for each correlation
#' @keywords bruker
#' @export
corlim_hits <- function( cd, corlim=seq(0,1,0.05), laglim = 0.6){

  #massage the raw cordata into a form we can work with:
  cld <- list()
  for(cc in 1:length(cd)){
    cld[[cc]] <- corlim_data(cd[[cc]],laglim)
  }

  for(ss in 1:length(cld))
    cld[[ss]]$cumscore = 0

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

  return(cld)
}
