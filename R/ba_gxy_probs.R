



#' Calculate the 'GXY' frame of a peptide, and use probabilities if hdroxylation for prolines at different positions
#' to estimate the probability of every possible level of hydroxylation
#'
#' @param seq a peptide sequence
#' @param xprob the probability of hydroxylation at the X position
#' @param yprob the probability of hydroxylation at the Y position
#' @param xyxprob the probability of hydroxylation at the X position if the Y position is also Proline
#' @param xyyprob the probability of hydroxylaiton at the Y position if the X position is also Proline
#' @param verbose set to TRUE for more verbose processing
#' @keywords "Mass spectrum",
#' @export
#' @examples
#' ba_gxy_probs
ba_gxy_probs <- function(seq, xprob = 0.05, yprob = 0.95, xyxprob = 0.75, xyyprob = 0.99, verbose = F){

  #ERROR CHECKING
  if(str_detect(seq,"G") == F){
    message("No G in this sequence - can't determine X-Y proline position\nReturning NA...")
    return (NA)
  }
  if(str_detect(seq,"P") == F){
    message("No P in this sequence - can't determine any X-Y proline positions\nReturning NA...")
    return (NA)
  }




  pclass <- c("X","Y")

  if(verbose)
    message(sprintf("Processing seq %s",seq))

  gs <- as.data.frame(str_locate_all(seq,"G"))
  #figure out the gs position
  gs <- as.integer(unique(gs$start%%3))

  message(sprintf("Gs at position %d",gs))

  #create a data frame of all the prolines
  ps <- as.data.frame(str_locate_all(seq,"P"))

  if(nrow(ps)==0)
    return(NA)

  ps$prob<- 0

  PPmotif <- F
  for(pp in 1:nrow(ps)){

    if(PPmotif){
      PPmotif <- F
    }
    else{
      ppos <- ps$start[pp]
      pmod <- ppos%%3
      pidx <- NA

      if(gs == 0){#pmod can be 1 or 2
        if(pmod == 1)pidx <- 1 #"X"
        if(pmod == 2)pidx <- 2 #"Y"
      }
      if(gs == 1){#pmod can be 0 or 2
        if(pmod == 2)pidx <- 1 #"X"
        if(pmod == 0)pidx <- 2 #"Y"
      }
      if(gs == 2){#pmod can be 0 or 1
        if(pmod == 0)pidx <- 1 #"X"
        if(pmod == 1)pidx <- 2 #"Y
      }

      if(pidx == 1)
        ps$prob[pp]<-xprob
      if(pidx == 2)
        ps$prob[pp]<-yprob

      message(sprintf("P found at %d, pmod is %d, pclass is %s, prob is %f",ppos,pmod,pclass[pidx],ps$prob[pp]))

      #SEE IF THERE'S A PROLINE AT THE OTHER POSITION
      if(pidx == 1){#X position
        #message("ppos is 0")
        if(ppos < str_length(seq)){
          #message("Y posn available")
          if(str_sub(seq,ppos+1,ppos+1) == "P"){
            message(sprintf("pp = %d: GPP motif found at position %d",pp,ppos))
            PPmotif <- T

            ps$prob[pp]<-0.5
            ps$prob[(pp+1)]<-0.95

          }
        }
      }
    }
  }

  #TODO: print ps if(verbose)
  #print(ps)

  nhyd <- ba_nhyd(ps$prob,verbose = F)
  #print(nhyd)

  return(nhyd)

}
