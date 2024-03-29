#' Calculate the 'GXY' frame of a peptide, and use probabilities of hdroxylation for prolines at different positions
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
ba_gxy_probs <- function(seq, xprob = 0.05, yprob = 0.95, xyxprob = 0.5, xyyprob = 0.99, verbose = F){
  #message("DEBUGGING!")
  if(verbose)
    message("inside ba_gxy_probs")

  nohyd <- data.frame(prob=1,nhyd=0)

  #ERROR CHECKING
  if(str_detect(seq,"G") == F){
    if(verbose)
      message("No G in this sequence - can't determine X-Y proline position\nReturning 0 hydroxylation prob of 1")
    return (nohyd)
  }
  if(str_detect(seq,"P") == F){
    if(verbose)
      message("No P in this sequence - can't determine any X-Y proline positions\nReturning 0 hydroxylation prob of 1")
    return (nohyd)
  }

  pclass <- c("X","Y")

  if(verbose)
    message(sprintf("Processing seq %s",seq))

  #figure out the G position - store it in gpos (NB - may need a strategy for the nonmodal G positions, possibly only process GXYG motifs..?)
  gs <- as.data.frame(str_locate_all(seq,"G"))
  # get the position value
  gm3 <- gs$start%%3
  # get the unique positions (should be 0 1 or 2)
  ux <- as.integer(unique(gm3))
  # get the modal position:
  gpos <- ux[which.max(tabulate(match(gm3,ux)))]

  if(verbose)
    message(sprintf("Gs at position %d",gpos))

  #create a data frame of all the prolines
  ps <- as.data.frame(str_locate_all(seq,"P"))

  if(nrow(ps)==0)# return a probability of 1 for zero hydroxylations
    return(nohyd)

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

      if(verbose)message(sprintf("Gpos is %d, pmod is %d",gpos,pmod))
      if(gpos == pmod){
        if(verbose)
          message(sprintf("Warning! P at position %d is out of sync with GXY",ps$start[pp]))
      }
      else{
        if(gpos == 0){#pmod can be 1 or 2 - CAN'T BE 0!
          if(pmod == 1)pidx <- 1 #"X"
          if(pmod == 2)pidx <- 2 #"Y"
        }
        if(gpos == 1){#pmod can be 0 or 2
          if(pmod == 2)pidx <- 1 #"X"
          if(pmod == 0)pidx <- 2 #"Y"
        }
        if(gpos == 2){#pmod can be 0 or 1
          if(pmod == 0)pidx <- 1 #"X"
          if(pmod == 1)pidx <- 2 #"Y
        }

        if(pidx == 1)
          ps$prob[pp]<-xprob
        if(pidx == 2)
          ps$prob[pp]<-yprob

        if(verbose)
          message(sprintf("P found at %d, pmod is %d, pclass is %s, prob is %f",ppos,pmod,pclass[pidx],ps$prob[pp]))

        #SEE IF THERE'S A PROLINE AT THE OTHER POSITION
        if(pidx == 1){#X position
          if(ppos < str_length(seq)){
            if(verbose)
              message(sprintf("Y posn available, n ppos = %d, sequence is %s",length(ppos),seq))

            if(str_sub(seq,ppos+1,ppos+1) == "P"){
              if(verbose)
                message(sprintf("pp = %d: GPP motif found at position %d",pp,ppos))

              PPmotif <- T
              ps$prob[pp]<-xyxprob
              ps$prob[(pp+1)]<-xyyprob

            }
          }
        }
      }
    }
  }

  #TODO: print ps if(verbose)
  #print(ps)

  # Now we have assigned hydroxylation probabilities to each proline in the sequence, we can calculate the
  # probability of each hydroxylation level
  nhyd <- ba_nhyd(ps$prob,verbose = F)

  return(nhyd)

}
