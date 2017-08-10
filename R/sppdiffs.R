

#' get the differences in spectra for two species.
#' NB the sequences are sorted by sequence position *within* this function
#'
#' @param mcs1 the mammalian collagen sequence for species 1
#' @param mcs2 the mammalian collagen sequence for species 2
#' @param name1 (optional) the name of species 1
#' @param name2 (optional) the name of species 2
#' @param verbose (optional) whether to print messages while processing
#' @param halt (optional) whether to pause until <return> is hit while processing
#' @return a vector of sequence positions for the position where the differences in the spectrum are
#' @export
#' @examples
#' sh <- load.mcs("sheep")
#' hu <- load.mcs()
#' diffs <- spp.diffs(sh,hu,"sheep","human")
spp.diffs <- function(mcs1, mcs2, name1="spp01", name2="spp02", verbose = F, halt = F){


  mcs1 <- mcs1[order(mcs1$seqpos),]
  mcs2 <- mcs2[order(mcs2$seqpos),]

  helixoffset <- 17

  if(verbose){
    if(halt)
      readline(sprintf("There are %d entries for spp1 and %d entries for spp2. We'll go through based on sequence position.  Hit <return>",nrow(mcs1),nrow(mcs2)))
    else
      message(sprintf("There are %d entries for spp1 and %d entries for spp2. We'll go through based on sequence position.",nrow(mcs1),nrow(mcs2)))
  }

  #We need to deal with the repetitions due to deam and hyd:
  ih=1
  is=1

  ndiffs = 0

  diff_pos <- NA

  #while(ih < nrow(mcs1) && is <nrow(mcs2)){
  #if testing, use:
  for(ii in 1:100){#nrow(mcs1)){

    #message(sprintf("in while loop, ih = %d, is = %d, s1$seqpos = %d, s1$seqpos = %d",ih,is,mcs1$seqpos[ih],mcs2$seqpos[is]))

    if(mcs1$seqpos[ih] == mcs2$seqpos[is]){

      #message("difference found")

      if(! (mcs1$seq[ih] == mcs2$seq[is])){

        indel <- F

        #NB - should use Smith-Waterman alignment really!
        # see the 'Biostrings' package - but see if we can find something quicker
        if(nchar(mcs1$seq[ih]) != nchar(mcs2$seq[is])){
          message("INDEL")
          indel <- T
        }
        else{
          message("SUBSITUTION")
        }

        #check the codes:
        validseq(mcs1$seq[ih])
        validseq(mcs2$seq[is])

        #if(verbose){
        message(sprintf("%s: helixpos %03d %s",name1,mcs1$seqpos[ih]-helixoffset,mcs1$seq[ih]))
        #}

        if(!indel){
          substr <- mcs1$seq[ih]
          for(cc in 1:nchar(substr)){
            ch <- substr(mcs1$seq[ih], cc, cc)
            cs <- substr(mcs2$seq[is], cc, cc)
            if(ch == cs)
              str_sub(substr,cc,cc) <- "."
            else
              str_sub(substr,cc,cc) <- "*"
          }
          message(sprintf("                   %s",substr))
        }

        #message(sprintf("%03d %s: %03d %s",is,name2,mcs2$seqpos[is],mcs2$seq[is]))
        message(sprintf("%s: helixpos %03d %s",name2,mcs2$seqpos[is]-helixoffset,mcs2$seq[is]))
        #readline("difference found. Hit <return>")
        ndiffs = ndiffs + 1
        diff_pos[ndiffs]<-mcs1$seqpos[ih]
      }
    }

    #iterate to the next sequence position
    ih_new = ih + 1
    while(mcs1$seqpos[ih] == mcs1$seqpos[ih_new] && ih_new < nrow(mcs1))
      ih_new = ih_new + 1
    ih = ih_new

    is_new = is + 1
    while(mcs2$seqpos[is] == mcs2$seqpos[is_new] && is_new < nrow(mcs2))
      is_new = is_new + 1
    is = is_new

    if(verbose){
      message(sprintf("ih is now %d, is is now %d, nrow 1 is %d, nrow 2 is %d",ih,is,nrow(mcs1),nrow(mcs2) ))
    }
  }

  message(sprintf("Found %d differences",ndiffs))
  return(diff_pos)

}
