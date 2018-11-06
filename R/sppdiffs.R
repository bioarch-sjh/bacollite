
#' @name spp.diffs
#' get the differences in sequence between two species.
#' NB the sequences are sorted by sequence position *within* this function,
#' so there's no need to sort them before calling it.
#'
#' @param mcs1 the mammalian collagen sequence for species 1
#' @param mcs2 the mammalian collagen sequence for species 2
#' @param name1 (optional) the name of species 1
#' @param name2 (optional) the name of species 2
#' @param verbose (optional) whether to print messages while processing
#' @param halt (optional) whether to pause until <return> is hit while processing
#' @return Returns a vector of sequence positions for the position where the differences in the spectrum are
#' @export
#' @examples
#' sh <- load.mcs("sheep")
#' hu <- load.mcs()
#' diffs <- spp.diffs(sh,hu,"sheep","human")
spp.diffs <- function(mcs1, mcs2, name1="spp01", name2="spp02", verbose = F, halt = F){

    #SORT THE PEPTIDES BY SEQUENCE POSITION
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
    idx1=1
    idx2=1

    ndiffs = 0

    diff_pos <- NA

    while(idx1 < nrow(mcs1) && idx2 <nrow(mcs2)){
    #if testing, use:
    #for(ii in 1:100){#nrow(mcs1)){

      if(verbose){
        message(sprintf("in while loop, idx1 = %d, idx2 = %d, s1$seqpos = %d, s1$seqpos = %d",idx1,idx2,mcs1$seqpos[idx1],mcs2$seqpos[idx2]))
        if(halt)
          readline("hit <return> to continue")
      }

      #if helix positions match - or they aren't present in the other seq (deals with INDELS)
      if(mcs1$seqpos[idx1] == mcs2$seqpos[idx2] ||
         ( (!mcs1$seqpos[idx1] %in% mcs2$seqpos)
           &&
           (!mcs2$seqpos[idx2] %in% mcs1$seqpos)
         )

      ){

        #message("difference found")

        if(! (mcs1$seq[idx1] == mcs2$seq[idx2])){

          indel <- F

          #NB - should use Smith-Waterman alignment really!
          # see the 'Biostrings' package - but see if we can find something quicker
          if(nchar(mcs1$seq[idx1]) != nchar(mcs2$seq[idx2])){
            message("INDEL")
            indel <- T
          }
          else{
            if(verbose)message("SUBSITUTION")
          }

          #check the codes:
          validseq(mcs1$seq[idx1])
          validseq(mcs2$seq[idx2])

          if(verbose){
            message(sprintf("%s: helixpos %03d\n%s",name1,mcs1$seqpos[idx1]-helixoffset,mcs1$seq[idx1]))
          #if(halt)
          #  readline("press <return> to continue")
          }

          if(!indel){
            substr <- mcs1$seq[idx1]
            for(cc in 1:nchar(substr)){
              ch <- substr(mcs1$seq[idx1], cc, cc)
              cs <- substr(mcs2$seq[idx2], cc, cc)
              if(ch == cs)
                str_sub(substr,cc,cc) <- "."
              else
                str_sub(substr,cc,cc) <- "*"
            }
            #floor (log10 (abs (x))) + 1 gets the number of digits - we'll need these to do the alignment properly
            if(verbose)message(sprintf("%s",substr))
          }

          #message(sprintf("%03d %s: %03d %s",idx2,name2,mcs2$seqpos[idx2],mcs2$seq[idx2]))
          if(verbose)message(sprintf("%s:\n%s: helixpos %03d\n",mcs2$seq[idx2],name2,mcs2$seqpos[idx2]-helixoffset))
          #readline("difference found. Hit <return>")
          ndiffs = ndiffs + 1
          diff_pos[ndiffs]<-mcs1$seqpos[idx1]
        }
      }





      #iterate to the next sequence position
      idx1_new = idx1 + 1
      while(mcs1$seqpos[idx1] == mcs1$seqpos[idx1_new] && idx1_new < nrow(mcs1))
        idx1_new = idx1_new + 1
      idx1 = idx1_new

      idx2_new = idx2 + 1
      while(mcs2$seqpos[idx2] == mcs2$seqpos[idx2_new] && idx2_new < nrow(mcs2))
        idx2_new = idx2_new + 1
      idx2 = idx2_new

      if(verbose){
        message(sprintf("idx1 idx2 now %d, idx2 is now %d, nrow 1 idx2 %d, nrow 2 is %d",idx1,idx2,nrow(mcs1),nrow(mcs2) ))
      }
    }

    if(verbose)message(sprintf("Found %d differences",ndiffs))
    return(diff_pos)

}
