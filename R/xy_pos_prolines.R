
#  TYPICAL USAGE
#    xy.pos.prolines("IGQPGAVGPAGIR")


#' count the prolines in the x and y position of a collagen peptide
#' @param seq the amino acid sequence
#' @param verbose default is False
#' @keywords isotopes mass sequence deamidation hyroxylation
#' @export
xy.ypos.prolines <- function(seq,verbose = F){

  gcount <- vector(length = 3)
  pcount <- vector(length = 3)

  if(verbose)
    message(sprintf("Sequence is %d long",length(seq)))

  for(ii in 1:nchar(seq)){

    pos <- (ii%%3)+1

    if(substr(seq,ii,ii) == "G")
      gcount[pos] = gcount[pos] + 1#str_count(seq[ii],"G")

    if(substr(seq,ii,ii) == "P")
      pcount[pos] = pcount[pos] + 1#str_count(seq[ii],"P")

    #if(verbose){
    #  message(sprintf("ii = %d, pos %d, letter = %s, gcount = %d %d %d, pcount = %d %d %d",ii,pos,
    #                  substr(seq,ii,ii),
    #                  gcount[1],gcount[2],gcount[3],
    #                  pcount[1],pcount[2],pcount[3]))
    #}
  }

  gmax = which.max(gcount)
  ypos = (gmax+1)%%3+1

  if(verbose){
    message(sprintf("sequence: %s",seq))
    message(sprintf("G count positions: %d %d %d, max %d",gcount[1],gcount[2],gcount[3],gmax))
    message(sprintf("P count positions: %d %d %d, yct %d",pcount[1],pcount[2],pcount[3],ypos))
  }

  return(pcount[ypos])
}


