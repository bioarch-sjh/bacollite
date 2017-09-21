
#' Wrapper for the 'R_iso_seq' q2e c function, to get theoretical peaks from a sequence
#' @param seqeunce the amino acid sequence
#' @keywords test
#' @export
#' @examples
#' ms_tpeaks("IGQPGAVGPAGIR")
ms_tpeaks <- function(sequence,verbose=F,ndeamidations=0,nhydroxylations=0){

  if(verbose){
    message(sprintf("Calculating isotope distributions for the peptide %s",sequence))
  }

  #CHECK THE STRING VERY CAREFULLY - THE C CODE IS FRAGILE!
  if(grepl('^[A-Z]+$', sequence)){

    result<-cppIso(sequence)

    #TODO: Check for a failed flag in result


    if(ndeamidations>0){
      max_deam <- stringr::str_count(sequence,"Q") + str_count(sequence,"N")
      if(ndeamidations>max_deam){
        message(sprintf("ERROR: %d deamidations are not possible for sequence %s",
                        ndeamidations,sequence))
        message(sprintf("       max number of deamidations is %d (count of 'Q' and 'N' in sequence)",max_deam))
        return (NA)
      }
      else{
        result$mass <- result$mass + (ndeamidations * 0.984015)
      }
    }

    if(nhydroxylations>0){
      max_hyd <- stringr::str_count(sequence,"P")
      if(nhydroxylations>max_hyd){
        message(sprintf("ERROR: %d hydroxylations are not possible for sequence %s",
                        nhydroxylations,sequence))
        message(sprintf("       max number of hydroxylations is %d (count of 'P' in sequence)",max_hyd))

        return (NA)
      }
      else{
        result$mass <- result$mass + (nhydroxylations * 16)
      }
    }

    return (result)

  }
  else{
    message(sprintf("ms_tpeaks: error processing string %s\n Input string must contain only uppercase alphabeticial characters",sequence))
    return(NA)
  }
}
