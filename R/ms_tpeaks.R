
#' Wrapper for the 'R_iso_seq' q2e c function, to get theoretical peaks from a sequence
#' @param seqeunce the amino acid sequence
#' @keywords test
#' @export
#' @examples
#' ms_tpeaks("IGQPGAVGPAGIR")
ms_tpeaks <- function(sequence,verbose=F){

  if(verbose){
    message(sprintf("Calculating isotope distributions for peptide %s",sequence))
  }

  #CHECK THE STRING VERY CAREFULLY - THE C CODE IS FRAGILE!
  if(grepl('^[A-Z]+$', sequence)){

    result<-cppIso(sequence)

    #TODO: Check for a failed flag in result

    return (result)

  }
  else{
    message("ms_tpeaks: Input string must contain only uppercase alphabeticial characters")
    return(NA)
  }
}
