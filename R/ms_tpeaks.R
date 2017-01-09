
#' Wrapper for the 'R_iso_seq' q2e c function, to get theoretical peaks from a sequence
#' @param seqeunce the amino acid sequence
#' @keywords test
#' @export
#' @examples
#' ms_tpeaks("IGQPGAVGPAGIR")
ms_tpeaks <- function(sequence){

  message(sprintf("Calculating isotope distributions for peptide %s",sequence))

  #CHECK THE STRING VERY CAREFULLY - THE C CODE IS FRAGILE!
  if(grepl('^[A-Z]+$', sequence)){

    failedflag <- 0
    #result <- .C("R_iso_seq", infile=as.character(sequence),mass=as.double(1:5),prob=as.double(1:5),errflag=as.integer(failedflag))

    R_iso_seq(rseq=as.character(sequence), xmass=as.double(1:5), xprob=as.double(1:5), failed=as.integer(failedflag))

    message(xmass)


    #TODO: Check the failed flag

    #return (result)

  }
  else{
    message("Input string must contain only uppercase alphabeticial characters")
    return(1)

  }
}
