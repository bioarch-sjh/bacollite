





#' This checks if the letters in a sequence are all valid amino acid codes
#'
#' @param seq character string of an amino acid sequence
#' @param halt (optional) whether to pause until <return> is hit while processing
#' @return a logical value: T if valid, F if not
#' @export
validseq <- function(seq, halt = F){
  valid <- T
  for(ii in 1:nchar(seq)){
    aa <- substr(seq, ii, ii)
    # codes from here:
    # http://www.biochem.ucl.ac.uk/bsm/dbbrowser/c32/aacode.html
    if(!str_detect("GPAVLIMCFYWHKRQNEDST",aa)){
      valid <- F
      if(halt)
        readline(sprintf("Found %s at position %d - not an amino acid code! hit <return>",aa,ii))
      else
        message(sprintf("Found %s at position %d - not an amino acid code! hit <return>",aa,ii))
    }
  }
  return(valid)
}
