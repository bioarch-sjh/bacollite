

#' Process a MALDI plate of replicates
#'
#' @param ts theoretical spectrum peaks for a particular peptide
#' @param pd dataframe containing the following fields:
#'   \item{froot}{The common part of the file name, e.g. "20211017_Filey_" as part of "20211017_Filey_A1.txt"}
#'   \item{sampleID}{A name for the sample, e.g. "A420"}
#'   \item{spot1}{The location of spot1, e.g. "A1"}
#'   \item{spot2}{The location of spot2, e.g. "A4"}
#'   \item{spot3}{The location of spot3, e.g. "A7"}
#' @param FUN the processing function, e.g. "averageMSD"
#' @param path the relative path to the data from where the function is called
#' @export
pplate <- function(pd,FUN=NULL,path="./",...){

  result <- list()
  for(rr in 1:nrow(pd)){
    fr <- sprintf("%s/%s",path,pd$froot[rr])
    spots <- c(pd$spot1[rr],pd$spot2[rr],pd$spot3[rr])
    reps<- load.sample(fr,pd$sampleID[rr],spots)
    result[[rr]] <- FUN(reps)
  }
  return(result)
}

