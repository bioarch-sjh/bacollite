

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
pplate <- function(path,pd,FUN=NULL,...){

  if(is.null(FUN)){
    message("No worker function provided to pplate")
    message("please supply a pplate compatible function name")
  }

  result <- list()
  for(rr in 1:nrow(pd)){
    fr <- sprintf("%s/%s",path,pd$froot[rr])

    #check file for spot 1 exists:
    fn1 = sprintf("%s%s.txt",fr,pd$spot1[rr])
    if(!file.exists(fn1)){
      message(sprintf("Entry %d, spot 1: Couldn't find file %s, skipping this sample",rr,fn1))
      next
    }
    else{
      #(sprintf("Sample 1 file is %s",fn1))
    }

    #check file for spot 2 exists:
    fn1 = sprintf("%s%s.txt",fr,pd$spot2[rr])
    if(!file.exists(fn1)){
      message(sprintf("Entry %d, spot 2: Couldn't find file %s, skipping this sample",rr,fn1))
      next
    }
    else{
      #message(sprintf("Sample 2 file is %s",fn1))
    }

    #check file for spot 3 exists:
    fn1 = sprintf("%s%s.txt",fr,pd$spot3[rr])
    if(!file.exists(fn1)){
      message(sprintf("Entry %d, spot 3: Couldn't find file %s, skipping this sample",rr,fn1))
      next
    }
    else{
      #message(sprintf("Sample 3 file is %s",fn1))
    }

    #spots <- c(pd$spot1[rr],pd$spot2[rr],pd$spot3[rr])
    #reps<- load.sample(fr,pd$sampleID[rr],spots)
    pdr <- pd[rr,]
    result[[rr]] <- FUN(path,pdr,...)
  }
  return(result)
}

