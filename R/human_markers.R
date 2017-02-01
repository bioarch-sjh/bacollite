


#' Load a set of peptide markers for human
#'
#' @param fn full path of the raw gpm text file
#' @export
#' @examples
#' hm<-load.human.markers()
load.human.markers <- function(fn=NA){

  if(is.na(fn)){
    fn<-system.file("extdata", "Human_Peptide_Markers.dat", package = "bacollite")
  }
  data <- read.table(fn,sep = ",",header=T)


  #convert the seq column to string format:
  data$seq <- sapply(data$seq, as.character)

  return(data)
}
