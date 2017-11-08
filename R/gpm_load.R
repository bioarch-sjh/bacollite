


#TODO: get the info from google sheets....
#' Load a set of peptides from the Global Protein Machine
#'
#' @param fn full path of the raw gpm text file
#' @export
load.gpm.raw <- function(fn=NA,colno=NA){


  gpm <- read.table(fn,
                    fill=T,header=T,sep="\t",comment.char="")

  gpm2 <- data.frame(seq=as.character(gpm$sequence),
                     nhyd = stringr::str_count(gpm$modifications,"P"),
                     nglut = (stringr::str_count(gpm$modifications,"Q")+stringr::str_count(gpm$modifications,"N")),
                     mass1=0, prob=0 )

  gpm2 <- unique(gpm2)

  for(i in 1:nrow(gpm2)){

    #get the mass data...
    #cd1 <- ms_tpeaks(gpm2$seq[i])
    #gpm2$mass1[i] = cd1$mass[1] +  (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)

    cd1<- ms_iso(gpm2$seq[i],ndeamidations=gpm2$nglut[i],nhydroxylations = gpm2$nhyd[i])

    #TODO: Need a better way to catch bad ndeamidations and nhydroxylations...
    #TODO: The whole prob field needs rethinking for load.gpm()
    if(is.data.frame(cd1)){
      gpm2$prob[i] = cd1$prob[1]
      gpm2$mass1[i] = cd1$mass[1]
    }
    else{
      gpm2$prob[i] = 0
      gpm2$mass1[i] = 0
    }



  }

  #remove outliers
  gpm2 <- gpm2[which(gpm2$mass1>795 & gpm2$mass1 <3500),]

  #Sort by mass
  gpm2 <- gpm2[order(gpm2$mass1),]

  #convert the seq column to string format:
  gpm2$seq <- sapply(gpm2$seq, as.character)

  #add the collagen number
  if(!is.na(colno)){
    gpm2$col<-colno
  }

  return(gpm2)

}




#TODO: get the info from google sheets....including the observation counts
#' Load a set of peptides from the Global Protein Machine
#'
#' @param col1fn full path of the raw gpm text file for collagen 1
#' @param col2fn full path of the raw gpm text file for collagen 1
#' @export
#' @examples
#' gpm<-load.gpm()
load.gpm <- function(col1fn=NA,col2fn=NA){

  if(is.na(col1fn))
    col1fn<-system.file("extdata", "gpm_human_collagen_COL1A1.dat", package = "bacollite")
  col1seqs<-load.gpm.raw(col1fn,colno=1)

  if(is.na(col2fn))
    col2fn<-system.file("extdata", "gpm_human_collagen_COL1A2.dat", package = "bacollite")
  col2seqs<-load.gpm.raw(col2fn,colno=2)

  gpm <- rbind(col1seqs,col2seqs)

  #Sort by mass
  gpm <- gpm[order(gpm$mass1),]

  return(gpm)

}

