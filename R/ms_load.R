
#' Load mass spec data from a 2-column space-delimited text file
#'
#' @param species the name of the species whose peptides will be loaded.
#' @export
load.sample <- function(froot,name,spots){

  # TODO: These are example values:
  # froot  <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
  # spots  <- c("G7","G10","G13")
  # sample <- "C1"

  #TODO: we need error checking on this!
  s1 <- read.table(sprintf("%s%s.txt",froot,spots[1]))
  s2 <- read.table(sprintf("%s%s.txt",froot,spots[2]))
  s3 <- read.table(sprintf("%s%s.txt",froot,spots[3]))

  colnames(s1) <- c("mass","intensity")
  colnames(s2) <- c("mass","intensity")
  colnames(s3) <- c("mass","intensity")

  sampleobject <- list("name" = name, "spot" = spots, "s1" = s1, "s2" = s2, "s3" = s3)

  return(sampleobject)

}



#' Load peptides from the mammalian collagen sequences googlesheet
#'
#' @param species the name of the species whose peptides will be loaded.
#' @examples
#' hcs <- load.mcs("human")
load.mcs<-function(spp, fromgs = T){

  sheet <- NA
  if(fromgs){
    message("Fetching the sheet of mammalian collagen sequences")
    #TODO: make sure this is public, OR store it as data..
    mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
    sheet <- gs_read(mcs)
  }
  else{
    message("ERROR: can't currently load from anywhere but google sheets\nuse fromgs=T")
    return(NA)
  }

  spidx<-ts_index(sheet,spp)
  if(spidx<0){
    message(sprintf("ERROR: cannot find sequence for %s",spp))
    return(NA)
  }

  message("Calculating sequences now")
  endcol<-ncol(sheet)
  start<-4
  count<-1


  for(j in start:endcol){
    #when we hit a cut point:
    if(grepl("K|R",sheet[spidx,j])){

      count<-count+1
      message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
      end=j

      sequence <- paste0(sheet[spidx,start:end],collapse="")

      nhyd <- str_count(sequence,"P")
      nglut <- str_count(sequence,"Q") + str_count(sequence,"N")

      phydp <- get_hydroxylations(sheet,start,end)

      #create the basics of the row
      hdat <- data.frame(sequence,start-4,end-start,nhyd)

      #add the probabilities
      if(nrow(phydp) >0 ){
        pnh <- ba_nhyd(phydp$prob)
        #print(pnh)
        message(sprintf("nhyd\tprob"))
        message(sprintf("%d\t%0.4f\n",pnh$nhyd,pnh$prob))

        #TODO: find a better way do create out data for the nhyds file

        hdat <- data.frame(hdat,t(pnh$prob))

        #readline("Calculated Nhyd prob\nhit<return> to continue")
      }
      else{
        message("No hydroxylation probabilities for this sequence")
        #TODO: later on - check if there any P's with no probs available

        hdat <- data.frame(hdat,1)
      }



    }
    else{

    }

  }

}
