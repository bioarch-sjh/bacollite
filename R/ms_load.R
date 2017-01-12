
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





#' Get the species index in the mammalian_collagen_sequences sheet
#'
#' @param species the name of the species whose peptides will be loaded.
#' @export
#' @examples
#' hcd <- ts_index(sheet,"human")
ts_index <- function(sheet,spp){

  #TODO: There is surely a more efficient way of doing this...
  found <- F
  spidx <- 0
  for(i in 1:nrow(sheet)){
    if(grepl(spp,sheet[i,1],ignore.case=TRUE)){
      if(!found){
        message(sprintf("FOUND: index is %d, search term is \"%s\"",i,sheet[i,1]))
        found <- T
        spidx <- i
      }
      else{
        message(sprintf("Further match found at  index %d, search term is \"%s\" - this entry will be ignored",i,sheet[i,1]))
      }
    }
  }
  if(!found){
    message("Match not found, exiting")
    return (-1)
  }
  return (spidx)
}





#we know start and end, so we can calculate n_hyds in this range
#phydp <- get_hydroxylations(sheet,start,end)
get_hydroxylations <- function(sheet,start,end,dopause=T){


  hidx<-ts_index(sheet,"hydroxylation")
  if(hidx<0){
    readline("Couldn't find the hydroxylation row - check the sheet!\nhit <return> to continue")
  }

  #message("pos\tprob")
  #TODO: this is redundant if the method below works...
  #for(i in start:end){
  #	message(sprintf("%d\t%f",i,sheet[hidx,i]))
  #}

  #TODO: Check this is a better method:
  d <- as.numeric(sheet[hidx,start:end])
  i <- (start-4):(end-4)

  hoffset <- -16
  output<-data.frame(helixpos=i[which(d>0)]+hoffset,pos=i[which(d>0)],prob=d[which(d>0)])

  message(sprintf("start=%d; end=%d\n",start,end))
  if(nrow(output)>0)
    print(output)
  else
    message("No ",appendLF=F)

  message("Hydroxylation probs found")
  #\nhit <return> to continue")

  return(output)

}







#' Load peptides from the mammalian collagen sequences googlesheet
#'
#' @param species the name of the species whose peptides will be loaded.
#' @export
#' @examples
#' hcs <- load.mcs("human")
load.mcs<-function(spp, sheet=bioarch_mammal_sequences,verbose=F){

#  sheet <- NA
#  if(fromgs){
#    message("Fetching the sheet of mammalian collagen sequences")
#    #TODO: make sure this is public, OR store it as data..
#    #mcs <- gs_title("Mammalian Collagen Sequences v0.0.1")
#    mcs <- gs_title("Mammalian Collagen Sequences v0.0.2")
#    sheet <- gs_read(mcs)
#  }
#  else{
#    message("ERROR: can't currently load from anywhere but google sheets\nuse fromgs=T")
#    return(NA)
#  }

  spidx<-ts_index(sheet,spp)
  if(spidx<0){
    message(sprintf("ERROR: cannot find sequence for %s",spp))
    return(NA)
  }

  message("Calculating sequences now")
  endcol<-ncol(sheet)
  start<-4
  count<-1

  sdata <- data.frame(seq=as.character(),nhyd=as.integer(),nglut=as.integer(),mass1=as.numeric(),prob1=as.numeric)


  for(j in start:endcol){
    #when we hit a cut point:
    if(grepl("K|R",sheet[spidx,j])){



      count<-count+1
      message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
      end=j

      sequence <- paste0(sheet[spidx,start:end],collapse="")

      #readline(sprintf("Found sequence %s\nhit <return> to process",sequence))

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

        if(verbose){
          message("Calculated Nhyd prob")
        }
      }
      else{
        message("No hydroxylation probabilities for this sequence")
        #TODO: later on - check if there any P's with no probs available

        hdat <- data.frame(hdat,1)
      }


      #readline("hit <return> to continue...\n")
      start = j+1

    }
    else{

    }
  }
  return(hdat)
}



