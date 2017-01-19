
#' Load mass spec data from a 2-column space-delimited text file
#'
#' @param species the name of the species whose peptides will be loaded.
#' @export
load.sample <- function(froot,name="Sample",spots){

  # TODO: These are example values:
  # froot  <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
  # spots  <- c("G7","G10","G13")
  # sample <- "C1"

  #check that froot exists
  if(is.na(froot)){
    message("Froot is not defined")
    return(NA)
  }


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
ts_index <- function(sheet,spp,verbose=F){

  #TODO: There is surely a more efficient way of doing this...
  found <- F
  spidx <- 0
  for(i in 1:nrow(sheet)){
    if(grepl(spp,sheet[i,1],ignore.case=TRUE)){
      if(!found){
        if(verbose){
          message(sprintf("FOUND: index is %d, search term is \"%s\"",i,sheet[i,1]))
        }
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
get_hydroxylations <- function(sheet,start,end,dopause=T,verbose=F){


  hidx<-ts_index(sheet,"hydroxylation")
  if(hidx<0){
    readline("Couldn't find the hydroxylation row - check the sheet!\nhit <return> to continue")
  }

  #TODO: Check this is a better method:
  #suppressWarnings because we don't mind that empty cells become NAs by coercion
  d <- suppressWarnings(as.numeric(sheet[hidx,start:end]))
  i <- (start-4):(end-4)

  hoffset <- -16
  output<-data.frame(helixpos=i[which(d>0)]+hoffset,pos=i[which(d>0)],prob=d[which(d>0)])

  if(verbose){
    message(sprintf("start=%d; end=%d\n",start,end))
    if(nrow(output)>0)
      print(output)
    else
      message("No ",appendLF=F)

    message("Hydroxylation probs found")
    #\nhit <return> to continue")
  }

  return(output)

}







#' Load peptides from the mammalian collagen sequences googlesheet
#'
#' @param species the name of the species whose peptides will be loaded.
#' @export
#' @examples
#' hcs <- load.mcs()
#' hcs <- load.mcs("goat")
load.mcs<-function(spp="human", sheet=bioarch_mammal_sequences,massmin=800,massmax=3500,verbose=F){

  spidx<-ts_index(sheet,spp)
  if(spidx<0){
    message(sprintf("ERROR: cannot find sequence for %s",spp))
    return(NA)
  }

  message(sprintf("Calculating sequences for %s now",spp))
  endcol<-ncol(sheet)
  start<-4
  count<-1

  sdata <- data.frame(seq=as.character(),nhyd=as.integer(),nglut=as.integer(),mass1=as.numeric(),prob=as.numeric())


  for(j in start:endcol){

    #when we hit a cut point, we can process the new peptide:
    if(grepl("K|R",sheet[spidx,j])){

      count<-count+1
      if(verbose){
        message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
      }
      end=j

      sequence <- paste0(sheet[spidx,start:end],collapse="")


      masses <- ms_tpeaks(sequence)

      if(verbose){
        message("Masses:")
        print(masses)
      }
      if(masses$mass[1]>massmin && masses$mass[1] < massmax){
        #readline(sprintf("Found sequence %s\nhit <return> to process",sequence))

        nhyd <- str_count(sequence,"P")
        nglut <- str_count(sequence,"Q") + str_count(sequence,"N")

        phydp <- get_hydroxylations(sheet,start,end)

        #create the basics of the row
        #hdat <- data.frame(sequence,start-4,end-start,nhyd)

        #add the probabilities
        foundhprobs<-F
        if(nrow(phydp) >0 ){
          pnh <- ba_nhyd(phydp$prob)
          foundhprobs<-T
          if(verbose){
            message("Calculated Nhyd prob:")
            print(pnh)
          }
        }
        else{
          if(verbose){
            message("No hydroxylation probabilities for this sequence")
            #TODO: later on - check if there any P's with no probs available
          }
          #hdat <- data.frame(hdat,1)
        }



        #New strategy is to have a row for each deam/hyd level - matches gpm dataframe.
        for(d in 0:nglut){
          if(foundhprobs){
            for(h in 1:nrow(pnh)){
              if(verbose){
                message(sprintf("Calculating hyd level %d of %d",h,nrow(pnh)))
              }
              #sdata <- data.frame(
              #           seq=as.character(),
              #           nhyd=as.integer(),
              #           nglut=as.integer(),
              #           mass1=as.numeric(),
              #           prob=as.numeric())

              newrow <- data.frame(
                           seq=as.character(sequence)
                          ,nhyd=pnh$nhyd[h]
                          ,nglut=d
                          ,mass1 = masses$mass[1] + (d*0.984015)+(pnh$nhyd[h]*16)
                          ,prob =  pnh$prob[h]
                        )
              sdata <- rbind(sdata,newrow)

              if(verbose){
                message(sprintf("newrow has %d rows and looks like:",nrow(sdata)))
                print(newrow)
                message(sprintf("Sdata now has %d rows and looks like:",nrow(sdata)))
                print(sdata[nrow(sdata),])
                answer<-readline("is hyd working?")
                if (substr(answer, 1, 1) == "n")
                  return (sdata)
              }
            }
          }
          else{#There are no hydoxylation probabilities

            #TODO: This is almost a direct copy of the above!
            newrow <- data.frame(
              seq=as.character(sequence)
              ,nhyd=nhyd
              ,nglut=d
              ,mass1 = masses$mass[1] + (d*0.984015)+(nhyd*16)
              ,prob =  1 #TODO: check this is always the case!
            )
            sdata <- rbind(sdata,newrow)

            if(verbose){
              message(sprintf("newrow has %d rows and looks like:",nrow(sdata)))
              print(newrow)

              message(sprintf("Sdata now has %d rows and the last row looks like:",nrow(sdata)))
              print(sdata[nrow(sdata),])

              answer<-readline("is nohyd working?")
              if (substr(answer, 1, 1) == "n")
                return (sdata)
            }
          }

        }

      }
      #readline("hit <return> to continue...\n")
      start = j+1

    }
    else{

    }
  }

  #Finally, let's sort the data by mass
  sdata<-sdata[order(sdata$mass1),]

  sdata$seq <- sapply(sdata$seq, as.character)

  return(sdata)
}



