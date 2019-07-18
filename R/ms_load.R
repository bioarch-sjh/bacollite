
#' Load mass spec data from a 2-column space-delimited text file
#'
#' @param froot the file path to the directory containing the data, plus the part of the filename that stays the same
#' @param name the name of the sample
#' @param spots a vector of spot names - there *must* be three of these
#' @param fext the file extension. Default is ".txt"
#' @export
#' @examples
#' froot  <- "~/tmp/bioarch_keri/20160909_Keri13/20160909_Keri13_0_"
#' spots  <- c("G7","G10","G13")
#' name <- "C1"
load.sample <- function(froot,name="Sample",spots,fext=".txt"){

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
  s1 <- read.table(sprintf("%s%s%s",froot,spots[1],fext))
  s2 <- read.table(sprintf("%s%s%s",froot,spots[2],fext))
  s3 <- read.table(sprintf("%s%s%s",froot,spots[3],fext))

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
get_hydroxylations <- function(sheet,start,end,spp="human",dopause=T,verbose=F){


  hidx<-ts_index(sheet,"hydroxylation")
  if(hidx<0){
    readline("Couldn't find the hydroxylation row - check the sheet!\nhit <return> to continue")
  }


  #suppressWarnings because we don't mind that empty cells become NAs by coercion
  d <- suppressWarnings(as.numeric(sheet[hidx,start:end]))

  # Need to remove probabilities where the seqeunce letter isn't 'P'
  # This can happen because the P may not be present in the species under consideration
  spidx<-ts_index(sheet,spp)
  seq <- sheet[spidx,start:end]
  isnt_p <- which(seq!="P")
  d[isnt_p] <- 0

  #handle indels
  #TODO: has to be a more efficient way of doing this
  if(anyNA(seq)){
    for(i in 1:length(seq)){
      if(is.na(seq[i])) d[i] = 0
    }
  }

  #create a helixpos index
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


#' Load sequences from the mammalian collagen sequences googlesheet
#'
#' @param species the name of the species whose peptides will be loaded.
#' @param sheet the spreadsheet to load the data from. Must be the same format as bioarch_mammal_sequences
#' @param verbose verbose processing with more detail. Useful for debugging
#' @param col3 whether to load the sequence from column3, or from the remainder of the spreadsheet. For debugging.
#' @export
#' @examples
#' hcs <- load.sequence()
#' hcs <- load.sequence("goat")
load.sequence<-function(spp="human",  sheet=bioarch_mammal_sequences, verbose = F, col3=F){

  spidx<-ts_index(sheet,spp)
  if(spidx<0){
    message(sprintf("ERROR: cannot find sequence for %s",spp))
    return(NA)
  }

  if(col3){
    sequence <- sheet[spidx,3]
  }
  else{
    startcol <-4

    message("Reading sequence from mcs data columns")

    endcol<-ncol(sheet)

    seqraw <- as.character(sheet[spidx,startcol:endcol])
    seqraw <- seqraw[!is.na(seqraw)]
    sequence <- paste0(seqraw,collapse="")
  }

  return(sequence)

}



#' Load peptides from the mammalian collagen sequences googlesheet and use the proline hydroxylation
#' probabilites contained in the sheet to estimate the probability of each number of hydroxyprolines in each peptide
#'
#' @param spp the name of the species whose peptides will be loaded.
#' @param sheet the spreadsheet to load the data from. Must be the same format as bioarch_mammal_sequences, which is the latest version of the Bioarch mammal sequence dataset
#' @param massmin the minimum mass for a sequence to be returned. Defaults to 800
#' @param massmax the maximum mass for a sequence to be returned. Defaults to 3500
#' @param verbose verbose processing with more detail. Useful for debugging
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

  if(verbose){
    message(sprintf("Calculating sequences for %s now",spp))
  }

  endcol<-ncol(sheet)
  shoff<-4 #sheet-based offset (column that the sequence starts in)

  start<-shoff
  count<-1

  sdata <- data.frame(seq=as.character(),nhyd=as.integer(),nglut=as.integer(),mass1=as.numeric(),prob=as.numeric(),seqpos=as.integer())


  runtoend <- F

  for(j in start:endcol){

    #when we hit a cut point, we can process the new peptide:
    if(grepl("K|R",sheet[spidx,j])){

      count<-count+1
      if(verbose){
        message(sprintf("%s\n%d:\t",sheet[spidx,j],count),appendLF=F)
      }
      end=j


      #TODO: <NA> is coerced into "NA" - not what we want! - trying this:
      seqraw <- as.character(sheet[spidx,start:end])
      seqraw <- seqraw[!is.na(seqraw)]
      sequence <- paste0(seqraw,collapse="")


      #TODO: Now we are calculating mass by *atom count*, we need to do the calculation afresh for each hyd/deam combination
      #      in the for loop below - this will take a little longer but will be more precise.
      #masses <- ms_tpeaks(sequence)
      masses <- ms_iso(sequence)

      if(verbose){
        message("Masses:")
        print(masses)
      }
      if(masses$mass[1]>massmin && masses$mass[1] < massmax){
        #readline(sprintf("Found sequence %s\nhit <return> to process",sequence))

        nhyd <- str_count(sequence,"P")
        nglut <- str_count(sequence,"Q") + str_count(sequence,"N")

        phydp <- get_hydroxylations(sheet,start,end,spp)

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
        else{#No "P" in the sequence
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


              masses <- ms_iso(sequence,ndeamidations=d,nhydroxylations=h-1)

              newrow <- data.frame(
                           seq=as.character(sequence)
                          ,nhyd=pnh$nhyd[h]
                          ,nglut=d
                                                   # no need to do this now we are using ms_iso
                          ,mass1 = masses$mass[1]  # + (d*0.984015)+(pnh$nhyd[h]*16)
                          ,prob =  pnh$prob[h]
                          ,seqpos = start - shoff + 1
                        )
              sdata <- rbind(sdata,newrow)

              if(verbose){
                message(sprintf("newrow has %d rows and looks like:",nrow(sdata)))
                print(newrow)
                message(sprintf("Sdata now has %d rows and looks like:",nrow(sdata)))
                print(sdata[nrow(sdata),])
                if(!runtoend){
                  answer<-readline("is hyd working? (y/n/r = yes/no/run to end)")
                  if (substr(answer, 1, 1) == "n")
                    return (sdata)
                  if (substr(answer, 1, 1) == "r")
                    runtoend = T
                }
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
              ,seqpos = start - shoff + 1
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


  #Generate the collagen labelling based on sequence position
  #
  # currently anything with seqpos>1040 is collagen 2 - that's the number of peptides from the "QLSYGY"
  # motif at the beginning of the per-peptide cells in the spreadsheet (column D)
  sdata$col <-1
  sdata$col[sdata$seqpos>1040] <- 2


  #Finally, let's sort the data by mass
  sdata<-sdata[order(sdata$mass1),]

  sdata$seq <- sapply(sdata$seq, as.character)

  return(sdata)
}



