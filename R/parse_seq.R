

#' Parse an amino acid sequence into a set of peptides, incorporating hydroxylation and deamidation modifications, and calculating mass
#'
#' @param seqeunce a character string holding the sequence, e.g. "GPPGAPGPPGPP"
#' @param cuts a regular expression listing the amino acids that the sequence should be cut at. Defaults to Trypsing cut sites K & R
#' @param skip any codes that mean the peptide is invalid. Defaults to "X"
#' @param massmin the minimum mass that a peptide should be. Defaults to 800
#' @param massmax the maximum mass that a peptide should be. Defaults to 3500
#' @param verbose whether to print messages to console during processing. Defaults to FALSE
#' @param max.missed.cleaves the number of missed cleaves to incorporate in the set of peptides.
#' @export
#' @examples
#' parse.seq("GPPGQKGPPGPQGPRGPPGPPGPM")
parse.seq <- function(sequence,cuts="K|R",skip="X", massmin = 800, massmax = 3500,verbose = F,max.missed.cleaves=0, cutbefore = F){

  #initialise some variables
  peptides <- NA
  pepidx <- 1
  len <- str_length(sequence)
  start <-1
  nextpos <- start
  num.mc <- 0
  tooheavy <- F
  badchar <- F
  attheend<-F
  pos <-1
  cut.offset <- 0

  #set the cut offset
  if(cutbefore){
    cut.offset <- 1
    if(verbose)message("Cutting before!")
  }

  #for(pos in 1:len){
  while(!attheend){

    #get the current amino acid
    aa <- str_sub(sequence, pos+cut.offset, pos+cut.offset)

    # if we are at a cut point, we can process the peptide, but note that
    # we have to keep track of the start position for the *next* peptide
    # which can get complicated as we deal with missed cleaves etc.
    if (str_detect(aa, cuts) || pos == len) {

      #record the position of the next start sequence
      nextpos[num.mc + 1] = pos

      #increment the number of missed cleaves
      num.mc <- num.mc + 1

      #get the current peptide sequence.
      ss <- str_sub(sequence,start,pos)

      #if there are any characters that should be skipped
      if(str_detect(ss,skip)){

        #do nothing
        if(verbose)
          message(sprintf("  can't use %s",ss))
        badchar <- T
      }

      #if there are no characters that should be skipped
      else{

        if(verbose)message(sprintf("found sequence %s",ss))

        #calculate the masses without PTMs
        masses <- ms_iso(ss)
        tooheavy <- F

        #if the masses are in range
        if(masses$mass[1]>massmin && masses$mass[1] < massmax){

          #get the number of possible hydroxylations
          nhyd <- str_count(ss,"P")

          #get the number of possible deamidations
          nglut <- str_count(ss,"Q") + str_count(ss,"N")

          #Go through each hydroxylation / deamidation combination
          for(hh in 0:nhyd){
            for(dd in 0:nglut){

              #calculate the mass for this PTM
              masses <- ms_iso(ss,ndeamidations=dd,nhydroxylations=hh)

              #create the new entry in the peptide list
              newrow <- data.frame(
                seq=as.character(ss)
                ,nhyd=hh
                ,nglut=dd
                #,mass1 = masses$mass[1] + (dd*0.984015)+(hh*16)
                ,mass1 = masses$mass[1]
                ,seqpos = start
                ,missed.cleaves = num.mc -1 # we have to subtract 1 because we've incremented num.mc already
              )

              #add the row to the list
              if(pepidx==1){
                peptides <- newrow
              }
              else{
                peptides <- rbind(peptides,newrow)
              }
              pepidx <- pepidx+1


            }
          }

        }
        else{
          if(verbose)message(sprintf("  mass %0.2f out of range for sequence %s",masses$mass[1],ss))
          if(masses$mass[1] > massmax)
            tooheavy <- T
        }
      }
      #if we are at the limit of mass or num.missed.cleaves, move the start pos on
      # or if too heavy
      # or if there's a bad character
      if(tooheavy || badchar || num.mc > max.missed.cleaves){
        start <- nextpos[1]+1
        tooheavy <- F
        badchar <- F
        #rewind the position if we've missed a cleave
        if(num.mc > 0)
          pos <- nextpos[1]#this will be incremented below!

        if(verbose)
          message(sprintf("nmc = %0.0f, start =%0.0f, pos = %0.0f\ncleave positions are %0.0f %0.0f %0.0f",num.mc,start,pos,nextpos[1],nextpos[2],nextpos[3]))

        num.mc <- 0
      }
    }

    pos <- pos + 1

    if(pos>len)
      break

  }

  #TODO: demistify why we have to use as.character for sequences...
  peptides$seq <- as.character(peptides$seq)
  return(peptides)

}
