
#' Calculate the atoms in an amino acid sequence
#' post-translation modifications (deamidation and hydroxylation)
#' @param seqeunce the amino acid sequence
#' @param ndeamidations the number of deamidations
#' @param nhydroxylations the number of hydroxylations
#' @keywords isotopes mass sequence deamidation hyroxylation
#' @return integer vector of atoms in the order C H N O S
#' @export
#' @examples
#' seq_to_atoms("IGQPGAVGPAGIR")
seq_to_atoms <- function(sequence,ndeamidations,nhydroxylations,verbose=F){


  if(verbose){
    message(sprintf("Calculating isotope distributions for the peptide %s",sequence))
  }

  #CHECK THE STRING VERY CAREFULLY - THE C CODE IS FRAGILE!
  if(grepl('^[A-Z]+$', sequence)){

    #result<-cppIso(sequence)
    atoms <- aa_seq_to_atoms(sequence)


    #REMEMBER: Order of atoms is C H N O S
    #                            1 2 3 4 5
    # From Kristine Korzow Richter:
    # hydroxylation or oxidation - gain one oxygen (O)
    # deamidation - gain one oxygen (O), loose one nitrogen (N) and one hydrogen (H)




    #Add the deamidations
    if(ndeamidations>0){
      max_deam <- stringr::str_count(sequence,"Q") + str_count(sequence,"N")
      if(ndeamidations>max_deam){
        message(sprintf("ERROR: %d deamidations are not possible for sequence %s",
                        ndeamidations,sequence))
        message(sprintf("       max number of deamidations is %d (count of 'Q' and 'N' in sequence)",max_deam))
        return (NA)
      }
      else{
        #gain an oxygen
        atoms[4] <- atoms[4] + ndeamidations
        #lose a nitrogen
        atoms[3] <- atoms[3] - ndeamidations
        #lose a hydrogen
        atoms[2] <- atoms[2] - ndeamidations

        if(atoms[3]<0 || atoms[2]<0)
          message(sprintf("ERROR - can't have a negative number of atoms! deamidation at level %d can't happen\n nNitrogen = %d, nHydrogen =%d",
                          ndeamidations,
                          atoms[3],
                          atoms[2]))
      }
    }

    if(nhydroxylations>0){
      max_hyd <- stringr::str_count(sequence,"P")
      if(nhydroxylations>max_hyd){
        message(sprintf("ERROR: %d hydroxylations are not possible for sequence %s",
                        nhydroxylations,sequence))
        message(sprintf("       max number of hydroxylations is %d (count of 'P' in sequence)",max_hyd))

        return (NA)
      }
      else{
        #gain an oxygen
        atoms[4] <- atoms[4] + nhydroxylations
      }
    }

    #Now we've done the mods, we can calculate the isotopes based on the atoms

    return (result)

  }
  else{
    message(sprintf("ms_iso: error processing string %s\n Input string must contain ONLY uppercase alphabeticial characters",sequence))
    return(NA)
  }
}
