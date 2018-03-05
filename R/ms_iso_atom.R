
#' Obtain the isotope distribution for a molecule consisting of atoms of C,H,N,O and S
#' post-translation modifications (deamitation and hydroxylation)
#' @param atoms vector containing the counts of atoms in the order C,H,N,O and S
#' @keywords isotopes mass sequence deamidation hyroxylation
#' @export
#' @examples
#' ms_iso_atom(atoms)
#' ms_iso_atom(c(47,69,15,12,0))
ms_iso_atom <- function(atoms,verbose=F){

  if(length(atoms)< 5){
    message("ms_iso_atom error: Vector with five entries needed\nreturning NA")
    return (NA)
    }
  if(is.numeric(atoms) == F){
    message("ms_iso_atom error: Vector with five integers needed\nreturning NA")
    return (NA)
    }
    
  
  
  if(verbose){
    message(sprintf("Calculating isotope distributions"))
  }

  result <- cppIsoAtom_p(atoms)


  return (result)

}
