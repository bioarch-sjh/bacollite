#' Get all spectrum data between a mass range
#'
#' @param ms the mass spectrum
#' @param lbl the lower bound limit
#' @param ubl the upper bound limit
#' @return the subrage
#' @example
#' sub <- ms_subrange(spectrum,1589.5,1595.5)
#' @export
ms_subrange <- function(ms,lbl,ubl){
  subms <- ms[
    ms[,1] <= ubl & ms[,1] > lbl
    ,]
  #TODO: find out why we *didn't* need the following and what the effects
  return(subms)
}
