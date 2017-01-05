
ms_subrange <- function(ms,lbl,ubl){
  subms <- ms[
    ms[,1] <= ubl & ms[,1] > lbl
    ,]
  #TODO: find out why we *didn't* need the following and what the effects
  return(subms)
}
