

#' Take the output of ms_fit and convert to hits for a range of correlation thresholds
#'
#' @param fn_hits the dataframe returned by ms_fit
#' @param fn_laglim the max permissable lag
#' @param fn_corlim the set of correlation thresholds
#' @return a list of dataframes, one for each replicate, each containing the following fields
#'   \item{cl}{The correlation limit}
#'   \item{nh}{The number of hits above the limit}
#'   \item{sc}{The sum of ion counts for hits}
#' @export
corlim_data <- function (fn_hits,fn_laglim,fn_corlim = seq(0,1,0.05)){

  ###########################
  ## Function arguments:
  #fn_hits <- s_hits2[[cc]]
  #fn_corlim <- seq(0,1,0.05)
  #fn_laglim <- laglim
  #
  ###########################
  # This should be  function:
  fn_result <- data.frame(cl = fn_corlim, nh = NA, sc = NA)


  for(ii in 1:length(fn_corlim)){

    lim <- fn_corlim[ii]

    #remove ions from matches that don't pass the limits
    fn_hits$ion1[(fn_hits$cor1<lim | abs(fn_hits$lag1)>fn_laglim)] <- 0.
    fn_hits$ion2[(fn_hits$cor2<lim | abs(fn_hits$lag2)>fn_laglim)] <- 0.
    fn_hits$ion3[(fn_hits$cor3<lim | abs(fn_hits$lag3)>fn_laglim)] <- 0.
    #message(fn_hits$ion1)

    #re-parse the hit data:
    fn_hits$hit1 <- T
    fn_hits$hit2 <- T
    fn_hits$hit3 <- T

    #remove hits that don't pass the limits:
    fn_hits$hit1[(fn_hits$cor1<lim | abs(fn_hits$lag1)>fn_laglim)] <- F
    fn_hits$hit2[(fn_hits$cor2<lim | abs(fn_hits$lag2)>fn_laglim)] <- F
    fn_hits$hit3[(fn_hits$cor3<lim | abs(fn_hits$lag3)>fn_laglim)] <- F
    #message(fn_hits$hit1)

    #Group the hits into a vector:
    rep_hits <- c(fn_hits$hit1, fn_hits$hit2, fn_hits$hit3)

    #Add the count of correct matches:
    fn_numhits <- length(rep_hits[rep_hits == T])
    fn_score = sum(fn_hits$ion1)+sum(fn_hits$ion2)+sum(fn_hits$ion3)
    #return

    fn_result$nh[ii] <- fn_numhits
    fn_result$sc[ii] <- fn_score

  }
  return(fn_result)
}
