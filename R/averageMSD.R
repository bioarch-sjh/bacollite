
#' Generate an average of three MALDI spectrum samples
#'
#' @param reps a data structure holding three objects, S1, S2 and S3, each of which has a column called `mass` and `intensity`
#' @param stats approx ccf ksmooth lm predict
#' @importFrom stats approx ccf ksmooth lm predict
#' @return a dataframe holding the following fields:
#'   \item{mass}{The interpolated mass values}
#'   \item{intensity}{The average intensity values}
#' @export
averageMSD<-function(reps,resolution = 0.01){
  minmass <- floor(min(c(reps$s1$mass,reps$s2$mass,reps$s3$mass)))
  maxmass <- ceiling(max(c(reps$s1$mass,reps$s2$mass,reps$s3$mass)))
  xout = seq(from = minmass, to = maxmass, by = resolution)

  yi1 <- approx(x=reps$s1$mass, y=reps$s1$intensity, xout=xout, method="linear", rule = 2)
  yi2 <- approx(x=reps$s2$mass, y=reps$s2$intensity, xout=xout, method="linear", rule = 2)
  yi3 <- approx(x=reps$s3$mass, y=reps$s3$intensity, xout=xout, method="linear", rule = 2)


  yavg <- (yi1$y + yi2$y + yi3$y)/3

  yavg = data.frame(mass = yi1$x,intensity = yavg,stringsAsFactors = F)

  #colnames(yavg) = c("mass","intensity")

  return(yavg)
}
