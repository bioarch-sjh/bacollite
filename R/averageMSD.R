averageMSD<-function(reps){
  minmass <- floor(min(c(reps$s1$mass,reps$s2$mass,reps$s3$mass)))
  maxmass <- ceiling(max(c(reps$s1$mass,reps$s2$mass,reps$s3$mass)))
  stepby <- 0.01
  xout = seq(from = minmass, to = maxmass, by = stepby)
  yi1 <- approx(x=reps$s1$mass, y=reps$s1$intensity, xout=xout, method="linear", rule = 2)
  yi2 <- approx(x=reps$s2$mass, y=reps$s2$intensity, xout=xout, method="linear", rule = 2)
  yi3 <- approx(x=reps$s3$mass, y=reps$s3$intensity, xout=xout, method="linear", rule = 2)

  yavg <- yi1
  yavg$y <- (yi1$y + yi2$y + yi3$y)/3

  return(yavg)
}
