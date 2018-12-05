

normwindow <- function(datax,txlim,myby){
  #Generate the x values
  xout = seq(from = txlim[1], to = txlim[2], by = myby)
  #Resample against xout.
  yout <- approx(x=datax[,1], y=datax[,2], xout=xout, method="linear", rule = 2)
  #renormalise this segment
  yout$y = yout$y/max(yout$y)

  return(yout)
}

###############################################################################
#' Locally align a pair of  mass specs
#'
#' @param data1 the first MS
#' @param data2 the secon MS
#' @param mass the mass around which to do the alighment
#' @param alim the window around the mass, defaluts to [-1,6]
#' @param laglim the largest permissable lag
#' @param by the resampling rate, default value is 0.005
#' @param doplot whether to generate a plot of the alignment
#' @param verbose whether to write messages whilst processing
#' @return a dataframe holding the following fields:
#'   \item{lag}{The lag for sample 1}
#'   \item{cor}{The correlation coefficient for sample 1}
#' @export
msms_align <- function(data1,data2,mass,alim=c(-1,6),laglim = 0.3,by=0.005,doplot=T, verbose=T){

  if(doplot){
    par(mar=c(0.9,2.3,2.9,.3), mfrow = c(3,1), oma=c(5,0,2,0))
  }

  txlim = mass+alim


  #Now we need to generate the correlations data by normalising the sampling frequency

  #STEPSIZE & MAX LAG
  myby <- by
  #mylagmax gives the 'reverse scaling' of the stepsize - useful when comparing etc.
  mylagmax <- 1/myby



  #create an interpolation (isodists is accurate to 2 decimal places)


  yii <- normwindow(data1,txlim,myby)
  yjj <- normwindow(data2,txlim,myby)

  #TODO: Pass this data out so we can plot it elsewhere
  if(doplot){
    plot(yii,type="l",col="red",xlim=c(txlim[1]-2,txlim[2]+2))
    lines(yjj)
    segments(x0=mass,y0=0,y1=1,col="purple")
    title(sprintf("Data resampled to resolution %0.3f Da, mylagmax = %f",myby,mylagmax), line = "-2")
  }

  ######################################################################################
  #Now we can do the cross-correlation:                                 4*mylagmax
  ccd <- ccf(yjj$y,yii$y,ylim=c(-0.1,1.0),plot=doplot,axes=F, lag.max = mylagmax, main = "")
  #message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
  #plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )


  cor = ccd$acf[,,1]
  lag = ccd$lag[,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res$cor),]


  out <- data.frame(
    cor = res_max$cor,
    lag = -res_max$lag * myby
  )

  #let's set a limit on the permissable lag:


  lag_range = c(-laglim,laglim)
  lagset = res[ (res$lag > -(laglim/myby)) & (res$lag < (laglim/myby)),]

  res_lrmax = lagset[which.max(lagset$cor),]

  td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f\n max inrange c = %.2f at lag %.3f",res_max$cor,res_max$lag*myby,res_lrmax$cor,res_lrmax$lag * myby)

  if(verbose){
    message("\nComparing max correlation with within-range correlation:")
    message(sprintf("  cor = %0.2f, lag = %0.2f\nlrcor = %0.2f, lrlag = %0.2f\n",out$cor,out$lag,res_lrmax$cor,res_lrmax$lag * myby))

    #Here's where we can do a more detailed analysis
    message(sprintf("max cor = %0.2f at lag %0.2f",out$cor, out$lag))

    message(sprintf("There are %d points in the correlation from %0.2f to %0.2f (scaled to %0.2f to %0.2f)",nrow(res),res$lag[1],res$lag[nrow(res)],res$lag[1]*myby,res$lag[nrow(res)]*myby))

    message(td)
    #readline(sprintf("Hit <return> for %s",td))
  }

  #Add the points of max alignment to the plot:
  if(doplot){
    points(x=res_max$lag,y=res_max$cor,pch=19,col="red")
    points(x=res_lrmax$lag,y=res_lrmax$cor,pch=10,col="green",cex = 3)
  }

  labelvals = c(-1,-laglim,0,laglim,1)

  if(doplot){
    title(td, line = "-2")
    #axis(1, at = c(-mylagmax/2,0,mylagmax/2), labels = c(-0.5,0,0.5))
    axis(1, at = mylagmax * labelvals, labels = labelvals)
    axis(2)
    #text(0,-0.4,td)
  }


  #Plot the calculated peaks as probabilities
  #until this is in the library, we have to load it:
  #source("plotseqpeaks.R")
  if(doplot){

    #plot the peaks from the sequence
    ###ba_plotseqpeaks(ts,txlim)
    #####plotseqpeaks(ocow,myxlim)
    plot(yii,type="l",col="red",ylim=c(0,1.1))

    if(verbose)
      message(sprintf("Lag is %0.3f",res_max$lag*myby))



    mycol="red"
    if(res_max$cor > 0.1){
      if(abs(res_max$lag*myby) < 0.4){
        mycol="green"
      }
      else{
        message(sprintf("Lag too great: %0.3f",res_max$lag*myby))
      }
    }
    else{
      message("Weak correlation - ignore")
    }

    lines(x=data2[,1],    y = 0.05 + data2[,2]/max(data2[data2[,1]>txlim[1] & data2[,1] < txlim[2],2]),col="grey50")
    lines(x=data2[,1]+out$lag,y = 0.1 + data2[,2]/max(data2[data2[,1]>txlim[1] & data2[,1] < txlim[2],2]),col=mycol)

    segments(x0=mass,y0=0,y1=1,col="purple")
  }

  return(out)

}
