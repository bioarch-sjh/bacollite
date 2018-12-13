


###############################################################################
# This replaces 'ms_offset_peaklineplot' in bioarch/dev folder
###############################################################################
###############################################################################
ms_offset_peaklineplot <- function(ms,offset,mycol){
  lines(x=(ms[,1]+offset),ms[,2]/max(ms[,2]),col=mycol)
}



###############################################################################
# This replaces 'ba_ms_align' in bioarch/dev folder
###############################################################################
###############################################################################
#' Align a theoretical peptide with a mass spec
#'
#' @param ts theoretical spectrum peaks for a particular peptide
#' @param data the MALDI
#' @param txlim the range of masses over which to carry out the alignment
#' @param gauss the level of gaussian smoothing. defaults to NA (no smoothing)
#' @param normlim the minimum upper value that the intensities should be normalised to. Defaults to NA (no minimum)
#' @param doplot whether to generate a plot of the alignment
#' @param verbose whether to write messages whilst processing
#' @param ccylim range of y axis in cross-correlation plot (defaults to [-0.1,0.5])
#' @return a dataframe holding the following fields:
#'   \item{lag1}{The lag for sample 1}
#'   \item{lag2}{The lag for sample 2}
#'   \item{lag3}{The lag for sample 3}
#'   \item{cor1}{The correlation coefficient for sample 1}
#'   \item{cor2}{The correlation coefficient for sample 2}
#'   \item{cor3}{The correlation coefficient for sample 3}
#'   \item{ion1}{The proportion of total ions for this peptide in sample 1}
#'   \item{ion2}{The proportion of total ions for this peptide in sample 2}
#'   \item{ion3}{The proportion of total ions for this peptide in sample 3}
#' @export
ms_align <- function(ts,data,txlim,gauss=NA,normlim=NA,cctitle=NA,doplot=F, verbose=F,ccylim=c(-0.1,0.5)){

  if(doplot){
    #save the old par
    #op <- par(no.readonly = TRUE)
    #set up a new plot window to show the alignment
    #dev.new()
    #set new par
    par(mar=c(4,4,0.5,0), mfrow = c(3,1), oma=c(0,0,0,0))
  }



  #Now we need to generate the correlations data by normalising the sampling frequency

  #STEPSIZE & MAX LAG
  myby <- 0.005 #125
  #mylagmax gives the 'reverse scaling' of the stepsize - useful when comparing etc.
  mylagmax <- 1/myby



  #create an interpolation (isodists is accurate to 2 decimal places)

  #Generate the x values
  xout = seq(from = txlim[1], to = txlim[2], by = myby)
  #Resample against xout.
  yii <- approx(x=data[,1], y=data[,2], xout=xout, method="linear", rule = 2)
  #renormalise this segment
  #yii$y = yii$y/max(yii$y)

  nmax = max(yii$y)
  if(!is.na(normlim)){
    nmax = max(nmax,normlim)
  }

  yii$y = (yii$y-min(yii$y))/(nmax-min(yii$y))

  #TODO: Pass this data out so we can plot it elsewhere
  if(doplot){
    plot(yii,type="l",col="red",xlab="Mass (Da)",ylab="Intensity")
  }

  #Now resample the theoretical data:
  yri <- approx(x=ts$mass,y=ts$prob,xout=xout, method="linear", rule = 2)
  #set yvals to zero
  yri$y[] <-0
  #go through each peak
  for(i in 1:length(ts$prob)){
    idx <- which.min(abs(yri$x-ts$mass[i]))
    yri$y[idx] <- ts$prob[i]
  }

  #Apply gaussian smoothing if set
  if(!is.na(gauss)){
    yrii <- ksmooth(yri$x,yri$y,"normal",bandwidth = gauss)
    yrii$y <- yrii$y / max(yrii$y)
    yri$y <- yrii$y
  }


  #TODO: Pass this data out so we can plot it elsewhere
  if(doplot){
    lines(yri)
    title(sprintf("Data/model resampled to %0.2f Da",myby), line = "-2")
  }

  ######################################################################################
  #Now we can do the cross-correlation:                                 4*mylagmax
  ccd <- ccf(yri$y,yii$y,ylim=ccylim,plot=doplot,axes=F, lag.max = mylagmax, main = "")
  #message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
  #plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )


  cor = ccd$acf[,,1]
  lag = ccd$lag[,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res$cor),]


  out <- data.frame(
    cor = res_max$cor,
    lag = res_max$lag * myby
  )

  #let's set a limit on the permissable lag:
  laglim = 0.3

  lag_range = c(-laglim,laglim)
  lagset = res[ (res$lag > -(laglim/myby)) & (res$lag < (laglim/myby)),]

  res_lrmax = lagset[which.max(lagset$cor),]

  if(is.na(cctitle)){
    cctitle <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f\n max inrange c = %.2f at lag %.3f",res_max$cor,res_max$lag*myby,res_lrmax$cor,res_lrmax$lag * myby)
  }
  
  if(verbose){
    message("\nComparing max correlation with within-range correlation:")
    message(sprintf("  cor = %0.2f, lag = %0.2f\nlrcor = %0.2f, lrlag = %0.2f\n",out$cor,out$lag,res_lrmax$cor,res_lrmax$lag * myby))

    #Here's where we can do a more detailed analysis
    message(sprintf("max cor = %0.2f at lag %0.2f",out$cor, out$lag))

    message(sprintf("There are %d points in the correlation from %0.2f to %0.2f (scaled to %0.2f to %0.2f)",nrow(res),res$lag[1],res$lag[nrow(res)],res$lag[1]*myby,res$lag[nrow(res)]*myby))

    message(cctitle)
    #readline(sprintf("Hit <return> for %s",cctitle))
  }

  if(doplot){
    points(x=res_max$lag,y=res_max$cor,pch=19,col="red")
    points(x=res_lrmax$lag,y=res_lrmax$cor,pch=10,col="green",cex = 3)
  }

  labelvals = c(-1,-laglim,0,laglim,1)

  if(doplot){
    title(cctitle, line = "-2")
    #axis(1, at = c(-mylagmax/2,0,mylagmax/2), labels = c(-0.5,0,0.5))
    axis(1, at = mylagmax * labelvals, labels = labelvals)
    axis(2)
    #text(0,-0.4,cctitle)
  }


  #Plot the calculated peaks as probabilities
  #until this is in the library, we have to load it:
  #source("plotseqpeaks.R")
  if(doplot){

    #plot the peaks from the sequence
    ba_plotseqpeaks(ts,txlim)
    #####plotseqpeaks(ocow,myxlim)

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

    lines(x=data[,1],    y = data[,2]/max(data[,2]),col="grey50")
    lines(x=data[,1]+res_max$lag*myby,y = data[,2]/max(data[,2]),col=mycol)

    title("Alignment",line = -2)
    
    #readline("hit <return> to close the plot window and carry on")
    #dev.off()
    #plot.new()
    #dev.set(dev.prev()) # go back to first
    #reset the par
    #op
    #par(op)
  }
  
  
  
  mtext(text="Mass (Da)",side=1,line=2,outer=TRUE)

  return(out)

}
