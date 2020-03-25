

#' Match a set of peptides with a Mass Spec
#'
#' @param peptides the set of peptides, in bacollite format
#' @param sample the MS sample data, in bacollite format - three replicates required
#' @param doplot boolean to say if the fit should be plotted
#' @param force boolean to say if *all* fits should be calculated, instead of just the good matches
#' @param vlevel verbose level. 0 means no comments, 3 means full comments
#' @param corlim threshold for correlation scores
#' @param laglim threshold for lag scores
#' @param gauss the level of gaussian smoothing. defaults to NA (no smoothing)
#' @param ionlim the fraction of the highest intensity that is used to calculate the scaling limit for ms_align
#' @param ignore_warnings whether to ignore the check for the input and calculated mass
#' @importFrom graphics axis legend lines mtext par plot points segments text title
#' @return a dataframe holding the following fields:
#'   \item{hit}{Whether a match was found for this peptide}
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
ms_fit<-function(peptides,sample,doplot=T,force=F,vlevel=0,corlim=0.0,laglim=0.6,gauss=NA,ionlim=NA,ignore_warnings=F,use_ms_iso=T){

  #initialise:
  plotno <- 0
  moff <- 1.5
  count <- 0

  hits <- logical(nrow(peptides))
  hits[] <- FALSE

  lags <- matrix(nrow=nrow(peptides),ncol=3)
  lags[] <- 0

  cors <- matrix(nrow=nrow(peptides),ncol=3)
  cors[] < -0

  ions <- matrix(nrow=nrow(peptides),ncol=3)
  ions[] <- 0

  sumions= c(sum(sample$s1$intensity),sum(sample$s2$intensity),sum(sample$s3$intensity))

  for(i in 1:nrow(peptides)){

    #get the mass data...
    #TODO: hasn't this been done?
    #cd1 <- q2e::q2e_tpeaks(peptides$seq[i])

    if(use_ms_iso){
    	cd1 <- ms_iso(peptides$seq[i],ndeamidations=peptides$nglut[i],nhydroxylations=peptides$nhyd[i])
    }
    else{
  		cd1 <- ms_tpeaks(peptides$seq[i])
  		cd1$mass <- cd1$mass + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)
  	}
    warn_mass_error<-F

    #TEST
    if(abs(cd1$mass[1]-peptides$mass1[i]) > 0.05){
      if(!ignore_warnings){
        readline(sprintf("ERROR: Peptide %d mass calculation problem:\ncd1$mass[1]=%f\npeptides$mass1 = %f\nnhyd = %d, nglut = %d, seq = %s",
                         i,cd1$mass[1],peptides$mass1[i],
                         peptides$nhyd[i],peptides$nglut[i],peptides$seq[i]))
        return (NA)
      }
      else{
        message(sprintf("ERROR: Peptide %d mass calculation problem:\ncd1$mass[1]=%f\npeptides$mass1 = %f\nnhyd = %d, nglut = %d, seq = %s",
                         i,cd1$mass[1],peptides$mass1[i],
                         peptides$nhyd[i],peptides$nglut[i],peptides$seq[i]))

      }
    }



    #readline(sprintf("Mass of sequence %s is %0.2f, nhyd = %d, nglut = %d, hit <return>",peptides$seq[i], cd1$mass[1], peptides$nhyd[i], peptides$nglut[i]))

    #TODO: consolidate with the _rams_ function...
    ###########################################################
    #RESTRICTION: Only do this for the range we have data for #
    if((max(cd1$mass) > 800 && min(cd1$mass) < 3500)| force) {

      if(vlevel>3){
        message(sprintf("Sequence is %s\nThere are %d prolines in this segment",
                        peptides$seq[i],peptides$nhyd[i]))
        #message(sprintf("There are %d glutamines  in this segment",nglut))
        message(sprintf("Mass range is %f to %f",min(cd1$mass),max(cd1$mass) ))
        #readline("hit <return> to continue")
      }
      count = count +1

      lbl <- min(cd1$mass) - moff# + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)
      ubl <- max(cd1$mass) + moff# + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)

      subms1 <- ms_subrange(sample$s1,lbl,ubl)
      subms2 <- ms_subrange(sample$s2,lbl,ubl)
      subms3 <- ms_subrange(sample$s3,lbl,ubl)

      mir <- 0.05

      enough_ions<-F
      if(max(subms1[,2]) > (mir*max(sample$s1[,2])) ){
        enough_ions<-T
      }
      if(max(subms2[,2]) > (mir*max(sample$s2[,2])) ){
        enough_ions<-T
      }
      if(max(subms3[,2]) > (mir*max(sample$s3[,2])) ){
        enough_ions<-T
      }


      if(enough_ions | force){
        if(vlevel > 0){
          message(sprintf("Max intensity  ratio sufficient in this segment (%f > %f) ", max(subms1[,2]), mir*max(sample$s1[,2]) ))
        }
        cdshift <-cd1
        #cdshift$mass <- cd1$mass + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)

        myxlim = c(lbl,ubl)

        #set the normalisation limits if needed
        normlim=c(NA,NA,NA)
        if(!is.na(ionlim)){
          normlim[1] <- max(sample$s1$intensity)*ionlim
          normlim[2] <- max(sample$s2$intensity)*ionlim
          normlim[3] <- max(sample$s3$intensity)*ionlim
        }

        align1 <- ms_align(cdshift,subms1,myxlim,gauss,normlim[1])
        align2 <- ms_align(cdshift,subms2,myxlim,gauss,normlim[2])
        align3 <- ms_align(cdshift,subms3,myxlim,gauss,normlim[3])

        if(vlevel > 0){
          message(sprintf("align1$lag is %0.3f, cor is %0.3f",align1$lag,align1$cor))
          message(sprintf("align2$lag is %0.3f, cor is %0.3f",align2$lag,align2$cor))
          message(sprintf("align3$lag is %0.3f, cor is %0.3f",align3$lag,align3$cor))
        }
        lags[i,1]<-align1$lag
        lags[i,2]<-align2$lag
        lags[i,3]<-align3$lag

        cors[i,1]<-align1$cor
        cors[i,2]<-align2$cor
        cors[i,3]<-align3$cor

        ions[i,1]<- sum( sample$s1$intensity[ (sample$s1$mass > (cd1$mass[1] + align1$lag - 0.5)) &  (sample$s1$mass < (cd1$mass[1] + align1$lag + 4.5)) ]  )/sumions[1]
        ions[i,2]<- sum( sample$s2$intensity[ (sample$s2$mass > (cd1$mass[1] + align2$lag - 0.5)) &  (sample$s2$mass < (cd1$mass[1] + align2$lag + 4.5)) ]  )/sumions[2]
        ions[i,3]<- sum( sample$s3$intensity[ (sample$s3$mass > (cd1$mass[1] + align3$lag - 0.5)) &  (sample$s3$mass < (cd1$mass[1] + align3$lag + 4.5)) ]  )/sumions[3]

        hitplot=F
        if(abs(align1$lag) < laglim & align1$cor > corlim)
        { hitplot=T}
        if(abs(align2$lag) < laglim & align2$cor > corlim)
        { hitplot=T}
        if(abs(align3$lag) < laglim & align3$cor > corlim)
        { hitplot=T}


        if(hitplot | force){

          plotno = plotno+1

          hits[i]<-T


          if( doplot ){

            if(vlevel > 0){
              message(sprintf("\nPlot number %d\nSegment at row %d of %d",
                            plotno,i,nrow(peptides)))
            }

            mymain <- sprintf(
              "plot %d, entry %d, mass %0.3f\npos = %d %s\nndeam = %d, nhyd = %d lag: %0.2f,%0.2f,%0.2f\ncor: %0.2f,%0.2f,%0.2f",
              plotno,i,peptides$mass1[i],
              peptides$seqpos[i],peptides$seq[i],
              peptides$nglut[i],peptides$nhyd[i],
              align1$lag,align2$lag,align3$lag,
              align1$cor,align2$cor,align3$cor)


            if(warn_mass_error){
              mymain <- sprintf("%s%s","(MASS ERROR) ",mymain)
            }

            myxxlim <- c(lbl-1,ubl+1)
            par(mar=c(3, 3, 6.5, 1))
            plot(1, type="n", xlab="Mass", ylab = "Probability",
                 xlim=myxxlim, ylim=c(0, 1), main=mymain)

            x <- cd1$mass # + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)
            y <- cd1$prob

            segments(x0=x, y0=y, y1=0, col="black")
            points(x, y, pch=21, col="black", bg="red")

            #####################

            ms_offset_peaklineplot(sample$s1,0,"grey")
            ms_offset_peaklineplot(sample$s2,0,"grey")
            ms_offset_peaklineplot(sample$s3,0,"grey")

            ms_offset_peaklineplot(sample$s1,align1$lag,"red")
            ms_offset_peaklineplot(sample$s2,align2$lag,"green")
            ms_offset_peaklineplot(sample$s3,align3$lag,"blue")

          }
        }
      }

    }
    else{
      if(vlevel > 0){
        message(sprintf("Mass %0.2f out of range",cd1$mass[1]))
      }
    }


  }

  if(vlevel > 0){
    message(sprintf("%d theoretical peptides in range",count))
  }

  alignments <- data.frame(hit=hits,lag=lags,cor=cors, ion = ions)
  colnames(alignments)<- c("hit","lag1","lag2","lag3","cor1","cor2","cor3","ion1","ion2","ion3")

  return(alignments)
}
