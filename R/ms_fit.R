

#' Match a set of peptides with a Mass Spec
#'
#' @param peptides the set of peptides
#' @param sample the MS sample data
#' @return a logical vector stating which peptides match the data
#' @export
ms_fit<-function(peptides,sample,doplot=T,vlevel=0,corlim=0.0,laglim=0.6){

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


  for(i in 1:nrow(peptides)){

    #get the mass data...
    #TODO: hasn't this been done?
    #cd1 <- q2e::q2e_tpeaks(peptides$seq[i])
    cd1 <- ms_tpeaks(peptides$seq[i])

    cd1$mass <- cd1$mass + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)

    #TEST
    if(abs(cd1$mass[1]-peptides$mass1[i]) > 0.01){
      readline(sprintf("ERROR: mass calculation problem:\ncd1$mass[1]=%f\npeptides$mass1 = %f",cd1$mass[1],peptides$mass1))
      return (NA)
    }



    #readline(sprintf("Mass of sequence %s is %0.2f, nhyd = %d, nglut = %d, hit <return>",peptides$seq[i], cd1$mass[1], peptides$nhyd[i], peptides$nglut[i]))

    #TODO: consolidate with the _rams_ function...
    ###########################################################
    #RESTRICTION: Only do this for the range we have data for #
    if(max(cd1$mass) > 800 && min(cd1$mass) < 3500){

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


      if(enough_ions){
        message(sprintf("Max intensity  ratio sufficient in this segment (%f > %f) ", max(subms1[,2]), mir*max(sample$s1[,2]) ))

        cdshift <-cd1
        #cdshift$mass <- cd1$mass + (peptides$nglut[i]*0.984015)+(peptides$nhyd[i]*16)

        myxlim = c(lbl,ubl)

        align1 <- ms_align(cdshift,subms1,myxlim)
        align2 <- ms_align(cdshift,subms2,myxlim)
        align3 <- ms_align(cdshift,subms3,myxlim)

        message(sprintf("align1$lag is %0.3f, cor is %0.3f",align1$lag,align1$cor))
        message(sprintf("align2$lag is %0.3f, cor is %0.3f",align2$lag,align2$cor))
        message(sprintf("align3$lag is %0.3f, cor is %0.3f",align3$lag,align3$cor))

        lags[i,1]<-align1$lag
        lags[i,2]<-align2$lag
        lags[i,3]<-align3$lag

        cors[i,1]<-align1$cor
        cors[i,2]<-align2$cor
        cors[i,3]<-align3$cor

        hitplot=F
        if(abs(align1$lag) < laglim & align1$cor > corlim)
        { hitplot=T}
        if(abs(align2$lag) < laglim & align2$cor > corlim)
        { hitplot=T}
        if(abs(align3$lag) < laglim & align3$cor > corlim)
        { hitplot=T}

        if(hitplot){

          plotno = plotno+1

          hits[i]<-T


          if(doplot){

            message(sprintf("\nPlot number %d\nSegment at row %d of %d",
                            plotno,i,nrow(peptides)))
            #0.984015 - if Q changes to E - add this much....

            #TODO: start-4 is used a few times - so it needs setting as a variable
            #if(nrow(phydp)>0){
            #	mymain <- sprintf(
            #	"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=%0.4f (max = %0.4f)",
            #	plotno,start-4,sequence,e,nglut,p,nhyd,pnh$prob[which(pnh$nhyd == p)],max(pnh$prob))
            #}
            #else{
            #	mymain <- sprintf(
            #	"plot %d, seqpos %d\n%s\nnglut = %d/%d, nhyd = %d/%d, hp=UNKNOWN",
            #	plotno,start-4,sequence,e,nglut,p,nhyd)
            #}

            mymain <- sprintf(
              "plot %d, entry %d, mass %0.3f\n%s\nndean = %d, nhyd = %d lag: %0.2f,%0.2f,%0.2f cor: %0.2f,%0.2f,%0.2f",
              plotno,i,peptides$mass1[i],
              peptides$seq[i],peptides$nglut[i],peptides$nhyd[i],
              align1$lag,align2$lag,align3$lag,
              align1$cor,align2$cor,align3$cor)

            myxxlim <- c(lbl-1,ubl+1)
            plot(1, type="n", xlab="Mass", ylab = "Probability",
                 xlim=myxxlim, ylim=c(0, 1), main=mymain)

            #legend('topright',spots,'pre-align'), lty = c(1,1,1,1),
            #	   col=c('red','green','blue','grey'),ncol=1,bty ="n")

            #ba_plotseqpeaks(cd1,myxlim)

            #cc = cc+1;
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

            #readline("hit <return> to continue")
          }
        }
      }

    }
    else{
      message(sprintf("Mass %0.2f out of range",cd1$mass[1]))
    }


  }

  message(sprintf("%d theoretical peptides in range",count))

  alignments <- data.frame(hit=hits,lag=lags,cor=cors)
  colnames(alignments)<- c("hit","lag1","lag2","lag3","cor1","cor2","cor3")

  return(alignments)
}
