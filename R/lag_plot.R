



#' Plot the lag correction and fit a line to it for a peptide sample
#'
#' @param fn full path of the raw gpm text fi
#' @export
#' @examples
#' lag_plot(markers,hits)
lag_plot <- function(markers,hits,no=1,ylim=c(-1.4,1.0),mypoly=2,title="Lag Plot",corlim=0.1,laglim=0.5){

  #TODO: check that the data coming in is sane

  #Plot the raw data:

  hits <- hits[[no]]

  th1 <- which(abs(hits$lag1)< laglim & hits$cor1 > corlim)
  hh1 <- hits[   th1,]
  mm1 <- markers[th1,]

  th2 <- which(abs(hits$lag2)< laglim & hits$cor2 > corlim)
  hh2 <- hits[   th2,]
  mm2 <- markers[th2,]

  th3 <- which(abs(hits$lag3)< laglim & hits$cor3 > corlim)
  hh3 <- hits[   th3,]
  mm3 <- markers[th3,]


  mypch=20
  mycex = 0.6

  plot(x=markers$mass1,y=hits$lag1,ylim=ylim,pch=mypch,cex=mycex,col="red",xlab="Mass (Da)", ylab= "Lag (Da)",main=title)
  points(x=markers$mass1,y=hits$lag2,pch=mypch,cex=mycex,col="green")
  points(x=markers$mass1,y=hits$lag3,pch=mypch,cex=mycex,col="blue")

  text(x=mm1$mass1,y=-1.5,labels = mm1$seq,cex=0.4,srt=90,adj=0,col= adjustcolor("red",alpha.f=0.2))
  text(x=mm2$mass1,y=-1.5,labels = mm2$seq,cex=0.4,srt=90,adj=0,col=adjustcolor("green",alpha.f=0.2))
  text(x=mm3$mass1,y=-1.5,labels = mm3$seq,cex=0.4,srt=90,adj=0,col=adjustcolor("blue",alpha.f=0.2))


  model1 <- lm(hh1$lag1 ~ poly(mm1$mass1,mypoly))
  pi1 <- predict(model1,data.frame(mm1$mass1,interval="confidence",level=0.99))
  lines(mm1$mass1,pi1,col="red")

  model2 <- lm(hh2$lag2 ~ poly(mm2$mass1,mypoly))
  pi2 <- predict(model2,data.frame(mm2$mass1,interval="confidence",level=0.99))
  lines(mm2$mass1,pi2,col="green")

  model3 <- lm(hh3$lag3 ~ poly(mm3$mass1,mypoly))
  pi3 <- predict(model3,data.frame(mm3$mass1,interval="confidence",level=0.99))
  lines(mm3$mass1,pi3,col="blue")

}
