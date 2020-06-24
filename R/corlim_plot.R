





#TODO: Document variables properly
#TODO: Pass in corlim (see line 44,45)

#' get the hits for a range of correlation limit values
#'
#' @param psample outputs from the ms_fit function for a set of peptides
#' @param corlim a sequence of correlation limit values for thresholding the data
#' @param laglim maximum acceptable value of lag for each correlation
#' @keywords bruker
#' @export
corlim_plot <- function(psample,sarea="",pnpep,pcld,presult,pnames,manname="unknown",pcols=c("#e2ba5e","#3597c6","#e3645f","#793787"),warn=T){


  calcid = presult$id[presult$score == max(presult$score)]
  if(length(calcid)>1)
    calcid = "unknown"


  title = sprintf("Sample '%s' %s\nmanual ID: '%s'; Calc ID: '%s'\nscores ",psample$name,sarea,manname,calcid)

  for(ss in 1:nrow(result)){
    title = sprintf("%s %s = %0.3f",title,presult$id[ss],presult$score[ss])
  }

  par(mar = c(4,4,5,4))
  plot(NA,xlab="Correlation Threshold",ylab = "Number of Hits",ylim=c(0,15),xlim = c(0,1))

  maincol="black"

  if(warn)
    if(manname != presult$id[presult$score == max(presult$score)])
      maincol="red"

  title(main = title,col.main=maincol)

  for(ss in 1:length(cld)){

    points(x=corlim,y=pcld[[ss]]$nh,col=pcols[ss])
    lines(x=corlim,y=pcld[[ss]]$nh,col=pcols[ss])

  }

  legend("topright",legend = corlab,col=pcols,lty = 1,pch=1)

}
