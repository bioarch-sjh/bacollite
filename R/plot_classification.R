
plot.classification <- function(result,samplename,manualID,corlim,cld,labs){

  if(max(result$score)>0.001){
    title = sprintf("Sample '%s': manual ID: '%s'; Calc ID: '%s'\nscores ",
                    samplename,manualID,result$id[result$score == max(result$score)])
  }
  else{
    title = sprintf("Sample '%s': manual ID: '%s'; Calc ID: 'unknown'\nscores ",samplename,manualID)
  }

  for(ss in 1:nrow(result)){
    title = sprintf("%s %s = %0.3f",title,result$id[ss],result$score[ss])
  }

  par(mar = c(4,4,5,4))
  plot(NA,xlab="Correlation Threshold",ylab = "Number of Hits",ylim=c(0,3*nrow(pepdata[[1]])),xlim = c(0,1), main = title)

  cols = c("orange","blue","red","black")

  for(pp in 1:length(cld)){
    points(x=corlim,y=cld[[pp]]$nh,col=cols[pp])
    #points(x=corlim,y=cld[[2]]$nh,col="blue")
    #points(x=corlim,y=cld[[3]]$nh,col="red")
    #points(x=corlim,y=cld[[4]]$nh,col="black")

    lines(x=corlim,y=cld[[pp]]$nh,col=cols[pp])
    #lines(x=corlim,y=cld[[2]]$nh,col="blue")
    #lines(x=corlim,y=cld[[3]]$nh,col="red")
    #lines(x=corlim,y=cld[[4]]$nh,col="black")
  }

  legend("topright",legend = labs,col=c("orange","blue","red","black"),lty = 1,pch=1)
}
