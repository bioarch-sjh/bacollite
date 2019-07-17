
#' plot the peaks of a mass spec
#'
#' @param data a data strucyture in MALDI format
#' @param myxlim c(min,max) for the range of the data you want to plot
#' @keywords bruker
#' @export
#' @examples
#' plotseqpeaks(cd1,myxlim)
ba_plotseqpeaks <- function(data,myxlim){

  x <- data$mass
  y <- data$prob

  plot(x, y, xlim = myxlim, ylim = c(0,1),xlab = "Mass (Da)", ylab="Intensity")

  #This is the function that draws the lines down:
  segments(x0=x, y0=y, y1=0, col=8)
  points(x, y, pch=21, col=1, bg=2)

}
