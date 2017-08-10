

#' Get the ratio of ion counts for two peptides
#'
#' @param sample the MALDI spectrum - three replictes, as loaded by load.sample()
#' @param peptide1 the first peptide, which forms the numerator of the ratio calculation
#' @param peptide2 the second peptide, which forms the denominator of the ratio calculation
#' @export
#' @examples
#' sample <- load.sample("~/samples","01",c("_1","_2","_3"))
#' peptides <- load.mcs()
#' ratios <- peptide.ratio(sample, peptides[1,], peptides[2,])
peptide.ratio <- function(sample, peptide1, peptide2){

  fita <- ms_fit(peptide1,sample,doplot=F,force = T)
  fitb <- ms_fit(peptide2,sample,doplot=F,force = T)

  r1 <- fita$ion1 / fitb$ion1
  r2 <- fita$ion2 / fitb$ion2
  r3 <- fita$ion3 / fitb$ion3

  return (c(r1,r2,r3))

}


##Create the space to hold the outputs:
#NEn$ratio1 = 0
#NEn$ratio2 = 0
#NEn$ratio3 = 0
#NEn$sample = NA

#NEsample = list()
#NEfita = list()
#NEfitb = list()
#for(cc in 1:nrow(NEn)){
  #message(sprintf("Processing entry %d",cc))

#  NEsample[[cc]] <- load.sample(NEn$root[cc],NEn$name[cc],
#                                c(
#                                  as.character(NEn$s1[cc]),
#                                  as.character(NEn$s2[cc]),
#                                  as.character(NEn$s3[cc])))
#  NEfita[[cc]] <- ms_fit(pda,NEsample[[cc]],doplot=F,force = T)
#  NEfitb[[cc]] <- ms_fit(pdb,NEsample[[cc]],doplot=F,force = T)

#  NEn$ratio1[cc] = NEfita[[cc]]$ion1 / NEfitb[[cc]]$ion1
#  NEn$ratio2[cc] = NEfita[[cc]]$ion2 / NEfitb[[cc]]$ion2
#  NEn$ratio3[cc] = NEfita[[cc]]$ion3 / NEfitb[[cc]]$ion3

#}

#NIt$ratio1 = 0
#NIt$ratio2 = 0
#NIt$ratio3 = 0
#NIt$sample = NA
#NIsample = list()
#NIfita = list()
#NIfitb = list()
#for(cc in 1:nrow(NIt)){
#  #message(sprintf("Processing entry %d, spot",cc))

#  NIsample[[cc]] <- load.sample(NIt$root[cc],NIt$name[cc],
#                                c(
#                                  as.character(NIt$s1[cc]),
#                                  as.character(NIt$s2[cc]),
#                                  as.character(NIt$s3[cc])))
#  NIfita[[cc]] <- ms_fit(pda,NIsample[[cc]],doplot=F,force = T)
#  NIfitb[[cc]] <- ms_fit(pdb,NIsample[[cc]],doplot=F,force = T)

#  NIt$ratio1[cc] = NIfita[[cc]]$ion1 / NIfitb[[cc]]$ion1
#  NIt$ratio2[cc] = NIfita[[cc]]$ion2 / NIfitb[[cc]]$ion2
#  NIt$ratio3[cc] = NIfita[[cc]]$ion3 / NIfitb[[cc]]$ion3
