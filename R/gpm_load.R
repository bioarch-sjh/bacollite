



#TODO: get the info from google sheets....
load.gpm <- function(fn="~/tmp/bioarch_keri/161208meeting/gmp_human_collagen.dat"){


  gpm <- read.table(fn,
                    fill=T,header=T,sep="\t",comment.char="")

  gpm2 <- data.frame(seq=as.character(gpm$sequence),
                     nhyd = stringr::str_count(gpm$modifications,"P"),
                     nglut = (stringr::str_count(gpm$modifications,"Q")+stringr::str_count(gpm$modifications,"N")),
                     mass1=0, prob1=0 )

  gpm2 <- unique(gpm2)

  for(i in 1:nrow(gpm2)){

    #get the mass data...

    cd1 <- q2e::q2e_tpeaks(gpm2$seq[i])
    gpm2$mass1[i] = cd1$mass[1] +  (gpm2$nglut[i]*0.984015)+(gpm2$nhyd[i]*16)
    gpm2$prob1[i] = cd1$prob[1]
  }

  #remove outliers
  gpm2 <- gpm2[which(gpm2$mass1>795 & gpm2$mass1 <3500),]

  #Sort by mass
  gpm2 <- gpm2[order(gpm2$mass1),]

  #convert the seq column to string format:
  gpm2$seq <- sapply(gpm2$seq, as.character)

  return(gpm2)

}
