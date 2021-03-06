

#' String of valid amino acid codes
#'
"amino_acid_codes"



#' Collagen sequences for a range of mammals generated by M.Collins / BioArch
#'
"bioarch_mammal_sequences"


#' EXPERIMENTAL "versioned" copy of bioarch_mammal_sequences
#'
"bms0.1.2"


#' Descriminatory peptide markers for cow
#'
#' @format A data frame with 5 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#' }
"dm_cow"


#' Descriminatory peptide markers for goat
#'
#' @format A data frame with 5 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#' }
"dm_goat"


#' Descriminatory peptide markers for sheep
#'
#' @format A data frame with 5 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#' }
"dm_sheep"


#' All peptides for cow
#'
#' @format A data frame with 905 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#' }
"pepcow"


#' Descriminatory peptide markers for goat
#'
#' @format A data frame with 909 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#' }
"pepgoat"


#' Descriminatory peptide markers for sheep
#'
#' @format A data frame with 5 rows and 7 variables:
#' \describe{
#'   \item{seq}{The peptide sequence}
#'   \item{nhyd}{The number of hydroxylations}
#'   \item{nglut}{The number of deamidations of glutamine}
#'   \item{mass1}{The mass of the first isotope of the peptide}
#'   \item{seqpos}{The position in the collagen sequence}
#'   \item{missed.cleaves}{The number of missed cleaves in the peptide}
#'   \item{collagen}{The collagen number that the sequence came from}
#'   \item{prob}{Estimate of the probability of the nhyd value}
#' }
"pepsheep"
