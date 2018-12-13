#include <Rcpp.h>
using namespace Rcpp;

#include "riso.h"


// [[Rcpp::export]]
DataFrame cppIso(CharacterVector cvseq){

	std::string cppseq = Rcpp::as<std::string>(cvseq);

	double	resultmass[NUMISOTOPES];
	double  resultprob[NUMISOTOPES];
	int failed;


	//void R_iso_seq(char * pseq,  double *resultmass, double *resultprob, int *failed){
	R_iso_seq(cppseq.c_str(),&(resultmass[0]), &resultprob[0], &failed);

  	NumericVector mass = NumericVector::create();
  	NumericVector prob = NumericVector::create();

  	mass.assign(resultmass,resultmass+NUMISOTOPES);
  	prob.assign(resultprob,resultprob+NUMISOTOPES);

	return Rcpp::DataFrame::create( Named("mass")= mass, Named("prob") = prob);
}




/*** R
x<-cppIso("GPPGPPGPQ")
*/



// [[Rcpp::export]]
DataFrame cppIsoAtom(int nC, int nH, int  nN, int  nO, int  nS ){

  double	resultmass[NUMISOTOPES];
  double  resultprob[NUMISOTOPES];
  int failed;

  R_iso_atom(nC, nH, nN, nO, nS, &(resultmass[0]), &resultprob[0], &failed);

  NumericVector mass = NumericVector::create();
  NumericVector prob = NumericVector::create();

  mass.assign(resultmass,resultmass+NUMISOTOPES);
  prob.assign(resultprob,resultprob+NUMISOTOPES);

  return Rcpp::DataFrame::create( Named("mass")= mass, Named("prob") = prob);
}




// [[Rcpp::export]]
DataFrame cppIsoAtom_p(IntegerVector Ratoms){

  double	resultmass[NUMISOTOPES];
  double  resultprob[NUMISOTOPES];
  int failed;

  //TODO: there's probably a neater way of doing this, but it's only 5 numbers..
  int atoms[NUMELEMENTS];
  for(int i=0;i<NUMELEMENTS;i++)
    atoms[i] = Ratoms[i];

  R_iso_atom_p(atoms , &(resultmass[0]), &resultprob[0], &failed);

  NumericVector mass = NumericVector::create();
  NumericVector prob = NumericVector::create();

  mass.assign(resultmass,resultmass+NUMISOTOPES);
  prob.assign(resultprob,resultprob+NUMISOTOPES);

  return Rcpp::DataFrame::create( Named("mass")= mass, Named("prob") = prob);
}



// [[Rcpp::export]]
IntegerVector aa_seq_to_atoms(CharacterVector cvseq){

  std::string cppseq = Rcpp::as<std::string>(cvseq);

  int numatoms[NUMELEMENTS];


  atomcount_from_peptide(cppseq.c_str(),numatoms);

  IntegerVector atomcount = IntegerVector::create();
  atomcount.assign(numatoms,numatoms+NUMELEMENTS);

  return atomcount;


}



