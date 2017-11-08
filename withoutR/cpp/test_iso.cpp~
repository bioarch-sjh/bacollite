
#include <stdlib.h>
#include <stdio.h>



#include "riso.h"



int main(){

	//void R_iso_atom(int nC, int nH, int  nN, int  nO, int  nS, 
	//				  double *resultmass, double *resultprob, int *failed)

	double	resultmass[NUMISOTOPES];
	double  resultprob[NUMISOTOPES];
	int numatoms[NUMELEMENTS];
	int failed;
  
    const char sequence[] = "GPPGPPGPQ";
  
  
  	printf("Testing a simple sequence first\n");
  	R_iso_seq(sequence,resultmass,resultprob,&failed);
  
  
  
  	printf("\n\n---\nNow looking at the \"atomic\" approach\n\n");
  
  
  
  
	R_iso_atom(36,54,10,11,0,resultmass,resultprob,&failed);



  	printf("\n\n---\nGetting atoms from seq and using that:\n\n");



    atomcount_from_peptide(sequence, numatoms);

	R_iso_atom_p(numatoms,resultmass,resultprob,&failed);


}



