#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"
#include "DP_RE.h"


void DP_single_population(int* ptrCol,int* ptrdimAfa,double* ptrAlphaPro,int nAlphaPro,int ENS,double* ptrResults){


	double Sum=0.0;
	int nCol_P=0;
	int iter=0;

	for (int j = 0; j < ptrdimAfa[1]; j++)	{
		Sum+=ptrCol[j];
		if (ptrCol[j] != 0){nCol_P++;}//Amount of non zero values in ptrCol
	}

	double* ptrCol_P; ptrCol_P=dvector(0,nCol_P-1);


	for (int j = 0; j < ptrdimAfa[1]; j++)	{
		if (ptrCol[j] != 0){
			ptrCol_P[iter]=ptrCol[j]/Sum;
			iter++;}
	}	


	for (int i = 0; i < nAlphaPro; i++){
		ptrResults[i]=DP_RE(ptrCol_P,nCol_P,ptrAlphaPro[i],ENS);
	}

	free_dvector(ptrCol_P,0,nCol_P-1);
}
