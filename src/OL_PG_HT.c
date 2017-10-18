#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_PG_HT(int* ptrIcol,int* ptrJcol,int* ptrdimAfa,double Alpha,double Beta){

	int sumI=0;
	int sumJ=0;

	int dpI=0;
	int dpJ=0;


	double* ptrPaI; ptrPaI=dvector(0,ptrdimAfa[1]-1);
	double* ptrPaJ; ptrPaJ=dvector(0,ptrdimAfa[1]-1);

	double* ptrAbcI; ptrAbcI=dvector(0,ptrdimAfa[1]-1);
	double* ptrAbcJ; ptrAbcJ=dvector(0,ptrdimAfa[1]-1);

	double numerator=0.0;
	double denominator=0.0;


	for (int j = 0; j < ptrdimAfa[1]; j++){
		sumI+=ptrIcol[j];
		sumJ+=ptrJcol[j];

		if (ptrIcol[j]==1){dpI++;}
		if (ptrJcol[j]==1){dpJ++;}
	}

	if (dpI==sumI){dpI--;}
	if (dpJ==sumJ){dpJ--;}//This might be erronous in Maciej's code, because in his code sumJ is sumI


	for (int j = 0; j < ptrdimAfa[1]; j++){
		ptrPaI[j]=((double)ptrIcol[j]/sumI)*(1.0-(double)dpI/sumI);
		ptrPaJ[j]=((double)ptrJcol[j]/sumJ)*(1.0-(double)dpJ/sumJ);

		ptrAbcI[j]=1.0-pow(1-ptrPaI[j],sumI);
		ptrAbcJ[j]=1.0-pow(1-ptrPaJ[j],sumJ);

		if (ptrAbcI[j]==0){ptrAbcI[j]=1.0;}
		if (ptrAbcJ[j]==0){ptrAbcJ[j]=1.0;}


		numerator	+= (pow(ptrPaI[j],Alpha) * pow(ptrPaJ[j],Beta))    /   (ptrAbcI[j]*ptrAbcJ[j]);
		denominator	+= pow(ptrPaI[j],2*Alpha)/ptrAbcI[j]   +   pow(ptrPaJ[j],2*Beta)/ptrAbcJ[j];
	}


	free_dvector(ptrPaI,0,ptrdimAfa[1]-1);
	free_dvector(ptrPaJ,0,ptrdimAfa[1]-1);
	free_dvector(ptrAbcI,0,ptrdimAfa[1]-1);
	free_dvector(ptrAbcJ,0,ptrdimAfa[1]-1);

	return 2*numerator/denominator;
}
