#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_PG(int* ptrIcol,int* ptrJcol,int* ptrdimAfa,double Alpha,double Beta){

	int sumI=0;
	int sumJ=0;

	double term1=0.0;
	double term2=0.0;
	double term3=0.0;
	double term4=0.0;

	for (int j = 0; j < ptrdimAfa[1]; j++){
		sumI+=ptrIcol[j];
		sumJ+=ptrJcol[j];
	}

	for (int j = 0; j < ptrdimAfa[1]; j++){
		term1+=pow((double) ptrIcol[j]/sumI,Alpha) * pow((double) ptrJcol[j]/sumJ,Beta);
		term2+=pow((double) ptrJcol[j]/sumJ,Alpha) * pow((double) ptrIcol[j]/sumI,Beta);
		term3+=pow((double) ptrIcol[j]/sumI,Alpha+Beta);
		term4+=pow((double) ptrJcol[j]/sumJ,Alpha+Beta);
	}

	return (term1+term2)/(term3+term4);
}

