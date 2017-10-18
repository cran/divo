#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_RD(int* ptrIcol,int* ptrJcol,int* ptrdimAfa,double Alpha){

	int sumI=0;
	int sumJ=0;

	double term1=0.0;

	double temp=0.0;

	double rd_out;

	for (int j = 0; j < ptrdimAfa[1]; j++){
		sumI+=ptrIcol[j];
		sumJ+=ptrJcol[j];
	}

	if (Alpha != 1){
	

		for (int j = 0; j < ptrdimAfa[1]; j++){
			term1+=pow((double) ptrIcol[j]/sumI,Alpha) * pow((double) ptrJcol[j]/sumJ,1.0-Alpha);
		}

		rd_out = (1.0/(Alpha-1.0)) * log(term1);
	}
	else {

		for (int j = 0; j < ptrdimAfa[1]; j++){
			temp += (double) ptrIcol[j]/sumI * log(((double) ptrIcol[j]/sumI)/((double) ptrJcol[j]/sumJ));
		}
		rd_out=temp;
	}


	return rd_out;
}
