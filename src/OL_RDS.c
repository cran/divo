#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_RDS(int* ptrIcol,int* ptrJcol,int* ptrdimAfa,double Alpha){

	int sumI=0;
	int sumJ=0;

	double term1=0.0;
	double term2=0.0;

	double rd_out;

	if (Alpha != 1){
		
		for (int j = 0; j < ptrdimAfa[1]; j++){
			sumI+=ptrIcol[j];
			sumJ+=ptrJcol[j];
		}

		for (int j = 0; j < ptrdimAfa[1]; j++){
			term1+=pow((double) ptrIcol[j]/sumI,Alpha) * pow((double) ptrJcol[j]/sumJ,1.0-Alpha);
			term2+=pow((double) ptrJcol[j]/sumJ,Alpha) * pow((double) ptrIcol[j]/sumI,1.0-Alpha);
		}

		rd_out= 0.5 * (1.0/(Alpha-1.0)) * (log(term1)+log(term2));
	}
	else {
		rd_out=0.0;
	}


	return rd_out;
}


