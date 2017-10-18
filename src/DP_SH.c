#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double DP_SH(double* ptrFreq,int nFreq){

	double Out=0.0;

	for (int i = 0; i < nFreq; i++){
		if (ptrFreq[i]==0){ptrFreq[i]=0.00000000001;}
		Out+=-log(ptrFreq[i])*ptrFreq[i];
	}

	return Out;
}
