#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double DP_HT(double* ptrFreq,int nFreq,double Alpha,int ENS,double Sum){

	double HtOut=0.0;
	double Temp1=0.0;
	double Temp2=0.0;

	if (Alpha==1){Alpha=1-0.00000000001;}
	
	for (int i = 0; i < nFreq; i++){
		Temp1+=pow(ptrFreq[i],Alpha)/(1-pow(1-ptrFreq[i],Sum));
		Temp2+=			  ptrFreq[i]/(1-pow(1-ptrFreq[i],Sum));
	}

	HtOut=log(Temp1)/(1-Alpha)-log(Temp2)/(1-Alpha);

	if (ENS==1){HtOut=exp(HtOut);}

	return HtOut;
}
