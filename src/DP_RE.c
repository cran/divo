#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"
#include "DP_SH.h"

double DP_RE(double* ptrFreq,int nFreq,double Alpha,int ENS){

	double Temp=0;
	double Re;

	if (Alpha==1){
		Re=DP_SH(ptrFreq,nFreq);
	}
	else{

		for (int i = 0; i < nFreq; i++){
			if (ptrFreq[i]==0){ptrFreq[i]=0.00000000001;}
			Temp+=pow(ptrFreq[i],Alpha);}

		Re=log(Temp)/(1-Alpha);
	}

	if (ENS==1){Re=exp(Re);}


	return Re;
}
