#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_SH2(int* ptrVec,int nVec,double Sum_dat){	

	double out=0.0;

	for (int i = 0; i < nVec; i++)
	{
		if (ptrVec[i]!=0){
			out+= -log(ptrVec[i]/Sum_dat)*ptrVec[i]/Sum_dat;
		}	
	}
	return out;
}
