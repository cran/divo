#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_SH(double* ptrVec,int nVec){	

	double out=0.0;

	for (int i = 0; i < nVec; i++)
	{
		if (ptrVec[i]==0){ptrVec[i]=0.00000000001;}

		out+= -log(ptrVec[i])*ptrVec[i];

	}

	return out;
}
