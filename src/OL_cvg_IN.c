#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"


double OL_cvg_IN(int* ptrAfa1D,int nAfa1D){

	int sum=0;
	int iter=0;

	for (int j = 0; j < nAfa1D; j++){

		sum+=ptrAfa1D[j];

		if (ptrAfa1D[j]==1){iter++;}
	}

	return (1 - ((double) iter)/sum);
}
