#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"


double OL_cvg(int* ptrCol,int* ptrdimAfa){

	int sumCol=0;
	int iter=0;

	for (int j = 0; j < ptrdimAfa[1]; j++){

		sumCol+=ptrCol[j];

		if (ptrCol[j]==1){iter++;}
	}

	return (1.0 - ((double) iter)/sumCol);
}
