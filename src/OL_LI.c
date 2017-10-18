#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_LI(int* ptrIcol,int* ptrJcol,int* ptrdimAfa){

	int Temp1=0;
	int Temp2=0;

	for (int j = 0; j < ptrdimAfa[1]; j++){

		if (ptrIcol[j]!=0 && ptrJcol[j]!=0){Temp1++;}

		if (ptrIcol[j]!=0){Temp2++;}
		if (ptrJcol[j]!=0){Temp2++;}
	}
	
	return 2*Temp1/(double)Temp2;
}
