#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_JI(int* ptrIcol,int* ptrJcol,int* ptrdimAfa){

	int Val_11=0;
	int Val_01=0;
	int	Val_10=0;

	for (int j = 0; j < ptrdimAfa[1]; j++){

		if (ptrIcol[j]!=0 && ptrJcol[j]!=0){Val_11++;}

		if (ptrIcol[j]==0 && ptrJcol[j]!=0){Val_01++;}

		if (ptrIcol[j]!=0 && ptrJcol[j]==0){Val_10++;}
	}

	return Val_11/(double)(Val_11+Val_01+Val_10);
}
