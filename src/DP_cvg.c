#include <R.h>
#include <Rmath.h>
#include <stdio.h>

double DP_cvg(int* ptrCol,int* ptrdimAfa){

	double Sum=0.0;
	int Count=1;//Value you want to count
	int nCount=0;

	double valOutCvg;

	for (int j = 0; j < ptrdimAfa[1]; j++) {
		Sum+=ptrCol[j];
		if (ptrCol[j] == Count){nCount++;}
	}

	valOutCvg=1.0-nCount/Sum;

	return valOutCvg;
}
