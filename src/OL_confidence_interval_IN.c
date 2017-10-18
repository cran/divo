#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

void OL_confidence_interval_IN(double* ptrResults1D,int nResults1D,double valCI,
	double* ptrOutMean,double* ptrOutMin,double* ptrOutMax){

	//////////////////////////////////////////////////////////////////////////
	/* 									OUTMEAN 							*/
	//////////////////////////////////////////////////////////////////////////

	double sum=0;

	for (int j = 0; j < nResults1D; j++){
		sum+=ptrResults1D[j];
	}

	ptrOutMean[0]=sum/nResults1D;

	//////////////////////////////////////////////////////////////////////////
	/* 						OUTMIN  	AND		OUTMAX  					*/
	//////////////////////////////////////////////////////////////////////////

	int nSorted=nResults1D;
	double* ptrSorted; ptrSorted=dvector(0,nSorted-1);

	for (int j = 0; j < nResults1D; j++){
		ptrSorted[j]=ptrResults1D[j];
	}


	//Sorting ptrSorted---- start
	double swap;
	int i, j;

	for (i = 0; i < nSorted; i++){
	    for ( j = i +1; j < nSorted; ++j){

	        if (ptrSorted[i] > ptrSorted[j]){
	            swap = ptrSorted[i];
	            ptrSorted[i] = ptrSorted[j];
	            ptrSorted[j] = swap;
	        }
	    }
	}
	//----------------------- end

	valCI=(1-valCI)/2;

	int tmp; tmp=((int)(nSorted*valCI));

	if (tmp>0){
		ptrOutMin[0]=ptrSorted[tmp];
		ptrOutMax[0]=ptrSorted[nSorted-1-tmp];
	}
	else{
		ptrOutMin[0]=ptrSorted[0];
		ptrOutMax[0]=ptrSorted[nSorted-1];
	}

	free_dvector(ptrSorted,0,nSorted-1);
}





