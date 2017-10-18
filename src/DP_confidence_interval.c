#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

void DP_confidence_interval(double** ptrResults,int* ptrdimResults,double valCI,
	double* ptrOutMean,int* ptrdimOutMean,double* ptrOutMin,int* ptrdimOutMin,double* ptrOutMax,int* ptrdimOutMax){

	//////////////////////////////////////////////////////////////////////////
	/* 									OUTMEAN 							*/
	//////////////////////////////////////////////////////////////////////////

	double sumResult=0.0;

	for (int j = 0; j < ptrdimResults[1]; j++)
	{
		for (int i = 0; i < ptrdimResults[0]; i++)
		{
			sumResult+=ptrResults[i][j];
		}

		ptrOutMean[j]=sumResult/ptrdimResults[0];

		sumResult=0.0;
	}

	
	//////////////////////////////////////////////////////////////////////////
	/* 						OUTMIN  	AND		OUTMAX  					*/
	//////////////////////////////////////////////////////////////////////////

	int iter=0;
	int nSortVec=ptrdimResults[0];
	double* ptrSortVec; ptrSortVec=dvector(0,nSortVec-1);

	int Temp1;
	int Temp2;
	double swap;
	int is, js;

	Temp1= (int) ptrdimResults[0]*(1-valCI)/2;

	for (int j = 0; j < ptrdimResults[1]; j++){

		for (int i = 0; i < ptrdimResults[0]; i++){
			ptrSortVec[i]=ptrResults[i][iter];
		}

		//Sort CODE start
		for (is = 0; is < nSortVec; is++){
		    for (js = is +1; js < nSortVec; ++js){
		        if (ptrSortVec[is] > ptrSortVec[js]){
		            swap = ptrSortVec[is];
		            ptrSortVec[is] = ptrSortVec[js];
		            ptrSortVec[js] = swap;
		        }
		    }
		}
		//Sort CODE end

		Temp2=ptrdimResults[0]-1-Temp1;//index for ptrOutMax,

		ptrOutMin[j]=ptrOutMean[j]-ptrSortVec[Temp1];
		ptrOutMax[j]=ptrSortVec[Temp2]-ptrOutMean[j];


		iter++;
	}

	free_dvector(ptrSortVec,0,nSortVec-1);
}







