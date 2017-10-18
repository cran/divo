#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"


void OL_confidence_interval(double** ptrResults,int* ptrdimResults,double valCI,
	double* ptrOutMean,int* ptrdimOutMean,double* ptrOutMin,int* ptrdimOutMin,double* ptrOutMax,int* ptrdimOutMax){
	
	//////////////////////////////////////////////////////////////////////////
	/* 									OUTMEAN 							*/
	//////////////////////////////////////////////////////////////////////////

	double sumResult=0.0;
	int os1=1;//offset1
	int os2=1;//offset2


	for (int j = 0; j < ptrdimResults[1]; j++)
	{
		for (int i = 0; i < ptrdimResults[0]; i++)
		{
			sumResult+=ptrResults[i][j];
		}

		//Filling OutMean matrix--------------------
		if ((j+os2)/(ptrdimOutMean[0])==(j+os2)/((double)ptrdimOutMean[0])){
			os1++;
			os2+=os1;
		}//printf("%i ",j+os2);

		ptrOutMean[j+os2]=sumResult/ptrdimResults[0];

		//Resetting sumResult
		sumResult=0.0;
	}

	//////////////////////////////////////////////////////////////////////////
	/* 						OUTMIN  	AND		OUTMAX  					*/
	//////////////////////////////////////////////////////////////////////////

	double normResult=0.0;
	int nnormVec=ptrdimResults[0];
	double* ptrnormVec1; ptrnormVec1=dvector(0,nnormVec-1);
	double* ptrnormVec2; ptrnormVec2=dvector(0,nnormVec-1);

	valCI=(1-valCI)/2;

	for (int i = 0; i < ptrdimResults[0]; i++){
		for (int j = 0; j < ptrdimResults[1]; j++){

			normResult+=pow(ptrResults[i][j],2);
		}
		ptrnormVec1[i]=sqrt(normResult);
		ptrnormVec2[i]=sqrt(normResult);

		//Resetting normResult
		normResult=0;
	}

	//Sorting ptrnormVec2---- start
	double swap;
	int i, j;

	for (i = 0; i < nnormVec; i++){
	    for ( j = i +1; j < nnormVec; ++j){

	        if (ptrnormVec2[i] > ptrnormVec2[j]){
	            swap = ptrnormVec2[i];
	            ptrnormVec2[i] = ptrnormVec2[j];
	            ptrnormVec2[j] = swap;
	        }
	    }
	}
	//----------------------- end

	int minSample; minSample=(int)(nnormVec*valCI);
	int maxSample; maxSample=nnormVec-1-minSample;


	for (int i = 0; i < nnormVec; i++)
	{

		if (ptrnormVec2[minSample]==ptrnormVec1[i]){
			//printf("Min: %f",ptrnormVec1[i]); //Delete?
			os1=1; os2=1; //offsets
			for (int j = 0; j < ptrdimResults[1]; j++){
				//Filling OutMin matrix--------------------
				if ((j+os2)/(ptrdimOutMin[0])==(j+os2)/((double)ptrdimOutMin[0])){
					os1++;
					os2+=os1;}

				ptrOutMin[j+os2]=ptrResults[i][j];
			}
		}


		if (ptrnormVec2[maxSample]==ptrnormVec1[i]){
			//printf("Max: %f",ptrnormVec1[i]); //Delete?
			os1=1; os2=1; //offsets	
			for (int j = 0; j < ptrdimResults[1]; j++){
				//Filling OutMax matrix--------------------
				if ((j+os2)/(ptrdimOutMax[0])==(j+os2)/((double)ptrdimOutMax[0])){
					os1++;
					os2+=os1;}
				
				ptrOutMax[j+os2]=ptrResults[i][j];
			}
		}

	}

	free_dvector(ptrnormVec1,0,nnormVec-1);
	free_dvector(ptrnormVec2,0,nnormVec-1);
}

