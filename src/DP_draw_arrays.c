#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"


void DP_draw_arrays(int* ptrX,int nAfa1D,double* ptrSizee,int* ptrAfa1D){

	double sum=0;
	double sum_dat=0;
	int sum_dat2;

	//Making  ptrProb
	double* ptrProb; ptrProb=dvector(0,nAfa1D-1);//(1)

	//Calculating sum
	for (int i=0; i<nAfa1D ;i++){
		sum = sum + ptrX[i];
		//sum += ptrX_[i]
	}

	//multiply by ptrSizee;
	sum_dat=sum * ptrSizee[0];
	//Should I round instead???, currently truncating!!!
	sum_dat2 = (int)(sum_dat);

	//Calculating probabilities
	for (int i=0; i<nAfa1D;i++){
		ptrProb[i] = ptrX[i]/sum;
	}

	
	GetRNGstate();//generates random number seed
	rmultinom(sum_dat2,ptrProb,nAfa1D,ptrAfa1D);
	PutRNGstate();


	free_dvector(ptrProb,0,nAfa1D-1);//(1)
}

