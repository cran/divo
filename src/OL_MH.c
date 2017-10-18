#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

double OL_MH(int* ptrIcol,int* ptrJcol,int* ptrdimAfa){

	long long int Temp1=0;
	int Temp2=0;
	int Temp3=0;
	long long int IcolSumSq=0;
	long long int JcolSumSq=0;
	long double Result=0;

	for (int j = 0; j < ptrdimAfa[1]; j++){
		Temp1+=(long long int)ptrIcol[j]*(long long int)ptrJcol[j];
		Temp2+= ptrIcol[j];
		Temp3+= ptrJcol[j];
		IcolSumSq+=(long long int)ptrIcol[j]*(long long int)ptrIcol[j];
		JcolSumSq+=(long long int)ptrJcol[j]*(long long int)ptrJcol[j];
	}

	//Get rid of this if you know that all columns contain at least one non zero values
	if(Temp2==0 || Temp3==0){error("MH.c Error: One of the columns is filled with zeros. Remove column.");}

	//Previous implementation, changed since it resulted in ver big numbers!
	//
	//long long int Num=0;
	//long long int Den=0;
	//Num= 2* Temp1* Temp2* Temp3;
	//Den= IcolSumSq*Temp3*Temp3 + JcolSumSq*Temp2*Temp2;
	

	//To look at individual variables:
	//
	//Result=(long double)Num/(long double)Den;
	//printf("Temp1: %llu \n",Temp1);
	//printf("Temp2: %llu \n",Temp2);
	//printf("Temp3: %llu \n",Temp3);
	//printf("IcolSumSq: %llu \n",IcolSumSq);
	//printf("JcolSumSq: %llu \n",JcolSumSq);
	//printf("Num: %llu \n",Num);
	//printf("Den: %llu \n",Den);
	//printf("Next \n");


	//In case numbers come too big, you will have to look into bignum implementation, for arbitrary number size handling!
	Result = 2 * (double long) Temp1 / (IcolSumSq* (double long)Temp3/Temp2 +JcolSumSq* (double long)Temp2/Temp3);


	return (double)Result;
}


