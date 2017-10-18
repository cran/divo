#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"
#include "OL_SH.h"
#include "OL_SH2.h"


double OL_I_Index(int* ptrIcol,int* ptrJcol,int* ptrdimAfa,double Alpha){
	
	double* ptrCols; ptrCols=dvector(0,1); 
	double*	ptrRows; ptrRows=dvector(0,ptrdimAfa[1]-1);

	double Sum_dat=0;
	double I_ind;
	double Temp1=0.0;
	double Temp2=0.0;
	double Fa;
	double P1;



	//Evaluating Sum_dat
	for (int j = 0; j < ptrdimAfa[1]; j++){Sum_dat+=ptrIcol[j]+ptrJcol[j];}

	//Making ptrCols and ptrRows
	ptrCols[0]=0.0;ptrCols[1]=0.0;
	for (int j = 0; j < ptrdimAfa[1]; j++){
		ptrCols[0]+=ptrIcol[j];
		ptrCols[1]+=ptrJcol[j];
		ptrRows[j]=(ptrIcol[j]+ptrJcol[j])/Sum_dat;
	}
	ptrCols[0]/=Sum_dat;ptrCols[1]/=Sum_dat;


	//-----------------------------------------------Execution-------------------------------------------

	if (Alpha==1)
	{
		I_ind=1-((OL_SH(ptrRows,ptrdimAfa[1])+OL_SH(ptrCols,2)-OL_SH2(ptrIcol,ptrdimAfa[1],Sum_dat)-OL_SH2(ptrJcol,ptrdimAfa[1],Sum_dat))/OL_SH(ptrCols,2));
	}
	else{

		for (int j = 0; j < ptrdimAfa[1]; j++){
			Temp2+=
			(pow(ptrIcol[j]/Sum_dat,Alpha))*pow(ptrRows[j]*ptrCols[0],(1-Alpha))+ 
			(pow(ptrJcol[j]/Sum_dat,Alpha))*pow(ptrRows[j]*ptrCols[1],(1-Alpha));
		}

		Temp1=pow(ptrCols[0],(2-Alpha))+pow(ptrCols[1],(2-Alpha));
		
		Fa=(1/(Alpha-1))*log(Temp2);
		P1=Fa/((1/(Alpha-1))*log(Temp1));
		I_ind=1-P1;

	}

	free_dvector(ptrCols,0,1);
	free_dvector(ptrRows,0,ptrdimAfa[1]-1);

	return I_ind;
}

//for (int ii = 0; ii < ptrdimAfa[1]; ii++)
//{
//	printf("%i ",ptrJcol[ii]);
//}
