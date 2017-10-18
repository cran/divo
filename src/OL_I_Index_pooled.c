#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"
#include "OL_SH.h"


double OL_I_Index_pooled(int* ptrAfa1D,int nAfa1D,int** ptrAfa,int* ptrdimAfa,double Alpha){
	
	//summing along columns
	//summing along rows

	double* ptrCols; ptrCols=dvector(0,ptrdimAfa[0]-1); 
	double*	ptrRows; ptrRows=dvector(0,ptrdimAfa[1]-1);
	double* ptrTable1; ptrTable1=dvector(0,nAfa1D-1);
	double* ptrTable2; ptrTable2=dvector(0,nAfa1D-1);

	double I_ind;

	//--------------------------------------------Initialization-------------------------------------------

	//
	double Sum_dat=0;for (int i=0; i<nAfa1D ;i++){Sum_dat += ptrAfa1D[i];}

	//Making Cols
	for (int i = 0; i < ptrdimAfa[0]; i++)
	{
		ptrCols[i]=0.0;

		for (int j = 0; j < ptrdimAfa[1]; j++)
		{
			ptrCols[i]+=ptrAfa[i][j];
			ptrTable1[i*ptrdimAfa[1]+j]=ptrAfa[i][j]/Sum_dat;
		}
		ptrCols[i]/=Sum_dat;
	}

	//Making Rows
	for (int j = 0; j < ptrdimAfa[1]; j++)
	{
		ptrRows[j]=0.0;

		for (int i = 0; i < ptrdimAfa[0]; i++)
		{
			ptrRows[j]+=ptrAfa[i][j];
		}
		ptrRows[j]/=Sum_dat;
	}

	//-----------------------------------------------Execution-------------------------------------------

	if (Alpha==1)
	{
		I_ind=1-((OL_SH(ptrRows,ptrdimAfa[1])+OL_SH(ptrCols,ptrdimAfa[0])-OL_SH(ptrTable1,nAfa1D))/OL_SH(ptrCols,ptrdimAfa[0]));
	}
	else{

		double Fa;
		double P1;
		double Temp1=0.0;
		double Temp2=0.0;

		//Asigning values to Table2
		for (int i = 0; i < ptrdimAfa[0]; i++){
			for (int j = 0; j < ptrdimAfa[1]; j++){
				ptrTable2[i*ptrdimAfa[1]+j]=ptrRows[j]*ptrCols[i];
			}
			Temp1+=pow(ptrCols[i],2-Alpha);
		}

		for (int i = 0; i < nAfa1D; i++){Temp2+=pow(ptrTable1[i],Alpha) * pow(ptrTable2[i],1-Alpha);}

		Fa=(1/(Alpha-1))*log(Temp2);
		P1=Fa/((1/(Alpha-1))*log(Temp1));

		I_ind=1-P1;
	}

	free_dvector(ptrCols,0,ptrdimAfa[0]-1); 
	free_dvector(ptrRows,0,ptrdimAfa[1]-1);
	free_dvector(ptrTable1,0,nAfa1D-1);
	free_dvector(ptrTable2,0,nAfa1D-1);

	return I_ind;
}

