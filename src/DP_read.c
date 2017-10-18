#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"
#include "DP_single_population.h"
#include "DP_single_population_cvg.h"
#include "DP_single_population_HT.h"
#include "DP_single_population_HT_cvg.h"
#include "DP_draw_arrays.h"
#include "DP_confidence_interval.h"
#include "DP_saveBootstrap.h"
#include "DP_cvg.h"



int DP_read(SEXP X, SEXP Test, SEXP Alpha, SEXP Resample, SEXP CI, 
		  SEXP Faaa, SEXP PlugIn, SEXP Sizee, SEXP AlphaPro, SEXP CVG, SEXP BS, 
		  SEXP OutMean, SEXP OutMin, SEXP OutMax, SEXP OutCvg) {

	//Declaring an integer to keep track of how many Protected variables there are
	int nprot=0;

	//////////////////////////////////////////////////////////////////////////
	/* 							READING IN VARIABLES						*/
	//////////////////////////////////////////////////////////////////////////

	//Reading in X variable---------------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(X = AS_INTEGER(X)); nprot++;
	//Defining pointers and assigning them to variables
	int* ptrX; ptrX = INTEGER(X);
	//Defining pointers and assigning them to variables
	SEXP dimX; PROTECT(dimX = allocVector(INTSXP,2)); nprot++;
	int* ptrdimX; ptrdimX=INTEGER(dimX);
	ptrdimX[0]=INTEGER(GET_DIM(X))[1]; ptrdimX[1]=INTEGER(GET_DIM(X))[0];
	//length of X array
	int nX; nX = ptrdimX[0] * ptrdimX[1];

	//Reading in test variable------------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(Test = AS_CHARACTER(Test)); nprot++;
	//Defining pointer
	char *Testptr[1];
	Testptr[0] = R_alloc(strlen(CHAR(STRING_ELT(Test, 0))),sizeof(char));
	//Assigning it to variable
	strcpy(Testptr[0], CHAR(STRING_ELT(Test, 0)));

	//Reading in Alpha variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(Alpha = AS_NUMERIC(Alpha)); nprot++;    
	//Defining pointers and assigning them to variables
	double *ptrAlpha; ptrAlpha = REAL(Alpha);

	//Reading in Resample variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(Resample = AS_INTEGER(Resample)); nprot++;    
	//Defining pointers and assigning them to variables
	int* ptrResample; ptrResample = INTEGER(Resample);

	//Reading in CI variable-----------------------------------------------
	PROTECT(CI = AS_NUMERIC(CI)); nprot++;    
	double *ptrCI; ptrCI = REAL(CI);

	//Reading in Faaa variable------------------------------------------------
	PROTECT(Faaa = AS_CHARACTER(Faaa)); nprot++;
	char *ptrFaaa[1];
	ptrFaaa[0] = R_alloc(strlen(CHAR(STRING_ELT(Faaa, 0))),sizeof(char));
	strcpy(ptrFaaa[0], CHAR(STRING_ELT(Faaa, 0)));

	//Reading in PlugIn variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(PlugIn = AS_LOGICAL(PlugIn)); nprot++;    
	//Defining pointers and assigning them to variables
	int *ptrPlugIn; ptrPlugIn = LOGICAL(PlugIn);

	//Reading in Sizee variable-----------------------------------------------
	PROTECT(Sizee = AS_NUMERIC(Sizee)); nprot++;    
	double *ptrSizee; ptrSizee = REAL(Sizee);

	//Reading in AlphaPro variable-----------------------------------------------
	PROTECT(AlphaPro = AS_NUMERIC(AlphaPro)); nprot++;
	double *ptrAlphaPro; ptrAlphaPro = REAL(AlphaPro);
	int nAlphaPro; nAlphaPro=length(AlphaPro);

	//Reading in CVG variable-------------------------------------------------
	PROTECT(CVG = AS_LOGICAL(CVG)); nprot++;    
	int *ptrCVG; ptrCVG = LOGICAL(CVG);

	//Reading in BS variable--------------------------------------------------
	PROTECT(BS = AS_LOGICAL(BS)); nprot++;    
	int *ptrBS; ptrBS = LOGICAL(BS);

	//Reading in OutMean variable---------------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(OutMean = AS_NUMERIC(OutMean)); nprot++;
	//Defining pointers and assigning them to variables
	double* ptrOutMean; ptrOutMean = REAL(OutMean);
	//Defining pointers and assigning them to variables
	SEXP dimOutMean; PROTECT(dimOutMean = allocVector(INTSXP,2)); nprot++;
	int* ptrdimOutMean; ptrdimOutMean=INTEGER(dimOutMean);
	ptrdimOutMean[0]=INTEGER(GET_DIM(OutMean))[1]; ptrdimOutMean[1]=INTEGER(GET_DIM(OutMean))[0];
	//length of X array
	int nOutMean; nOutMean = ptrdimOutMean[0] * ptrdimOutMean[1];

	//Reading in OutMin variable---------------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(OutMin = AS_NUMERIC(OutMin)); nprot++;
	//Defining pointers and assigning them to variables
	double* ptrOutMin; ptrOutMin = REAL(OutMin);
	//Defining pointers and assigning them to variables
	SEXP dimOutMin; PROTECT(dimOutMin = allocVector(INTSXP,2)); nprot++;
	int* ptrdimOutMin; ptrdimOutMin=INTEGER(dimOutMin);
	ptrdimOutMin[0]=INTEGER(GET_DIM(OutMin))[1]; ptrdimOutMin[1]=INTEGER(GET_DIM(OutMin))[0];
	//length of X array
	int nOutMin; nOutMin = ptrdimOutMin[0] * ptrdimOutMin[1];

	//Reading in OutMax variable---------------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(OutMax = AS_NUMERIC(OutMax)); nprot++;
	//Defining pointers and assigning them to variables
	double* ptrOutMax; ptrOutMax = REAL(OutMax);
	//Defining pointers and assigning them to variables
	SEXP dimOutMax; PROTECT(dimOutMax = allocVector(INTSXP,2)); nprot++;
	int* ptrdimOutMax; ptrdimOutMax=INTEGER(dimOutMax);
	ptrdimOutMax[0]=INTEGER(GET_DIM(OutMax))[1]; ptrdimOutMax[1]=INTEGER(GET_DIM(OutMax))[0];
	//length of X array
	int nOutMax; nOutMax = ptrdimOutMax[0] * ptrdimOutMax[1];

	//Reading in OutCvg variable-----------------------------------------------
	PROTECT(OutCvg = AS_NUMERIC(OutCvg)); nprot++;
	double *ptrOutCvg; ptrOutCvg = REAL(OutCvg);
	int nOutCvg; nOutCvg=length(OutCvg);

	//REASIGNING VECTOR-----------------------------------------------------
	//----------------------------------------------------------------------

	//Making Afa1D----------------------------------------------------------
	int nAfa1D=nX;
	
	SEXP Afa1D;
	PROTECT(Afa1D = allocMatrix(INTSXP, ptrdimX[1], ptrdimX[0])); nprot++;
	int* ptrAfa1D; ptrAfa1D = INTEGER(Afa1D);
	for (int i = 0; i <nAfa1D; i++){ptrAfa1D[i]=ptrX[i];}

	//Making Afa------------------------------------------------------------
	int* ptrdimAfa; ptrdimAfa=ptrdimX;

	int** ptrAfa;
	ptrAfa=(int**) malloc((unsigned) ptrdimAfa[0]*sizeof(int*));
	for (int i = 0; i < ptrdimAfa[0]; i++){ptrAfa[i]=(ptrAfa1D+i*ptrdimAfa[1]);}


	//////////////////////////////////////////////////////////////////////////
	/* 							OUTPUT ALLOCATION							*/
	//////////////////////////////////////////////////////////////////////////

	int nResults1D;	
	if (*ptrPlugIn){nResults1D=nOutMean;}
	else{nResults1D=nOutMean*ptrResample[0];}

	SEXP Results1D;
	PROTECT(Results1D = allocVector(REALSXP,nResults1D)); nprot++;
	double *ptrResults1D; ptrResults1D = REAL(Results1D);

	//Making Results------------------------------------------------------------
	SEXP dimResults; PROTECT(dimResults = allocVector(INTSXP,3)); nprot++;
	int* ptrdimResults; ptrdimResults=INTEGER(dimResults);

	//Getting correct ptrdimResults for Plugin, TRUE False
	if (*ptrPlugIn)  {ptrdimResults[0]=1;ptrdimResults[1]=nOutMean;}
	else{ptrdimResults[0]=ptrResample[0];ptrdimResults[1]=nOutMean;}

	double** ptrResults;
	ptrResults=(double**) malloc((unsigned) ptrdimResults[0]*sizeof(double*));
	for (int i = 0; i < ptrdimResults[0]; i++){ptrResults[i]=(ptrResults1D+i*ptrdimResults[1]);}

	double* iterptrResults=ptrResults1D;

	//////////////////////////////////////////////////////////////////////////
	/* 								EXECUTION								*/
	//////////////////////////////////////////////////////////////////////////


	// PLUGIN: TRUE----------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	if (*ptrPlugIn) {

		//Doesn't matter if X or XT is chosen, sum is the same
		double Sum_dat=0;for (int i=0; i<nX ;i++){Sum_dat += ptrX[i];}

		//METHOD: DP RE----------------------------------------------------------------
		//-----------------------------------------------------------------------------
		if ( strcmp(Testptr[0],"DP") == 0 && strcmp(ptrFaaa[0],"RE") == 0) {
			int ENS=0;
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int j = 0; j < ptrdimAfa[0]; j++){
					ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
					DP_single_population_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
			else{
				for (int j = 0; j < ptrdimAfa[0]; j++){
					DP_single_population(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
		}
		//METHOD: ENS RE----------------------------------------------------------------
		//------------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"ENS") == 0 && strcmp(ptrFaaa[0],"RE") == 0) {
			int ENS=1;
			//CVG: ---------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int j = 0; j < ptrdimAfa[0]; j++){
					ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
					DP_single_population_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
			else{
				for (int j = 0; j < ptrdimAfa[0]; j++){
					DP_single_population(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
		}
		//METHOD: DP HT----------------------------------------------------------------
		//-----------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"DP") == 0 && strcmp(ptrFaaa[0],"HT") == 0) {
			int ENS=0;
			//CVG: ---------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int j = 0; j < ptrdimAfa[0]; j++){
					ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
					DP_single_population_HT_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
			else{
				for (int j = 0; j < ptrdimAfa[0]; j++){
					DP_single_population_HT(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
		}
		//METHOD: ENS HT----------------------------------------------------------------
		//------------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"ENS") == 0 && strcmp(ptrFaaa[0],"HT") == 0) {
			int ENS=1;
			//CVG: ---------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int j = 0; j < ptrdimAfa[0]; j++){
					ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
					DP_single_population_HT_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
			else{
				for (int j = 0; j < ptrdimAfa[0]; j++){
					DP_single_population_HT(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
					iterptrResults+=nAlphaPro;
				}
			}
		}
		//----------------------------------------------------------------
		else {error("Method String didn't pass properly");}

	}
	// PLUGIN: FALSE---------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	else {

		//METHOD: DP RE----------------------------------------------------------------
		//-----------------------------------------------------------------------------
		if ( strcmp(Testptr[0],"DP") == 0 && strcmp(ptrFaaa[0],"RE") == 0) {
			for (int i = 0; i < ptrResample[0]; i++){
				DP_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){DP_saveBootstrap(ptrAfa,ptrdimAfa);}
				int ENS=0;
				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {
					for (int j = 0; j < ptrdimAfa[0]; j++){
						ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
						DP_single_population_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
				else{
					for (int j = 0; j < ptrdimAfa[0]; j++){
						DP_single_population(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
			}
		}
		//METHOD: ENS RE----------------------------------------------------------------
		//------------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"ENS") == 0 && strcmp(ptrFaaa[0],"RE") == 0) {
			for (int i = 0; i < ptrResample[0]; i++){
				DP_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){DP_saveBootstrap(ptrAfa,ptrdimAfa);}
				int ENS=1;
				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {
					for (int j = 0; j < ptrdimAfa[0]; j++){
						ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
						DP_single_population_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
				else{
					for (int j = 0; j < ptrdimAfa[0]; j++){
						DP_single_population(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
			}
		}
		//METHOD: DP HT----------------------------------------------------------------
		//-----------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"DP") == 0 && strcmp(ptrFaaa[0],"HT") == 0) {
			for (int i = 0; i < ptrResample[0]; i++){
				DP_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){DP_saveBootstrap(ptrAfa,ptrdimAfa);}
				int ENS=0;
				//CVG: ---------------------------------------------------------------------
				if (ptrCVG[0]) {
					for (int j = 0; j < ptrdimAfa[0]; j++){
						ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
						DP_single_population_HT_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
				else{
					for (int j = 0; j < ptrdimAfa[0]; j++){
						DP_single_population_HT(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}			
			}
		}
		//METHOD: ENS HT----------------------------------------------------------------
		//------------------------------------------------------------------------------
		else if (strcmp(Testptr[0],"ENS") == 0 && strcmp(ptrFaaa[0],"HT") == 0) {
			for (int i = 0; i < ptrResample[0]; i++){
				DP_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){DP_saveBootstrap(ptrAfa,ptrdimAfa);}		
				int ENS=1;
				//CVG: ---------------------------------------------------------------------
				if (ptrCVG[0]) {
					for (int j = 0; j < ptrdimAfa[0]; j++){
						ptrOutCvg[j]=DP_cvg(ptrAfa[j],ptrdimAfa);
						DP_single_population_HT_cvg(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,ptrOutCvg[j],iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}
				else{
					for (int j = 0; j < ptrdimAfa[0]; j++){
						DP_single_population_HT(ptrAfa[j],ptrdimAfa,ptrAlphaPro,nAlphaPro,ENS,iterptrResults);
						iterptrResults+=nAlphaPro;
					}
				}	
			}
		}
		//METHOD: Error-----------------------------------------------------------------
		//------------------------------------------------------------------------------
		else {error("Method String didn't pass properly");}


	}


	//////////////////////////////////////////////////////////////////////////
	/* 							RESULTS ANALYSIS							*/
	//////////////////////////////////////////////////////////////////////////	


	if (*ptrPlugIn){
		for (int ij = 0; ij < nOutMean; ij++){ptrOutMean[ij]=ptrResults1D[ij];}
	}
	
	else{
		
		double valCI=ptrCI[0];
		DP_confidence_interval(ptrResults,ptrdimResults,valCI,
	ptrOutMean,ptrdimOutMean,ptrOutMin,ptrdimOutMin,ptrOutMax,ptrdimOutMax);

	}

	//////////////////////////////////////////////////////////////////////////
	/* 							FREE MEMORY 								*/
	//////////////////////////////////////////////////////////////////////////	


	free(ptrAfa);
	free(ptrResults);

	return(nprot);

}





