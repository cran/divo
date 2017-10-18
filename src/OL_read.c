#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"
#include "OL_I_Index_pooled.h"
#include "OL_I_Index.h"
#include "OL_MH.h"
#include "OL_PG.h"
#include "OL_PG_HT.h"
#include "OL_LI.h"
#include "OL_JI.h"
#include "OL_RDS.h"
#include "OL_RD.h"
#include "OL_draw_arrays.h"
#include "OL_saveBootstrap.h"

#include "OL_confidence_interval_IN.h"
#include "OL_cvg_IN.h"
#include "OL_confidence_interval.h"
#include "OL_cvg.h"

int OL_read(SEXP X, SEXP Test, SEXP Alpha, SEXP Resample, SEXP CI, 
		  SEXP Faaa, SEXP PlugIn, SEXP Sizee, SEXP Beta, SEXP CVG, SEXP BS, 
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
	//length of Alpha 
	//int nAlpha = length(Alpha);

	//Reading in Resample variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(Resample = AS_INTEGER(Resample)); nprot++;    
	//Defining pointers and assigning them to variables
	int* ptrResample; ptrResample = INTEGER(Resample);
	//length of Resample
	//int nResample = length(Resample);

	//Reading in CI variable-----------------------------------------------
	PROTECT(CI = AS_NUMERIC(CI)); nprot++;    
	double *ptrCI; ptrCI = REAL(CI);
	//int nCI = length(CI);

	//Reading in Faaa variable------------------------------------------------
	PROTECT(Faaa = AS_CHARACTER(Faaa)); nprot++;
	char *Faaaptr[1];
	Faaaptr[0] = R_alloc(strlen(CHAR(STRING_ELT(Faaa, 0))),sizeof(char));
	strcpy(Faaaptr[0], CHAR(STRING_ELT(Faaa, 0)));

	//Reading in PlugIn variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(PlugIn = AS_LOGICAL(PlugIn)); nprot++;    
	//Defining pointers and assigning them to variables
	int *ptrPlugIn; ptrPlugIn = LOGICAL(PlugIn);
	//length of PlugIn 
	//int nPlugIn = length(PlugIn);

	//Reading in Sizee variable-----------------------------------------------
	PROTECT(Sizee = AS_NUMERIC(Sizee)); nprot++;    
	double *ptrSizee; ptrSizee = REAL(Sizee);
	//int nSizee = length(Sizee);

	//Reading in Beta variable------------------------------------------------
	//PROTECT(Beta = AS_CHARACTER(Beta)); nprot++;
	//char *Betaptr[1];
	//Betaptr[0] = R_alloc(strlen(CHAR(STRING_ELT(Beta, 0))),sizeof(char));
	//strcpy(Betaptr[0], CHAR(STRING_ELT(Beta, 0)));
	//Reading in Alpha variable-----------------------------------------------
	//Protect variables so that R doesn't delete them
	PROTECT(Beta = AS_NUMERIC(Beta)); nprot++;    
	//Defining pointers and assigning them to variables
	double *ptrBeta; ptrBeta = REAL(Beta);

	//Reading in CVG variable-------------------------------------------------
	PROTECT(CVG = AS_LOGICAL(CVG)); nprot++;    
	int *ptrCVG; ptrCVG = LOGICAL(CVG);
	//int nCVG = length(CVG);

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

	//////////////////////////////////////////////////////////////////////////
	/* 							OUTPUT ALLOCATION							*/
	//////////////////////////////////////////////////////////////////////////

	//Results------------------------------------------------------------------
	//-------------------------------------------------------------------------
	int nResults1D;
	if ( strcmp(Testptr[0],"IN") == 0) {
		if (*ptrPlugIn){nResults1D=1;}
		else{nResults1D=ptrResample[0];}
	}
	else{
		if (*ptrPlugIn){nResults1D=(ptrdimX[0]*(ptrdimX[0]-1))/2;}
		else{nResults1D=((ptrdimX[0]*(ptrdimX[0]-1))/2)*ptrResample[0];}
	}
	
	SEXP Results1D;
	PROTECT(Results1D = allocVector(REALSXP,nResults1D)); nprot++;
	double *ptrResults1D; ptrResults1D = REAL(Results1D);

	//Making Results-----------------------------------------------------------
	SEXP dimResults; PROTECT(dimResults = allocVector(INTSXP,2)); nprot++;
	int* ptrdimResults; ptrdimResults=INTEGER(dimResults); 

	//Getting correct ptrdimResults for corresponding Tests and Plugins
	if ( strcmp(Testptr[0],"IN") == 0) {
		if (*ptrPlugIn){ptrdimResults[0]=1;ptrdimResults[1]=1;}
		else{ptrdimResults[0]=ptrResample[0];ptrdimResults[1]=1;}
	}
	else{
		if (*ptrPlugIn){ptrdimResults[0]=1;ptrdimResults[1]=(ptrdimX[0]*(ptrdimX[0]-1))/2;}
		else{ptrdimResults[0]=ptrResample[0];ptrdimResults[1]=(ptrdimX[0]*(ptrdimX[0]-1))/2;}
	}

	double** ptrResults;
	ptrResults=(double**) malloc((unsigned) ptrdimResults[0]*sizeof(double*));
	for (int i = 0; i < ptrdimResults[0]; i++){ptrResults[i]=(ptrResults1D+i*ptrdimResults[1]);}

	//ResCvg-------------------------------------------------------------------
	//-------------------------------------------------------------------------
	
	//Note ResCvg is not necessary for "IN" function but is constructed anyway

	int nResCvg1D;
	if ( strcmp(Testptr[0],"IN") == 0) {
		if (*ptrPlugIn){nResCvg1D=1;}
		else{nResCvg1D=ptrResample[0];}
	}
	else{
		if (*ptrPlugIn){nResCvg1D=ptrdimX[0];}
		else{nResCvg1D=ptrdimX[0]*ptrResample[0];}
	}


	SEXP ResCvg1D;
	PROTECT(ResCvg1D = allocVector(REALSXP,nResCvg1D)); nprot++;
	double *ptrResCvg1D; ptrResCvg1D = REAL(ResCvg1D);

	//Making Results-----------------------------------------------------------
	SEXP dimResCvg; PROTECT(dimResCvg = allocVector(INTSXP,2)); nprot++;
	int* ptrdimResCvg; ptrdimResCvg=INTEGER(dimResCvg); 

	//Getting correct ptrdimCvg for corresponding Tests and Plugins
	if ( strcmp(Testptr[0],"IN") == 0) {	
		if (*ptrPlugIn){ptrdimResCvg[0]=1;ptrdimResCvg[1]=1;}
		else{ptrdimResCvg[0]=ptrResample[0];ptrdimResCvg[1]=1;}
	}
	else{
		if (*ptrPlugIn){ptrdimResCvg[0]=1;ptrdimResCvg[1]=ptrdimX[0];}
		else{ptrdimResCvg[0]=ptrResample[0];ptrdimResCvg[1]=ptrdimX[0];}		
	}


	double** ptrResCvg;
	ptrResCvg=(double**) malloc((unsigned) ptrdimResCvg[0]*sizeof(double*));
	for (int i = 0; i < ptrdimResCvg[0]; i++){ptrResCvg[i]=(ptrResCvg1D+i*ptrdimResCvg[1]);}


	//////////////////////////////////////////////////////////////////////////
	/* 								EXECUTION								*/
	//////////////////////////////////////////////////////////////////////////
	
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
	
	if (*ptrPlugIn) {

		//METHOD: IN-------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_29, everything working
		if (strcmp(Testptr[0],"IN") == 0) {
			int iter=0;
			int iterCvg=0;

			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {//Coverage

				ptrResCvg1D[iterCvg]=OL_cvg_IN(ptrAfa1D,nAfa1D);
				iterCvg++;

				double valAlpha; valAlpha=ptrResCvg1D[0];
				ptrResults1D[iter]=OL_I_Index_pooled(ptrAfa1D,nAfa1D,ptrAfa,ptrdimAfa,valAlpha);
				iter++;
			}
			else{//No Coverage
				ptrResults1D[iter]=OL_I_Index_pooled(ptrAfa1D,nAfa1D,ptrAfa,ptrdimAfa,ptrAlpha[0]);
				iter++;
			}
		}
		//METHOD: INP------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, everything working except of CVG, which if TRUE fucks up everything
		else if (strcmp(Testptr[0],"INP") == 0) {
			int iter=0;
			int iterCvg=0;			
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {//Coverage
				for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
					ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
					iterCvg++;
				}


				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						double valAlpha; valAlpha=(ptrResCvg[0][ic]+ptrResCvg[0][jc])/2;
						ptrResults1D[iter]=OL_I_Index(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
						iter++;
					}
				}
			}
			else {//No Coverage
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_I_Index(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
						iter++;
					}
				}				
			}
		}
		//METHOD: MH-------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, only works with very small data set, otherwise not		
		else if (strcmp(Testptr[0],"MH") == 0) {
			//For CVG=FALSE
			int iter=0;
			for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
				for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
					ptrResults1D[iter]=OL_MH(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
					iter++;
				}
			}
		}
		//METHOD: PG-------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, everything working except of CVG, which if TRUE fucks up everything		
		else if (strcmp(Testptr[0],"PG") == 0) {
			int iter=0;
			int iterCvg=0;
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
					ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
					iterCvg++;
				}
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_PG(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrResCvg[0][ic],ptrResCvg[0][jc]);
						iter++;
					}
				}
			}
			else {
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_PG(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0],ptrBeta[0]);
						iter++;
					}
				}				
			}
		}
		//METHOD: PG_HT----------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, everything working except of CVG, which if TRUE fucks up everything		
		else if (strcmp(Testptr[0],"PG_HT") == 0) {
			int iter=0;
			int iterCvg=0;
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
					ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
					iterCvg++;
				}			
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_PG_HT(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrResCvg[0][ic],ptrResCvg[0][jc]);
						iter++;
					}
				}
			}
			else{
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_PG_HT(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0],ptrBeta[0]);
						iter++;
					}
				}				
			}
		}
		//METHOD: RD-------------------------------------------------------------------
		//-----------------------------------------------------------------------------	
		//Checked:2017_07_30, This function gives me an error, even with small data set. Compare with RDS as this one is working	
		else if (strcmp(Testptr[0],"RD") == 0) {
			int iter=0;
			int iterCvg=0;
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
					ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
					iterCvg++;
				}									
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						double valAlpha; valAlpha=(ptrResCvg[0][ic]+ptrResCvg[0][jc])/2;
						ptrResults1D[iter]=OL_RD(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
						iter++;
					}
				}
			}
			else{
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_RD(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
						iter++;
					}
				}
			}
		}
		//METHOD: RDS------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, everything working except of CVG, which if TRUE fucks up everything			
		else if (strcmp(Testptr[0],"RDS") == 0) {
			int iter=0;
			int iterCvg=0;
			//CVG: --------------------------------------------------------------------
			if (ptrCVG[0]) {
				for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
					ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
					iterCvg++;
				}						
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						double valAlpha; valAlpha=(ptrResCvg[0][ic]+ptrResCvg[0][jc])/2;
						ptrResults1D[iter]=OL_RDS(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
						iter++;
					}
				}
			}
			else{
				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_RDS(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
						iter++;
					}
				}
			}
		}
		//METHOD: LI-------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_30, everything is working!		
		else if (strcmp(Testptr[0],"LI") == 0) {
			int iter=0; 
			for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
				for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
					ptrResults1D[iter]=OL_LI(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
					iter++;
				}
			}
		}
		//METHOD: JI-------------------------------------------------------------------
		//-----------------------------------------------------------------------------	
		//Checked:2017_07_30, everything is working!	
		else if (strcmp(Testptr[0],"JI") == 0) {
			int iter=0;
			for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
				for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
					ptrResults1D[iter]=OL_JI(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
					iter++;
				}
			}
		}
		else {error("Method String didn't pass properly");}
	}
	// PLUGIN: FALSE---------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------
	else {
		//METHOD: IN-------------------------------------------------------------------
		//-----------------------------------------------------------------------------
		//Checked:2017_07_29
		if ( strcmp(Testptr[0],"IN") == 0) {
			int iter=0;
			int iterCvg=0;
			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);} //Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage

					ptrResCvg1D[iterCvg]=OL_cvg_IN(ptrAfa1D,nAfa1D);
					iterCvg++;
	
					double valAlpha; valAlpha=ptrResCvg1D[0];
					ptrResults1D[iter]=OL_I_Index_pooled(ptrAfa1D,nAfa1D,ptrAfa,ptrdimAfa,valAlpha);
					iter++;
				}
				else{//No Coverage
					ptrResults1D[iter]=OL_I_Index_pooled(ptrAfa1D,nAfa1D,ptrAfa,ptrdimAfa,ptrAlpha[0]);
					iter++;
				}
			}
		}
		//METHOD: INP------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"INP") == 0) {
			int iter=0;
			int iterCvg=0;
			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);} //Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage
					for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
						ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
						iterCvg++;
					}
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							double valAlpha; valAlpha=(ptrResCvg[i][ic]+ptrResCvg[i][jc])/2;							
							ptrResults1D[iter]=OL_I_Index(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
							iter++;
						}
					}
				}
				else {//No Coverage
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_I_Index(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
							iter++;
						}
					}				
				}
			}
		}
		//METHOD: MH-------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"MH") == 0) {
			int iter=0;
			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);} //Bootstrap

				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_MH(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
						iter++;
					}
				}

			}
		}
		//METHOD: PG-------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"PG") == 0) {
			int iter=0;
			int iterCvg=0;

			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage
					for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
						ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
						iterCvg++;
					}
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_PG(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrResCvg[i][ic],ptrResCvg[i][jc]);
							iter++;
						}
					}
				}
				else {//No Coverage
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_PG(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0],ptrBeta[0]);
							iter++;
						}
					}				
				}
			}
		}
		//METHOD: PG_HT----------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"PG_HT") == 0) {
			int iter=0;
			int iterCvg=0;

			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage
					for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
						ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
						iterCvg++;
					}			
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_PG_HT(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrResCvg[i][ic],ptrResCvg[i][jc]);
							iter++;
						}
					}
				}
				else{//No Coverage
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_PG_HT(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0],ptrBeta[0]);
							iter++;
						}
					}				
				}
			}
		}
		//METHOD: RD-------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"RD") == 0) {
			int iter=0;
			int iterCvg=0;

			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage
					for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
						ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
						iterCvg++;
					}									
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							double valAlpha; valAlpha=(ptrResCvg[i][ic]+ptrResCvg[i][jc])/2;
							ptrResults1D[iter]=OL_RD(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
							iter++;
						}
					}
				}
				else{//No Coverage
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_RD(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
							iter++;
						}
					}
				}
			}
		}
		//METHOD: RDS------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"RDS") == 0) {
			int iter=0;
			int iterCvg=0;

			//Resample: --------------------------------------------------------------------			
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				//CVG: --------------------------------------------------------------------
				if (ptrCVG[0]) {//Coverage
					for (int ii = 0; ii < ptrdimAfa[0]; ii++){//Note that index ends before ptrdimAfa[0] unlike the other loops that compare two columes, which end before ptrdimAfa[0]-1 and ptrdimAfa[0]
						ptrResCvg1D[iterCvg]=OL_cvg(ptrAfa[ii],ptrdimAfa);
						iterCvg++;
					}						
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							double valAlpha; valAlpha=(ptrResCvg[i][ic]+ptrResCvg[i][jc])/2;
							ptrResults1D[iter]=OL_RDS(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,valAlpha);
							iter++;
						}
					}
				}
				else{//No Coverage
					for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
						for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
							ptrResults1D[iter]=OL_RDS(ptrAfa[ic],ptrAfa[jc],ptrdimAfa,ptrAlpha[0]);
							iter++;
						}
					}
				}
			}
		}
		//METHOD: LI-------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"LI") == 0) {
			int iter=0;

			//Resample: --------------------------------------------------------------------			
			for (int i = 0; i < ptrResample[0]; i++){

				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_LI(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
						iter++;
					}
				}
			}
		}
		//METHOD: JI-------------------------------------------------------------------
		//-----------------------------------------------------------------------------		
		else if (strcmp(Testptr[0],"JI") == 0) {
			int iter=0;

			//Resample: --------------------------------------------------------------------
			for (int i = 0; i < ptrResample[0]; i++){
			
				OL_draw_arrays(ptrX,nAfa1D,ptrSizee,ptrAfa1D);
				if (*ptrBS){OL_saveBootstrap(ptrAfa,ptrdimAfa);}//Bootstrap

				for (int ic = 0; ic < ptrdimAfa[0]-1; ic++){
					for (int jc = ic+1; jc < ptrdimAfa[0]; jc++){
						ptrResults1D[iter]=OL_JI(ptrAfa[ic],ptrAfa[jc],ptrdimAfa);
						iter++;
					}
				}

			}
		}
		else {error("Method String didn't pass properly");}

	}

	//////////////////////////////////////////////////////////////////////////
	/* 							RESULTS ANALYSIS							*/
	//////////////////////////////////////////////////////////////////////////	

	if ( strcmp(Testptr[0],"IN") == 0) {
		double valCI=ptrCI[0];
		OL_confidence_interval_IN(ptrResults1D,nResults1D,valCI,
			ptrOutMean,ptrOutMin,ptrOutMax);

		//If ptrOutCvg is wanted for IN , then fill it with values here from ptrResCvg1D 
	}
	else{

		if (ptrCVG[0]) {
			if (strcmp(Testptr[0],"INP") == 0 || strcmp(Testptr[0],"PG") == 0 || 
				strcmp(Testptr[0],"PG_HT") == 0 || strcmp(Testptr[0],"RD") == 0 || 
				strcmp(Testptr[0],"RDS") == 0){
				for (int j = 0; j < ptrdimResCvg[1]; j++){
					for (int i = 0; i < ptrdimResCvg[0]; i++){
						if (i==0){ptrOutCvg[j]=0.0;}//Making sure ptrOutCvg is only filled with zeros at the start
						ptrOutCvg[j]+=ptrResCvg[i][j]/ptrdimResCvg[0];
					}
				}
			}
		}

		double valCI=ptrCI[0];
		OL_confidence_interval(ptrResults,ptrdimResults,valCI,
			ptrOutMean,ptrdimOutMean,ptrOutMin,ptrdimOutMin,ptrOutMax,ptrdimOutMax);
	}

	free(ptrResults);
	free(ptrResCvg);
	free(ptrAfa);
	
	return(nprot);	
}


