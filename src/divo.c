#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//For registration purposes:
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
//-------------------------

#include "nrutil.h"

#include "DP_read.h"
#include "OL_read.h"


SEXP read(SEXP DPorOL,

		  SEXP DP_X, SEXP DP_Test, SEXP DP_Alpha, SEXP DP_Resample, SEXP DP_CI, 
		  SEXP DP_Faaa, SEXP DP_PlugIn, SEXP DP_Sizee, SEXP DP_AlphaPro, SEXP DP_CVG, 
		  SEXP DP_BS, SEXP DP_OutMean, SEXP DP_OutMin, SEXP DP_OutMax, SEXP DP_OutCvg,

		  SEXP OL_X, SEXP OL_Test, SEXP OL_Alpha, SEXP OL_Resample, SEXP OL_CI, 
		  SEXP OL_Faaa, SEXP OL_PlugIn, SEXP OL_Sizee, SEXP OL_Beta, SEXP OL_CVG, 
		  SEXP OL_BS, SEXP OL_OutMean, SEXP OL_OutMin, SEXP OL_OutMax, SEXP OL_OutCvg
		  ) {

	//////////////////////////////////////////////////////////////////////////
	/* 								Initializing							*/
	//////////////////////////////////////////////////////////////////////////

	//Declaring an integer to keep track of Protected variables---------------
	int nprot=0;
	int DP_nprot=0;
	int OL_nprot=0;

	//Reading in DPorOL variable----------------------------------------------
	PROTECT(DPorOL = AS_CHARACTER(DPorOL)); nprot++;
	char *ptrDPorOL[1];//Defining pointer
	ptrDPorOL[0] = R_alloc(strlen(CHAR(STRING_ELT(DPorOL, 0))),sizeof(char));
	strcpy(ptrDPorOL[0], CHAR(STRING_ELT(DPorOL, 0)));//Assigning to variable

	//////////////////////////////////////////////////////////////////////////
	/* 							Choosing DP or OL							*/
	//////////////////////////////////////////////////////////////////////////

	if (strcmp(ptrDPorOL[0],"DP") == 0) {//Calling Diversity Profile functions

		DP_nprot = DP_read(DP_X, DP_Test, DP_Alpha, DP_Resample, DP_CI, 
			 			   DP_Faaa,DP_PlugIn, DP_Sizee, DP_AlphaPro, DP_CVG, 
			 			   DP_BS, DP_OutMean, DP_OutMin, DP_OutMax, DP_OutCvg);

	}
	else if (strcmp(ptrDPorOL[0],"OL") == 0) {//Calling Overlap functions

		OL_nprot = OL_read(OL_X, OL_Test, OL_Alpha, OL_Resample, OL_CI, 
		  				   OL_Faaa, OL_PlugIn, OL_Sizee, OL_Beta, OL_CVG, 
		  				   OL_BS, OL_OutMean, OL_OutMin, OL_OutMax, OL_OutCvg);

	}
	else {error("You neither selected DP nor OL");}


	//////////////////////////////////////////////////////////////////////////
	/* 							Unprotect and Return						*/
	//////////////////////////////////////////////////////////////////////////

	UNPROTECT(nprot+DP_nprot+OL_nprot);
	return(DPorOL);
}


	//////////////////////////////////////////////////////////////////////////
	/* 						For registration purposes						*/
	//////////////////////////////////////////////////////////////////////////


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n} //This defines {"myCall, (DL_FUNC) &myCall, 31"}

//Note that callMethods = R_CallDef, I think
static const R_CallMethodDef callMethods[]  = { 
	CALLDEF(read,31), //alternativley one could use: {"read", (DL_FUNC) &read, 31}
  	{NULL, NULL, 0}
};


// Register the .Call routines. No .C() or .Fortran() or .External() routines,so pass those arrays as NULL.
// Note that the function below is named after the file.c name 
void R_init_BiodivoTools(DllInfo *info){
	R_registerRoutines(info,NULL, callMethods,NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}






