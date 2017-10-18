//http://stackoverflow.com/questions/4638568/write-2d-array-to-a-file-in-c
#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "nrutil.h"

void DP_saveBootstrap(int** ptrAfa,int* ptrdimAfa){

	FILE *fp;

	fp = fopen("Matricies.txt","a");

    for (int j=0;j<ptrdimAfa[1];j++) {
		for (int i=0;i<ptrdimAfa[0];i++) {
	    	fprintf(fp,"%d ",ptrAfa[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fclose(fp);

}
