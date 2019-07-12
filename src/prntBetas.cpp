#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "SingAnal.h"

void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar) { 
    double TwoSpeciesLink(double *Params, int NPar);
    void varcov(double xxLike(double *p, int N), double *Params, int NPar, double **covar);
    void varcov_mt(double xxLike(double *p, int N), double *Params, int NPar, double **covar);
    extern TSingSpec SingSpec; FILE *g=SingSpec.g; //time_t t1,t2; float t3=0;
    int ii,jj,iii,jjj,imat,flg=0;    double xa,va; 

	printf("\nNumber of parameters           = %d\n", NPar);
    printf("Number of function calls           = %1.0f\n",Work[2]);
    fprintf(g,"\nNumber of parameters           = %d\n", NPar);
    fprintf(g,"Number of function calls           = %1.0f\n",Work[2]);

    if (Work[3]<SingSpec.LikeNRSig) {
        fprintf(g,"\n**** Numerical convergence may not have been reached.\n");
        fprintf(g,"     Parameter estimates converged to approximately \n");
        fprintf(g,"     %4.2f significant digits.\n\n",Work[3]);
    }
    fprintf(g,"-2log(likelihood)              = %6.4f\n",2.0*MaxLL);
    fprintf(g,"AIC                            = %6.4f\n",2.0*(MaxLL+NPar));

    if (SingSpec.maxfn>0 && SingSpec.novar>1) {
        if (SingSpec.nthreads<2) 
			varcov(xxLike,Params,NPar,covar);
		else
			varcov_mt(xxLike,Params,NPar,covar);
        for (ii=0,flg=1; ii<NPar; ii++) flg=flg*(covar[ii][ii]>0);
        if (flg<1) fprintf(g,"\n***** neg. std.err(s) in VC matrix\n");
    }    
    char **CovNames=new char*[NPar];  //time(&t1);
    fprintf(g,"\nUntransformed Estimates of coefficients for covariates (Beta's)\n");
    fprintf(g,"======================================================================\n");
    fprintf(g,"                                          estimate    std.error\n");
    for (imat=jj=jjj=0; imat<=5; imat++) {
        for (ii=iii=0; ii<SingSpec.NParKK[imat]; ii++,jj++) {
			xa=SingSpec.BetaFixed[imat][ii]; va=0;
            CovNames[jj]=new char[64]; sprintf(CovNames[jj],"%c%d",'A'+imat,ii+1);
            if (xa>1.0e44) { xa=Params[jjj]; va=covar[jjj][jjj]; }
            fprintf(g,"%-4s %-32s : %10.6f  %10.6f\n",
                    CovNames[jj],SingSpec.Betaname[imat][ii],xa,sqrt(va));
            if (SingSpec.BetaFixed[imat][ii]>1.0e44) { iii++; jjj++; }
        }
    }

    if (SingSpec.novar==3 || SingSpec.novar==5 || SingSpec.novar==6) {
        fprintf(g,"\nVariance-Covariance Matrix of Untransformed estimates (Beta's):\n       ");
        for (ii=0; ii<NPar; ii++) fprintf(g,"  %7s  ",CovNames[ii]); fprintf(g,"\n");
        for (ii=0; ii<NPar; ii++) { fprintf(g,"%7s ",CovNames[ii]);
            for (jj=0; jj<NPar; jj++) fprintf(g," %10.6f",covar[ii][jj]); fprintf(g,"\n");
        }
        fprintf(g,"------------------------------\n");
    }
    for (jj=0; jj<NPar; jj++) delete [] CovNames[jj]; delete [] CovNames;
	//time(&t2); t3=t2-t1; printf("time to print betas: %f secs.\n",t3);
}
