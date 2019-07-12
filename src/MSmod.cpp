#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif
	
void MSMod(int NBoot, int LoFBoot, int NPar, double *Params){
    // Data[][] contains detection/nondetection info
    // Missing[][] contains whether observations were missing
    // SiteCov[][] contains site-specific covariates
    // SampCov[][][] contains sampling occasion covariates
    // N is number of sites
    // T is number of sampling occasions
    // NSiteCov is number of site covariates
    // NSampCov is number of sampling covariates
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    double MS1Like(double *Pms, int NPar); 
    void CMGradient(double *Params, double *Grad, int NPar);
    bool InvertMatrix(double **a, int n);
    double ran1(long *idum);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void MSModCorr(int NBoot, int LofBoot,int NPar, double *Params);
	void chisq_test(void);
    void DoAmoeba( double Params, int NPar, double (*funk)(double [], int), double xincr);	
	void try_optmz_randiv(double (*LinkFn)(double *Params, int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);	
 	
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;

    int ii, jj, TruSt, maxstate=0, N=SingSpec.N, T=SingSpec.T;
    double **covar, *Hess,*Work,MaxLL;
    
    fprintf(g,"\nMulti-state model\n-----------------\n"); 
	ii=0; 
	if (SingSpec.NrowsDM[0]>2) 
	    if (strncmp(SingSpec.Realname[0][3],"th1.01",6)==0) ii=1;
	if (ii) MSModCorr(NBoot,LoFBoot,NPar,Params);
	else {
      for (ii=0; ii<N; ii++) for (jj=0; jj<T; jj++) 
		if (SingSpec.Data[ii][jj]>maxstate) maxstate = SingSpec.Data[ii][jj];
      SingSpec.Nstates=maxstate+1;
      if (maxstate>2) {
        printf("\n**** error: max of 2 states for this model parameterization ***\n");
        exit(1);
      }
      if (SingSpec.Verbose) printf("npar=%d nstates=%d\n",NPar,SingSpec.Nstates); 
      /*////////////////////////////////////////////////////////////////////////////
        Each row of Params corresponds to the estimated parameters.
        Parameters are in the order;  psi1,psi2 covariates, p1,p2 covariates, dlta covars
      */////////////////////////////////////////////////////////////////////////////
      if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, MS1Like, 1.0);
      // set up covariance matrix of beta parameters
      covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii] = new double[NPar];
      Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2];
	  try_optmz_randiv(MS1Like, NPar, Params, Work, Hess, &MaxLL);
		
      //optmiz(MS1Like, NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
    
      // output results
      prntBetas(MS1Like,NPar, Params, Work, MaxLL, covar);
    
      // calc average probability of occupancy (PsiBar) , variances
      fprintf(g,"\n============================================================\n");
	  for (TruSt=1; TruSt<SingSpec.Nstates; TruSt++) { 
		print_individual_estimates2(g,TruSt-1,0,TruSt-1,1,covar,Params,SingSpec.LinkFn,0,0); // PsiX
	  }
	  print_individual_estimates2(g,SingSpec.Nstates-1,1,0,T,covar,Params,SingSpec.LinkFn,0,0); // P1
	  print_individual_estimates2(g,SingSpec.Nstates+T-1,1,T,T,covar,Params,SingSpec.LinkFn,0,0); // P2
	  print_individual_estimates2(g,SingSpec.Nstates+T+T-1,2,0,T,covar,Params,SingSpec.LinkFn,0,0);  // delta
	  varcov_real(0,0,Params,NPar,covar);

	  //chisq_test();
      // delete dynamic variables
      for (ii=0; ii<NPar; ii++) delete[] covar[ii]; delete[] covar; 
	  delete[] Hess;  delete[] Work;
	}
}
