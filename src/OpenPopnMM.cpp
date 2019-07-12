#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "rngs.h"
// #include <afx.h>
#include "SingAnal.h"

#define MISSNG -1
void OpenPopnMM(int NBoot, char *s, int NPar, double *Params, int nboot2) {
    
    void readfixed(FILE *f, char *s, double *Params); 
    // declare functions used here
    double OpenModLinkMM(double *Params, int NPar);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	double gof_open(double *Params, int prnt);
	void chisq_test(void);
	int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **cov);	
	void try_optmz_randiv(double (*LinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);	
    
    extern TSingSpec SingSpec;  FILE *g=SingSpec.g; 
    int site,seasn,ipar, ngrps, PrmyPeriods=SingSpec.PrmyPeriods, N=SingSpec.N;   SingSpec.Verbose=1;
		
    if (SingSpec.Verbose>0) printf("OpenPopnMM...(model=%d)\n",SingSpec.Model); 
	PrmyPeriods=SingSpec.PrmyPeriods=SingSpec.NrowsDM[2]+1;
    
    ngrps=SingSpec.NrowsDM[4]+1; if (ngrps<1) ngrps=1;

    SingSpec.psibar=new double**[PrmyPeriods+1];
    for (seasn=0; seasn<=PrmyPeriods; seasn++) {
        SingSpec.psibar[seasn]=new double*[N];
        for (site=0; site<N; site++) SingSpec.psibar[seasn][site]=new double[2];
    }
    
    fprintf(g,"\nMulti-season multi-method model - ");
    fprintf(g,"psi,theta(),gam(),esp(),p() parameterization\n");
	fprintf(g,"=======================================================================\n");
    if (N==0) { printf("\n\nERROR: No data has been entered!\n"); exit(1); }
    
    fprintf(g,"\n%d Primary periods\n  %d Methods\nSecondary periods:",PrmyPeriods,SingSpec.NMethods);
    if (PrmyPeriods<1) { fprintf(g,"\nError: invalid number of primary periods (%d)\n",PrmyPeriods); exit(1); }
    if (SingSpec.NMethods<1) { fprintf(g,"\nError: invalid number of methods (%d)\n",SingSpec.NMethods); exit(1); }
    for (seasn=0; seasn<PrmyPeriods; seasn++) {
        fprintf(g," %d",SingSpec.SecPeriods[seasn]); 
        if (SingSpec.SecPeriods[seasn]<1) {
		  fprintf(g,"\nError: invalid number of secondary periods (%d)\n",SingSpec.SecPeriods[seasn]);
		  exit(1);
        }
    } fprintf(g,"\n");
    
    fprintf(g,"\nModel(%d):%s\n\n",SingSpec.Model,SingSpec.modname);
    
    // all the data has been read in, now set things up for the minimisation
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, OpenModLinkMM, 0.3);
    // set up covariance matrix of beta parameters
    double **covar, *Grad1;    covar = new double*[NPar];    Grad1 = new double[NPar+2];
    for (ipar=0; ipar<NPar; ipar++) covar[ipar] = new double[NPar];
    
    double *Hess,*Work,MaxLL;  Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2];
    //           call optimization routine
    fprintf(g,"\n\nOpen Population Multi-method Model:\n\n");  // output results
    fprintf(g,"Number of primary sampling periods = %d\n", PrmyPeriods);
    fprintf(g,"Number of methods/survey           = %d\n", SingSpec.NMethods);
    
	try_optmz_randiv(OpenModLinkMM, NPar, Params, Work, Hess, &MaxLL);

    // output results
    fprintf(g,"\nModel has been fit using the logistic link.\n\n");
    prntBetas(OpenModLinkMM,NPar, Params, Work, MaxLL, covar);
    
    char starg[12]; strcpy(starg,"Psi"); seasn=PrmyPeriods; // if psi(i), parameterization...
    if(SingSpec.Model<2 || SingSpec.Model>=4) seasn=1;               // else # psi parms to print=1
	
	print_real_estimates(OpenModLinkMM, NPar, Params, covar);

    // delete dynamic variables
    delete[] Grad1; delete[] Hess; delete[] Work;
    for (seasn=0; seasn<=PrmyPeriods; seasn++) delete [] SingSpec.psibar[seasn];  delete[] SingSpec.psibar; 
    for (ipar=0; ipar<NPar; ipar++) delete [] covar[ipar];  delete[] covar; 
}
