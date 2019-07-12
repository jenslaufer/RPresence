#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "rngs.h"
// #include <afx.h>
#include "SingAnal.h"

void StagEntModl(int NBoot, char *s, int NPar, double *Params) {
    
    // declare function used here
    double StaggEntryModLink(double *Params, int NPar);
    int optmiz(double xxlike(double *parms, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params,
               double *Hess, double *Grad1, double *MaxLL, double *Work);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void chisq_test(void);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);	
	void try_optmz_randiv(double (*xxFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);
    
    extern TSingSpec SingSpec;
    FILE *g;  g=SingSpec.g;
    int ii, N=SingSpec.N, T=SingSpec.T, Prmy=SingSpec.PrmyPeriods;
    int dptr,eptr,igam,ieps;
    
    fprintf(g,"\nStaggered-Entry model\n---------------------\n");
    if (N==0) { printf("\n\nERROR: No data has been entered!\n"); exit(1); }
    if (SingSpec.NrowsDM[3]<1) {SingSpec.PrmyPeriods=1; SingSpec.SecPeriods[0]=T;}
    Prmy=SingSpec.PrmyPeriods;
    fprintf(g,"\n%d Primary periods\n",SingSpec.PrmyPeriods);
    fprintf(g,"\nModel(%d):%s\n\n",SingSpec.Model,SingSpec.modname);
    
    double *y;  y = new double[NPar+1]; 
    
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, StaggEntryModLink, 0.3);
    // set up covariance matrix of beta parameters
    double **covar, *Grad1;    covar = new double*[NPar];    Grad1 = new double[NPar+2];
    for (ii=0; ii<NPar; ii++) covar[ii] = new double[NPar];
    
    double *Hess,*Work,MaxLL;  Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2];
    //           call optimization routine
    fprintf(g,"\n\nSingle-season Staggered Entry Population Model:\n\n");  // output results
	
	try_optmz_randiv(StaggEntryModLink, NPar, Params, Work, Hess, &MaxLL);
    // output results
    prntBetas(StaggEntryModLink,NPar, Params, Work, MaxLL, covar);

    igam=1+T; ieps=igam+SingSpec.PrmyPeriods-1;
    eptr=ieps+SingSpec.PrmyPeriods-1; dptr=eptr+T;

    print_individual_estimates2(g,0,0,0,1,covar,Params,SingSpec.LinkFn,0,0); // Psi
 
    print_individual_estimates2(g,eptr,4,0,T-SingSpec.PrmyPeriods,covar,Params,SingSpec.LinkFn,0,0); // e

    print_individual_estimates2(g,dptr,5,0,T-SingSpec.PrmyPeriods,covar,Params,SingSpec.LinkFn,0,0); // d

    print_individual_estimates2(g,1,1,0,T,covar,Params,SingSpec.LinkFn,0,0);   // p
    
    if (Prmy>1) {
        print_individual_estimates2(g,igam,2,0,Prmy-1,covar,Params,SingSpec.LinkFn,0,0);   // gam
        print_individual_estimates2(g,ieps,3,0,Prmy-1,covar,Params,SingSpec.LinkFn,0,0);   // eps
    }

    varcov_real(0,0,Params,NPar,covar);

	chisq_test();
    // delete dynamic variables
    for (ii=0; ii<NPar; ii++) { delete [] covar[ii]; } delete[] covar;
    delete[] Grad1;    delete[] Hess;    delete[] Work;    delete[] y;
}
