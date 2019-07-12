#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SingAnal.h"

void TwoSpCorrDet(int NBoot, char *s, int NPar, double *Params){ 

    double TwoSpCorrDetLink(double *Params, int NPar);
    int optmiz(double xxlike(double *parms, int NPars),int NPar, int Sig, int MaxFN, int iopt, double *Params,
               double *Hess, double *Grad, double *MaxLL, double *Work);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);	

    extern TSingSpec SingSpec;
    int ii, jj, N=SingSpec.N, T=SingSpec.T;
	double *Hess, *Grad, *Work, MaxLL, **covar; FILE *g;

    g=SingSpec.g;
    fprintf(g,"Two Species Corr. Detections, Single Season Model:\n---------------------------------\n");
    printf("alt_parm=%s\n",SingSpec.Realname[0][2]);
    if (strstr(SingSpec.Realname[0][2],"psiBa")>(void*)NULL || strstr(SingSpec.Realname[0][2],"psiB2")>(void*)NULL) {
        SingSpec.Alt_Param_Checked=1;
        fprintf(g,"-- alternate parameterization requested: psiA,psiBA,psiBa,pA,pB,rA,rBA,rBa\n");
    }
    else 
        if (strstr(SingSpec.Realname[0][2],"nu")>(void*)NULL) {
            SingSpec.Alt_Param_Checked=2;
            fprintf(g,"-- alternate parameterization requested: psiA,psiBa,nu,pA,pB,rA,rBa,rho\n");
        }
        else {
            SingSpec.Alt_Param_Checked=0;
            fprintf(g,"-- standard parameterization requested: psiA,psiB,phi,pA,pB,rA,rB,delta\n");
        }
    if (SingSpec.Nstates==2) { 
        SingSpec.N/=2; N/=2; SingSpec.NMiss/=2;
        fprintf(g,"-- stacked input: 1st %d records=species A, last %d records=species B\n",N,N);
    }
    else  fprintf(g,"-- compressed input detected: 1=species A, 2=species B, 3=both\n");
    printf("alt_parm_checked=%d\n",SingSpec.Alt_Param_Checked);
        
    // try amoeba routine to get starting values
    
    fprintf(g,"\nModel:%s NPar=%d s=%s\n",SingSpec.modname,NPar,s);

    if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, TwoSpCorrDetLink, 1.0);
    Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2]; Grad = new double[NPar+2];
    covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii]=new double[NPar];

    optmiz(TwoSpCorrDetLink,NPar, 7, SingSpec.maxfn, 0, Params-1, Hess, Grad-1, &MaxLL, Work);
	
	if (SingSpec.Verbose>1) {
		for (ii=0; ii<SingSpec.N; ii++) {
	      for (jj=0; jj<SingSpec.T; jj++) fprintf(g,"%d",SingSpec.Data[ii][jj]);
	      fprintf(g," %12.6f %12.6f %12.6f\n",SingSpec.det_hist_frq[ii],SingSpec.expval[ii],SingSpec.det_hist_frq[ii]-SingSpec.expval[ii]);
	    }
	}
    // output results
    prntBetas(TwoSpCorrDetLink,NPar, Params, Work, MaxLL, covar);
    
    //   call link function to print re-parameterized estimates
    //SingSpec.ifn=-2; fx=TwoSpeciesLink(Params,NPar);
    
    //  print_individual_estimates(FILE *g, char *lbl, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn);
    fprintf(g,"\n============================================================\n");
    print_individual_estimates2(g,0,0,0,1,covar,Params,SingSpec.LinkFn,0,0); // PsiA
    print_individual_estimates2(g,1,0,1,1,covar,Params,SingSpec.LinkFn,0,0);  // psiBA
    print_individual_estimates2(g,2,0,2,1,covar,Params,SingSpec.LinkFn,0,0); // psiBa
    
    print_individual_estimates2(g,3,0,3,T,covar,Params,SingSpec.LinkFn,0,0); // thetaA
    print_individual_estimates2(g,3+T,0,3+T,T,covar,Params,SingSpec.LinkFn,0,0); // thetaA'
    print_individual_estimates2(g,3+2*T,0,3+2*T,T,covar,Params,SingSpec.LinkFn,0,0);  // thetaBA
    print_individual_estimates2(g,3+3*T,0,3+3*T,T,covar,Params,SingSpec.LinkFn,0,0);  // thetaBA'
    print_individual_estimates2(g,3+4*T,0,3+4*T,T,covar,Params,SingSpec.LinkFn,0,0);  // thetaBa
    print_individual_estimates2(g,3+5*T,0,3+5*T,T,covar,Params,SingSpec.LinkFn,0,0);  // thetaBa'
    
    print_individual_estimates2(g,3+6*T,1,0,T,covar,Params,SingSpec.LinkFn,0,0);  // pA
    print_individual_estimates2(g,3+7*T,1,T,T,covar,Params,SingSpec.LinkFn,0,0);   // pB
    print_individual_estimates2(g,3+8*T,1,2*T,T,covar,Params,SingSpec.LinkFn,0,0);  // rA
    print_individual_estimates2(g,3+9*+T,1,3*T,T,covar,Params,SingSpec.LinkFn,0,0);
    print_individual_estimates2(g,3+10*T,1,4*T,T,covar,Params,SingSpec.LinkFn,0,0); // rBa
    
    print_individual_estimates2(g,3+11*T,2,0,3,covar,Params,SingSpec.LinkFn,0,0); // th0pi
    // delete dynamic variables
    delete[] Hess;    delete[] Work;    delete[] Grad; 
    for (ii=0; ii<NPar; ii++) delete[] covar[ii];    delete[] covar;
}
