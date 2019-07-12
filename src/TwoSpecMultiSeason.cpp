#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SingAnal.h"

void TwoSpecMultiSeason(int NBoot, char *s, int NPar, double *Params){ 
    double TwoSpeciesMSLink(double *Params, int NPar);
    void print_individual_estimates2(FILE *g, int realparmno, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);	
	void try_optmz_randiv(double (*xxFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);

    extern TSingSpec SingSpec;
    int K, ii, jj, PrmyPeriods=SingSpec.PrmyPeriods, N=SingSpec.N, T=SingSpec.T, tmp_linkfn;
    double *Hess, *Grad, *Work, MaxLL, **covar;
    int GAMBAAdifferent,ipar,jpar; FILE *g;
    
    g=SingSpec.g;
    SingSpec.psibar=new double**[PrmyPeriods+1]; 
    for (ii=0; ii<=PrmyPeriods; ii++) {SingSpec.psibar[ii]=new double*[N]; 
	  for (jj=0; jj<N; jj++) SingSpec.psibar[ii][jj]=new double[2];
	}
    fprintf(g,"Two Species, Multi-Season Model:\n---------------------------------\n");
    printf("Two Species, Multi-Season Model:\n---------------------------------\n");
    printf("alt_parm=%s\n",SingSpec.Realname[0][2]);
    if (strstr(SingSpec.Realname[0][2],"psiBa") != (void*)NULL || strstr(SingSpec.Realname[0][2],"psiB2") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=1;
        fprintf(g,"-- alternate parameterization requested: psiA,psiBA,psiBa,pA,pB,rA,rBA,rBa\n");
    }
    if (strstr(SingSpec.Realname[0][2],"nu") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=2;
        fprintf(g,"-- alternate parameterization requested: psiA,psiBa,nu,pA,pB,rA,rBa,rho\n");
    }
    if (strstr(SingSpec.Realname[0][2],"tau") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=3;
        fprintf(g,"-- alternate parameterization requested: psiAb,psiBa,tau,pA,pB,rAb,rBa,ups\n");
    }
    if (SingSpec.Alt_Param_Checked==0) 
            fprintf(g,"-- standard parameterization requested: psiA,psiB,phi,pA,pB,rA,rB,delta\n");
    if (SingSpec.Nstates==2) { 
        SingSpec.N/=2; N/=2; SingSpec.NMiss/=2; 
        fprintf(g,"-- stacked input: 1st %d records=species A, last %d records=species B\n",N,N);
    }
    else  fprintf(g,"-- compressed input detected: 1=species A, 2=species B, 3=both\n");
    printf("alt_parm_checked=%d\n",SingSpec.Alt_Param_Checked);
    
    fprintf(g,"\n%d Primary periods\n  Secondary periods:",PrmyPeriods);
    if (PrmyPeriods<1) {
        fprintf(g,"\nError: invalid number of primary periods (%d)\n",PrmyPeriods);
        exit(1);
    }
    for (ii=0; ii<PrmyPeriods; ii++) {
        fprintf(g," %d",SingSpec.SecPeriods[ii]); 
        if (SingSpec.SecPeriods[ii]<1) {
            fprintf(g,"\nError: invalid number of secondary periods (%d)\n",SingSpec.SecPeriods[ii]);
            exit(1);
        }
    } fprintf(g,"\n");
    
    fprintf(g,"Number of sites                = %d\n", SingSpec.N);
    fprintf(g,"Number of sampling occasions   = %d\n", T);
    fprintf(g,"Number of missing observations = %d\n\n", SingSpec.NMiss);
    
    // try amoeba routine to get starting values
    
    fprintf(g,"\nModel:%s Par=%d s=%s\n",SingSpec.modname,NPar,s);

	if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, TwoSpeciesMSLink, 1.0);
 
    Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2]; Grad = new double[NPar+2];
    covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii]=new double[NPar];
    printf("calling optmiz...\n");
	try_optmz_randiv(TwoSpeciesMSLink, NPar, Params, Work, Hess, &MaxLL);
    prntBetas(TwoSpeciesMSLink,NPar, Params, Work, MaxLL, covar); K=PrmyPeriods-1;
    
    fprintf(g,"\n============================================================\n");
    tmp_linkfn=SingSpec.LinkFn;
    print_individual_estimates2(g,0,0,0,1,covar,Params,tmp_linkfn,0,0);  // psiA
    print_individual_estimates2(g,1,0,1,1,covar,Params,tmp_linkfn,0,0);  // psiB or psiBA
    if (SingSpec.Alt_Param_Checked!=PSIBA_PSIBa) tmp_linkfn=expLnk;
    print_individual_estimates2(g,2,0,2,1,covar,Params,tmp_linkfn,0,0);  // phi or psiBa or nu
    ipar=jpar=3; 
	for (ii=1; ii<3; ii++) {                                //  ii=1 -> gam...,   ii=2 -> eps...
		tmp_linkfn=SingSpec.LinkFn; 
		print_individual_estimates2(g,ipar,ii,0,K,covar,Params,tmp_linkfn,0,0);  // gamAB
		ipar+=K;
		print_individual_estimates2(g,ipar,ii,K,K,covar,Params,tmp_linkfn,0,0);  // gamAb
		ipar+=K;
		if (SingSpec.Alt_Param_Checked==PSIBA_ODDS_RATIO) tmp_linkfn=expLnk;
		print_individual_estimates2(g,ipar,ii,2*K,K,covar,Params,tmp_linkfn,0,0);  // gamBAA
		ipar+=K; jpar=3*K;
		GAMBAAdifferent=(SingSpec.NrowsDM[1]>((PrmyPeriods-1)*4+1));
		if (GAMBAAdifferent) {
			print_individual_estimates2(g,ipar,ii,jpar,K,covar,Params,tmp_linkfn,0,0);  // gamBAa
			ipar+=K; jpar+=K;
		}
		print_individual_estimates2(g,ipar,ii,jpar,K,covar,Params,tmp_linkfn,0,0);  // gamBaA
		ipar+=K; jpar+=K;
		if (GAMBAAdifferent) {
			print_individual_estimates2(g,ipar,ii,jpar,K,covar,Params,tmp_linkfn,0,0);  // gamBaa
			ipar+=K; jpar+=K;
		}
	}
    tmp_linkfn=SingSpec.LinkFn;
    print_individual_estimates2(g,ipar  ,  3,0,T,covar,Params,tmp_linkfn,0,0); // pA
    print_individual_estimates2(g,ipar+T,  3,T,T,covar,Params,tmp_linkfn,0,0);  // pB
    print_individual_estimates2(g,ipar+2*T,3,2*T,T,covar,Params,tmp_linkfn,0,0);  // rA
    print_individual_estimates2(g,ipar+3*T,3,3*T,T,covar,Params,tmp_linkfn,0,0);  // rB or rBA
    if (SingSpec.Alt_Param_Checked!=PSIBA_PSIBa) tmp_linkfn=expLnk;           
    print_individual_estimates2(g,ipar+4*T,3,4*T,T,covar,Params,tmp_linkfn,0,0);  //  dlta or rBa
    // delete dynamic variables
    delete[] Hess;    delete[] Work;    delete[] Grad; 
    for (ii=0; ii<NPar; ii++) delete[] covar[ii];    delete[] covar;
	printf("exiting TwoSpecMultiSeason.cpp...\n");
}
