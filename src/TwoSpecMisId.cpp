#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SingAnal.h"

void TwoSpecMisId(int NPar, double *Params) { 
  FILE *g;
  void TwoSpCorrDet(int NBoot, char *s, int NPar, double *Params);
  double TwoSpMisIdLink(double *Params, int NPar);
  double TwoSpMisIdMSLink(double *Params, int NPar);
  void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
  void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
  void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
  double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
  void try_optmz_randiv(double (*xxLinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);
  int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **covar) ;

  extern TSingSpec SingSpec;
  int ii;
  double *Hess, *Work, MaxLL, **covar;
  double (*xxLike)(double *p, int N);
  g=SingSpec.g;
  xxLike=&TwoSpMisIdLink;
  if (SingSpec.NrowsDM[2]>0) {
	  xxLike=TwoSpMisIdMSLink;
	  fprintf(g,"Two Species, Multi-Season Model with mis-identification:\n---------------------------------\n");
  }
  else 
    fprintf(g,"Two Species, Single Season Model with mis-identification:\n---------------------------------\n");
  SingSpec.Alt_Param_Checked=0;
  if (strstr(SingSpec.Realname[0][2],"psiBa") != (void*)NULL || strstr(SingSpec.Realname[0][2],"psiB2") != (void*)NULL) {
      SingSpec.Alt_Param_Checked=1;
      fprintf(g,"-- alternate parameterization requested: psiA,psiBA,psiBa, pA,pB,rA,rBA,rBa\n");
  }
  if (SingSpec.Alt_Param_Checked==0) {
          fprintf(g,"-- standard parameterization requested: psiA,psiB,phi, pA,pB,rA,rB,delta\n");
  }
  
  // try amoeba routine to get starting values
    
  fprintf(g,"\nModel:%s  NPar=%d\n",SingSpec.modname,NPar);

  if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, xxLike, 1.0);
  Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2]; 
  covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii]=new double[NPar];

  try_optmz_randiv(xxLike, NPar, Params, Work, Hess, &MaxLL);

  // output results
  prntBetas(xxLike,NPar, Params, Work, MaxLL, covar);
  fprintf(g,"\n============================================================\n");
  print_real_estimates(xxLike, NPar, Params, covar);
	
  // delete dynamic variables
  delete[] Hess;    delete[] Work; 
  for (ii=0; ii<NPar; ii++) delete[] covar[ii];    delete[] covar;
}
