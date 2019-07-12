#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

double get_mlogit_parm(double *Params, int parm_grp,int iparm, int istart, int ilen,int incr, int site, int srvy) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    int i; double rslt=0,sum=0,psi;
    //printf("get_mlogit_parm(grp=%d iparm=%d istart=%d ilen=%d incr=%d)\n",parm_grp,iparm,istart,ilen,incr);
    for (i=0; i<ilen; i+=incr) {
        //printf("calling getParam(i=%d grp -1 site 0)\n",i);
        psi = getParam(i+istart,parm_grp,-1,site,srvy,Params,0,expLnk); 
        //printf("get_mlogit psi=%f\n",psi);
        if (SingSpec.fixed[i+istart]<-998) sum+=psi;
    }
    if (SingSpec.fixed[iparm]<-998) 
        rslt=getParam(iparm,parm_grp,-1,site,srvy,Params,0,expLnk)/(1+sum);
    else rslt=SingSpec.fixed[iparm];
    //printf("rslt=%f\n",rslt);
    return(rslt);
}

void MSOpenPop(int NBoot, char *ss, int NPar, double *Params) {
    
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int iLinkFn, int fsite, int nsites);
    double MSOpenLink(double Params[], int npar);
    double MSOpenLink1a(double Params[], int npar);
    double (*MSOpenLinkFn)(double Params[], int npar);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void chisq_test(void);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);		
    int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **covar);
	double Random(void);
	void try_optmz_randiv(double (*MSOpenLinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);
	
    int i,j,nocc,N,ii, maxstate=0, site, nseasns, NrealPar=0;
    double **parmEst, **covar, *Hess,*Work, MaxLL;
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    
    nocc = SingSpec.T; 	N = SingSpec.N; MSOpenLinkFn=&MSOpenLink; 
    nseasns=SingSpec.PrmyPeriods;
    if (SingSpec.Nstates>0) maxstate=SingSpec.Nstates-1;
	if (SingSpec.Alt_Param_Checked==2) maxstate=3;  //  if integrated-habitat model, make sure 3 states
    SingSpec.Nstates=maxstate+1;
    
	//double *prd, *dpsi_dth; prd=new double[NPar]; dpsi_dth=new double[NPar];
    if (SingSpec.Verbose>1) printf("altParam:%d\n",SingSpec.Alt_Param_Checked);
    fprintf(g,"\nMulti-state-multi-season Model:\n");   // output results
    if (SingSpec.Alt_Param_Checked==1) fprintf(g,"  Psi-R parameterization selected\n");
    if (SingSpec.Alt_Param_Checked==2) fprintf(g,"  2-habitat integrated model parameterization selected\n");
	if (SingSpec.T==nseasns) {
		MSOpenLinkFn=&MSOpenLink1a; SingSpec.Alt_Param_Checked=4;
		fprintf(g,"  1srvy/seasn function used.\n");
	}
    
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, MSOpenLinkFn, 1.0);
    
    // set up covariance matrix of beta parameters
    covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii] = new double[NPar];
	Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2];
	
	try_optmz_randiv(MSOpenLinkFn, NPar, Params, Work, Hess, &MaxLL);

    // output results
    prntBetas(MSOpenLinkFn,NPar, Params, Work, MaxLL, covar);
	for (i=0; i<6; i++) NrealPar+=SingSpec.NrowsDM[i];
    parmEst=new double*[NrealPar];  
	for (i=0; i<NrealPar; i++) {
	    parmEst[i]=new double[N];
	    for (site=0; site<N; site++) parmEst[i][site]=SingSpec.realParmEst[i][site];
	} 	
	if (SingSpec.novar>3) { //     print estimates
      switch (SingSpec.Alt_Param_Checked) {
      case 1: 
          for (i=0; i<SingSpec.NrowsDM[0]; i++)
              print_individual_estimates2(g,i,0,i,1,covar,Params,SingSpec.LinkFn,0,0); 
          for (j=0; j<SingSpec.NrowsDM[4]; j++)
              print_individual_estimates2(g,i+j,4,j,1,covar,Params,SingSpec.LinkFn,0,0); 
          break;
      case 2:
          i=nseasns-1;
          print_individual_estimates2(g,0,0,0,1,covar,Params,SingSpec.LinkFn,0,0); ii=1;// pi
          print_individual_estimates2(g,1,0,1,1,covar,Params,SingSpec.LinkFn,0,0); ii=2;//Ps iA
          print_individual_estimates2(g,2,0,2,1,covar,Params,SingSpec.LinkFn,0,0); ii=3;//Ps iB
          print_individual_estimates2(g,3    ,1,  0,i,covar,Params,SingSpec.LinkFn,0,0); ii=i;//eta0 AA
          print_individual_estimates2(g,3+i  ,1,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eta1 AA
          print_individual_estimates2(g,3+i*2,1,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eta0 BA
          print_individual_estimates2(g,3+i*3,1,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eta1 BA
          print_individual_estimates2(g,3+i*4,2,0,i,covar,Params,SingSpec.LinkFn,0,0); ii=i;//gam AA
          print_individual_estimates2(g,3+i*5,2,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//gam AB
          print_individual_estimates2(g,3+i*6,2,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//gam BA
          print_individual_estimates2(g,3+i*7,2,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//gam BB
          print_individual_estimates2(g,3+i*8,3,0,i,covar,Params,SingSpec.LinkFn,0,0); ii=i;//eps AA
          print_individual_estimates2(g,3+i*9,3,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eps AB
          print_individual_estimates2(g,3+i*10,3,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eps BA
          print_individual_estimates2(g,3+i*11,3,ii,i,covar,Params,SingSpec.LinkFn,0,0); ii+=i;//eps BB
          print_individual_estimates2(g,3+i*12,4,0,nocc,covar,Params,SingSpec.LinkFn,0,0); ii=nocc;//p1   
          print_individual_estimates2(g,3+i*12+nocc,4,ii,nocc,covar,Params,SingSpec.LinkFn,0,0);//p2   
          break;
      case 3:
          site=0;  
          for (i=0; i<SingSpec.NrowsDM[0]; i++)
              print_individual_estimates2(g,i,0,i,1,covar,Params,SingSpec.LinkFn,0,0); 
          for (j=0; j<SingSpec.NrowsDM[4]; j++)
              print_individual_estimates2(g,i+j,4,j,1,covar,Params,SingSpec.LinkFn,0,0); 
          fprintf(g,"\n");
      default:
	  	for (i=0; i<NrealPar; i++) for (site=0; site<N; site++) SingSpec.realParmEst[i][site]=parmEst[i][site];
	    i=print_real_estimates(MSOpenLinkFn, NPar, Params, covar);
  
      }  // end switch
	}  // end if (SingSpec.novar>3)
    chisq_test(); 	
    // delete dynamic variables
    for (i=0; i<NrealPar; i++) free(parmEst[i]); free(parmEst);	 
}
