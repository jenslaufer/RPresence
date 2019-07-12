#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

double get_mlogit_parm2(double *Params, int parm_grp, int ii, int i1, int i2, int i3, int site, int srvy) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;   double rslt=0,sum=0,psi[4];
	//if (site==0) printf("get_mlogit_parm2(%f %f %f   %d %d %d %d)",Params[0],Params[1],Params[2],ii,i1,i2,i3);
	psi[0] = getParam(i1,parm_grp,-1,site,srvy,Params,0,expLnk); 
    psi[1] = getParam(i2,parm_grp,-1,site,srvy,Params,0,expLnk); 
    psi[2] = 0; if (i3>=0) psi[2]=getParam(i3,parm_grp,-1,site,srvy,Params,0,expLnk); 
	sum=1+psi[0]+psi[1]+psi[2]; 
    if (SingSpec.fixed[ii]<-998) rslt=psi[ii]/sum;
    else rslt=SingSpec.fixed[ii]; 
	//if(site==0) printf(" psi=%f %f %f   rslt=%f\n",psi[0],psi[1],psi[2],rslt);
    return(rslt);
}

void print_mlogit_parms( int parm_grp, int parm_num, int ii, int i1, int i2, int i3, int NPar, double *Params, double **covar) {
    
	double psi=0,psi2,se,ci1,ci2,*prd,*dpsi_dth,eps=.1e-10,var,lastpsi;
	extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
	char **lineout,s[1024]; bool alleq=1;
	int  i,j,k,site, N=SingSpec.N;
	
    site=0; prd=new double[NPar]; dpsi_dth=new double[NPar]; 
    fprintf(g,"\n   parameter     site   estimate    std.error       95%% conf. interval\n");
    lineout=new char*[N]; 
    for (i=0; i<N; i++) lineout[i]=new char[200];
    for (site=0; site<N; site++) {
        lastpsi=psi; psi=get_mlogit_parm2(Params,parm_grp, ii, i1,i2,i3, site,0); 
		if (site>0) if (psi!=lastpsi) alleq=0;
        if (SingSpec.fixed[ii]<-998) {
            for (j=0; j<NPar; j++) prd[j]=0;
            for (j=0; j<NPar; j++) {
                Params[j]+=eps; 
                psi2=get_mlogit_parm2(Params,parm_grp, ii, i1,i2,i3, site,0);
                Params[j]-=eps; 
                dpsi_dth[j]=(psi2-psi)/eps;
                for (k=0; k<NPar; k++) prd[k]+=dpsi_dth[j]*covar[j][k];
            }
            for (j=0,var=0; j<NPar; j++) var+=prd[j]*dpsi_dth[j];
            se=sqrt(var); ci1=psi-1.96*se; ci2=psi+1.96*se;
            j=sprintf(lineout[site],"%12s %4d %12.5f %12.5f  %12.5f - %8.5f\n",
			          SingSpec.Realname[parm_grp][parm_num],site+1,psi,se,ci1,ci2);
        }
        else
            if (site==0) j=sprintf(lineout[0],"%12s **** %12.5f %s*** fixed\n",
			        SingSpec.Realname[parm_grp][parm_num],psi,s);
    }
    j=N; if (alleq) j=1; for (site=0; site<j; site++) fprintf(g,"%s",lineout[site]);

}

void MSModCorr(int NBoot, int LoFBoot, int NPar, double *Params)  {
    // Data[][] contains detection/nondetection info
    // Missing[][] contains whether observations were missing
    // SiteCov[][] contains site-specific covariates
    // SampCov[][][] contains sampling occasion covariates
    // N is number of sites
    // T is number of sampling occasions
    // NSiteCov is number of site covariates
    // NSampCov is number of sampling covariates
	void DoAmoeba(double *Params, int NPar, double (*funk)(double [], int), double xincr);
	double MS1CorrLike(double *Params, int NPar);
    void CMGradient(double *Params, double *Grad, int NPar);
    bool InvertMatrix(double **a, int n);
    double ran1(long *idum);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void chisq_test(void);
	void try_optmz_randiv(double (*LinkFn)(double *Params, int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);
	
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;

    int i, ii, jj, maxstate=0, N=SingSpec.N, T=SingSpec.T;
    double **covar, *Hess,*Work,MaxLL;
	
    
    fprintf(g,"\nMulti-state model w/ correlated detections\n-----------------\n");
    for (ii=0; ii<N; ii++) for (jj=0; jj<T; jj++) 
		if (SingSpec.Data[ii][jj]>maxstate) maxstate = SingSpec.Data[ii][jj];
    SingSpec.Nstates=maxstate+1;
    if (SingSpec.Verbose) printf("npar=%d nstates=%d\n",NPar,SingSpec.Nstates); 
    if (maxstate>3) {
        printf("\n**** error: max of 3 states for this model parameterization ***\n");
        exit(1);
    }
/*////////////////////////////////////////////////////////////////////////////
        Each row of Params corresponds to the estimated parameters.
        Parameters are in the order;  psi1,psi2 covariates, p1,p2 covariates, dlta covars
*/////////////////////////////////////////////////////////////////////////////
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, MS1CorrLike, 1.0);
		
    // set up covariance matrix of beta parameters
    covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii] = new double[NPar]; 
    Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2];
    printf("calling optmiz...npar=%d\n",NPar);
    try_optmz_randiv(MS1CorrLike, NPar, Params, Work, Hess, &MaxLL);
    //optmiz(MS1CorrLike, NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
    
    // output results
    prntBetas(MS1CorrLike,NPar, Params, Work, MaxLL, covar);
    
    // calc average probability of occupancy (PsiBar) , variances
    fprintf(g,"\n============================================================\n"); printf("=====\n");
	
	print_mlogit_parms(0,0,0,0,1,2, NPar,Params, covar);
	print_mlogit_parms(0,1,1,0,1,2, NPar,Params, covar);
	print_mlogit_parms(0,2,2,0,1,2, NPar,Params, covar);
	i=SingSpec.Nstates-1;
	for (ii=0; ii<2; ii++, i+=T) {
		print_individual_estimates2(g,i,0,i,T,covar,Params,SingSpec.LinkFn,0,0); // th1.11
	}
	print_mlogit_parms(0,i  ,0,i,i+T,-1, NPar,Params, covar); //th2.01
	print_mlogit_parms(0,i+T,1,i,i+T,-1, NPar,Params, covar); //th2.02
	i+=(2*T); 
	print_mlogit_parms(0,i  ,0,i,i+T,-1, NPar,Params, covar); //th2.11
	print_mlogit_parms(0,i+T,1,i,i+T,-1, NPar,Params, covar); //th2.12
	i+=(2*T); 
	print_mlogit_parms(0,i  ,0,i,i+T,-1, NPar,Params, covar); //th2.21
	print_mlogit_parms(0,i+T,1,i,i+T,-1, NPar,Params, covar); //th2.22

	i+=(2*T); 
	print_mlogit_parms(0,i    ,0,i,i+T,i+T+T, NPar,Params, covar); //th3.01
	print_mlogit_parms(0,i+T  ,1,i,i+T,i+T+T, NPar,Params, covar); //th3.02
	print_mlogit_parms(0,i+T+T,2,i,i+T,i+T+T, NPar,Params, covar); //th3.03
	i+=(3*T); 
	print_mlogit_parms(0,i    ,0,i,i+T,i+T+T, NPar,Params, covar); //th3.11
	print_mlogit_parms(0,i+T  ,1,i,i+T,i+T+T, NPar,Params, covar); //th3.12
	print_mlogit_parms(0,i+T+T,2,i,i+T,i+T+T, NPar,Params, covar); //th3.13
	i+=(3*T); 
	print_mlogit_parms(0,i    ,0,i,i+T,i+T+T, NPar,Params, covar); //th3.21
	print_mlogit_parms(0,i+T  ,1,i,i+T,i+T+T, NPar,Params, covar); //th3.22
	print_mlogit_parms(0,i+T+T,2,i,i+T,i+T+T, NPar,Params, covar); //th3.23
	i+=(3*T); 
	print_mlogit_parms(0,i    ,0,i,i+T,i+T+T, NPar,Params, covar); //th3.31
	print_mlogit_parms(0,i+T  ,1,i,i+T,i+T+T, NPar,Params, covar); //th3.32
	print_mlogit_parms(0,i+T+T,2,i,i+T,i+T+T, NPar,Params, covar); //th3.33

	i=26*T+3;
	print_individual_estimates2(g,i,2,0,1,covar,Params,SingSpec.LinkFn,0,0); // pi12,pi22
	
	i++; 
	print_mlogit_parms(2,1,0,i,i+1,-1, NPar,Params, covar);
	print_mlogit_parms(2,2,1,i,i+1,-1, NPar,Params, covar);
	i+=2; 
	print_mlogit_parms(2,3,0,i,i+1,i+2, NPar,Params, covar);
	print_mlogit_parms(2,4,1,i,i+1,i+2, NPar,Params, covar);
	print_mlogit_parms(2,5,2,i,i+1,i+2, NPar,Params, covar);
	
	chisq_test();
	//varcov_real(0,0,3,0,   0,0,Params,NPar,covar);

    // delete dynamic variables
    for (ii=0; ii<NPar; ii++) delete[] covar[ii]; delete[] covar; 
	delete[] Hess;  delete[] Work;
}
