#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "rngs.h"
// #include <afx.h>
#include "SingAnal.h"
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif

void OpenPopSD(int NPar, double *Params, int nboot2) {    
    // declare functions used here
    double OpenModLinkSD(double *Params, int NPar);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	void prnt_th0pi(int NPar, double *Params, double **covar, int ii);
	void chisq_test(void);
	int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **cov);
	void try_optmz_randiv(double (*LinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);	
    void gofSD(int NPar, double *Params, int LoFBoot);
	
    extern TSingSpec SingSpec;  FILE *g=SingSpec.g; 
    int i,ii=0,k, site,seasn,ipar, ngrps, pDM, th0pistrt, pstart,
          PrmyPeriods=SingSpec.PrmyPeriods, N=SingSpec.N, T=SingSpec.T, lmt;
    double sum=0;
    
    if (SingSpec.Verbose>0) printf("OpenPopn...(model=%d, apc=%d)\n",SingSpec.Model,SingSpec.Alt_Param_Checked); 
    ngrps=SingSpec.NrowsDM[4]+1; if (ngrps<1) ngrps=1;

    SingSpec.psibar=new double**[PrmyPeriods+1];
    for (seasn=0; seasn<=PrmyPeriods; seasn++) {
        SingSpec.psibar[seasn]=new double*[N];
        for (site=0; site<N; site++) SingSpec.psibar[seasn][site]=new double[2];
    }
    if (SingSpec.NrowsDM[2]>0) {
		fprintf(g,"\nMulti-season-Correlated-Detections model - \n=============================================\n");
		if (SingSpec.Verbose) printf("\nMulti-season-Correlated-Detections model - \n=============================================\n");
	}
	else {
		fprintf(g,"\nSingle-season-Correlated-Detections model - \n=============================================\n");
		if (SingSpec.Verbose) printf("\nSingle-season-Correlated-Detections model - \n=============================================\n");
	}
    if (N==0) { printf("\n\nERROR: No data has been entered!\n"); exit(1); }
	SingSpec.PrmyPeriods=PrmyPeriods=SingSpec.NrowsDM[1]+1;
	if (PrmyPeriods<2) SingSpec.SecPeriods[0]=SingSpec.T;
    
    fprintf(g,"\n%d Primary periods\n  Secondary periods:",PrmyPeriods);
    if (PrmyPeriods<1) {
        fprintf(g,"\nError: invalid number of primary periods (%d)\n",PrmyPeriods);
		printf("\nError: invalid number of primary periods (%d)\n",PrmyPeriods);
        exit(1);
    }
    for (seasn=0; seasn<PrmyPeriods; seasn++) {
        fprintf(g," %d",SingSpec.SecPeriods[seasn]); 
        if (SingSpec.SecPeriods[seasn]<1) {
            fprintf(g,"\nError: invalid number of secondary periods (%d)\n",SingSpec.SecPeriods[seasn]);
			printf("\nError: invalid number of secondary periods (%d)\n",SingSpec.SecPeriods[seasn]);
            exit(1);
        }
    } 
    fprintf(g,"\n\nModel(%d):%s\n\n",SingSpec.Model,SingSpec.modname);
    
    // all the data has been read in, now set things up for the minimisation
    if(SingSpec.UseAmoeba) DoAmoeba(Params,NPar,OpenModLinkSD,0.3);
    // set up covariance matrix of beta parameters
    double **covar, *Grad1;    covar = new double*[NPar];    Grad1 = new double[NPar+2];
    for (ipar=0; ipar<NPar; ipar++) covar[ipar] = new double[NPar];
    
    double *Hess,*Work,MaxLL;  Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2];
    //           call optimization routine
    fprintf(g,"\n\nOpen Population Model:\n\n");  // output results
	if (SingSpec.Verbose) printf("\n\nOpen Population Model:\n   calling optmiz...\n");
    
	try_optmz_randiv(OpenModLinkSD, NPar, Params, Work, Hess, &MaxLL); if (SingSpec.Verbose) printf("   optmiz done\n");
    // output results
    fprintf(g,"\nModel has been fit using the logistic link.\n\n");
    prntBetas(OpenModLinkSD,NPar, Params, Work, MaxLL, covar);
    
    char starg[12]; strcpy(starg,"Psi"); seasn=PrmyPeriods; // if psi(i), parameterization...
    print_individual_estimates2(g,0,0,0,1,covar,Params,SingSpec.LinkFn,0,0);  // psi
    for (seasn=i=0; seasn<PrmyPeriods; i+=SingSpec.SecPeriods[seasn++]) 
		print_individual_estimates2(g,i+1,0,i+1,SingSpec.SecPeriods[seasn],covar,Params,SingSpec.LinkFn,0,0);  // th0
	for (seasn=0; seasn<PrmyPeriods; i+=SingSpec.SecPeriods[seasn++])   
		print_individual_estimates2(g,i+1,0,i+1,SingSpec.SecPeriods[seasn],covar,Params,SingSpec.LinkFn,0,0);  // th1
	if (PrmyPeriods>1) {
		print_individual_estimates2(g,2*T+1,1,0,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // gamma
        print_individual_estimates2(g,2*T+PrmyPeriods,2,0,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // epsilon
    }
	print_individual_estimates2(g,2*T+PrmyPeriods*2-1,3,0,SingSpec.NrowsDM[3],covar,Params,SingSpec.LinkFn,0,0); // p
    if (SingSpec.NrowsDM[5]>0) {
		ii=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3]+SingSpec.NrowsDM[4];
        print_individual_estimates2(g,ii,5,0,SingSpec.NrowsDM[5],covar,Params,SingSpec.LinkFn,0,0); // pi (for het. model)
	}
	pstart=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]; pDM=3;
	if (SingSpec.NrowsDM[2]>0)
		if (SingSpec.Realname[2][0][0]=='t') { pstart=SingSpec.NrowsDM[0]; pDM=1; }
	th0pistrt=pstart+SingSpec.NrowsDM[pDM];
	
    if (SingSpec.fixed[th0pistrt]>1)  //  if th0pi is fixed and > 1, 
        prnt_th0pi(NPar, Params, covar, 0);      //     th0pi computed as equlib value
    else {
        print_individual_estimates2(g,th0pistrt,4,0,SingSpec.NrowsDM[4],covar,Params,SingSpec.LinkFn,0,0); // th0pi
		prnt_th0pi(NPar, Params, covar, 1);      //     th0 for 1st segment
	}
    ii=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3]+SingSpec.NrowsDM[4];
		
	//if (SingSpec.novar>3) varcov_real(0,0,Params,NPar,covar);
 /////////////
	double *lastpsi,*lastlam,*grad,*grad2,*prd,**lamout,**lamoutse,psi0,gam,eps,psiI,lam,xlam=0,xlam2=0,xpsi=0,xpsi2=0;
	char lbl[80]="psi"; double xpsise,vc,ci1,ci2,verysmallnumber=.1e-10; int j; lmt=SingSpec.lmt;
	if (lmt>SingSpec.N) lmt=SingSpec.N;
	if (PrmyPeriods>1) {
		fprintf(g,"\n\n DERIVED parameters - %s2,%s3,%s4,...\n",lbl,lbl,lbl);
		fprintf(g,"\n        Site                        %s(t)  Std.err     95%% conf. interval\n",lbl);
		lastpsi=new double[PrmyPeriods]; lastlam=new double[PrmyPeriods]; 
		for (seasn=0; seasn<PrmyPeriods; seasn++) lastlam[seasn]=lastpsi[seasn]=-1;
		grad=new double[NPar]; grad2=new double [NPar]; prd=new double[NPar];
		lamout=new double*[lmt]; lamoutse=new double*[lmt]; 
		for (site=0; site<lmt; site++) {
			lamout[site]=new double[PrmyPeriods]; 
			lamoutse[site]=new double[PrmyPeriods];
			lamout[site][0]= lamoutse[site][0]=-1;
		}	
		for (site=0; site<lmt; site++) { 
			for (k=1; k<PrmyPeriods; k++) {            //  compute psi(i) for each primary period, site
				for (i=0; i<NPar; i++) {
					psi0 = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); 
					SingSpec.psibar[0][0][0]=psiI=psi0;
					for (seasn=1; seasn<PrmyPeriods; seasn++) {
						gam = getParam(2*T+seasn,1,seasn-1,site,seasn-1,Params,0,1);
						eps = getParam(2*T+PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
						if (SingSpec.Alt_Param_Checked) eps=1-eps;
						lam=psiI; psiI= psiI*(1-eps)+(1-psiI)*gam; 	lam=psiI/lam; 
						if (seasn==k) { xpsi=psiI; xlam=lam;}
						SingSpec.psibar[0][0][0]=psiI; 
					}
					Params[i]-=verysmallnumber;
					psi0 = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); 
					SingSpec.psibar[0][0][0]=psiI=psi0;
					for (seasn=1; seasn<PrmyPeriods; seasn++) {
						gam = getParam(2*T+seasn,1,seasn-1,site,seasn-1,Params,0,1);
						eps = getParam(2*T+PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
						if (SingSpec.Alt_Param_Checked) eps=1-eps;
						lam=psiI; psiI= psiI*(1-eps)+(1-psiI)*gam; lam=psiI/lam;
						if (seasn==k) { xpsi2=psiI; xlam2=lam;}
						SingSpec.psibar[0][0][0]=psiI;
					}
					grad[i]=(xpsi2-xpsi)/verysmallnumber; grad2[i]=(xlam2-xlam)/verysmallnumber;
					Params[i]+=verysmallnumber;
				}
				//        Print derived psi's (or eps or gam)
				//           matrix multiply grad * varcov(betas)
				for (i=0; i<NPar; prd[i++]=sum)for (j=0,sum=0; j<NPar; j++) sum+=grad[j]*covar[j][i];
				//           matrix multiply {grad * varcov(betas)} * grad'
				for (i=0,vc=0; i<NPar; i++) vc+=prd[i]*grad[i]; 
				if (SingSpec.Model!=4) 
					if (lastpsi[k]!=xpsi) {
						xpsise=sqrt(vc); ci1=xpsi-1.96*xpsise; ci2=xpsi+1.96*xpsise;
						fprintf(g,"%s(%2d)  %4d %8s      :     %8.4f%8.4f   %8.4f -%7.4f \n",lbl,k+1,site+1,
								SingSpec.sitename[site],xpsi,xpsise,ci1,ci2);
					}
				lastpsi[k]=xpsi; SingSpec.psibar[k-1][0][0]=xpsi;
				//        Save derived lambda's for printing later
				//           matrix multiply grad * varcov(betas)
				for (i=0; i<NPar; prd[i++]=sum)for (j=0,sum=0; j<NPar; j++) sum+=grad2[j]*covar[j][i];
				//           matrix multiply {grad * varcov(betas)} * grad'
				for (i=0,vc=0; i<NPar; i++) vc+=prd[i]*grad2[i];
				lamout[site][k]=xlam; lamoutse[site][k]=sqrt(vc);
			}
		}
		fprintf(g,"\n\n DERIVED parameters - lam2,lam3,lam4,...\n"); strcpy(lbl,"lam");
		fprintf(g,"\n        Site                        %s(t)  Std.err     95%% conf. interval\n",lbl);
		for (site=0; site<lmt; site++) {                // 
			for (k=1; k<PrmyPeriods; k++) {            //  compute psi(i) for each primary period, site
				ci1=lamout[site][k]-1.96*lamoutse[site][k]; ci2=lamout[site][k]+1.96*lamoutse[site][k];
				j=0; if (site>0) if (lamout[site][k]!=lamout[site-1][k]) j=1; if (site==0) j=1;
				if (j>0)
					fprintf(g,"%s(%2d)  %4d %8s      :     %8.4f%8.4f   %8.4f -%7.4f\n",lbl,k+1,site+1,
							SingSpec.sitename[site],lamout[site][k],lamoutse[site][k],ci1,ci2);
			}  //  end for k=1 to PrmyPeriods-1
		}   // end for site=0 to N
		
		delete [] lastpsi; delete [] lastlam; delete [] prd; delete [] grad; delete [] grad2;
		for (site=0; site<lmt; site++) { delete [] lamout[site]; delete [] lamoutse[site]; }
		delete [] lamout; delete [] lamoutse;		
	}  //  end if PrmyPeriods>1

	chisq_test();
	if (nboot2>0 && PrmyPeriods<2) gofSD(NPar, Params, nboot2) ;

 /////////////
    
    // delete dynamic variables
    for (ipar=0; ipar<NPar; ipar++) delete [] covar[ipar];  delete[] covar; 
    for (seasn=0; seasn<=PrmyPeriods; seasn++) delete [] SingSpec.psibar[seasn];  delete[] SingSpec.psibar; 
    delete[] Grad1; delete[] Hess; delete[] Work;
}

void prnt_th0pi(int NPar, double *Params, double **covar, int ii) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	extern TSingSpec SingSpec; FILE *g=SingSpec.g;
	double *tmpprd,sum,xpsise,xpsi=-1,ci1,ci2,th0,th1,prd,verysmallnumber=.000000000001,psi0,vc,*grad1,prev_psi,prev_se,th0pi; 
	int i,j,seasn,site,PrmyPeriods=SingSpec.PrmyPeriods,th0strt,T=SingSpec.T,N=SingSpec.N,th0pistrt; 
    char lbl[2][20]={"th0pi","th0(1)"},lbl1[20];

    fprintf(g,"\n\n DERIVED parameters\n"); 
    if (ii==1) fprintf(g,"      th0(1) = th0pi*th0 + (1-th0pi)*th1 = Pr(1st segment is used)\n"); 
    if (ii==0) fprintf(g,"      th0pi = th0(1) = (th0(1) / (1+th0(1)-th1(1)) = Pr(segment 0 is used)\n"); 

    fprintf(g,"\n        Site                     %s  Std.err     95%% conf. interval\n",lbl[ii]);
	tmpprd=new double[NPar+1]; grad1=new double[NPar+1]; int irow0=1, irow1=T+1;
	th0strt=0; xpsi=xpsise=-1; 
	th0pistrt=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3];
    for (seasn=0; seasn<PrmyPeriods; seasn++) {
		i=sprintf(lbl1,"%s(%d)",lbl[ii],seasn+1);
		for (site=0; site<N; site++) { 
		    for (i=0; i<NPar; i++) {
		        th0=getParam(th0strt+irow0,0,irow0,site,irow0-1,Params,0,SingSpec.LinkFn);
			    th1=getParam(th0strt+irow1,0,irow1,site,irow0-1,Params,0,SingSpec.LinkFn);
				prd=1+th0-th1; xpsi=th0/prd; 
				if (ii>0) {
					th0pi=getParam(th0pistrt+seasn,4,seasn,site,seasn,Params,0,SingSpec.LinkFn);
					xpsi=th0pi*th1+(1-th0pi)*th0;
				}
				Params[i]+=verysmallnumber;
		        th0=getParam(th0strt+irow0,0,irow0,site,irow0-1,Params,0,SingSpec.LinkFn);
		        th1=getParam(th0strt+irow1,0,irow1,site,irow0-1,Params,0,SingSpec.LinkFn);
			    prd=1+th0-th1; psi0=th0/prd;
				if (ii>0) {
					th0pi=getParam(th0pistrt+seasn,4,seasn,site,seasn,Params,0,SingSpec.LinkFn);
					psi0=th0pi*th1+(1-th0pi)*th0;
				}
				grad1[i]=(psi0-xpsi)/verysmallnumber;
				Params[i]-=verysmallnumber;
		    }
		    //  matrix multiply grad * varcov(betas)
			for (i=0; i<NPar; tmpprd[i++]=sum)for (j=0,sum=0; j<NPar; j++) sum+=grad1[j]*covar[j][i];
	        //           matrix multiply {grad * varcov(betas)} * grad'
		    for (i=0,vc=0; i<NPar; i++) vc+=tmpprd[i]*grad1[i];
		    prev_psi=xpsi; prev_se=xpsise; xpsise=sqrt(vc); 
			if (fabs(xpsi-prev_psi)>.00001 || fabs(xpsise-prev_se)>.00001) { // || SingSpec.dups) {
				ci1=xpsi-1.96*xpsise; ci2=xpsi+1.96*xpsise;
				//        Print derived th0pi's 
				fprintf(g,"%-12s %4d %-16s:%8.4f %8.4f   %8.4f -%7.4f\n",
					lbl[ii],site+1,SingSpec.sitename[site],xpsi,xpsise,(ci1>0 ? ci1:0),(ci2<1 ? ci2:1));
			}
		}
		irow0+=SingSpec.SecPeriods[seasn]; irow1+=SingSpec.SecPeriods[seasn];
	}
	delete [] tmpprd; delete [] grad1;
}
