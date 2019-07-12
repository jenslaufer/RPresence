#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "rngs.h"
// #include <afx.h>
#include "SingAnal.h"

#define MISSNG -1
void OpenPopn(int NBoot, char *s, int NPar, double *Params, int nboot2) {

    void readfixed(FILE *f, char *s, double *Params);
    // declare functions used here
    double OpenModLink(double *Params, int NPar);
	double OpenModLink1(double *Params, int NPar);
	double (*fnptr)(double *Params, int NPar);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void varcov_real(int site, int srvy, double *Params, int NPar, double **covar);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    void Open_bootstrap(int NBoot, int NPar, double *Params, double **OrigPsi, double **OrigGam, double **OrigEps, double **OrigP);
	double gof_open(double *Params, int prnt);
	void chisq_test(void);
    void OpenPopnMM(int NBoot, char *s, int NPar, double *Params, int nboot2);
	int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **covar);
	double psi_i(int site, int seasn, double psi0, double *Params);
	void try_optmz_randiv(double (*LinkFn)(double *Params, int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);
	int optmiz(double xxLike(double *p, int N), int n, int nsig, int maxfn, int iopt, double *x, double *h, double *g, double *rf, double *w);

    int i,j,k, site,srvy,seasn,ipar, ngrps, PrmyPeriods, N, T,maxsite=1;
    extern TSingSpec SingSpec;  FILE *g=SingSpec.g;
	PrmyPeriods=SingSpec.PrmyPeriods, N=SingSpec.N, T=SingSpec.T;

    if (SingSpec.Verbose>0) printf("OpenPopn...(model=%d)\n",SingSpec.Model);
    ngrps=SingSpec.NrowsDM[4]+1; if (ngrps<1) ngrps=1;
    
    SingSpec.psibar=new double**[PrmyPeriods+1];
    for (seasn=0; seasn<=PrmyPeriods; seasn++) {
        SingSpec.psibar[seasn]=new double*[N];
        for (site=0; site<N; site++) SingSpec.psibar[seasn][site]=new double[2];
    }
	fnptr=OpenModLink; 
    fprintf(g,"\nMulti-season model - ");
    switch(SingSpec.Model) {
      case 2:      fprintf(g,"psi(),gam(),p() parameterization\n"); break;
      case 3:      fprintf(g,"psi(),eps(),p() parameterization\n"); break;
      case 4:      fprintf(g,"psi,gam(),p(), eps=1-gam parameterization\n"); break;
      case 5:      fprintf(g,"psi,gam(),eps(),p() autologistic parameterization\n"); break;
	  case 6:      fprintf(g,"psi,gam(),eps(),p() \n"); fnptr=OpenModLink1; break;
      default:     fprintf(g,"psi,gam(),esp(),p() parameterization\n");
	}   // end switch
	fprintf(g,"=======================================================================\n");
    if (N==0) { printf("\n\nERROR: No data has been entered!\n"); exit(1); }

    fprintf(g,"\n%d Primary periods\n  Secondary periods:",PrmyPeriods);
    if (PrmyPeriods<1) {
        fprintf(g,"\nError: invalid number of primary periods (%d)\n",PrmyPeriods);
        exit(1);
    }
    for (seasn=0; seasn<PrmyPeriods; seasn++) {
        fprintf(g," %d",SingSpec.SecPeriods[seasn]);
        if (SingSpec.SecPeriods[seasn]<1) {
            fprintf(g,"\nError: invalid number of secondary periods (%d)\n",SingSpec.SecPeriods[seasn]);
            exit(1);
        }
    } fprintf(g,"\n");

    if (SingSpec.Verbose>0) printf("\nModel(%d):%s\n\n",SingSpec.Model,SingSpec.modname);
	fprintf(g,"\nModel(%d):%s\n\n",SingSpec.Model,SingSpec.modname);

    // all the data has been read in, now set things up for the minimisation
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, fnptr, 0.3);
    // set up covariance matrix of beta parameters
    double **covar, *Grad1;    covar = new double*[NPar];    Grad1 = new double[NPar+2];
    for (ipar=0; ipar<NPar; ipar++) covar[ipar] = new double[NPar];

    double *Hess,*Work,MaxLL;  Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2];

    if (SingSpec.Verbose>0) printf("\n\nOpen Population Model:\n\nNumber of primary sampling periods = %d\n", PrmyPeriods);
	fprintf(g,"\n\nOpen Population Model:\n\nNumber of primary sampling periods = %d\n", PrmyPeriods);
    //           call optimization routine
	//try_optmz_randiv(OpenModLink, NPar, Params, Work, Hess, &MaxLL);
	try_optmz_randiv(fnptr, NPar, Params, Work, Hess, &MaxLL);
	if (SingSpec.Verbose>1) printf("optimize done... maxll=%f\n",MaxLL);


    // output results
    fprintf(g,"\nModel has been fit using the logistic link.\n\n");
    prntBetas(fnptr,NPar, Params, Work, MaxLL, covar);

	i=print_real_estimates(fnptr, NPar, Params, covar);

    char starg[12]; strcpy(starg,"Psi"); seasn=PrmyPeriods; // if psi(i), parameterization...
	if (1==2) {
    if(SingSpec.Model<2 || SingSpec.Model>=4) seasn=1;               // else # psi parms to print=1
    print_individual_estimates2(g,0,0,0,seasn,covar,Params,SingSpec.LinkFn,0,0);  // psi
    //if (SingSpec.Verbose) {
    //    fprintf(g,"avg psi for neighboring sitesr:\n");
    //    for (site=0; site<N; site++) fprintf(g,"site %d avg psi*n=%f\n",i,SingSpec.psibar[0][site][0]); fprintf(g,"\n");
    //}
    if (SingSpec.NrowsDM[0]>PrmyPeriods) {   //  print gam if in combined 1st design matrix
        print_individual_estimates2(g,seasn,0,seasn,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // gam/eps
        seasn+=PrmyPeriods-1;
    }
    if (SingSpec.NrowsDM[1]>0) {
        print_individual_estimates2(g,seasn,1,0,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // gamma
        seasn+=PrmyPeriods-1;
	}
    if (SingSpec.NrowsDM[0]>seasn) {   //  print gam/eps if in combined 1st design matrix
        if (SingSpec.Realname[0][seasn][0]=='e') {
            print_individual_estimates2(g,seasn,0,seasn,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // gam/eps
            seasn+=PrmyPeriods-1;
        }
    }
    if (SingSpec.NrowsDM[2]>0) {
        print_individual_estimates2(g,seasn,2,0,PrmyPeriods-1,covar,Params,SingSpec.LinkFn,0,0); // eps
        seasn+=PrmyPeriods-1;
    }
    if (SingSpec.NrowsDM[0]>seasn) {   //  print p if in combined 1st design matrix
        i=SingSpec.NrowsDM[0]-seasn;
        print_individual_estimates2(g,seasn,0,seasn,i,covar,Params,SingSpec.LinkFn,0,0); // gam/eps
    }
    print_individual_estimates2(g,SingSpec.PrmyPeriods*2+1,3,0,SingSpec.NrowsDM[3],covar,Params,SingSpec.LinkFn,0,0); // p

    if (SingSpec.NrowsDM[4]>0)   //  print Pledger mixture parms...
        print_individual_estimates2(g,SingSpec.PrmyPeriods*2+1+T,4,0,SingSpec.NrowsDM[4],covar,Params,SingSpec.LinkFn,0,0); // p
    }
	//if (SingSpec.novar>3) varcov_real(0,0,Params,NPar,covar);	
    if (NBoot>0)  {
	    // store Original real parameters Psi,gamma,eps,p
		double **OrigPsi,**OrigGam,**OrigEps,**OrigP; SingSpec.simulating=1;
		OrigPsi = new double*[N]; OrigGam = new double*[N]; OrigEps = new double*[N]; OrigP = new double *[N];
		for (site=0; site<N; site++) {
			OrigPsi[site]=new double[SingSpec.PrmyPeriods+1]; OrigP[site]=new double[T];
			OrigGam[site]=new double[SingSpec.PrmyPeriods]; OrigEps[site]=new double[SingSpec.PrmyPeriods];
			//double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params,int transf, int LinkFn) {			
			OrigPsi[site][0]=getParam(0,0,0,site,0,Params,0,1);   //  Phi[site][0][0][0];
			for (k=0; k<SingSpec.PrmyPeriods; k++) {
				OrigGam[site][k]=getParam(1+k,1,k,site,k,Params,0,1);   // =Phi[site][1][0][k+1];
				OrigEps[site][k]=getParam(k+SingSpec.PrmyPeriods,2,k,site,k,Params,0,1);   // =1-Phi[site][0][0][k+1];
				OrigPsi[site][k+1]=OrigPsi[site][k]*(1.-OrigEps[site][k])+(1-OrigPsi[site][k])*OrigGam[site][k];
			}
			for (k=0; k<T; k++) OrigP[site][k]=getParam(2*SingSpec.PrmyPeriods-1+k,3,k,site,k,Params,0,1);   //=SingSpec.closedp[site][k][0];
		}		
		Open_bootstrap(NBoot, NPar, Params, OrigPsi, OrigGam, OrigEps, OrigP);
		SingSpec.simulating=0;
		for (site=0; site<N; site++) {
			delete [] OrigPsi[site]; delete [] OrigGam[site]; delete [] OrigEps[site]; delete [] OrigP[site];
		}
		delete [] OrigPsi; delete [] OrigGam; delete [] OrigEps; delete [] OrigP;		
	}
    char lbl[12];
    double psi0,psiI,xpsi=0,xpsi2=-1,lam,xlam=0,xlam2=0,ci1,ci2,eps=0,sum,vc,*prd,*grad,*grad2,
          verysmallnumber=.000000000001, **psiout, **psioutse, **lamout,**lamoutse,eps2;

    if (SingSpec.Model<5 && SingSpec.novar>3) {
        //   print derived psi's (psi2,psi3,...)
        fprintf(g,"\n"); for (i=0; i<80; i++) fprintf(g,"="); strcpy(lbl,"lam");
        grad=new double[NPar]; grad2=new double [NPar]; prd=new double[NPar];
        psiout=new double*[N]; psioutse=new double*[N]; lamout=new double*[N]; lamoutse=new double*[N];
        for (site=0; site<N; site++) {
            psiout[site]=new double[PrmyPeriods];  lamout[site]=new double[PrmyPeriods];
            psioutse[site]=new double[PrmyPeriods];lamoutse[site]=new double[PrmyPeriods];
            psiout[site][0]=psioutse[site][0]=lamout[site][0]=lamoutse[site][0]=-1;
        }
        for (site=0; site<N; site++) {
            for (k=1; k<PrmyPeriods; k++) {            //  compute psi(i) for each primary period, site
                for (i=0; i<NPar; i++) {
                    psi0 = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
                    SingSpec.psibar[0][0][0]=psi0;
                    for (seasn=1; seasn<PrmyPeriods; seasn++) {
                        psiI=SingSpec.psibar[0][0][0]= psi_i(site,seasn,psi0,Params); lam=psiI/psi0; psi0=psiI;
                        if (seasn==k) { xpsi=psiI; xlam=lam;}
                    }
                    Params[i]-=verysmallnumber;
                    psi0=SingSpec.psibar[0][0][0]=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
                    for (seasn=1; seasn<PrmyPeriods; seasn++) {
                        psiI=SingSpec.psibar[0][0][0]=psi_i(site,seasn,psi0,Params); lam=psiI/psi0; psi0=psiI;
                        if (seasn==k) { xpsi2=psiI; xlam2=lam;}
                    }
                    grad[i]=(xpsi2-xpsi)/verysmallnumber; grad2[i]=(xlam2-xlam)/verysmallnumber;
                    Params[i]+=verysmallnumber;
                }
                //           matrix multiply grad * varcov(betas)
                for (i=0; i<NPar; prd[i++]=sum) for (j=0,sum=0; j<NPar; j++) sum+=grad[j]*covar[j][i];
                //           matrix multiply {grad * varcov(betas)} * grad'
                for (i=0,vc=0; i<NPar; i++) vc+=prd[i]*grad[i];
				psiout[site][k]=xpsi; psioutse[site][k]=sqrt(vc);
                //           matrix multiply grad2 * varcov(betas) to get var(lam)
                for (i=0; i<NPar; prd[i++]=sum)for (j=0,sum=0; j<NPar; j++) sum+=grad2[j]*covar[j][i];
                //           matrix multiply {grad2 * varcov(betas)} * grad2'
                for (i=0,vc=0; i<NPar; i++) vc+=prd[i]*grad2[i];
                lamout[site][k]=xlam; lamoutse[site][k]=sqrt(vc);
            }
        }
		if (SingSpec.NrowsDM[0]<2)  {
			fprintf(g,"\n\n DERIVED parameters - psi2,psi3,psi4,...\n"); strcpy(lbl,"psi"); sum=0;
			fprintf(g,"\n        Site                        %s(t)  Std.err     95%% conf. interval\n",lbl);
			for (k=1; k<PrmyPeriods; k++) {
				maxsite=1;
				for (site=0; site<N; site++) if ((psiout[site][k]-psiout[0][k])>.0001) maxsite=N;
				for (site=0; site<maxsite; site++) {                //
					ci1=psiout[site][k]-1.96*psioutse[site][k]; ci2=psiout[site][k]+1.96*psioutse[site][k];
                    fprintf(g,"%s(%2d)  %4d %8s      :     %8.4f%8.4f   %8.4f -%7.4f\n",lbl,k+1,site+1,
                            SingSpec.sitename[site],psiout[site][k],psioutse[site][k],ci1,ci2);
				}
			}
		}
		        //      gam                        eps
		if (SingSpec.NrowsDM[1]<1 || SingSpec.NrowsDM[2]<1) { seasn=1;
			strcpy(lbl,"eps"); if (SingSpec.NrowsDM[1]<1) strcpy(lbl,"gam");
			fprintf(g,"\n\n DERIVED parameters - %s1,%s2,%s3,...\n",lbl,lbl,lbl);
			fprintf(g,"\n        Site                        %s(t)  Std.err     95%% conf. interval\n",lbl);
			double psi_ip1,psi_i,gam_i;
			for (site=0; site<N; site++) {                             //  compute either gam or eps from eq.
				for (k=1; k<PrmyPeriods; k++) {                        //  psi(i+1)=psi(i)*(1-eps)+(1-psi(i))*gam
					for (i=0; i<NPar; i++) {
						psi_i = getParam(k-1,0,k-1,site,k-1,Params,0,1);         // gam_i is either gam or eps
						psi_ip1 = getParam(k,0,k,site,k,Params,0,1);             //  eps is either eps or gam
						gam_i = getParam(PrmyPeriods+k-1,1+(SingSpec.NrowsDM[1]<1),k-1,site,k-1,Params,0,1);
						if (SingSpec.NrowsDM[1]<1) eps=(psi_ip1-psi_i*(1-gam_i))/(1-psi_i+verysmallnumber);
						else 					   eps=1-(psi_ip1-(1-psi_i)*gam_i)/(psi_i+verysmallnumber);
						Params[i]-=verysmallnumber;
						psi_i = getParam(k-1,0,k-1,site,k-1,Params,0,1);
						psi_ip1 = getParam(k,0,k,site,k,Params,0,1);
						gam_i = getParam(PrmyPeriods+k-1,1+(SingSpec.NrowsDM[1]<1),k-1,site,k-1,Params,0,1);
						if (SingSpec.NrowsDM[1]<1) eps2=(psi_ip1-psi_i*(1-gam_i))/(1-psi_i+verysmallnumber);
						else 					   eps2=1-(psi_ip1-(1-psi_i)*gam_i)/(psi_i+verysmallnumber);
						grad[i]=(eps2-eps)/verysmallnumber;
						Params[i]+=verysmallnumber;
					}
					//           matrix multiply grad * varcov(betas)
					for (i=0; i<NPar; prd[i++]=sum) for (j=0,sum=0; j<NPar; j++) sum+=grad[j]*covar[j][i];
					//           matrix multiply {grad * varcov(betas)} * grad'
					for (i=0,vc=0; i<NPar; i++) vc+=prd[i]*grad[i];
					psiout[site][k]=eps; psioutse[site][k]=sqrt(vc);
				}
			}
			for (k=1; k<PrmyPeriods; k++) {
				maxsite=1;
				for (site=0; site<N; site++) if ((psiout[site][k]-psiout[0][k])>.0001) maxsite=N;
				for (site=0; site<maxsite; site++) {                //
					ci1=psiout[site][k]-1.96*psioutse[site][k]; ci2=psiout[site][k]+1.96*psioutse[site][k];
                    fprintf(g,"%s(%2d)  %4d %8s      :     %8.4f%8.4f   %8.4f -%7.4f\n",lbl,k,site+1,
                            SingSpec.sitename[site],psiout[site][k],psioutse[site][k],ci1,ci2);
				}
			}
		}
        fprintf(g,"\n\n DERIVED parameters - lam2,lam3,lam4,...\n"); strcpy(lbl,"lam"); sum=0;
        fprintf(g,"\n        Site                        %s(t)  Std.err     95%% conf. interval\n",lbl);
        for (k=1; k<PrmyPeriods; k++) {            //  compute psi(i) for each primary period, site
			maxsite=1;
			for (site=0; site<N; site++) if ((lamout[site][k]-lamout[0][k])>.0001) maxsite=N;
			for (site=0; site<maxsite; site++) {                //
                ci1=lamout[site][k]-1.96*lamoutse[site][k]; ci2=lamout[site][k]+1.96*lamoutse[site][k];
                j=0; if (site>0) if (lamout[site][k]!=lamout[site-1][k]) j=1; if (site==0) j=1;
                fprintf(g,"%s(%2d)  %4d %8s      :     %8.4f%8.4f   %8.4f -%7.4f\n",lbl,k+1,site+1,
                        SingSpec.sitename[site],lamout[site][k],lamoutse[site][k],ci1,ci2);
            }
			sum+=SingSpec.det_hist_frq[site];
        }
		chisq_test();

        delete [] prd; delete [] grad; delete [] grad2;
        for (site=0; site<N; site++) {
            delete [] psiout[site]; delete [] psioutse[site];  delete [] lamout[site]; delete [] lamoutse[site];
        }
        delete [] psiout; delete [] psioutse; delete [] lamout; delete [] lamoutse;
    }  //     end if singspec.model < 5

	//   *****************   GOODNESS OF FIT TEST **************************
    if (nboot2>0) {
        double estpsi,estp,estgam,esteps,TestStat2=0,TestStat2s=0,avgts=0,mints=.1e40,maxts=-9e9,*psave;
        int isim,occ,p_val=0, detect, **OrigData, **OrigOpenData; SingSpec.simulating=1;
        psave=new double[NPar]; for (i=0; i<NPar; i++) psave[i]=Params[i];
		//  copy orig data
        OrigData = new int*[N]; OrigOpenData = new int*[N];  // store Original Data and Parameter Estimates
        for (site=0; site<N; site++) {
            OrigData[site] = new int[T]; OrigOpenData[site]=new int[T];
            for (srvy=0; srvy<T; srvy++) {
				OrigData[site][srvy] = SingSpec.Data[site][srvy];
				if (srvy<SingSpec.PrmyPeriods) OrigOpenData[site][srvy] = SingSpec.OpenData[site][srvy];
			}
        }
        TestStat2=gof_open(psave,SingSpec.Verbose);
		if (SingSpec.Verbose>1) printf("teststat=%f\n",TestStat2);
        for (isim=0; isim<nboot2; isim++) {
            for (i=0; i<NPar; i++) Params[i]=psave[i];
            for (site=0; site<N; site++) {
                estpsi=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
                occ=(Random()<estpsi);
                for (seasn=srvy=0; seasn<PrmyPeriods; seasn++) {
					SingSpec.OpenData[site][seasn]=0; detect=0;
                    for (j=0; j<SingSpec.SecPeriods[seasn]; j++,srvy++) {
                        estp=getParam(2*PrmyPeriods-1,3,srvy,site,srvy,Params,0,SingSpec.LinkFn);
                        SingSpec.Data[site][srvy]=occ*(Random()<estp);
						if (OrigData[site][srvy]<0) SingSpec.Data[site][srvy]=-1;
						if (SingSpec.Data[site][srvy]>0) detect=1;
                    }
                    SingSpec.OpenData[site][seasn]=detect;
					if (OrigOpenData[site][srvy]<0) SingSpec.OpenData[site][srvy]=-1;
                    if (seasn<(PrmyPeriods-1)) {
                        estgam=getParam(1+seasn,1,seasn,site,seasn,Params,0,SingSpec.LinkFn);
                        esteps=getParam(PrmyPeriods+seasn,2,seasn,site,seasn,Params,0,SingSpec.LinkFn);
                        occ=occ*(Random()<(1.-esteps))+(1-occ)*(Random()<estgam);
                    }
                }
            }

            i=optmiz(fnptr,NPar, 7, SingSpec.maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work);
            TestStat2s=gof_open(Params,SingSpec.Verbose);
            if (SingSpec.Verbose>1) printf("  simulation %d teststat=%f maxll=%f\n",isim,TestStat2s,MaxLL);
            avgts+=TestStat2s; if (TestStat2s>=TestStat2) p_val++;
            if (TestStat2s<mints) mints=TestStat2s;
            if (TestStat2s>maxts) maxts=TestStat2s;
        }
		SingSpec.simulating=0;
        fprintf(g,"\n------------------------------------------------------\n");
        fprintf(g,"Goodness of fit results \n\n");
        avgts/=(double)nboot2;
        fprintf(g,"Test Statistic (data)    = %15.4f\n",TestStat2);
        fprintf(g,"\nFrom %d parametric bootstraps...\n",nboot2);
        fprintf(g,"Probability of test statistic >= observed = %6.4f\n",p_val/(nboot2+1.0));
        fprintf(g,"Lowest  simulated Test Stat = %15.4f\n",mints);
        fprintf(g,"Average simulated Test Stat = %15.4f\n",avgts);
        fprintf(g,"Highest simulated Test Stat = %15.4f\n",maxts);
        fprintf(g,"Estimate of c-hat = %6.4f   (=TestStat/AvgTestStat)\n",TestStat2/avgts);
        fprintf(g,"------------------------------------------------------\n");

		//  retrieve orig data
        for (site=0; site<N; site++) {
            for (srvy=0; srvy<T; srvy++) {
				SingSpec.Data[site][srvy] = OrigData[site][srvy];
				if (srvy<SingSpec.PrmyPeriods) SingSpec.OpenData[site][srvy] = OrigOpenData[site][srvy];
			}
        }
    }

    // delete dynamic variables
    for (ipar=0; ipar<NPar; ipar++) delete [] covar[ipar];  delete[] covar;
    for (seasn=0; seasn<=PrmyPeriods; seasn++) delete [] SingSpec.psibar[seasn];  delete[] SingSpec.psibar;
    delete[] Grad1; delete[] Hess; delete[] Work;
}

double psi_i(int site, int seasn, double psi0, double *Params) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	extern TSingSpec SingSpec;	int PrmyPeriods=SingSpec.PrmyPeriods; double psiI,gam,eps;

	switch (SingSpec.Model) {
	case 2: psiI = getParam(seasn,0,seasn,site,0,Params,0,SingSpec.LinkFn);
	        gam = getParam(PrmyPeriods+seasn-1,1,seasn-1,site,seasn-1,Params,0,1); // psi,gam
            eps = 1-gam+(gam-psiI)/psi0;
	        break;
	case 3: psiI = getParam(seasn,0,seasn,site,0,Params,0,SingSpec.LinkFn);
		    eps = getParam(PrmyPeriods+seasn-1,1,seasn-1,site,seasn-1,Params,0,1); // psi,eps
            gam = (psi0*eps-psi0+psiI)/(1-psi0);
	        break;
	case 4: gam = getParam(seasn,0,seasn,site,seasn-1,Params,0,1);  // gam=1-eps parmeritization
			eps=1-gam;
	        psiI= psi0*(1-eps)+(1-psi0)*gam;
	        break;
	default:
	        gam = getParam(seasn,1,seasn-1,site,seasn-1,Params,0,1);
	        eps = getParam(PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
			psiI= psi0*(1-eps)+(1-psi0)*gam;
	}
	return(psiI);
}

