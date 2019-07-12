#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "SingAnal.h"
#include "rngs.h"
#define PDL 0

void roylegof(double *Params) {
	extern TSingSpec SingSpec; 
	double cdtr(double x, int i, int *j);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	
	FILE *g; g=SingSpec.g;
	int i,j,k,t,ndf=0,T=SingSpec.T,lmt2,lmt1=SingSpec.rlmt; 
	double chi,sumchi=0,*f,**frq0,**expcnt,xlam,p,smexpcnt=0,smfrq0=0,eps=1e-12;
	f=new double[lmt1+1];	frq0=new double*[T]; expcnt=new double*[T];
	for (i=0; i<T; i++) {
		frq0[i]=new double[lmt1+1]; expcnt[i]=new double[lmt1+1];
		for (j=0; j<=lmt1; j++) frq0[i][j]=expcnt[i][j]=0;
	}
	fprintf(g,"\nAbundance GOF test\n"); fprintf(g,  "------------------\n");	
	for (i=0; i<SingSpec.N; i++) {
		xlam=getParam(0,0,0,i,0,Params,0,expLnk); f[0]=exp(-xlam); lmt2=lmt1;
		for (k=1; k<=lmt1; k++) {
			f[k]=f[k-1]*xlam/(double)k;  //f=dpois(j,lam)=lam^k*exp(-lam)/k!
			if (f[k]<eps && lmt2==lmt1) lmt2=k;			
		}
		for (t=0; t<SingSpec.T; t++) {
			p=getParam(1+t,1,t,i,t,Params,0,logitLnk);	
			if (SingSpec.Data[i][t]>=0 && SingSpec.Data[i][t]<=lmt1) frq0[t][SingSpec.Data[i][t]]++;			
			for (k=0; k<=lmt2; k++) {	//  for each possible site density...
				for (j=0; j<=k; j++) //  for each possible count...
					expcnt[t][j]+=f[k]*SingSpec.choose[k][j]*pow(p,j)*pow(1.-p,k-j);
			}
		}
	}
	for (j=0; j<T; j++) { fprintf(g,"survey %d\n",j+1);
   	  fprintf(g," count      Observed        Expected       chi-square\n");
	  smfrq0=smexpcnt=0;
	  for (i=0; i<lmt1; i++) {
		if (expcnt[j][i]>0.00005) {
			chi=0;
			if (expcnt[j][i]>=2.0) {
				chi=(frq0[j][i]-expcnt[j][i])*(frq0[j][i]-expcnt[j][i])/expcnt[j][i];
				sumchi+=chi; ndf++;
			}
			else { smfrq0+=frq0[j][i]; smexpcnt+=expcnt[j][i];}
			fprintf(g,"  %2d    %12.4f    %12.4f    %12.4f\n",i,frq0[j][i],expcnt[j][i],chi);
		}
	  }
	}
    if (smfrq0>0 || smexpcnt>0) {   //  add in pooled obs,exp,chi
		chi=(smfrq0-smexpcnt)*(smfrq0-smexpcnt)/smexpcnt;
		sumchi+=chi; ndf++;
		fprintf(g," pooled %12.4f    %12.4f    %12.4f\n",smfrq0,smexpcnt,chi);
	}
	delete []f;	delete []frq0; delete []expcnt; 
	fprintf(g,"\nTotal chi-square  :%f\nDegrees of freedom:%d\n",sumchi,ndf);
	fprintf(g,"Probability       :%f\n",1-cdtr(sumchi,ndf,&i));
}


void RoyleMod(int NBoot, int goftest, int NPar, double *Params) { 
    // Data[][] contains species counts,  N is number of sites,  T is number of sampling occasions

    double JARLike(double *Params, int NPar);    double JARLikeNB(double *Params, int NPar);
    double NHetLike(double *Params, int NPar);   double NHetLikeNB(double *Params, int NPar);
    double ran1(long *idum);   
    int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, 
		 int iopt, double *Params, double *Hess, double *Grad, double *MaxLL, double *Work);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	void roylegof(double *params);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
	
    int i, j, N, T, ii, srvy, k,site=0,ireal=0, altmodel=0;     // delcare some index counters
    extern TSingSpec SingSpec;  FILE *g; g=SingSpec.g;
    double xlam,psise,xp,xpsi,xn, xnse, xlam1,xpsi1,xpsi2,xn1,xn2,N2=0,a;
	N=SingSpec.N; T=SingSpec.T; 
	
	if (SingSpec.Realname[1][0][0]=='c') altmodel=1;
	if (SingSpec.Realname[0][0][0]=='u') altmodel+=2;
    printf("altmodel=%d\n",altmodel);
    //  Set up vectors for likelihood function
    
    SingSpec.choose = new double*[SingSpec.rlmt+1];  
    SingSpec.choose[0]=new double[1]; SingSpec.choose[0][0]=1;
    SingSpec.choose[1]=new double[2]; SingSpec.choose[1][0]=SingSpec.choose[1][1]=1;
    for (i=2; i<=SingSpec.rlmt; i++) {
        SingSpec.choose[i]=new double[i+1];
        SingSpec.choose[i][0]=SingSpec.choose[i][i]=1;
        SingSpec.choose[i][1]=SingSpec.choose[i][i-1]=i; 
        for (j=2; j<=(i/2); j++) 
            SingSpec.choose[i][j]=SingSpec.choose[i][i-j]=SingSpec.choose[i][j-1]/j*(i-j+1);
    }
    SingSpec.maxx = new int[N+1];
    for (i=0; i<N; i++) {
        for (j=0,SingSpec.maxx[i]=SingSpec.Data[i][0]; j<T; j++)  {
            if (SingSpec.Data[i][j]>SingSpec.maxx[i]) SingSpec.maxx[i]=SingSpec.Data[i][j];
        }   
		N2+=SingSpec.det_hist_frq[i];
	}
    fprintf(g,"rlmt=%d\n\n",SingSpec.rlmt);

    if(SingSpec.UseAmoeba) DoAmoeba(Params,NPar,JARLike,1.0);
   
             // set up covariance matrix of beta parameters
    double **covar, **Grad1, **vc; double **prd,xp2,xlam2,eps=1.e-10,sum;
	int Nreal=T+1; if (altmodel==3) Nreal=T+2;
    covar = new double*[NPar+1]; Grad1 = new double*[Nreal+1]; vc=new double*[Nreal+1];
    for (i=0; i<=Nreal; i++) { 
		Grad1[i] = new double[NPar+1]; for (j=0; j<=NPar; j++) Grad1[i][j]=0; 
		vc[i] = new double[Nreal+1]; for (j=0; j<=Nreal; j++) vc[i][j]=0; 
	}  
    for (i=0; i<=NPar; i++) { covar[i] = new double[NPar+1]; for (j=0; j<=NPar; j++) covar[i][j]=0; }  
    double *Hess,*Work,MaxLL;    Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2];
    if (SingSpec.Verbose>0) printf("optimize...\n");
    if(altmodel==0) { //  Royle Point-count model
		fprintf(g,"\n\nRoyle Count Model (Poisson)\n");
        optmiz(JARLike,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
		printf("JARLike\n");
		//if (goftest>0) roylegof
		printf("optmiz done\n");
        prntBetas(JARLike,NPar, Params, Work, MaxLL, covar);
    }
    if(altmodel==1) {//  Nichols heterogeneity model
		fprintf(g,"\n\nNichols Heterogeneity Model (Poisson)\n");
		printf("NHetLike\n");
        optmiz(NHetLike,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
		//if (goftest>0) roylegof();
        prntBetas(NHetLike,NPar, Params, Work, MaxLL, covar);
    }
    if(altmodel==2) { //  Royle Point-count model w/ neg. binomial prior dist.
		fprintf(g,"\n\nRoyle Count Model (Neg. Bin.)\n");
		printf("JARLikeNB\n");
        optmiz(JARLikeNB,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
		//if (goftest>0) roylegof();
        prntBetas(JARLikeNB,NPar, Params, Work, MaxLL, covar);
    }
    if(altmodel==3) {//  Nichols heterogeneity model w/ neg bin
		fprintf(g,"\n\nNichols Heterogeneity Model (Neg. Bin.)\n"); 
		printf("NHetLikeNB\n");
        optmiz(NHetLikeNB,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
		//if (goftest>0) roylegof();
        prntBetas(NHetLikeNB,NPar, Params, Work, MaxLL, covar);
    }	
	if (SingSpec.Verbose>0) printf("done... printing estimates...\n");
    fprintf(g,"\n============================================================\n");
	print_individual_estimates2(g,0,0,0,1,covar,Params,3,0,0); ii=1;// lambda
	if (SingSpec.NrowsDM[0]>1) {
		if (SingSpec.Realname[0][1][0]=='a') 
			print_individual_estimates2(g,ii++,0,1,1,covar,Params,3,0,0); // a
		else
			print_individual_estimates2(g,ii++,0,1,1,covar,Params,SingSpec.LinkFn,0,0); // psi
	}    
	fprintf(g,"\n============================================================\n");
	ii=1; print_individual_estimates2(g,ii,1,0,T,covar,Params,SingSpec.LinkFn,0,0); // r or c
    fprintf(g,"\n============================================================\n");
    
    if (SingSpec.novar>3) {
        if (SingSpec.Verbose>0) printf("computing var/cov matrix...\n");
        prd=new double*[Nreal]; 
        for (ii=0; ii<Nreal; ii++) prd[ii]=new double[NPar];
		for (ii=0; ii<NPar; ii++) {
			xlam =getParam(0,0,0,site,0,Params,1,1); Params[ii]+=eps;
			xlam2=getParam(0,0,0,site,0,Params,1,1); Params[ii]-=eps;
			Grad1[0][ii]=(xlam2-xlam)/eps; ireal=0;
			if (altmodel==2 || altmodel==3) {
				xlam =getParam(1,0,1,site,0,Params,1,1); Params[ii]+=eps;
				xlam2=getParam(1,0,1,site,0,Params,1,1); Params[ii]-=eps;
				Grad1[1][ii]=(xlam2-xlam)/eps; ireal=1;
			}
			for (srvy=0; srvy<T; srvy++) {
				xp =getParam(++ireal,1,srvy,site,srvy,Params,0,1); Params[ii]+=eps;
				xp2=getParam(  ireal,1,srvy,site,srvy,Params,0,1); Params[ii]-=eps;
				Grad1[ireal][ii]=(xp2-xp)/eps;
			}
		}
		//           matrix multiply grad * varcov(betas)
		for (i=0; i<Nreal; i++) for (j=0; j<NPar; prd[i][j++]=sum) for (k=0,sum=0; k<NPar; k++) sum+=Grad1[i][k]*covar[k][j];
		//           matrix multiply {grad * varcov(betas)} * grad'
		for (i=0; i<Nreal; i++) for (j=0; j<Nreal; vc[i][j++]=sum) for (k=0,sum=0; k<NPar; k++) sum+=prd[i][k]*Grad1[j][k];
		if (SingSpec.NrowsDM[0]>=0) {
			if (SingSpec.Verbose>0) printf("printing derived parameters...\n");
			fprintf(g,"\nDerived  parameter psi - occupancy estimate  std.err  95%% confidence interval\n");
			fprintf(g,"--------------------------           --------  -------  ------------------------\n");
			for (site=0; site<N; site++) 
				if (SingSpec.NParKK[0]>1 || site==0) {			
					xlam=getParam(0,0,0,site,0,Params,1,1); xpsi=1-exp(-xlam);
					if (SingSpec.Realname[0][0][0]=='u') {
						a=xlam*xlam/getParam(1,0,1,site,0,Params,1,1);
						xpsi=1-pow(a/(xlam+a),a);
					}
					for (i=0; i<SingSpec.NParKK[0]; i++) prd[0][i]=0;
					for (i=0; i<SingSpec.NParKK[0]; i++) {
						Params[i]+=eps;
						xlam=getParam(0,0,0,site,0,Params,1,1); xpsi1=1-exp(-xlam);
						if (SingSpec.Realname[0][0][0]=='u') {
							a=xlam*xlam/getParam(1,0,1,site,0,Params,1,1);
							xpsi1=1-pow(a/(xlam+a),a);
						}
						Grad1[0][i]=(xpsi-xpsi1)/eps; Params[i]-=eps; 
						for (j=0; j<SingSpec.NParKK[0]; j++) prd[0][j]+=Grad1[0][i]*covar[i][j];
					}			
					for (i=0,sum=0; i<SingSpec.NParKK[0]; i++) sum+=prd[0][i]*Grad1[0][i];
					psise=sqrt(sum); xpsi1=xpsi-psise*1.96; xpsi2=xpsi+psise*1.96;
					//Params=xpsi-1.96*sqrt(covar[0][0]); xlam1=getParam(0,0,0,site,0,Params,1,1);
					//Params=xpsi+1.96*sqrt(covar[0][0]); xlam2=getParam(0,0,0,site,0,Params,1,1);
					//Params=xpsi;
					//xpsi=1-exp(-xlam); psise=sqrt(exp(-xlam-xlam)*vc[0][0]);  // vc=var-cov matrix of real params
					//xpsi1=1-exp(-xlam1); xpsi2=1-exp(-xlam2);                 // covar=var-cov matrix of betas
					fprintf(g,"psi(%30s) : %6.4f   %6.4f %9.4f -%7.4f\n", 
						SingSpec.sitename[site],xpsi,psise,xpsi1,xpsi2);
				}
			fprintf(g,"\nDerived  parameter N - Tot. Abund. estimate  std.err  95%% confidence interval\n");
			fprintf(g,"--------------------------           --------  -------  ------------------------\n");
			for (site=0; site<N; site++) 
				if (SingSpec.NParKK[0]>1 || site==0) {			
					xlam=getParam(0,0,0,site,0,Params,1,1); xn=Params[0];
					Params[0]=xn-1.96*sqrt(covar[0][0]); xlam1=getParam(0,0,0,site,0,Params,1,1);
					Params[0]=xn+1.96*sqrt(covar[0][0]); xlam2=getParam(0,0,0,site,0,Params,1,1);
					Params[0]=xn;
					xn=xlam*N2; xnse=sqrt(N2*N2*vc[0][0]);
					xn1=xlam1*N2; xn2=xlam2*N2;
					fprintf(g,"N(%31s) : %6.2f   %6.2f %9.2f -%7.2f\n",SingSpec.sitename[site],xn,xnse,xn1,xn2);
				}
			}
	}
	//   new gof test
	if (SingSpec.NrowsDM[0]==1) roylegof(Params);
	
    if (1==2) {
		double *fest=(double*)malloc((1+SingSpec.rlmt)*sizeof(double));
		int Abund_Index_Model=SingSpec.TSpecificP; double lam,sumf,exlam;
		for (i=0; i<NPar; i++) { if(Params[i]<-22) Params[i]=-22; if(Params[i]>22) Params[i]=22;}
		for (i=0; i<N; i++) {      // print scaled f's	
			printf("area(%d)=%f\n",i,SingSpec.area[i]);
			lam=getParam(0,0,0,i,0,Params,1,1)*SingSpec.area[i]; 
			if (lam<.1e-8)lam=.1e-8;
			if(!Abund_Index_Model && lam>40) lam=40;  // avoid numerical errors if lam too large
			exlam=exp(-lam); fest[0]=exlam; sumf=exlam;
			for (N=1; N<=SingSpec.rlmt; sumf+=fest[N++]) fest[N]=fest[N-1]*lam/(double)N;//f=dpois(j,lam)=lam^N*exp(-lam)/N!
			printf("%4d %f \n",i,lam);
			// for (N=0; N<=SingSpec.rlmt; N++) printf(" %f",fest[N]/(sumf+.1e-8)); printf("\n");
		}
	}
    //fprintf(g,"\nVariance-Covariance Matrix\n");
    // Ouput the VC matrix
    //for (i=0; i<NPar; i++) {
    //for (j=0; j<NPar; j++)  fprintf(g,"%7.4f ",covar[i][j]); fprintf(g,"\n");
    //}
    // Use bootstrapping to estimate V-C matrix
    if (NBoot > 0) {
    } // end if(NBoot>0)
    
    fprintf(g,"------------------------------------------------------");
    printf("\ndone - check confirm window...\n");
    // delete dynamic variables
    for (i=0; i<NPar; i++) delete[] covar[i];
    for (i=0; i<Nreal; i++) delete[] Grad1[i];
    delete[] covar; delete[] Grad1;
    delete[] Hess;   delete[] Work;
}
