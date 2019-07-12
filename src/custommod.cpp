#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
//#include <search.h>

int dblcmp(const void *v1, const void *v2) { 
    if ((*(double *)v1 - *(double *)v2) < 0) return(-1);
    if ((*(double *)v1 - *(double *)v2) > 0) return(1);
    return(0);
}

void CustomMod(int NBoot, int LoFBoot, int NPar, double *Params, double poolcutoff)
{
    // Data[][] contains detection/nondetection infonew
    // Missing[][] contains whether observations were missing
    // SiteCov[][] contains site-specific covariates
    // SampCov[][][] contains sampling occasion covariates
    // N is number of sites
    // T is number of sampling occasions
    // NSiteCov is number of site covariates
    // NSampCov is number of sampling covariates
    // LoFTest indicates whether model is assessed for lack-of-fit
    
    double CMLike(double *Params, int NPar); 
    double gof(double *psi_hat, double **p_hat, double ***th_hat, double *pi_hat, double **p10_hat, double **b_hat, int prnt, double poolcutoff, int modtype);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    double ran1(long *idum);
    int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params, double *Hess, double *Grad, double *MaxLL, double *Work);
    bool InvertMatrix(double **a, int n); time_t time1,time2;
    void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
    void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
	double invlogit(double x); double logit(double x);
	void chisq_test(void);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
	void try_optmz_randiv(double (*LinkFn)(double *Params, int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);	
    
    // declare variables
    char starg[80], *tmpstr;
    int site, srvy, kk, ll, BS, i, j, k, N=SingSpec.N, T=SingSpec.T, ii, iBS=0, NN=0, iisite,
          **OrigData, nper=0, frstp, ndet, itemp,ier, gofmtype=0, NrealPar=0;
    long seed=SingSpec.seed;
    double *psi_hat, **p_hat, ***th_hat, *pi_hat, **p10_hat, **b_hat;
    double **covar, *Grad1, *Hess,*Work,MaxLL;
    double *OrigParams, AvTS, BootTestStat=0; 
    double *Bootpsi_hat, **Bootp_hat, ***Bootth_hat, *Bootpi_hat, **Bootp10_hat, **Bootb_hat, TestStat=0; 
    double **OrigSiteCov, ***OrigSampCov, *BTSsave, *Origexpval, *Origfrq;
    double **BootParams, *ObsParams, *AvParams, **BootCovar;
    double prd,*tmpprd,verysmallnumber=.000001,sum,vc,ci1,ci2;
    tmpstr=(char*)malloc(256);
	
    SingSpec.psibar=new double**[1]; SingSpec.psibar[0]=new double*[N];
    for (site=0; site<N; site++) SingSpec.psibar[0][site]=new double[2];
	
    fprintf(g,"\nCustom Model: %s\n\n",SingSpec.modname);   // output results
    if (SingSpec.FalsePos) fprintf(g,"  with 'false positive' detections\n\n");
	if (SingSpec.NrowsDM[0]==1) SingSpec.NMethods=1;	
    
    if (!SingSpec.FalsePos)
        for (site=0,sum=0; site<N; site++)                      //    loop for each site (site=0 to N)...
            for (srvy=0; srvy<T; srvy++)                              //   loop for each survey (srvy=0 to T)...
                if (SingSpec.Data[site][srvy]>1) SingSpec.Data[site][srvy]=1;          //   if data>2, data=1
    
    if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, CMLike, 1.0); 
    // set up covariance matrix of beta parameters
    
    covar = new double*[NPar]; Grad1 = new double[NPar];
    for (site=0; site<NPar; site++) {
        covar[site] = new double[NPar];
        for (srvy=0; srvy<NPar; srvy++) covar[site][srvy]=0;
    }
    Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2]; 
    if (SingSpec.Verbose>0) printf("calling optmiz...\n");
	
	try_optmz_randiv(CMLike, NPar, Params, Work, Hess, &MaxLL);
	
    //optmiz(CMLike,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess-1, Grad1-1, &MaxLL, Work-1);
    //optmiz(CMLike,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work);
    if (SingSpec.Verbose>0) printf("optmiz done\n");
    fprintf(g,"Number of significant digits   = %5.1f\n", Work[3]);
    if (Work[3]<SingSpec.LikeNRSig) {
        fprintf(g,"\n**** Numerical convergence may not have been reached.\n");
        fprintf(g,"     Parameter estimates converged to approximately \n");
        fprintf(g,"     %4.2f significant digits.\n\n",Work[3]);
    }
    if (SingSpec.LinkFn==1) fprintf(g,"\nModel has been fit using the logistic link.\n\n");
    else if (SingSpec.LinkFn==4) fprintf(g,"\nModel has been fit using the sin link.\n\n");
    else fprintf(g,"\nModel has been fit using the logistic link for p, cloglog for psi.\n\n");
    
    // output results
    prntBetas(CMLike,NPar, Params, Work, MaxLL, covar);
    
    fprintf(g,"\n============================================================\n");
    strcpy(starg,"Psi"); site=0; 
	print_individual_estimates2(g,0,0,site,1,covar,Params,SingSpec.LinkFn,0,SingSpec.N+SingSpec.nphantom);
    fprintf(g,"\n============================================================\n");
    if (SingSpec.FalsePos) {
        strcpy(starg,"p11"); site=0; print_individual_estimates2(g,1,1,site,0,covar,Params,SingSpec.LinkFn,0,N);
        strcpy(starg,"p10"); site=T; print_individual_estimates2(g,1+T,1,site,T,covar,Params,SingSpec.LinkFn,0,N);
        strcpy(starg,"b"); site+=T; print_individual_estimates2(g,1+T+T,1,site,T,covar,Params,SingSpec.LinkFn,0,N);
    }
    if (SingSpec.NMethods>1) {
        nper=T/SingSpec.NMethods;
        strcpy(starg,"theta"); site=1; 
        print_individual_estimates2(g,1,0,1,T/SingSpec.NMethods,covar,Params,SingSpec.LinkFn,0,N);
        fprintf(g,"\n============================================================\n");
    }
    site=0;
    if (!SingSpec.FalsePos) {
        strcpy(starg,"p"); 
        print_individual_estimates2(g,1+SingSpec.NMethods,1,site,0,covar,Params,SingSpec.LinkFn,0,N);
    }
    if (SingSpec.NrowsDM[2]>0) {
        strcpy(starg,"pi"); 
        print_individual_estimates2(g,1+SingSpec.NMethods+T,2,0,0,covar,Params,SingSpec.LinkFn,0,N);
    }
	if (SingSpec.prnt_derived) {
		fprintf(g,"\n============================================================\n");
		fprintf(g,"\n\n DERIVED parameter - Psi-conditional = [Pr(occ | detection history)]\n"); 
		fprintf(g,"\n        Site                     psi-cond  Std.err     95%% conf. interval\n");
		tmpprd=new double[NPar+1]; frstp=1; 
		int X0=0; if (SingSpec.FalsePos) X0=1;  //  if false pos, then check if at least 1 sure detection (2), otherwise (1)
		double tmpsave,p10,cond_psi,uncond_psi,cond_psi2,uncond_psi2,cond_psise,transformed_cond_psise; 
		TestStat = CMLike(Params, NPar); 
		for (site=0; site<N; site++) { ndet=0;
			for (srvy=0; srvy<T; srvy++) if (SingSpec.Data[site][srvy]>X0) ndet++;  //check for >= 1 (sure) detections			
			cond_psise=0; cond_psi=ci1=ci2=1;
			if (ndet<1) {
				uncond_psi=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); prd=1;
				if (SingSpec.FalsePos) 
					for (srvy=0; srvy<T; srvy++)  {
						p10=getParam(srvy+1+T,1,srvy+T,site,srvy,Params,0,SingSpec.LinkFn);
						if (SingSpec.Data[site][srvy]==0) prd*=(1-p10);
						else prd*=p10;
					}
				cond_psi=(SingSpec.expval[site]-(1-uncond_psi)*prd)/SingSpec.expval[site];
				for (i=0; i<NPar; i++) {              
					tmpsave=Params[i]; Params[i]+=verysmallnumber;
					TestStat = CMLike(Params, NPar); 
					uncond_psi2=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); prd=1;
					if (SingSpec.FalsePos) 
						for (srvy=0; srvy<T; srvy++)  {
							p10=getParam(srvy+1+T,1,srvy+T,site,srvy,Params,0,SingSpec.LinkFn);
							if (SingSpec.Data[site][srvy]==0) prd*=(1-p10);
							else prd*=p10;
						}
					cond_psi2=(SingSpec.expval[site]-(1-uncond_psi2)*prd)/SingSpec.expval[site];
					Grad1[i]=(cond_psi2-cond_psi)/verysmallnumber;
					Params[i]=tmpsave;
					TestStat = CMLike(Params, NPar); 
				}
				//  matrix multiply grad * varcov(betas)
				for (i=0; i<NPar; tmpprd[i++]=sum)for (j=0,sum=0; j<NPar; j++) sum+=Grad1[j]*covar[j][i];
				//           matrix multiply {grad * varcov(betas)} * grad'
				for (i=0,vc=0; i<NPar; i++) vc+=tmpprd[i]*Grad1[i];
				cond_psise=sqrt(vc);
				//   compute logit of std.err(cond.psi) = se(cond_psi) / (cond_psi * (1-cond_psi)) 
				transformed_cond_psise=cond_psise/cond_psi/(1-cond_psi); 
				ci1=invlogit(logit(cond_psi)-1.96*transformed_cond_psise); 
				ci2=invlogit(logit(cond_psi)+1.96*transformed_cond_psise);
			} 				//        Print derived psi's 
			fprintf(g,"psi-cond  %4d %8s        : %8.4f%8.4f   %8.4f -%7.4f \n",site+1,
					SingSpec.sitename[site],cond_psi,cond_psise,ci1,ci2);
			if (site==0) uncond_psi=cond_psi;
		} 
		delete [] tmpprd;
	}
	chisq_test();
    // Assess model fit here, but output later
    if (LoFBoot>0) { 
		if (SingSpec.NMethods>1)   gofmtype=1;  // multi-method
		if (SingSpec.Model==7)     gofmtype=2;  // spatial-correlation
		if (SingSpec.NrowsDM[1]<1) gofmtype=3;  // heterogeneous detections
		if (SingSpec.FalsePos)     gofmtype=4;  //  false-positives
        tmpstr = getenv("USERNAME");
        fprintf(g,"\n============================================================\n");
        fprintf(g,"\nAssessing Model Fit for Single-season model:\n");
        psi_hat = new double[N]; p_hat = new double*[N]; th_hat = new double**[N]; pi_hat = new double[N];
		p10_hat=new double*[N]; b_hat=new double*[N];
        frstp=1; if (SingSpec.NMethods>1) frstp=T/SingSpec.NMethods+1;
        for (site=0; site<N; site++) {
            p_hat[site] = new double[T]; p10_hat[site] = new double[T]; 
			b_hat[site] = new double[T]; th_hat[site]=new double*[T]; 
			for (j=0; j<T; j++) p10_hat[site][j]=b_hat[site][j]=0;
            psi_hat[site]=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);   
            if (LoFBoot<5 && site==0) fprintf(g,"\npsi_hat=%f\np_hat=",psi_hat[site]);
            for (srvy=0; srvy<T; srvy++) {
                p_hat[site][srvy]=getParam(srvy+frstp,1,srvy,site,srvy,Params,0,SingSpec.LinkFn); 
				if (SingSpec.FalsePos) {
					p10_hat[site][srvy]=getParam(frstp+T+srvy,1,T+srvy,site,srvy,Params,0,SingSpec.LinkFn);
					b_hat[site][srvy]=getParam(frstp+T+T+srvy,1,T+T+srvy,site,srvy,Params,0,SingSpec.LinkFn);
				}
                if (LoFBoot<5 && site==1) fprintf(g," %f",p_hat[site][srvy]);
				th_hat[site][srvy]=new double[3]; 
            }
            if (LoFBoot<5 && site==1) fprintf(g,"\nth_hat=");
            if (SingSpec.NMethods>1) {
                for (srvy=0; srvy<nper; srvy++) {
					th_hat[site][srvy][0]=getParam(srvy+1,0,srvy+1,site,srvy,Params,0,SingSpec.LinkFn); 
					if (LoFBoot<5 && site==1) fprintf(g," %f",th_hat[site][srvy][0]);
				}
				if (LoFBoot<5 && site==1) fprintf(g,"\n");
            }
        }
        TestStat = CMLike(Params, NPar);
        
        TestStat=gof(psi_hat, p_hat, th_hat, pi_hat, p10_hat, b_hat, 2, poolcutoff,gofmtype); //  4th arg=print flag
        if (tmpstr) if (strstr(tmpstr,"hines")!=NULL && SingSpec.Verbose>0) printf("lofboot=%d  DataTestStat = %f\n",LoFBoot,TestStat);
        if (SingSpec.Verbose>0) printf("\nlofboot=%d  DataTestStat = %f\n",LoFBoot,TestStat);
        
		OrigData=SingSpec.Data; Origfrq=SingSpec.det_hist_frq; //  save pointers to orig data...
		
		for (NN=site=0; site<N; site++) NN+=floor(0.5+SingSpec.det_hist_frq[site]);
		
		SingSpec.Data = new int*[NN]; SingSpec.OpenData = new int*[NN]; SingSpec.det_hist_frq=new double[NN]; // store Original Data and Parameter Estimates
        for (site=0; site<NN; site++) { 
			SingSpec.Data[site] = new int[T]; SingSpec.OpenData[site]=new int[T]; SingSpec.det_hist_frq[site]=1; 
			for (srvy=0; srvy<T; srvy++) SingSpec.Data[site][srvy]=-1;
		}
		//         redimension psibar for expanded det. histories
		for (i=0; i<N; i++) delete[] SingSpec.psibar[0][i]; 
		SingSpec.psibar[0]=new double*[NN];
		for (site=0; site<NN; site++) {
			SingSpec.psibar[0][site]=new double[2];
			SingSpec.psibar[0][site][0]=SingSpec.psibar[0][site][1]=0;
		}
		for (i=0; i<6; i++) NrealPar+=SingSpec.NrowsDM[i];
		
		for (i=0; SingSpec.realParmEst[i]!=NULL; i++) delete[] SingSpec.realParmEst[i]; delete SingSpec.realParmEst; 
		//printf("redim realParmEst...(%d to %d)\n",i,NrealPar);
		// redim realParmEst...
		SingSpec.realParmEst = new double*[NrealPar+1];      
		for (i=0; i<NrealPar; i++) SingSpec.realParmEst[i] = new double[NN]; SingSpec.realParmEst[NrealPar]=NULL;
		Origexpval = SingSpec.expval; SingSpec.expval=new double[NN];   //  redim expval...
	
        OrigParams = Params; Params = new double[NPar+1]; //  Save Params as OrigParams, allocate new Params vector...
        for (i=0; i<NPar; i++) Params[i] = (OrigParams[i]>4.0?4.0:OrigParams[i]);
		
        int maxfn=SingSpec.maxfn, p_value=0, seasn, nseasns, meth, dsum,nmiss_seasn,locc=1; AvTS = 0.0;  // add in bootstrap loop
        double LoBootTestStat=9.99999e44, HiBootTestStat=-LoBootTestStat;
        BTSsave=new double[LoFBoot];  fprintf(g,"\n"); 
		nseasns=1; if (SingSpec.NMethods>1) nseasns=T/SingSpec.NMethods;
        Bootpsi_hat = new double[NN]; Bootp_hat = new double*[NN]; Bootth_hat = new double**[NN]; 
		Bootpi_hat = new double[NN];  Bootp10_hat=new double*[NN]; Bootb_hat=new double*[NN];	
		for (site=0; site<NN; site++) {
            Bootp_hat[site] = new double[T]; Bootp10_hat[site] = new double[T]; 
			Bootb_hat[site] = new double[T]; Bootth_hat[site]=new double*[T]; 	
			for (j=0; j<T; j++) { Bootp10_hat[site][j]=Bootb_hat[site][j]=0; Bootth_hat[site][j]=new double[3]; }
		}
        for(BS=0; BS<LoFBoot; BS++) { time(&time1);       // generate Data, missing values are fixed
			SingSpec.simulating=1;
            if (SingSpec.Verbose>0) printf("                              %d/%d\r",BS+1,LoFBoot);
            for (site=iisite=0; site<N; site++) {
			  for (ii=0,locc=1; ii<floor(0.5+Origfrq[site]); ii++,iisite++) {
                if (ran1(&seed)<psi_hat[site]) {    // if site occupied
                    for (srvy=0; srvy<T; srvy++) {
						seasn=(int)srvy/SingSpec.NMethods;
						if (SingSpec.NMethods>1) {
							if ((srvy % SingSpec.NMethods)==0) locc=ran1(&seed)<th_hat[site][seasn][0]; // if first method of survey...// sim local site use 
						}
						if (gofmtype==2) locc=ran1(&seed)<th_hat[site][srvy][0];
                        if (OrigData[site][srvy]>=0) {
							SingSpec.Data[iisite][srvy] =(ran1(&seed) < (p_hat[site][srvy]*locc));
							if (SingSpec.FalsePos) SingSpec.Data[iisite][srvy]+=(ran1(&seed) < b_hat[site][srvy]);
						}
					}
                } else {    // if site vacant
                    for (srvy=0; srvy<T; srvy++) {
                        if (OrigData[site][srvy]>=0) SingSpec.Data[iisite][srvy] = 0; //  don't change missing data
						if (SingSpec.FalsePos) SingSpec.Data[iisite][srvy] = (ran1(&seed) < p10_hat[site][srvy]);
                    }
                }
				for (srvy=seasn=0; seasn<nseasns; seasn++) {             //  for multi-method... treat methods as secondary periods...
					for (meth=dsum=nmiss_seasn=0; meth<SingSpec.NMethods; srvy++,meth++) {
						if (SingSpec.Data[iisite][srvy]>0 && !SingSpec.FalsePos) SingSpec.Data[iisite][srvy]=1;
						if (SingSpec.Data[iisite][srvy]>MISSNG) dsum+=SingSpec.Data[iisite][srvy];
						else  nmiss_seasn++;
					}
					SingSpec.OpenData[iisite][seasn] = MISSNG;
					if (nmiss_seasn<SingSpec.NMethods) SingSpec.OpenData[iisite][seasn] = (dsum>0);
				}
			  }
            } 
            // data generated, now fit model                                 set init values for parms
            for (i=0; i<NPar; i++) Params[i] = OrigParams[i];  // this line moved from inside if stmt below
            if (LoFBoot<11 && SingSpec.Verbose>=0) {
                fprintf(g,"simulation number %d data:\n",BS+1);
                for (site=0; site<NN; site++) {
                    fprintf(g,"%2d:",site+1);
                    for (srvy=0; srvy<T; srvy++) fprintf(g,"%c",(SingSpec.Data[site][srvy]>=0?'0'+SingSpec.Data[site][srvy]:'.'));
                    fprintf(g," ");
                    for (srvy=0; srvy<nseasns; srvy++) fprintf(g,"%c",(SingSpec.OpenData[site][srvy]>=0?'0'+SingSpec.OpenData[site][srvy]:'.'));
					fprintf(g,"\n");
                }
            }
			SingSpec.N=NN;
			
            if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, CMLike, 0.2);  //*********** optimize **********
            ier=optmiz(CMLike,NPar, SingSpec.LikeNRSig, maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work);
			
            for (site=0; site<NN; site++) {
                Bootpsi_hat[site]=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
                for (srvy=0; srvy<T; srvy++) 
                    Bootp_hat[site][srvy]=getParam(srvy+frstp,1,srvy,site,srvy,Params,0,SingSpec.LinkFn);
				if (SingSpec.NMethods>1)
					for (srvy=0; srvy<nper; srvy++) 
						Bootth_hat[site][srvy][0]=getParam(srvy+1,0,srvy+1,site,srvy,Params,0,SingSpec.LinkFn);
            }            
            if (ier==131) {
                printf("iteration %d:Error - max function calls exceeded\n",BS+1);
                fprintf(g,"iteration %d:Error - max function calls exceeded\n",BS+1);
            }            
            else {
                BootTestStat = gof(Bootpsi_hat,Bootp_hat,Bootth_hat,Bootpi_hat,Bootp10_hat,Bootb_hat,1+(LoFBoot<5),poolcutoff,gofmtype); 
                if(BootTestStat >= TestStat) p_value++;
                AvTS += BootTestStat; BTSsave[iBS++]=BootTestStat;
                if (BootTestStat<LoBootTestStat) LoBootTestStat=BootTestStat;
                if (BootTestStat>HiBootTestStat) HiBootTestStat=BootTestStat;
            }
            time(&time2); 
            if (SingSpec.Verbose>0) {
                printf("          bootstrap %3d / %d  TS=%f c-hat = %f beta:%f\n",
                       BS+1,LoFBoot,BootTestStat,TestStat/(AvTS/(BS+1)),Params[0]);
            }
        }   // end bootstrap loop
        SingSpec.simulating=0;
        qsort((void*)BTSsave, (size_t)iBS,sizeof(double),dblcmp);
        AvTS /= (double)iBS;
        // Output results
        fprintf(g,"------------------------------------------------------\n");
        fprintf(g,"Test Statistic (data)    = %15.4f\n",TestStat);
        fprintf(g,"\nFrom %d parametric bootstraps...\n",iBS);
        fprintf(g,"Probability of test statistic >= observed = %6.4f\n",p_value/(iBS+1.0));
        fprintf(g,"Lowest  simulated Test Stat = %15.4f\n",LoBootTestStat);
        fprintf(g,"Average simulated Test Stat = %15.4f\n",AvTS);
        fprintf(g,"Median simulated Test Stat = %15.4f\n",BTSsave[(int)iBS/2]);
        fprintf(g,"Highest simulated Test Stat = %15.4f\n",HiBootTestStat);
        fprintf(g,"Estimate of c-hat = %6.4f   (=TestStat/AvgTestStat)\n",TestStat/AvTS);
        fprintf(g,"------------------------------------------------------\n");
        fprintf(g,"- Distribution of simulated test statistics ----------\n");
        fprintf(g," TestStat  Prop.\n");
        double cut=HiBootTestStat/20.,w; j=0; w=cut;
        for (i=0; i<iBS; i++) {
            if (BTSsave[i]>cut) { 
                fprintf(g,"%8.4f %8.4f:",cut+w/2.,j/(double)iBS); 
                for (k=0; k<=(int)(80*j/(double)iBS); k++) fprintf(g,"*"); fprintf(g,"\n");
                cut+=HiBootTestStat/20; j=0;
            }
            j++;
        }
        delete [] BTSsave;
        fprintf(g,"------------------------------------------------------\n");
        //printf("Test Statistic (data)    = %15.4f\n",TestStat);
        
        // copy orginal data and paramter estimates back into correct arrays;
        delete Params; Params=OrigParams; 
		delete SingSpec.expval; SingSpec.expval=Origexpval;
		delete SingSpec.det_hist_frq; SingSpec.det_hist_frq=Origfrq;
        for (site=0; site<NN; site++) delete [] SingSpec.Data[site]; delete SingSpec.Data; 
		SingSpec.Data=OrigData;    SingSpec.N=N;    
        // delete dynamic variables
        for (site=0; site<N; site++) { 
			delete [] p_hat[site]; delete [] p10_hat[site]; delete [] b_hat[site];  
			for (j=0; j<T; j++)  delete [] th_hat[site][j]; 
			delete [] th_hat[site]; 
		}
        for (site=0; site<NN; site++) { 
			delete [] Bootp_hat[site]; delete [] Bootp10_hat[site]; delete [] Bootb_hat[site];
			for (j=0; j<T; j++) delete [] Bootth_hat[site][j];
			delete [] Bootth_hat[site];
		}
        delete[] p_hat;    delete[] Bootp_hat; 
		delete [] p10_hat; delete [] Bootp10_hat;
		delete [] b_hat;   delete [] Bootb_hat;
		delete [] th_hat;  delete [] Bootth_hat;
		delete[] psi_hat;  delete[] Bootpsi_hat; 
        
    }/////  end of if(Lofboot) stmt ///////////////////////////////////////////////////////////
    
    if (NBoot > 0) {    // Use bootstrapping to estimate V-C matrix
        printf("bootstrap estimate v-c matrix...\n"); SingSpec.simulating=1;
        OrigData = new int*[N]; OrigSiteCov = new double*[N]; OrigSampCov = new double**[N];
        for (site=0; site<N; site++) {
            OrigData[site] = new int[T]; OrigSiteCov[site] = new double[SingSpec.NSiteCov];
            OrigSampCov[site] = new double*[T];
            for (srvy=0; srvy<T; srvy++) {
                OrigSampCov[site][srvy] = new double[SingSpec.NSampCov];
                OrigData[site][srvy] = SingSpec.Data[site][srvy];
                for (kk=0; kk<SingSpec.NSampCov; kk++) OrigSampCov[site][srvy][kk] = SingSpec.SampCov[kk][site][srvy];
            }
            for (srvy=0; srvy<SingSpec.NSiteCov; srvy++) OrigSiteCov[site][srvy] = SingSpec.SiteCov[site][srvy];
        }
        BootParams = new double*[NBoot];  ObsParams = new double[NPar];  AvParams = new double[NPar+1];
        for (site=0; site<NPar; site++) { ObsParams[site] = Params[site];  AvParams[site] = 0.0; }
        AvParams[NPar] = 0.0;     // this is where PsiBar is stored
        
        for (int bb=0; bb<NBoot; bb++) {
            BootParams[bb] = new double[NPar+1];
            // Bootstrap data
            for (site=0; site<N; site++) {
                itemp = (int) (N*ran1(&seed));
                if (SingSpec.Verbose>0) printf("%3d",itemp);
                for (srvy=0; srvy<T; srvy++) {
                    SingSpec.Data[site][srvy] = OrigData[itemp][srvy];
                    if (SingSpec.Verbose>0) printf(" %2d",SingSpec.Data[site][srvy]);
                    for (kk=0; kk<SingSpec.NSampCov; kk++) {
                        SingSpec.SampCov[kk][site][srvy] = OrigSampCov[itemp][srvy][kk];
                    }
                }
                if (SingSpec.Verbose>0) printf("\n");
                for (srvy=0; srvy<SingSpec.NSiteCov; srvy++) {
                    SingSpec.SiteCov[site][srvy] = OrigSiteCov[itemp][srvy];
                }
            }
            
            // setup and minimise
            for (srvy=0; srvy<NPar; srvy++) Params[srvy] = ObsParams[srvy];
            if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, CMLike, 0.2);
            optmiz(CMLike,NPar, SingSpec.BootNRSig, SingSpec.maxfn, 0, Params, Hess, Grad1-1, &MaxLL, Work);
    
            // store estimates
            for (site=0; site<NPar; site++) {
                BootParams[bb][site] = Params[site]; AvParams[site] += BootParams[bb][site];
            }
            // calc PsiBar for this bootstrap
            BootParams[bb][NPar] = 0.0;
            for (site=0; site<N; site++) BootParams[bb][NPar]+=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);     
            BootParams[bb][NPar] /= (double)N; AvParams[NPar] += BootParams[bb][NPar];
            if (SingSpec.Verbose>0) printf("%d %f %f %f %3.1f\n",bb+1,BootParams[bb][0],BootParams[bb][1],BootParams[bb][2],Work[3]);
        }// end bootstrap loop
        
        fprintf(g,"\navg param values for %d bootstrap iterations\n",NBoot);
        for (site=0; site<NPar+1; site++) {
            AvParams[site] /= NBoot; fprintf(g,"%d %f\n",site,AvParams[site]);
        }
        
        BootCovar = new double*[NPar+1];
        for (site=0; site<NPar+1; site++) BootCovar[site] = new double[NPar+1];
        for (site=0; site<NPar+1; site++) {
            for (srvy=site; srvy<NPar+1; srvy++) {
                BootCovar[site][srvy] = 0.0;
                for (int bb=0; bb<NBoot; bb++) 
                    BootCovar[site][srvy] += (BootParams[bb][site]-AvParams[site])*(BootParams[bb][srvy]-AvParams[srvy]);
                BootCovar[site][srvy] /= (NBoot-1.0);
                BootCovar[srvy][site] = BootCovar[site][srvy];
            }
        }
        
        // output VC matrix
        fprintf(g,"\nBootstrap Variance-Covariance Matrix\n");
        for (srvy=0; srvy<2; srvy++) 
            for (site=0; site<SingSpec.NParKK[srvy]; site++) {
                ll=(int)SingSpec.DMat[srvy][0][site];
                if(ll>0) fprintf(g,"%15s","intercept");
                if(ll<-1000) fprintf(g,"%15s",SingSpec.CovNames[-1001-ll]);
                if(ll<-2001) fprintf(g,"%15s",SingSpec.CovNames2[-2001-ll]);
            }
        fprintf(g,"%15s\n","Psi-Bar");
        
        for (site=0; site<NPar+1; site++) {
            for (srvy=0; srvy<NPar+1; srvy++) fprintf(g,"%15.4f",BootCovar[site][srvy]);
            fprintf(g,"\n");
        }
        fprintf(g,"\nBootstrap estimate of SE for overall proportion of sites occupied = %15.4f\n",sqrt(BootCovar[NPar][NPar]));
        
        // copy observed data back into global variables
        for (site=0; site<N; site++) {
            for (srvy=0; srvy<T; srvy++) {
                SingSpec.Data[site][srvy] = OrigData[site][srvy];
                for (kk=0; kk<SingSpec.NSampCov; kk++) 
                    SingSpec.SampCov[kk][site][srvy] = OrigSampCov[site][srvy][kk];
            }
            for (srvy=0; srvy<SingSpec.NSiteCov; srvy++) 
                SingSpec.SiteCov[site][srvy] = OrigSiteCov[site][srvy];
        }
        
        // tidy up dynamic variables
        for (site=0; site<N; site++) {
            delete[] OrigData[site]; delete[] OrigSiteCov[site];
            for (srvy=0; srvy<T; srvy++) delete[] OrigSampCov[site][srvy];
            delete[] OrigSampCov[site];
            
        }
        delete[] OrigData;  delete[] OrigSampCov;
        delete[] OrigSiteCov; delete[] ObsParams;  delete[] AvParams;
        for (site=0; site<NBoot; site++) delete[] BootParams[site];  delete[] BootParams;
        for (site=0; site<=NPar; site++) delete[] BootCovar[site];  delete[] BootCovar;
        
    }// end if(NBoot>0)
	SingSpec.simulating=0;
    if (SingSpec.Verbose>1) printf("delete dynamic variables... npar=%d\n",NPar);
    // delete dynamic variables 
    for (site=0; site<NPar; site++) delete[] covar[site]; delete[] covar;  
	delete[] Grad1; delete[] Hess;  delete[] Work;
}
