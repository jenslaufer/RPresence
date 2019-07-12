#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "rngs.h"
// #include <afx.h>
#include "SingAnal.h"
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif

void gofSD(int NPar, double *Params, int LoFBoot) {    
    // declare functions used here
    double OpenModLinkSD(double *Params, int NPar);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
	void chisq_test(void);
	double gof(double *psi, double **hat, double ***th, double *pi, double **p10, double **b, int prnt, double poolcutoff, int mtype);
	double ran1(long *idum);
	int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params, double *Hess, double *Grad, double *MaxLL, double *Work);
    int dblcmp(const void *v1, const void *v2);
	
    extern TSingSpec SingSpec;
	double *Grad1, *Hess,*Work,MaxLL;
	double *OrigParams, AvTS, BootTestStat=0; 
    double *Bootpsi_hat, **Bootp_hat, ***Bootth_hat, *Bootpi_hat, TestStat=0; 
    double *BTSsave;
    double *psi_hat,**p_hat,*pi_hat,***th_hat,poolcutoff=2;
	char *tmpstr;
	int i,j,k,ier,frstp=1,BS,gofmtype=2,T=SingSpec.T,N=SingSpec.N,srvy,site,**OrigData;
	FILE *g=SingSpec.g; time_t time1,time2;   long seed=SingSpec.seed;
	Work=new double[3*NPar+1]; Grad1=new double[NPar+1]; Hess=new double[NPar*(NPar+1)];

       tmpstr = getenv("USERNAME");
        fprintf(g,"\n============================================================\n");
        fprintf(g,"\nAssessing Model Fit for Single-season model:\n");
        psi_hat = new double[N]; p_hat = new double*[N]; th_hat = new double**[N]; pi_hat = new double[N];
        Bootpsi_hat = new double[N]; Bootp_hat = new double*[N]; Bootth_hat = new double**[N]; Bootpi_hat=new double[N];
        frstp=1+2*SingSpec.T;
        for (site=0; site<N; site++) {
            p_hat[site] = new double[T]; Bootp_hat[site] = new double[T]; 
			th_hat[site]=new double*[T]; Bootth_hat[site] = new double*[T];  
            psi_hat[site]=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);   // psi_hat[site]=get_psi(Params,NPar,site);
			pi_hat[site]=getParam(1+3*T,4,0,site,0,Params,0,SingSpec.LinkFn);  
            if (LoFBoot<5 && site==0) fprintf(g,"\npsi_hat=%f\np_hat=",psi_hat[site]);
            for (srvy=0; srvy<T; srvy++) {
                p_hat[site][srvy]=getParam(srvy+frstp,3,srvy,site,srvy,Params,0,SingSpec.LinkFn); // p_hat[site][srvy]=get_p(frstp+srvy,Params,NPar,site,srvy,dmnum,indxp);
                if (LoFBoot<5 && site==1) fprintf(g," %f",p_hat[site][srvy]);
				th_hat[site][srvy]=new double[3]; Bootth_hat[site][srvy]=new double[3];
            }
            for (srvy=0; srvy<T; srvy++) {
				th_hat[site][srvy][0]=getParam(srvy+1,0,srvy+1,site,srvy,Params,0,SingSpec.LinkFn); //get_th(Params,NPar,site,srvy);
				th_hat[site][srvy][1]=getParam(T+1+srvy,0,T+1+srvy,site,srvy,Params,0,SingSpec.LinkFn);
			}          
        }
        TestStat = OpenModLinkSD(Params, NPar);
        
        TestStat=gof(psi_hat, p_hat, th_hat, pi_hat, p_hat, p_hat, 2, poolcutoff, gofmtype); //  4th arg=print flag
        if (tmpstr) if (strstr(tmpstr,"hines")!=NULL && SingSpec.Verbose>0) printf("lofboot=%d  DataTestStat = %f\n",LoFBoot,TestStat);
        if (SingSpec.Verbose>0) printf("\nlofboot=%d  DataTestStat = %f\n",LoFBoot,TestStat);
        OrigData = new int*[N];   // store Original Data and Parameter Estimates
        for (site=0; site<N; site++) {
            OrigData[site] = new int[T];
            for (srvy=0; srvy<T; srvy++) OrigData[site][srvy] = SingSpec.Data[site][srvy];
        }
        OrigParams = new double[NPar]; 
        for (site=0; site<NPar; site++) OrigParams[site] = (OrigParams[site]>4.0?4.0:Params[site]);
        int maxfn=SingSpec.maxfn, p_value=0, nseasns, locc=1; AvTS = 0.0;  // add in bootstrap loop
        double LoBootTestStat=9.99999e44, HiBootTestStat=-LoBootTestStat;
        BTSsave=new double[LoFBoot];  int iBS=0; fprintf(g,"\n"); 
		nseasns=1; 
        for(BS=0; BS<LoFBoot; BS++) { time(&time1);       // generate Data, missing values are fixed
            if (SingSpec.Verbose>0) printf("                              %d/%d\r",BS+1,LoFBoot);
            for (site=0; site<N; site++) {
                if (ran1(&seed)<psi_hat[site]) {    // if site occupied
					locc = ran1(&seed)<pi_hat[site];
                    for (srvy=0; srvy<T; srvy++) {
						locc=ran1(&seed)<th_hat[site][srvy][locc]; // sim local site use
                        if (OrigData[site][srvy]>=0) SingSpec.Data[site][srvy] =(ran1(&seed) < (p_hat[site][srvy]*locc));
					}
                } else {    // if site vacant
                    for (srvy=0; srvy<T; srvy++) {
                        if (OrigData[site][srvy]>=0) SingSpec.Data[site][srvy] = 0; //  don't change missing data
                    }
                }
            } 
            // data generated, now fit model                                 set init values for parms
            for (srvy=0; srvy<NPar; srvy++) Params[srvy] = OrigParams[srvy];  // this line moved from inside if stmt below
            if (LoFBoot<11 && SingSpec.Verbose>1) {
                fprintf(g,"simulation number %d data:\n",BS+1);
                for (site=0; site<N; site++) {
                    fprintf(g,"%2d:",site+1);
                    for (srvy=0; srvy<T; srvy++) fprintf(g,"%c",(SingSpec.Data[site][srvy]>=0?'0'+SingSpec.Data[site][srvy]:'.'));
                    fprintf(g," ");
                    for (srvy=0; srvy<nseasns; srvy++) fprintf(g,"%c",(SingSpec.OpenData[site][srvy]>=0?'0'+SingSpec.OpenData[site][srvy]:'.'));
					fprintf(g,"\n");
                }
            }
            if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, OpenModLinkSD, 0.2);
            ier=optmiz(OpenModLinkSD,NPar, SingSpec.LikeNRSig, maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work);
            
            for (site=0; site<N; site++) {
                Bootpsi_hat[site]=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
				Bootpi_hat[site]=getParam(1+3*T,4,0,site,0,Params,0,SingSpec.LinkFn);  
                for (srvy=0; srvy<T; srvy++) {
					Bootp_hat[site][srvy]=getParam(srvy+frstp,3,srvy,site,srvy,Params,0,SingSpec.LinkFn);
					Bootth_hat[site][srvy][0]=getParam(srvy+1,0,srvy+1,site,srvy,Params,0,SingSpec.LinkFn);
					Bootth_hat[site][srvy][1]=getParam(T+srvy+1,0,T+srvy+1,site,srvy,Params,0,SingSpec.LinkFn);
				}
            }
			if (SingSpec.Verbose>1) {
				fprintf(g,"boot psi:%f ",Bootpsi_hat[0]);
				for (srvy=0; srvy<T; srvy++) fprintf(g,"%f ",Bootp_hat[0][srvy]);
				for (srvy=0; srvy<T; srvy++) fprintf(g,"%f ",Bootth_hat[0][srvy][0]);
				for (srvy=0; srvy<T; srvy++) fprintf(g,"%f ",Bootth_hat[0][srvy][1]);
				fprintf(g,"\n");
			}
            
            if (ier==131) {
                printf("iteration %d:Error - max function calls exceeded\n",BS+1);
                fprintf(g,"iteration %d:Error - max function calls exceeded\n",BS+1);
                if (SingSpec.Verbose>0) printf("iteration %d:Error - max function calls exceeded\n",BS+1);
            }            
            else {
                BootTestStat = gof(Bootpsi_hat,Bootp_hat,Bootth_hat,Bootpi_hat,Bootp_hat,Bootp_hat,1+(LoFBoot<5),poolcutoff,gofmtype); 
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
        for (site=0; site<NPar; site++) Params[site] = OrigParams[site];
        for (site=0; site<N; site++) for (srvy=0; srvy<T; srvy++) SingSpec.Data[site][srvy]=OrigData[site][srvy];
        
        // delete dynamic variables
        for (site=0; site<N; site++) {
			delete[] OrigData[site]; delete[] p_hat[site]; delete[] Bootp_hat[site];
			for (j=0; j<T; j++) {delete [] th_hat[site][j]; delete [] Bootth_hat[site][j]; }
			delete [] th_hat[site]; delete [] Bootth_hat[site];
		}
        delete[] psi_hat;  delete[] p_hat;  delete [] th_hat;  delete [] pi_hat;   
		delete[] Bootpsi_hat; delete[] Bootp_hat; delete [] Bootth_hat; delete [] Bootpi_hat;
        delete[] OrigData;  delete[] OrigParams;
		delete [] Hess; delete [] Grad1; delete [] Work;
}
