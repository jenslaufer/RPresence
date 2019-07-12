#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "SingAnal.h"
#include "rngs.h"
#define VARFLG 1

void PreDefMod(int NBoot, int LoFBoot) { 
    // Data[][] contains detection/nondetection info
    // Missing[][] contains whether observations were missing
    // N is number of sites
    // T is number of sampling occasions
    // Groups is number of groups for mixture models
    // TSpecificP indicates whether detection probabilities are time specific
    
    // declare PDamoeba and PDLike
    
    int dblcmp(const void *v1, const void *v2); 
    double PDLike(double *Params, int NPar);    // void PDGradient(double *Params, double *Grad, int NPar);
    double ran1(long *idum);	
    int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params, double *Hess, double *Grad, double *MaxLL, double *Work);
    void varcov(double xxLike(double *p, int N), double *Params, int NPar, double **covar);
    double gof(double *psi_hat, double **p_hat, double ***th_hat, double *pi_hat, double **p10, double **b, int prnt, double poolcutoff, int mtype);
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);	
    
    extern TSingSpec SingSpec; FILE *g;    
    int i,j,ii, jj, kk, bb, k, itemp, dbg=SingSpec.Verbose, gofmtype=0;     // delcare some indx counters
    int N=SingSpec.N, T=SingSpec.T;
    double sum,*Params, y,TestStat,poolcutoff=.1e-6,*Hess,*Work,MaxLL;
    double *psi_hat,*pi_hat,**x_hat,***p_hat,**covar, **Grad1, **Grad2;
    
    int NPar = 1 /*psi*/ + (SingSpec.Groups-1) + SingSpec.Groups*(1 + (SingSpec.TSpecificP*(SingSpec.T-1)));
    SingSpec.NrowsDM[0]=NPar; SingSpec.NParKK[0]=NPar;
    g=SingSpec.g; SingSpec.NMiss = 0; if (SingSpec.MissClass) NPar--;
    printf("N=%d ,T=%d, NPar=%d, TSpecificP=%d\n",SingSpec.N,SingSpec.T,NPar,SingSpec.TSpecificP);
    Params=new double[NPar+1]; for (i=0; i<=NPar; i++) Params[i]=0;
    psi_hat=new double[N]; pi_hat=new double[N]; p_hat=new double**[N]; x_hat=new double*[N]; 
    for (i=0; i<N; i++) { 
		p_hat[i]=new double*[T]; x_hat[i]=new double[T]; 
		for (j=0; j<T; j++) p_hat[i][j]=new double[2];
	}
    /*////////////////////////////////////////////////////////////////////////////
       Each row of Params corresponds to the estimated parameters.
       Parameters are in the order;
       psi
       G1, G2, ... G[Groups-1]
       p1, (p2, .., pT if required) for G1
       p1, (p2, .., pT if required) for G2
       .
       :
     */////////////////////////////////////////////////////////////////////////////
    
    Params[0]=.5; Params[1]=.5;
    if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, PDLike, 1.0);

    // set up covariance matrix of beta parameters
    covar = new double*[NPar+1]; Grad1 = new double*[NPar+1];  Grad2 = new double*[NPar+1];
    for (ii=0; ii<=NPar; ii++) {
        covar[ii] = new double[NPar+1];  Grad1[ii] = new double[NPar+1]; Grad2[ii] = new double[NPar+1];
        for (jj=0; jj<=NPar; jj++) Grad1[ii][jj]=0;
    }
    Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2];

    printf("calling optmiz...\n");
    optmiz(PDLike,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
    
	fprintf(g,"**************\n"); y = MaxLL;
    // Output results
    if (SingSpec.MissClass)
        fprintf(g,"\n\nMis-classification Model: Detection probabilities are ");
    else
        fprintf(g,"\n\nPredefined Model: Detection probabilities are ");
    if (!SingSpec.TSpecificP) fprintf(g,"NOT ");
    fprintf(g,"time-specific\n\n");
    fprintf(g,"Number of groups               = %d\n", SingSpec.Groups);
    fprintf(g,"Number of parameters           = %d\n", NPar);
    if (Work[3]<SingSpec.LikeNRSig) {
        fprintf(g,"**** Numerical convergence may not have been reached.\n");
        fprintf(g,"**** Parameter estimates converged to approximately %f\n",Work[3]);
        fprintf(g,"**** significant digits.\n");
    }
    fprintf(g,"Number of function calls           = %1.0f\n",Work[2]);
    fprintf(g,"-2log(likelihood)              = %f\n",2.0*y);
    fprintf(g,"AIC                            = %f\n",2.0*(y+NPar));

    if (dbg>0) {
        for (ii=0,kk=0; ii<NPar; ii++)
            for (jj=0; jj<=ii; jj++) {covar[ii][jj]=Hess[kk++]; covar[jj][ii]=covar[ii][jj];}
        fprintf(g,"Beta and Hessian matrix:\n");
        for (ii=0; ii<NPar; ii++) {
            fprintf(g,"%8.4f ",Params[ii]);
            for (jj=0; jj<NPar; jj++) fprintf(g," %12.4f",covar[ii][jj]);
            fprintf(g,"\n");
        } 
        fprintf(g,"\n");
    }

    varcov(PDLike,Params,NPar,covar);
    
    if (dbg>1) {
        fprintf(g,"Beta and Beta Var-cov matrix:\n");
        for (ii=0; ii<NPar; ii++) {
            fprintf(g,"%f ",Params[ii]);
            for (jj=0; jj<NPar; jj++) fprintf(g," %f",covar[ii][jj]);
            fprintf(g,"\n");
        } fprintf(g,"\n");
    }
    // find gradient of real parameters and store in Grad1
    Grad1[0][0] = cos(Params[0])/2.0; // grad. psi
    int indx = 1;    sum=1.0;
    if (SingSpec.Groups>2) {                     // grad. Group probs
        for (ii=0; ii<SingSpec.Groups-1; ii++) sum += exp(Params[ii+1]);
        for (ii=1; ii<SingSpec.Groups; ii++, indx++) {
            Grad1[ii][ii] = exp(Params[ii])*(sum - exp(Params[ii]))/(sum*sum);
            for (jj=ii+1; jj<SingSpec.Groups; jj++) {
                Grad1[jj][ii] = Grad1[ii][jj] = -exp(Params[ii])*exp(Params[jj])/(sum*sum);
            }
        }
    }
    if (SingSpec.Groups==2) { sum=1/(sin(Params[1])/2.+.5); Grad1[1][1]=cos(Params[1])/2.; indx++;}
    // everything else in Params uses the sin link
    for (ii=indx; ii<NPar; ii++) { Grad1[ii][ii] = cos(Params[ii])/2.0;}
    
    // convert beta VC to real VC
    // matrix multiply V*G, store result in Grad2
    for (ii=0; ii<NPar; ii++) {
        for (jj=0; jj<NPar; jj++) {
            Grad2[ii][jj] = 0.0;
            for (kk=0; kk<NPar; kk++) Grad2[ii][jj] += covar[ii][kk]*Grad1[kk][jj];
        }
    }
    // matrix multiply transpose(G)*Grad2, store result in V
    for (ii=0; ii<NPar; ii++) {
        for (jj=0; jj<NPar; jj++) {
            covar[ii][jj] = 0.0;
            for (kk=0; kk<NPar; kk++) covar[ii][jj] += Grad1[kk][ii]*Grad2[kk][jj];
        }
    }
    k=0;
	psi_hat[0]=(sin(Params[0])+1.)/2.; pi_hat[0]=1./sum;
	if (!SingSpec.MissClass) {
        fprintf(g,"\nProportion of sites occupied    (Psi)   = %6.4f", psi_hat[0]);
        if (covar[0][0] < 0.0 && covar[0][0]>-.00000000001) covar[0][0]=0;
        if (covar[0][0] < 0.0) fprintf(g," (XXXX)\n"); else fprintf(g," (%6.4f)\n",sqrt(covar[0][0]));
        if (SingSpec.Groups>1) 
            fprintf(g,"Probability of group membership (Theta) = %6.4f (%6.4f)",1.0/sum,sqrt(covar[1][1]));
    }
    else  {
        fprintf(g,"Probability of group classification (Theta) = %6.4f",1.0/sum); 
        k=-1;
    }
    // sum calculated above in VC section
    if (SingSpec.Groups!=2) 
        for (ii=0; ii<SingSpec.Groups-1; ii++) fprintf(g,", %6.4f",exp(Params[++k])/sum);
    else {fprintf(g,", %6.4f",1.-1./sum); k++;}
    fprintf(g,"\nDetection probabilities         (p):\n");
    fprintf(g,"  grp   srvy      p            se(p)\n");
    fprintf(g,"  ---   ----   ---------    -----------\n");
    for (ii=0; ii<SingSpec.Groups; ii++) { k++;;
        for (jj=0; jj<SingSpec.T; jj++) {
            if (SingSpec.TSpecificP && jj>0) k++;
            //if (SingSpec.TSpecificP || jj==0) {
                p_hat[0][jj][ii]=(sin(Params[k])+1.)/2.;
                fprintf(g,"%4d %6d   %9.6f ",ii+1,jj+1,p_hat[0][jj][ii]);
                if (covar[k][k]<0.) fprintf(g,"    (XXXXXXXXX)"); 
                else fprintf(g,"    (%9.6f)\n",sqrt(covar[k][k]));
            //}
        }
        fprintf(g,"\n");
    }
    fprintf(g,"Variance-Covariance Matrix\n");    // Ouput the VC matrix
    if (!SingSpec.MissClass) fprintf(g,"   psi  ");
    for (ii=0; ii<SingSpec.Groups-1; ii++) fprintf(g,"Group %d ",ii+2);
    if (SingSpec.TSpecificP) {
        for (ii=0; ii<SingSpec.Groups; ii++) 
            for (jj=0; jj<SingSpec.T; jj++) fprintf(g,"  p%d(G%d)",jj+1,ii+1);
    } 
    else 
        for (ii=0; ii<SingSpec.Groups; ii++) fprintf(g,"  p(G%d)",ii+1);
    fprintf(g,"\n");
    for (ii=0; ii<NPar; ii++) {
        for (jj=0; jj<NPar; jj++)  fprintf(g,"%7.4f ",covar[ii][jj]); fprintf(g,"\n");
    }
    for (ii=1; ii<SingSpec.N; ii++) {
        psi_hat[ii]=psi_hat[0]; pi_hat[ii]=pi_hat[0];
        for (jj=0; jj<T; jj++) {
            p_hat[ii][jj][0]=p_hat[0][jj][0]; p_hat[ii][jj][1]=p_hat[0][jj][1];
        }
    }   
    TestStat=gof(psi_hat, x_hat, p_hat, pi_hat, x_hat, x_hat, 22, poolcutoff, gofmtype); //  4th arg=print flag
    // Assess model fit here, but output later
	char *tmpstr; tmpstr=(char*)malloc(256);
	double TestStat2=0,*Bootpsi_hat,***Bootp_hat,*Bootpi_hat,BootTestStat=0,*OrigParams,AvTS,*BTSsave;
	int site,srvy,BS,**OrigData; 
	time_t time1,time2,site_occ,site_grp; long seed=SingSpec.seed; 
	Bootpsi_hat=new double[N]; Bootpi_hat=new double[N]; Bootp_hat=new double**[N];
	for (i=0; i<N; i++) {
	  Bootp_hat[i]=new double*[T];
	  for (j=0; j<T; j++) {
		  Bootp_hat[i][j]=new double[2];
	  }
	}
    if (LoFBoot>0) { 
        tmpstr = getenv("USERNAME");
        fprintf(g,"\n============================================================\n");
        fprintf(g,"\nAssessing Model Fit for Single-season predefined model:\n");

        if (strstr(tmpstr,"hines")!=NULL && SingSpec.Verbose>0) printf("lofboot=%d  DataTestStat = %f %f\n",LoFBoot,TestStat,TestStat2);
        if (SingSpec.Verbose>0) printf("\nlofboot=%d  DataTestStat = %f\n",LoFBoot,TestStat);
        OrigData = new int*[N];  // store Original Data and Parameter Estimates
        for (site=0; site<N; site++) {
            OrigData[site] = new int[T]; 
            for (srvy=0; srvy<T; srvy++) OrigData[site][srvy] = SingSpec.Data[site][srvy];
        }
        OrigParams = new double[NPar]; 
        for (site=0; site<NPar; site++) OrigParams[site] = (OrigParams[site]>4.0?4.0:Params[site]);
        int maxfn=SingSpec.maxfn, p_value=0, nseasns, ier; 
        double LoBootTestStat=9.99999e44, HiBootTestStat=-LoBootTestStat;
        BTSsave=new double[LoFBoot];  int iBS=0; AvTS = 0.0; nseasns=1; 
		fprintf(g,"\n"); 	
        for(BS=0; BS<LoFBoot; BS++) { time(&time1);       // generate Data, missing values are fixed
            if (SingSpec.Verbose>0) printf("                              %d/%d\r",BS+1,LoFBoot);
            for (site=0; site<N; site++) {
                site_occ=(ran1(&seed)<psi_hat[site]); 
				site_grp=(ran1(&seed)<pi_hat[site]);
                for (srvy=0; srvy<T; srvy++) 
                    if (OrigData[site][srvy]>=0) 
						SingSpec.Data[site][srvy]=site_occ*(ran1(&seed)<(p_hat[site][srvy][site_grp]));
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
			if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, PDLike, 0.2);

            ier=optmiz(PDLike,NPar, SingSpec.LikeNRSig, maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
            
            for (site=0; site<N; site++) {
                Bootpsi_hat[site]=(sin(Params[0])+1.)/2.;
				Bootpi_hat[site]=(sin(Params[1])+1.)/2.;
                for (srvy=0; srvy<T; srvy++) {
					Bootp_hat[site][srvy][0]=(sin(Params[srvy+2])+1.)/2.;
					Bootp_hat[site][srvy][1]=(sin(Params[T+srvy+2])+1.)/2.;
				}
            }
			if (SingSpec.Verbose>1) {
				fprintf(g,"boot psi:%f theta:%f ",Bootpsi_hat[0],Bootpi_hat[0]);
				for (srvy=0; srvy<T; srvy++) fprintf(g,"%f ",Bootp_hat[0][srvy][0]);
				fprintf(g,"\n");
			}
            
            if (ier==131) {
                printf("iteration %d:Error - max function calls exceeded\n",BS+1);
                fprintf(g,"iteration %d:Error - max function calls exceeded\n",BS+1);
                if (SingSpec.Verbose>0) printf("iteration %d:Error - max function calls exceeded\n",BS+1);
            }            
            else {
				BootTestStat = gof(Bootpsi_hat,x_hat,Bootp_hat,Bootpi_hat, x_hat, x_hat, -2,poolcutoff, gofmtype); 
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
        for (site=0; site<N; site++) {delete[] OrigData[site]; delete[] p_hat[site]; delete[] Bootp_hat[site];}
        delete[] OrigData;   delete[] p_hat;       delete[] Bootp_hat;
        delete[] psi_hat;    delete[] Bootpsi_hat; delete[] OrigParams;
        
    }/////  end of if(Lofboot) stmt ///////////////////////////////////////////////////////////

	
    // Use bootstrapping to estimate V-C matrix
    if (NBoot > 0) {  printf("bootstraping...\n");
		SelectStream(0); PlantSeeds(-1);   //  initialize random number generator (in rngs.c) 
        double **BootParams, *ObsParams, *AvParams; int **OrigData; 
        double **BootCovar;    BootCovar = new double*[NPar]; 
        for (ii=0; ii<NPar; ii++) { BootCovar[ii] = new double[NPar];  } 
        OrigData = new int*[SingSpec.N]; 
        for (ii=0; ii<SingSpec.N; ii++) { 
            OrigData[ii] = new int[SingSpec.T]; 
            for (jj=0; jj<SingSpec.T; jj++) OrigData[ii][jj] = SingSpec.Data[ii][jj];
        }
        BootParams = new double*[NBoot]; ObsParams = new double[NPar]; AvParams = new double[NPar]; 
        for (ii=0; ii<NPar; ii++) { ObsParams[ii] = Params[ii];  AvParams[ii] = 0.0; } 
        
        for (bb=0; bb<NBoot; bb++) { printf(".......... bootstrap iteration %d\r",bb+1);
            BootParams[bb] = new double[NPar]; 
            
            // Bootstrap data
            for (ii=0; ii<SingSpec.N; ii++) {
                itemp = (int)(SingSpec.N*Random());
                for (jj=0; jj<SingSpec.T; jj++) SingSpec.Data[ii][jj] = OrigData[itemp][jj];
            }
             
            // setup and minimise
            for (jj=0; jj<NPar; jj++) {  Params[jj] = ObsParams[jj];  }
			if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, PDLike, 0.2);
            
            //printf("optmiz(%d %d ..)\n",NPar,SingSpec.BootNRSig);
            optmiz(PDLike,NPar, SingSpec.BootNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1[0]-1, &MaxLL, Work);
            //printf("%f %f %f %f func evals\n========== ",Work[0],Work[1],Work[2],Work[3]);
            //for (ii=0; ii<NPar; ii++) printf(" %f",Params[ii]); printf("\n");
            // store estimates
            BootParams[bb][0] = (sin(Params[0])+1.0)/2.0; 
            indx = 1;     sum=1.0; 
            if (SingSpec.Groups>2) {                     // grad. Group probs
                for (ii=0; ii<SingSpec.Groups-1; ii++) sum += exp(Params[ii+1]);
                for (ii=1; ii<SingSpec.Groups-1; ii++, indx++) BootParams[bb][ii] = exp(Params[ii])/sum;
            }
            if (SingSpec.Groups==2) {BootParams[bb][1]=sin(Params[1])/2+.5; indx++;}
            // everything else in Params uses the sin link
            for (ii=indx; ii<NPar; ii++) BootParams[bb][ii] = (sin(Params[ii])+1.0)/2.0;
            for (ii=0; ii<NPar; ii++) AvParams[ii] += BootParams[bb][ii];  
        }  // end bootstrap loop
        
        for	(ii=0; ii<NPar; ii++) { 
            AvParams[ii] /= (double)NBoot; 
        } 
        for (ii=0; ii<NPar; ii++) { 
            for (jj=ii; jj<NPar; jj++) {
                BootCovar[ii][jj] = 0.0;
                for (int bb=0; bb<NBoot; bb++) 
                    BootCovar[ii][jj] += (BootParams[bb][ii]-AvParams[ii])*(BootParams[bb][jj]-AvParams[jj]);
                BootCovar[ii][jj] /= (double)(NBoot-1);
                BootCovar[jj][ii] = BootCovar[ii][jj];
            }
        }
        
        // output VC matrix
        fprintf(g,"Bootstrap SE & Variance-Covariance Matrix\n");
        if (!SingSpec.MissClass) fprintf(g,"   SE          psi "); 
        for (ii=0; ii<SingSpec.Groups-1; ii++) { fprintf(g,"Group %d ",ii+2); }
        if (SingSpec.TSpecificP) 
            for (ii=0; ii<SingSpec.Groups; ii++) 
                for (jj=0; jj<SingSpec.T; jj++) fprintf(g,"  p%d(G%d)",jj+1,ii+1);
        else 
            for (ii=0; ii<SingSpec.Groups; ii++) fprintf(g,"  p(G%d)",ii+1); 
        fprintf(g,"\n");
        
        for (ii=0; ii<NPar; ii++) { fprintf(g,"%9.6f ",sqrt(BootCovar[ii][ii]));
            for (jj=0; jj<NPar; jj++) fprintf(g,"%9.6f ",BootCovar[ii][jj]);  fprintf(g,"\n");
        }
        fprintf(g,"\n");
        fprintf(g,"Bootstrap estimate of SE for\n proportion of sites occupied = %f\n",sqrt(BootCovar[0][0]));
        
        // copy observed data back into Data and Missing
        for (ii=0; ii<SingSpec.N; ii++) for (jj=0; jj<SingSpec.T; jj++) SingSpec.Data[ii][jj] = OrigData[ii][jj];
        
        // tidy up dynamic variables
        for (ii=0; ii<SingSpec.N; ii++) { delete[] OrigData[ii]; }
        delete[] OrigData;  delete[] ObsParams;  delete[] AvParams;
        for (ii=0; ii<NBoot; ii++) { delete[] BootParams[ii]; } delete[] BootParams;
        for (ii=0; ii<NPar; ii++) { delete[] BootCovar[ii]; } delete[] BootCovar;
        
    } // end if(NBoot>0)
    
    fprintf(g,"------------------------------------------------------");
    // delete dynamic variables
    delete[] Params; 
    for (ii=0; ii<NPar; ii++) { delete[] covar[ii]; delete[] Grad1[ii]; delete[] Grad2[ii];  }
    delete[] covar; delete[] Grad1; delete[] Grad2;
    delete[] Hess;   delete[] Work;
}

