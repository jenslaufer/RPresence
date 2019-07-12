#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "SingAnal.h"
#include "rngs.h"


void Simulation(double PSI,double P,int NTot,int NInt,int TInt,int TOther,int Removal, int NSims,char *fname) {
    
	void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    double PDLike(double *Params, int NPar);
    double CMLike(double *Params, int NPar);
    // void PDGradient(double *Params, double *Grad, int NPar);
    double ran1(long *idum);   
    int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params, double *Hess, double *Grad, double *MaxLL, double *Work);
	void varcov(double xxLike(double *p, int N), double *Params, int NPar, double **covar);
    extern TSingSpec SingSpec; FILE *g=SingSpec.g;
    
    FILE *f;
    int ii, jj, kk, ll, n, NPar=2, ZeroCount=0, NotConv = 0, seen=0, temp, index, isw=0;
    double *psi, **p, sum, *Params, **covar, *Grad1, *Hess,*Work,MaxLL, AvPsi=0, SEPsi=0, AvSEPsi=0, AvNaive=0; 
    time_t t1; char timstr[40]; time(&t1); strcpy(timstr,ctime(&t1)); timstr[strlen(timstr)-1]='\0';
    
    printf("simulation(psi=%f, p=%f ntot=%d nint=%d tint=%d tother=%d removal=%d nsims=%d fname=%s)\n",
           PSI,P,NTot,NInt,TInt,TOther,Removal,NSims,fname);
    if (PSI<0) { isw=1; PSI=-PSI; }
    // get input values from file or screen
    SingSpec.Groups = 1; SingSpec.TSpecificP = 0;  SingSpec.Verbose=0;
    if (strlen(fname)>0) {
        printf("reading simulation input from %s\n",fname);
        if ((f=fopen(fname,"r"))==NULL) { printf("error opening simulation input file (%s)\n",fname); exit(1);}
        fscanf(f,"%d %d %d %d", &NTot , &NInt , &TInt , &TOther);
        printf("NTot=%d Nint=%d TInt=%d TOther=%d\n",NTot,NInt,TInt,TOther);
        psi = new double[NTot];  p = new double*[NTot];
        for (ii=0; ii<NTot; ii++) {
            fscanf(f,"%lf", &psi[ii]); p[ii] = new double[TInt]; printf("%d %f",ii,psi[ii]);
            for (jj=0; jj<TInt; jj++) { fscanf(f,"%lf", &p[ii][jj]); printf(" %f",p[ii][jj]);} printf("\n");
        }
        fscanf(f,"%d %d %d", &SingSpec.Groups, &Removal , &NSims);  fclose(f); 
        printf("Groups=%d Removal=%d NSims=%d\n",SingSpec.Groups,Removal,NSims);
    } 
    else {
        psi = new double[NTot];  p = new double*[NTot];
        for (ii=0; ii<NTot; ii++) {
            psi[ii] = PSI;  p[ii] = new double[TInt];
            for (jj=0; jj<TInt; jj++) p[ii][jj] = P;
        }
    }

    NPar = 2 + 2*(SingSpec.Groups-1); SingSpec.N = NTot; SingSpec.T = TInt; 
	
	SingSpec.expval = new double[NTot+1]; SingSpec.finalBetaEst = new double[NPar+1];  Params = new double[NPar+1];  
    
    covar = new double*[NPar]; Grad1 = new double[NPar]; 
    for (ii=0; ii<NPar; ii++) covar[ii] = new double[NPar]; 
    SingSpec.finalVCmat = new double*[NPar+1]; for (ii=0; ii<=NPar; ii++) SingSpec.finalVCmat[ii]=new double[NPar+1];	
	
	SingSpec.realParmEst = new double*[NPar+1];	SingSpec.finalBetaEst = new double[NPar]; SingSpec.fixed = new double[NPar];
	for (ii=0; ii<NPar; ii++) { SingSpec.realParmEst[ii] = new double[NTot]; SingSpec.fixed[ii]=-999; } SingSpec.realParmEst[NPar]=NULL;
	SingSpec.psibar=new double**[1]; SingSpec.psibar[0]=new double*[NTot];  for (ii=0; ii<NTot; ii++) SingSpec.psibar[0][ii]=new double[2];
	SingSpec.NrowsDM[0]=1; SingSpec.NrowsDM[1]=SingSpec.T; SingSpec.NParKK[0]=SingSpec.NParKK[1]=1;
	SingSpec.DMat[0]=new double*[1]; SingSpec.DMat[0][0]=new double[1]; SingSpec.DMat[0][0][0]=1; SingSpec.DMat[1]=new double*[SingSpec.T];
	SingSpec.DMat_ptr[0]=new int*[1]; SingSpec.DMat_ptr[0][0]=new int[1]; SingSpec.DMat_ptr[0][0][0]=0; SingSpec.DMat_ptr[1]=new int*[SingSpec.T];
	SingSpec.LnkFn[0]=new int[1]; SingSpec.LnkFn[0][0]=0;	SingSpec.LnkFn[1]=new int[SingSpec.T]; 
	SingSpec.BetaFixed[0]=new double[1]; SingSpec.BetaFixed[1]=new double[1]; SingSpec.BetaFixed[0][0]=SingSpec.BetaFixed[1][0]=1.1e44;
	
	for (ii=0; ii<SingSpec.T; ii++) { 
	    SingSpec.DMat[1][ii]=new double[1]; SingSpec.DMat[1][ii][0]=1; 
	    SingSpec.DMat_ptr[1][ii]=new int[1]; SingSpec.DMat_ptr[1][ii][0]=0; SingSpec.LnkFn[1][ii]=0;
	}
    
    Hess = new double[NPar*(NPar+1)/2+2];  Work = new double[3*NPar+2]; 

	SelectStream(0); PlantSeeds(-1);   //  initialize random number generator (in rngs.c)
	SingSpec.det_hist_frq=new double[NTot]; SingSpec.Data = new int*[NTot]; 
	for (ii=0; ii<NTot; ii++) { SingSpec.Data[ii] = new int[TInt];	SingSpec.det_hist_frq[ii]=1;}
    for (kk=0; kk<NSims; kk++) { 
        // generate Data
        for (ii=0,n=0; ii<NTot; ii++) {
            seen = 0; 
            if (Random() < psi[ii]) {    // species is present
                for (jj=0; jj<TInt; jj++) { 
                    SingSpec.Data[ii][jj] = 0; 
                    if (ii>=NInt && jj>=TOther) SingSpec.Data[ii][jj] = -1;
                    else { 
                        if (Random() < p[ii][jj]) { 
                            SingSpec.Data[ii][jj] = 1; 
                            if(jj==0) seen = 1; 
                            if (Removal>0) { 
                                for (ll=jj+1; ll<TInt; ll++) SingSpec.Data[ii][ll] = -1; 
                                break;
                            }
                        } 
                    }
                }
            } 
            else {    // species is absent
                for (jj=0; jj<TInt; jj++) {
                    if (ii>=NInt && jj>=TOther) SingSpec.Data[ii][jj] = -1; 
                    else SingSpec.Data[ii][jj] = 0;
                }
            }
            if (seen>0) n++;
        }   // end ii loop to generate data
        if (n==0) ZeroCount++;
        
        /*////////////////////////////////////////////////////////////////////////////
           Each row of Params corresponds to the estimated parameters.
           Parameters are in the order;
           psi
           G1, G2, ... G[Groups-1]
           p1, (p2, .., pT if required) for G1
           p1, (p2, .., pT if required) for G2
           .
           .
           .
         */////////////////////////////////////////////////////////////////////////////
        
        if (n>0) { 
			if (isw==0) DoAmoeba(Params, NPar, PDLike, 1.0); // intialise Params
			else DoAmoeba(Params, NPar, CMLike, 1.0); // intialise Params

            // set up covariance matrix of beta parameters
            if (isw==0) temp = optmiz(PDLike,NPar, 4, SingSpec.maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work);
            else temp = optmiz(CMLike,NPar, 4, SingSpec.maxfn, 0, Params-1, Hess, Grad1-1, &MaxLL, Work); 
			ii=temp;
            if (Work[3]<4) NotConv++;   // sufficient convergence not reached
            else { 
                if (isw==0) varcov(PDLike,Params,NPar,covar);
                else varcov(CMLike,Params,NPar,covar);
                // find gradient of real parameters and store in row 0 of Grad1
                if (isw==0) Grad1[0] = cos(Params[0])/2.0; // grad. psi
				else Grad1[0]=exp(-Params[0])/pow(1+exp(-Params[0]),2);
                index = 1; sum=1.0;
                if (SingSpec.Groups>1) {                     // grad. Group probs
                    for (ii=0; ii<SingSpec.Groups-1; ii++) sum += exp(Params[ii+1]);
                    for (ii=1; ii<SingSpec.Groups; ii++, index++) Grad1[ii] = exp(Params[ii])*(sum - exp(Params[ii]))/(sum*sum);
                }
                // everything else in Params uses the sin link
                for (ii=index; ii<NPar; ii++) {
					if (isw==0) Grad1[ii] = cos(Params[ii])/2.0;
				    else Grad1[ii]=exp(-Params[ii])/pow(1+exp(-Params[ii]),2);
				}

                // convert beta VC to real VC
                for (ii=0; ii<NPar; ii++)
                    for (jj=0; jj<NPar; jj++) covar[ii][jj] *= Grad1[ii]*Grad1[jj];

                if (isw==0) AvPsi += (sin(Params[0])+1)/2.0; else AvPsi+=1/(1+exp(-Params[0]));
                if (isw==0) SEPsi += pow((sin(Params[0])+1)/2.0,2); else SEPsi += pow((1/(1+exp(-Params[0]))),2);
                AvSEPsi += sqrt(covar[0][0]);
                AvNaive += ((double)n)/NTot;
            } // end if converged
        } // end if (n>0)
        printf("simulation %d               \r",kk+1);
    }   // end simulation loop
    printf("\n");

    AvPsi /= (NSims-ZeroCount-NotConv);
    AvSEPsi /= (NSims-ZeroCount-NotConv);
    SEPsi -= (NSims-ZeroCount-NotConv)*AvPsi*AvPsi;
    SEPsi /= (NSims-1.0-ZeroCount-NotConv);
    SEPsi = sqrt(SEPsi);
    AvNaive /= (NSims-ZeroCount-NotConv);

    fprintf(g,"\nSimulation Results\n\n");
    fprintf(g,"Total number of sites sampled                :%d \n",NTot);
    fprintf(g,"Number of sites sampled more intensively     :%d \n",NInt);
    fprintf(g,"Number of visits to intensively sampled sites:%d \n",TInt);
    fprintf(g,"Number of visits to other sites              :%d \n",TOther);
    if (Removal) {
        fprintf(g,"Sites NOT visited after first detection.\n");
    }
    fprintf(g,"\nNumber of simulations                               :%d \n",NSims);
    fprintf(g,  "Number of times species was not detected at any site:%d \n",ZeroCount);
    fprintf(g,  "Number of times convergence not achieved            :%d \n",NotConv);
    if (strlen(fname)>0) {
        fprintf(g,"\nRunning simulations from file   :%s \n",fname);
        fprintf(g,  "Number of groups in fitted model:%d \n",SingSpec.Groups);
    } 
    else {
        fprintf(g,"\nDetection probability            :%f \n",p[0][0]);
        fprintf(g,  "True proportion of sites occupied:%f \n",psi[0]);
    }
    fprintf(g,"\nAverage naive estimate from single visit   :%f \n",AvNaive);
    fprintf(g,  "Average estimate of occupancy probability  :%f \n",AvPsi);
    fprintf(g,  "Simulation based estimate of standard error:%f \n",SEPsi);
    fprintf(g,  "Average estimate of the standard error     :%f \n",AvSEPsi);
	
    fclose(g);
	exit(0);
}
