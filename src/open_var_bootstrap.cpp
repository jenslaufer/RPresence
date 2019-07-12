#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
//#include <search.h>

void Open_bootstrap(int NBoot, int NPar, double *Params, double **OrigPsi, double **OrigGam, double **OrigEps, double **OrigP){
    // declare functions used here
    double OpenModLink(double *Params, int NPar);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
    double ran1(long *idum);
	void try_optmz_randiv(double (*LinkFn)(double *Params, int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);	
    
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    double xran,xnboot=0,avg,se,MaxLL,psi,gam,eps,cp;
    int site,BS,ii,jj,k,kj,occ,sitelmt,N=SingSpec.N,T=SingSpec.T;
    long seed=SingSpec.seed; time_t time1,time2;
    
    //printf("\nNPar=%d nboot=%d prmypers=%d\n",NPar,NBoot,SingSpec.PrmyPeriods);
    if (NBoot>0) {
        fprintf(g,"\n============================================================\n");
        fprintf(g,"\nBootstrap estimation of parameter estimates/std. errors(seed=%ld):\n\n",seed);
		double **OrigPsi,**OrigGam,**OrigEps,**OrigP;
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
        double *OrigParams = new double[NPar+2]; printf("origparms=%p\n",OrigParams);
        for (ii=0; ii<NPar; ii++) OrigParams[ii]=Params[ii];
        int **OrigData = new int*[N];   // store Original Data and Parameter Estimates
        double **Boot_avgpsi = new double*[N];  double **Boot_sdpsi = new double*[N];
        double **Boot_avggam = new double*[N];  double **Boot_sdgam = new double*[N];
        double **Boot_avgeps = new double*[N];  double **Boot_sdeps = new double*[N];
        double **Boot_avgp = new double*[N];    double **Boot_sdp = new double*[N];
        for (ii=0; ii<N; ii++) { 
            Boot_avgpsi[ii] = new double[SingSpec.PrmyPeriods];
            Boot_avggam[ii] = new double[SingSpec.PrmyPeriods];
            Boot_avgeps[ii] = new double[SingSpec.PrmyPeriods];
            Boot_sdpsi[ii] = new double[SingSpec.PrmyPeriods];
            Boot_sdgam[ii] = new double[SingSpec.PrmyPeriods];
            Boot_sdeps[ii] = new double[SingSpec.PrmyPeriods];
            Boot_avgp[ii] = new double[T]; Boot_sdp[ii] = new double[T];
            OrigData[ii] = new int[T+1];
            for (k=0; k<SingSpec.PrmyPeriods; k++)  
                Boot_avgpsi[ii][k]=Boot_avggam[ii][k]=Boot_avgeps[ii][k]=
                  Boot_sdpsi[ii][k]=Boot_sdgam[ii][k]=Boot_sdeps[ii][k]=0;
            for (k=0; k<T; k++) {
                Boot_avgp[ii][k]=Boot_sdp[ii][k]=0;
                OrigData[ii][k] = SingSpec.Data[ii][k];
            }
        }
        double *Hess = new double[NPar*(NPar+1)/2+2]; double *Work = new double[3*NPar+2]; 
        
        fprintf(g,"psi=%f\ngam=",OrigPsi[0][0]);
        for (k=1; k<SingSpec.PrmyPeriods; k++) fprintf(g," %f",OrigGam[0][k-1]); fprintf(g,"\neps=");
        for (k=1; k<SingSpec.PrmyPeriods; k++) fprintf(g," %f",OrigEps[0][k-1]); fprintf(g,"\np  =");
        for (k=0; k<T; k++) fprintf(g," %f",OrigP[0][k]); fprintf(g,"\n\n");
        
        time(&time1); sitelmt=N; if (SingSpec.lmt==1) sitelmt=1;
        for(BS=0; BS<NBoot; BS++) {        // generate Data, missing values are fixed
            printf("                          bootstrap: %d/%d\r",BS+1,NBoot);
            for (ii=0; ii<N; ii++) { occ=0;
                if (ran1(&seed)<=OrigPsi[ii][0]) occ=1;   //  phi[ii][0][0][0] = initial psi
                for (k=kj=0; k<SingSpec.PrmyPeriods; k++) {
                    for (jj=0; jj<SingSpec.SecPeriods[k]; kj++,jj++) {
                        if (OrigData[ii][kj]>=0) {            
                            if (ran1(&seed) <= occ*OrigP[ii][kj]) SingSpec.Data[ii][kj] = 1;
                            else                                  SingSpec.Data[ii][kj] = 0;
                        }
                    }
                    xran=ran1(&seed); //OrigEps is actually 1-eps
                    occ=occ*(xran<=(1-OrigEps[ii][k]))+(1-occ)*(xran<=OrigGam[ii][k]);  //occ=occ*(1-eps)+(1-occ)*gamma
                }
                if(SingSpec.Verbose) {
                    fprintf(g,"%d>>",BS+1); 
                    for (k=0; k<T; k++) fprintf(g,"%d",SingSpec.Data[ii][k]); 
                    fprintf(g,"\n");
                }
            }
            // data generated, now fit model                              set init values for parms
            for (jj=0; jj<NPar; jj++) {
                xran=OrigParams[jj]; if (xran<-1.5) xran=-1.5;
                if(xran>1.5) xran=1.5;
                Params[jj] = xran; //OrigParams[jj]; 
            }
            
            if(SingSpec.UseAmoeba) DoAmoeba(Params, NPar, OpenModLink, 0.2);

            //printf("call optmiz...\n");
	        try_optmz_randiv(OpenModLink, NPar, Params, Work, Hess, &MaxLL);
            for (k=0; k<NPar; k++) fprintf(g," %f",Params[k]);
            fprintf(g," %f",OrigPsi[0][0]);
            for (k=1; k<SingSpec.PrmyPeriods; k++) fprintf(g," %f",OrigGam[0][k]);
            for (k=1; k<SingSpec.PrmyPeriods; k++) fprintf(g," %f",OrigEps[0][k]);
            fprintf(g," %f\n",OrigP[0][0]);
            if (Work[3]>2.) {
                xnboot++;
                for (ii=0; ii<sitelmt; ii++) { kj=0; psi=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
                    Boot_avgpsi[ii][0]+=psi;
                    Boot_sdpsi[ii][0]+=psi*psi;
                    for (k=0; k<SingSpec.PrmyPeriods; k++) {
                        if (k<(SingSpec.PrmyPeriods-1)) {
                            gam=getParam(1+k,1,k,site,k,Params,0,SingSpec.LinkFn);
							eps=getParam(k+SingSpec.PrmyPeriods,2,k,site,k,Params,0,SingSpec.LinkFn);
                            Boot_avggam[ii][k]+=gam;
                            Boot_sdgam[ii][k] +=gam*gam;
                            Boot_avgeps[ii][k]+=eps;
                            Boot_sdeps[ii][k] +=eps*eps;
                            psi=psi*(1-eps)+(1-psi)*gam;
                            Boot_avgpsi[ii][k+1]+=psi;
                            Boot_sdpsi[ii][k+1] +=psi*psi;
                        }
                        for (jj=0; jj<SingSpec.SecPeriods[k]; jj++,kj++) {
                            if (kj>=T) printf("\n******** subscript error *****\n");
							cp=getParam(2*SingSpec.PrmyPeriods-1+k,3,kj,site,kj,Params,0,SingSpec.LinkFn);
                            Boot_avgp[ii][kj]+=cp;
                            Boot_sdp[ii][kj]+=cp*cp;
                        }
                    }
                }
            }
        }   // end bootstrap loop
        time(&time2); 
        
        // Output results
        fprintf(g,"%d good reps out of %d bootstrap simulations\n\n",(int)xnboot,NBoot);
        fprintf(g,"Bootstrap estimates:\n");
        fprintf(g,"                    data         bootstrap\n");
        fprintf(g,"parm       site   estimate    estimate  std.err\n");
        fprintf(g,"----     ------     ----      -----    --------\n");
        for (k=0; k<SingSpec.PrmyPeriods; k++) {
            for (ii=0; ii<sitelmt; ii++) {
                avg=Boot_avgpsi[ii][k]; se=(Boot_sdpsi[ii][k]-avg*avg/xnboot)/(xnboot-1); 
                avg=avg/xnboot; se=sqrt(se);
                fprintf(g,"psi(%d) %8d  %9.4f %9.4f %9.4f\n",k+1,ii+1,OrigPsi[ii][k],avg,se);
            }
        }
        for (k=0; k<(SingSpec.PrmyPeriods-1); k++) {
            for (ii=0; ii<sitelmt; ii++) {
                avg=Boot_avggam[ii][k]; se=(Boot_sdgam[ii][k]-avg*avg/xnboot)/(xnboot-1); 
                avg=avg/xnboot; se=sqrt(se);
                fprintf(g,"gam(%d) %8d  %9.4f %9.4f %9.4f\n",k+1,ii+1,OrigGam[ii][k],avg,se);
            }
        }
        for (k=0; k<(SingSpec.PrmyPeriods-1); k++) {
            for (ii=0; ii<sitelmt; ii++) {
                avg=Boot_avgeps[ii][k]; se=(Boot_sdeps[ii][k]-avg*avg/xnboot)/(xnboot-1); 
                avg=avg/xnboot; se=sqrt(se);
                fprintf(g,"eps(%d) %8d  %9.4f %9.4f %9.4f\n",k+1,ii+1,OrigEps[ii][k],avg,se);
            }
        }
        for (k=0; k<T; k++) {
            for (ii=0; ii<sitelmt; ii++) {
                avg=Boot_avgp[ii][k]; se=(Boot_sdp[ii][k]-avg*avg/xnboot)/(xnboot-1); 
                avg=avg/xnboot; se=sqrt(se);
                fprintf(g,"p(%d)   %8d  %9.4f %9.4f %9.4f\n",k,ii+1,OrigP[ii][k],avg,se);
            }
        }
        fprintf(g,"------------------------------------------------------\n");
        printf("\ncpu time for bootstrap:%ld\n",(long)(time2-time1));
        fprintf(g,"cpu time for bootstrap:%ld\n",(long)(time2-time1));
        fprintf(g,"------------------------------------------------------\n");
        
        // copy orginal data and paramter estimates back into correct arrays;
        for (ii=0; ii<NPar; ii++) Params[ii] = OrigParams[ii];
        for (ii=0; ii<N; ii++) for (jj=0; jj<T; jj++) SingSpec.Data[ii][jj]=OrigData[ii][jj];
        // delete dynamic variables
		for (site=0; site<N; site++) {
			delete [] OrigPsi[site]; delete [] OrigGam[site]; delete [] OrigEps[site]; delete [] OrigP[site];
		}
		delete [] OrigPsi; delete [] OrigGam; delete [] OrigEps; delete [] OrigP;	
		
        for (ii=0; ii<N; ii++) {  
            delete OrigData[ii]; delete [] OrigPsi[ii]; 
            delete[] Boot_avgpsi[ii]; delete[] Boot_sdpsi[ii];
            delete[] Boot_avggam[ii]; delete[] Boot_sdgam[ii];
            delete[] Boot_avgeps[ii]; delete[] Boot_sdeps[ii];
            delete[] Boot_avgp[ii]; delete[] Boot_sdp[ii];
        }
        delete[] OrigData; delete[] OrigParams;
        delete[] Boot_avgpsi; delete[] Boot_avggam; delete[] Boot_avgeps; delete[] Boot_avgp;
        delete[] Boot_sdpsi; delete[] Boot_sdgam; delete[] Boot_sdeps; delete[] Boot_sdp;
        delete[] Hess; delete[] Work;
    }/////  end of if(NBoot) stmt /////////
}
