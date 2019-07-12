#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "SingAnal.h"

double xxLike(double *Params, int NPar);

#include <pthread.h>
struct vlltype { 
	double *p;                  // Params vector
	int npar,ii,jj;                   // NPar
	double *dlta;
	double (*pfunc)(double *p, int np);  //  Likelihood function 
	double rv;                  // return value
};	
void *vcdiag(void *arg) {
	double xxLike(double *Params, int NPar);
	double x,fx,f2H,fH,f2h,fh,dlta; int i; vlltype *llargs = (vlltype*)arg; 
	i=llargs->ii; x=llargs->p[i]; dlta=llargs->dlta[0]; fx=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->p[i]=x+2*dlta; f2H=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->p[i]=x+dlta;    fH=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->p[i]=x-2*dlta; f2h=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->p[i]=x-dlta;    fh=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->rv = -f2H + 16.0*fH - 30.0*fx + 16.0*fh - f2h;
	llargs->p[i]=x; 
	return(llargs);
}
void *vcoffdiag(void *arg) {
	double xxLike(double *Params, int NPar);
	double x,y, dlta_i, dlta_j,fHiHj,fHihj,fhiHj,fhihj; int i,j; vlltype *llargs = (vlltype*)arg; 
	i=llargs->ii;  dlta_i=llargs->dlta[0]; x=llargs->p[i];
	j=llargs->jj;	dlta_j=llargs->dlta[1]; y=llargs->p[j];
	llargs->p[i]=x+dlta_i; llargs->p[j]=y+dlta_j; fHiHj=llargs->pfunc(llargs->p, llargs->npar); 
	                        llargs->p[j]=y-dlta_j; fHihj=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->p[i]=x-dlta_i; llargs->p[j]=y+dlta_j; fhiHj=llargs->pfunc(llargs->p, llargs->npar); 
	                        llargs->p[j]=y-dlta_j; fhihj=llargs->pfunc(llargs->p, llargs->npar); 
	llargs->rv = fHiHj - fHihj - fhiHj + fhihj;
	llargs->p[i]=x; llargs->p[j]=y;
	return(llargs);
}

void varcov_mt(double xxLike(double *p, int N), double *Params, int NPar, double **covar) {
    /*      SUBROUTINE covar computes the var-covar matrix using 
       Finite-difference Approximations of Derivatives (Central difference approximations)*/
   
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
	int InvertMatrix(double **a, int n);
    int ii,jj,i,k,nthreads=SingSpec.nthreads,thnum,lasti,lastj; time_t t1,t2,t0; double ETA=.1e-12, eps; //double ETA=.1e-19, eps;
    float t3=0; double *h; h = new double[NPar+2];

	pthread_t pth[nthreads];

	vlltype *llargs; llargs=new vlltype[nthreads]; time(&t0);
	for (i=0; i<nthreads; i++) { 
		llargs[i].p=new double[NPar+1]; 
		llargs[i].npar=NPar; 
		llargs[i].ii=llargs[i].jj=0;
		llargs[i].dlta=new double[2];
		llargs[i].pfunc=xxLike;
	}	
    eps=pow(ETA, 1.0/3.0);  //  ETA=machine precision (smallest x, such that x+1 not = 1)
    eps=pow(pow(10.0, -SingSpec.LikeNRSig), 1.0/3.0);
    if (SingSpec.Verbose>1) fprintf(g,"LikeNRSig=%d eps=%g ETA=%g\n",SingSpec.LikeNRSig,eps,ETA);
	
	double *save_expval; save_expval=new double[SingSpec.N];
	for (ii=0; ii<SingSpec.N; ii++) save_expval[ii]=SingSpec.expval[ii];

    for (ii=0; ii<NPar; ii++) {
		h[ii] = eps*(1.0+fabs(Params[ii]));
        if (SingSpec.Verbose>1) fprintf(g,"h(%d)=%e parm(%d)=%e\n",ii,h[ii],ii,Params[ii]);
    }
    // set up covariance matrix of beta parameters
    time(&t1); if (SingSpec.Verbose>0) printf("Computing varcov matrix...\n");
	for (k=0; k<nthreads; k++) for (jj=0; jj<NPar; jj++) llargs[k].p[jj] = Params[jj];  //fx = xxLike(Par[0], NPar);, f2H, fH, f2h, fh
    if (SingSpec.novar>1) {
		for (ii=k=lasti=0; ii<NPar; ii++) {     //  for each diagonal element,...    (k=thread number)
			if (1) printf("computing diagonal %d/%d...(%d)\n",ii+1,NPar,SingSpec.ifn);
			llargs[k].ii=ii; llargs[k].dlta[0]=h[ii]; 
			pthread_create(&pth[k],NULL,vcdiag,(void*)&llargs[k]); SingSpec.ifn++;  
			if (++k>=nthreads) {  
				for (thnum=0; thnum<nthreads; thnum++) {
					pthread_join(pth[thnum],NULL);
					covar[lasti+thnum][lasti+thnum] = llargs[thnum].rv / (12.0*h[lasti+thnum]*h[lasti+thnum]);
				}
				k=0; lasti+=nthreads;
			}
		}
		for (thnum=0; thnum<k; thnum++) {
			pthread_join(pth[thnum],NULL); 
			covar[lasti+thnum][lasti+thnum] = llargs[thnum].rv / (12.0*h[lasti+thnum]*h[lasti+thnum]);
		}
		SingSpec.ifn+=k*5;
		for (ii=0; ii<NPar; ii++) {      //        for each row...
			time(&t2); 
			if (difftime(t2,t1)>2.0)  { 
				printf("computing row %d/%d...(%d) %f\n",ii+1,NPar,SingSpec.ifn,difftime(t2,t1)); t1=t2;
			}
			for (k=0,jj=lastj=ii+1; jj<NPar; jj++) {     //   for each column...
				llargs[k].ii=ii; llargs[k].dlta[0]=h[ii]; 
				llargs[k].jj=jj; llargs[k].dlta[1]=h[jj];  //   fHiHj ,fHihj fhiHj   fhihj
				pthread_create(&pth[k],NULL,vcoffdiag,(void*)&llargs[k]); SingSpec.ifn++; 
				if (++k>=nthreads) {  
					for (thnum=0; thnum<nthreads; thnum++) { 
						pthread_join(pth[thnum],NULL);			
						covar[lastj+thnum][ii]=covar[ii][lastj+thnum]= llargs[thnum].rv / (4*h[ii]*h[lastj+thnum]);
					}
					k=0; lastj+=nthreads;
				}
			}
			for (thnum=0; thnum<k; thnum++) {
				pthread_join(pth[thnum],NULL); 
				covar[lastj+thnum][ii]=covar[ii][lastj+thnum]= llargs[thnum].rv / (4*h[ii]*h[lastj+thnum]);
			}
			if (covar[ii][ii]==0.) covar[ii][ii]=eps;
		}
	}
    delete[] h; if (SingSpec.Verbose>0) printf("100%% done           \n");
    if (SingSpec.Verbose>1) {
        fprintf(g,"hessian from varcov:\n");
        for (ii=0; ii<NPar; ii++) {
            for (jj=0; jj<NPar; jj++) fprintf(g," %e",covar[ii][jj]);
            fprintf(g,"\n");
        }
    }
	for (ii=0; ii<SingSpec.N; ii++) SingSpec.expval[ii]=save_expval[ii];
	ii=InvertMatrix(covar,NPar);
	for (ii=0; ii<NPar; ii++) for (jj=0; jj<NPar; jj++) SingSpec.finalVCmat[ii][jj]=covar[ii][jj];
	for (i=0; i<nthreads; i++) { delete [] llargs[i].p; delete [] llargs[i].dlta;}
	delete llargs;	delete[] save_expval;
	time(&t2); t3=t2-t0;  fprintf(g,"CPU time to compute varcov matrix: %1.1f min.\n",t3/60.);
	printf("CPU time to compute varcov matrix: %1.1f min.\n",t3/60.);
	//SingSpec.finalVCmat=covar; 
}
	
void varcov(double xxLike(double *p, int N), double *Params, int NPar, double **covar) {
    /*      SUBROUTINE covar computes the var-covar matrix using 
       Finite-difference Approximations of Derivatives (Central difference approximations)*/
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
	int InvertMatrix(double **a, int n);
    int ii,jj ; time_t t1,t2,t0; float t3=0; double ETA=.1e-12, eps; //double ETA=.1e-19, eps;
    double fHiHj, fHihj, fhiHj, fhihj, fx, f2H, fH, f2h, fh;
    double *h,*Par;    h = new double[NPar+2]; Par=new double[NPar+1];
    eps=pow(ETA, 1.0/3.0);  //  ETA=machine precision (smallest x, such that x+1 not = 1)
    eps=pow(pow(10.0, -SingSpec.LikeNRSig), 1.0/3.0);
    if (SingSpec.Verbose>1) fprintf(g,"LikeNRSig=%d eps=%g ETA=%g\n",SingSpec.LikeNRSig,eps,ETA);
	
	double *save_expval; save_expval=new double[SingSpec.N]; time(&t0);
	for (ii=0; ii<SingSpec.N; ii++) save_expval[ii]=SingSpec.expval[ii];

    for (ii=0; ii<NPar; ii++) {
		Par[ii]=Params[ii]; h[ii] = eps*(1.0+fabs(Params[ii]));
        if (SingSpec.Verbose>1) fprintf(g,"h(%d)=%e parm(%d)=%e\n",ii,h[ii],ii,Params[ii]);
    }
    // set up covariance matrix of beta parameters
    time(&t1); if (SingSpec.Verbose>0) printf("Computing varcov matrix...\n");
    if (SingSpec.novar>1) {
      for (ii=0; ii<NPar; ii++) {  
	    if (SingSpec.Verbose>0) printf("computing row %d/%d...\n",ii+1,NPar);
        fx = xxLike(Params, NPar);
        Par[ii] = Params[ii] + 2*h[ii]; f2H = xxLike(Par, NPar);
        Par[ii] = Params[ii] + h[ii];    fH = xxLike(Par, NPar);
        Par[ii] = Params[ii] - 2*h[ii]; f2h = xxLike(Par, NPar);
        Par[ii] = Params[ii] - h[ii];    fh = xxLike(Par, NPar);
        covar[ii][ii] = (-f2H + 16.0*fH - 30.0*fx + 16.0*fh - f2h) / (12.0*h[ii]*h[ii]);
        fHiHj = ((fH-fx)/h[ii] - (fx-fh)/h[ii])/h[ii];
        Par[ii] = Params[ii];
        for (jj=ii+1; jj<NPar; jj++) {
            Par[ii] = Params[ii] + h[ii]; Par[jj] = Params[jj] + h[jj]; fHiHj = xxLike(Par, NPar);      
                                           Par[jj] = Params[jj] - h[jj]; fHihj = xxLike(Par, NPar); 
            Par[ii] = Params[ii] - h[ii];                                fhihj = xxLike(Par, NPar);										   
                                           Par[jj] = Params[jj] + h[jj]; fhiHj = xxLike(Par, NPar);      

            Par[ii] = Params[ii]; Par[jj] = Params[jj];
            covar[jj][ii]=covar[ii][jj]= (fHiHj - fHihj - fhiHj + fhihj) / (4*h[ii]*h[jj]);
        }
        if (covar[ii][ii]==0.) covar[ii][ii]=eps;
        time(&t2); 
        if (difftime(t2,t1)>2.0) { 
            if (SingSpec.Verbose>0) printf("%1.0f%% done (%d/%d)\r",ii*100/(double)NPar,ii,NPar); 
            t1=t2; 
        }
      }
	}
    delete[] h; if (SingSpec.Verbose>0) printf("100%% done           \n");
    if (SingSpec.Verbose>1) {
        fprintf(g,"hessian from varcov:\n");
        for (ii=0; ii<NPar; ii++) {
            for (jj=0; jj<NPar; jj++) fprintf(g," %e",covar[ii][jj]);
            fprintf(g,"\n");
        }
    }
	for (ii=0; ii<SingSpec.N; ii++) SingSpec.expval[ii]=save_expval[ii];
	printf("invert matrix...\n"); ii=InvertMatrix(covar,NPar); printf("done\n");
	for (ii=0; ii<NPar; ii++) for (jj=0; jj<NPar; jj++) SingSpec.finalVCmat[ii][jj]=covar[ii][jj];
	delete[] save_expval;	delete[] Par;
	time(&t2); t3=t2-t0; fprintf(g,"CPU time to compute varcov matrix: %1.1f min.\n",t3/60.);
	printf("CPU time to compute varcov matrix: %1.1f minutes\n",t3/60.);
	//SingSpec.finalVCmat=covar;
}
