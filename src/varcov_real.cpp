#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "SingAnal.h"

double xpsi(int i, double *Params) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    int ii=0,jj=0; double psi,gam,eps;
    if (i>=(SingSpec.PrmyPeriods-1)) i=i-SingSpec.PrmyPeriods+1; // eps(psi) instead of gam(psi)
    psi=getParam(0,0,0,0,0,Params,0,1); SingSpec.psibar[0][0][0]=psi;
    if (i<1) return(psi);
    for (jj=1; jj<SingSpec.PrmyPeriods; jj++) {
        gam=getParam(jj,1,jj-1,ii,jj-1,Params,0,1);
        eps=getParam(SingSpec.PrmyPeriods+jj-1,2,jj-1,ii,jj-1,Params,0,1);
        if (SingSpec.Alt_Param_Checked) eps=1-eps;  // alt-parm: use persistance instead of extinction
        psi=psi*(1.-eps)+(1.-psi)*gam; SingSpec.psibar[jj][0][0]=psi;
        if (jj==i) return(psi);
    }
    return(psi);
}

double realparm(int i, int *parmfixed, double *Params, double **bigdm, int **bigdm_ptr, int ncols, int site, int srvy, int *lnkfn) {
    extern TSingSpec SingSpec;
    double p1, sum=0, tmpdm; int j,tmpptr;
    p1=SingSpec.fixed[i]; *parmfixed=0;
    if (p1<=-998) {
        for (j=0; j<ncols; j++) {
            tmpdm=bigdm[i][j]; tmpptr=bigdm_ptr[i][j];
            if(tmpptr>3000)  tmpdm=xpsi(i-1,Params);
			else {
				if(tmpptr<0) tmpdm=SingSpec.SampCov[-tmpptr-1][site][srvy];
				if(tmpptr>0) tmpdm=SingSpec.SiteCov[site][tmpptr-1];
			}
            sum+=tmpdm*Params[j];
        }
        switch (lnkfn[i]) {
            case 0: 
            case logitLnk: p1 = (sum>-60 ? (sum<60 ? 1/(1+exp(-sum)) : 1) : 0); break;
            case loglogLnk: p1 = (sum>-60 ? (sum<60 ? 1-exp(-exp(sum)) : 1-1.0e-50) : 1.0e-50); break;
            case expLnk: p1 = (sum<50. ? exp(sum) : exp(50.)); break;
            case sinLnk: p1 = (sin(sum)+1.)/2.; break;
            case IDLnk: p1 = sum; break;
            default: p1 = (sum>-60 ? (sum<60 ? 1/(1+exp(-sum)) : 1) : 0);
		}
    }
    else *parmfixed=1; //printf("i=%d p1=%20.17f\n",i,p1);
    return(p1);
}

void varcov_real(int site, int srvy, double *Params, int NPar, double **covar) {
    //      SUBROUTINE covar computes the var-covar matrix of real parameters using delta method
    extern TSingSpec SingSpec; FILE *g=SingSpec.g; char **name;
    int i,j,k,nrows,ncols=0,dm,dm2,ii,jj,*lnkfn, parmfixed, **bigdm_ptr;
    double p1,p2,sum,**gradXcov,**bigdm,**grad1,se,eps=.0000000001,*p_real,**vc_real;
	float t3; time_t t1,t2; time(&t1);
    if (SingSpec.Verbose>0) printf("Computing varcov matrix of real parameters...\n");
    
    for (dm=0,nrows=0; dm<6; dm++)                 //  count no. rows,cols in big design matrix
        for (i=0; i<SingSpec.NrowsDM[dm]; i++,nrows++)
            for (dm2=0,ncols=0; dm2<6; dm2++)
                for (j=0; j<SingSpec.NParKK[dm2]; j++) ncols++;
    
    name=new char*[nrows]; for (i=0; i<nrows; i++) name[i]=new char [80];
    bigdm=new double*[nrows]; bigdm_ptr=new int*[nrows];  
	for (i=0; i<nrows; i++) { bigdm[i]=new double [ncols]; bigdm_ptr[i]=new int[ncols];}
    grad1=new double*[nrows];  for (i=0; i<nrows; i++) grad1[i]=new double [ncols];
    lnkfn=new int[nrows];
    gradXcov=new double*[nrows];  for (i=0; i<nrows; i++) gradXcov[i]=new double [ncols];
    vc_real=new double*[nrows]; p_real=new double[nrows]; 
    for (i=0; i<nrows; i++) vc_real[i]=new double [nrows];
    for (i=0; i<nrows; i++) for (j=0; j<ncols; j++) bigdm[i][j]=bigdm_ptr[i][j]=0;
    ii=jj=0;//   paste 6 small design matrices into big design matrix
    for (dm=0; dm<6; dm++) {
        for (i=0; i<SingSpec.NrowsDM[dm]; i++) 
            for (j=0; j<SingSpec.NParKK[dm]; j++) {
                bigdm[i+ii][j+jj]=SingSpec.DMat[dm][i][j]; 
				bigdm_ptr[i+ii][j+jj]=SingSpec.DMat_ptr[dm][i][j];
                strcpy(name[i+ii],SingSpec.Realname[dm][i]);
            }
        ii+=SingSpec.NrowsDM[dm]; jj+=SingSpec.NParKK[dm];
    }
#ifdef dbg
    printf("parms:\n"); for (i=0; i<ncols; i++) printf(" %20.17f",Params[i]); printf("\n");
    printf("------------\ndesign matrix:\n");
    for (i=0; i<nrows; i++) {
	    printf("%12s ",name[i]);
    	for (j=0; j<ncols; j++) 
			if (bigdm_ptr[i][j]==0) printf(" %f",bigdm[i][j]); 
			else 
				if (bigdm_ptr[i][j]>0) printf(" %8.8s",SingSpec.CovNames[bigdm_ptr[i][j]-1]); 
				else printf(" %8s",SingSpec.CovNames2[-bigdm_ptr[i][j]+1]); 
		printf("\n");
	}
    printf("------------\nbeta var-cov matrix:\n");
    for (i=0; i<ncols; i++) {for (j=0; j<ncols; j++) printf(" %f",covar[i][j]); printf("\n");}
#endif
    for (i=0; i<nrows; i++) {                //   compute gradient = dbeta/dreal ... dim=nrows,ncols
        p1=p2=p_real[i]=realparm(i,&parmfixed,Params,bigdm,bigdm_ptr,ncols,site,srvy,lnkfn); 
        for (k=0; k<ncols; k++) {
            Params[k]+=eps; //   if param not fixed...
            if (parmfixed==0) p2=realparm(i,&parmfixed,Params,bigdm,bigdm_ptr,ncols,site,srvy,lnkfn);
            grad1[i][k]=(p2-p1)/eps;
            Params[k]-=eps;
        }
		printf("%d/%d\n",i+1,nrows);
    }                                      //  matrix mult grad1*covar
#ifdef dbg
    fprintf(g,"grad1:\n");
    for (i=0; i<nrows; i++) {for (j=0; j<ncols; j++) fprintf(g," %f",grad1[i][j]); fprintf(g,"\n");}
#endif
    for (i=0; i<nrows; i++) {
        for (j=0; j<ncols; j++) {
            for (sum=0,k=0; k<ncols; k++) sum+=grad1[i][k]*covar[k][j]; gradXcov[i][j]=sum;
        }
    }                                      // matrix mult grad1*covar*grad1'  
    for (i=0; i<nrows; i++) {
        for (j=0; j<nrows; j++) {
            for (sum=0,k=0; k<ncols; k++) sum+=gradXcov[i][k]*grad1[j][k]; vc_real[i][j]=sum;
        }
    }  
    if (SingSpec.novar>3) {
        fprintf(g,"\nReal parameters: (computed using covariates from 1st site and 1st survey)\n\n");
        fprintf(g,"Real parameter :              estimate   SE(estimate)\n");
        for (i=0; i<nrows; i++) { 
            se=sqrt(vc_real[i][i]); //ci1=p_real[i]-1.96*se; ci2=p_real[i]+1.96*se;
            fprintf(g,"%6i %-20s %10.4f  %10.4f\n",i+1,name[i],p_real[i],se);
        }
        if (SingSpec.novar==6) {
			fprintf(g,"Variance-covariance matrix of real parameters:\n");
			for (i=0; i<nrows; i++) {
				for (j=0; j<nrows; j++) {
					fprintf(g," %10.6f",vc_real[i][j]);
					if ((j%8)==7) fprintf(g,"\n");
				}
				fprintf(g,"\n");
			}
		}
    }
    for (i=0; i<nrows; i++) { delete[] gradXcov[i]; delete[] grad1[i]; delete[] name[i]; }
    delete[] grad1; delete[] gradXcov; delete[] name; delete[] lnkfn;
	 time(&t2); t3=t2-t1; printf("\nVarCov(real):CPU time= %1.0f seconds (%1.2f min)\n",t3,t3/60.);
}
