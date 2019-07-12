#include <stdio.h>
#include <math.h>
#include "SingAnal.h"
#define NRANSI
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double amotry(double **p, double y[], double *psum, int ndim,
              double (*funk)(double [], int),  int ihi, double fac) {
    int j;	double fac1,fac2,ytry,*ptry;          // funk needs to be in the form:
    ptry = new double[ndim];                      //    double funk(double *Params, int ndim)
    fac1=(1.0-fac)/ndim; fac2=fac1-fac;
    for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry, ndim);
    if (ytry < y[ihi]) {
        y[ihi]=ytry;
        for (j=0;j<ndim;j++) { psum[j] += ptry[j]-p[ihi][j]; p[ihi][j]=ptry[j]; }
    }
    delete[] ptry; 
    //	printf("amotry:p="); for (j=0; j<ndim; j++) printf(" %f",p[ihi][j]); printf(" amotry=%f\n",ytry);
    return ytry;
}


void amoeba(double **p, double y[], int ndim, double (*funk)(double [], int), int *nfunk){

// funk needs to be in the form ... double funk(double *Params, int ndim)

    double amotry(double **p, double y[], double *psum, int ndim, double (*funk)(double [], int), int ihi, double fac);
	extern TSingSpec SingSpec; FILE *g=SingSpec.g;
    int i,ihi,ilo,inhi,j,mpts=ndim+1,NMAX=ndim*5000;   // maximum number of evaluations
    double rtol,sum=0,swap,ysave,ytry,ftol=.0001;
    
    double *psum = new double[ndim]; *nfunk=0;
    for (j=0;j<ndim;j++) {sum=0; for (i=0;i<mpts;i++) sum+=p[i][j]; psum[j]=sum;}
	
	if (SingSpec.Verbose>1) {
		printf("amoeba:p=\n"); for (j=0; j<ndim; j++) { for (i=0; i<mpts; i++) printf(" %f",p[i][j]); printf(" %f\n",psum[j]);}
	    printf(" amoeba:sum=%f\n",sum);
	}
    for (;;) {
        ilo=0; ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
        for (i=0;i<mpts;i++) {
            if (y[i] <= y[ilo]) ilo=i;
            if (y[i] > y[ihi]) { inhi=ihi; ihi=i;} 
            else if (y[i] > y[inhi] && i != ihi) inhi=i;
        }
        rtol= fabs(y[ihi]-y[ilo]); // convergence criteria
        if (rtol < ftol) {  
			SWAP(y[0],y[ilo])
            for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
            break;
        }
        if (*nfunk >= NMAX) { 
            printf("\n****Maximum function calls (NMAX) exceeded (%d) in amoeba\n",NMAX);  
            fprintf(g,"\n****Maximum function calls (NMAX) exceeded (%d) in amoeba\n",NMAX);  
            delete[] psum; 
            return; 
        }
        *nfunk += 2;
        ytry=amotry(p,y,psum,ndim,funk,ihi,(double) -1.);
        if (ytry <= y[ilo]) {
            i=0; ytry=amotry(p,y,psum,ndim,funk,ihi,(double)2.0);
        }
        else if (ytry >= y[inhi]) {
            ysave=y[ihi]; ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0;i<mpts;i++) {
                    if (i != ilo) {
                        for (j=0;j<ndim;j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=(*funk)(psum, ndim);
                    }
                }
                *nfunk += ndim;
                for (j=0;j<ndim;j++) {for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j]; psum[j]=sum;}
            }
        } else --(*nfunk);
    }
	delete[] psum;
}

void DoAmoeba(double *Params, int NPar, double (*funk)(double [], int), double xincr){

  void amoeba(double **p, double y[], int ndim, double (*funk)(double [], int), int *nfunk);
  int ii,jj,nfunk=0; double *y,**Par; 
  y=new double[NPar+1]; Par=new double*[NPar+1]; Par[0]=new double[NPar];
  for (jj=0; jj<NPar; jj++) Par[0][jj]=Params[jj]=0;
  
  y[0] =(*funk)(Params, NPar); 
  for (ii=1; ii<NPar+1; ii++) {
	  Par[ii]=new double[NPar]; for (jj=0; jj<NPar; jj++) Par[ii][jj] = Params[jj];
      Par[ii][ii-1] += xincr; 
      y[ii] = (*funk)(Params, NPar); 
  }
  amoeba(Par, y, NPar, funk, &nfunk); 
  delete [] y; for (ii=0; ii<=NPar; ii++) delete[] Par[ii]; delete[] Par;
}
#undef SWAP
#undef NRANSI
