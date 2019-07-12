#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "SingAnal.h"
#define MAX(a,b) (a>b?a:b)

double xxLike(double *Params, int NPar);

void optmz2(double *a, int n, double *z, double sig, double *w, int *irx, int mk, double eps) {
/*
C     SUBROUTINE OPTMZ2 IS CALLED BY OPTMIZ WHEN MAXIMIZING THE LOG-LIKELIHOOD FUNCTION.
 */
    int i,j,ii,jj,ij,mm=0;  double ti,v,tim=1/sig,al,r,b,gm,y;	//printf("optmz2...\n");

    if (n<=1) {
      a[1]+=sig*z[1]*z[1]; *irx=1; if (a[1]>0) return;
      a[1]=0; *irx=0; 
	  return;
	}
	if (sig==0 || *irx==0) return;
	if (sig<0) {

      ti=1/sig; jj=0;
      if (mk!=0) {       //                             L*W = Z ON INPUT
          for (j=1; j<=n; j++) {jj+=j; if (a[jj]!=0) ti+=(w[j]*w[j]/a[jj]); }
      }
      else {
          for (j=1; j<=n; j++) w[j]=z[j];
          for (j=1; j<=n; j++) {
              jj+=j; v=w[j];
              if (a[jj]<=0) w[j]=0;
              else {
                  ti+=(v*v)/a[jj];
                  if (j!=n) {
                      ij=jj;
                      for (i=j+1; i<=n; i++) { ij+=(i-1); w[i]-=v*a[ij]; }
                  }
              }
          }
      }
      //   SET TI, TIM AND W
	  if (*irx<=0) { 
	    ti=0; *irx=(-*irx-1); tim=ti; ii=jj; i=n; mm=1;
        for (j=1; j<=n; j++) {
          if (a[ii]!=0) tim=ti-(w[i]*w[i])/a[ii];
          w[i]=ti; ti=tim; ii-=i; i--;
        }
	  }
	  else 
	    if (ti>0) { 
          ti=eps/sig; if (eps==0) (*irx)--; tim=ti; ii=jj; i=n; mm=1;
          for (j=1; j<=n; j++) {
            if (a[ii]!=0) tim=ti-(w[i]*w[i])/a[ii];
            w[i]=ti; ti=tim; ii-=i; i--;
          }
	    }
        else 
		    if ((mk-1)<=0) { mm=0; tim=1/sig;}
	}

    jj=0;   //              UPDATE A   x140
    for (j=1; j<=n; j++) {
        jj+=j; ij=jj; v=z[j]; //       UPDATE A(J,J)
        if (a[jj]>0) {
		  al=v/a[jj]; ti=w[j]; if (mm==0) ti=tim+v*al;
          r=ti/tim; a[jj]*=r; if (r==0 || j==n) break;
          //                   UPDATE REMAINDER OF COL J
          b=al/ti;
          if (r<=4) for (i=j+1; i<=n; i++) { ij+=(i-1); z[i]-=v*a[ij]; a[ij]+=b*z[i];}
          else for (gm=tim/ti, i=j+1; i<=n; i++) { ij+=(i-1); y=a[ij]; a[ij]=b*z[i]+y*gm;	z[i]-=v*y;}
          tim=ti;
		} else {
          if (*irx>0 || sig<0 || v==0) ti=tim; 
		  else {
            *irx=1-*irx; a[jj]=v*v/tim;
            if (j==n) return;
            for (i=j+1; i<=n; i++) { ij+=(i-1); a[ij]=z[i]/v; }
            return;
		  }
		}
    }
    *irx=abs(*irx);
}

#include <pthread.h>

struct lltype { 
	double *p;                  // Params vector
	double dlta;
	int ii, idiff, npar; 
	double (*pfunc)(double *p, int np);  //  Likelihood function 
	double rv;                  // return value
};

void *xxLike2(void *arg) {
	double xxLike(double *Params, int NPar);
	int i; double x,dlta,y1,y2=0; lltype *llargs = (lltype*)arg; 
	i=llargs->ii; dlta=llargs->dlta; x=llargs->p[i]; 
	llargs->p[i]=x+dlta; y1=llargs->pfunc(llargs->p, llargs->npar); 
	if (llargs->idiff==2) { llargs->p[i]=x-dlta; y2=llargs->pfunc(llargs->p, llargs->npar);} 
	llargs->rv=y1-y2; llargs->p[i]=x;
	return(llargs);
}

int optmiz_mt(double xxLike(double *p, int N), int n, int nsig, int maxfn, int iopt, double *x, double *h, double *g, double *rf, double *w) {
	void optmz2(double *a, int n, double *z, double sig, double *w, int *irx, int mk, double eps);
    double hh,eps,hjj,v,df,relx=0,gs0=0,diff,aeps=0,alpha=0,ff=0,tot,f=0,z=0,gys,dgs,sig,zz=0,gnrm;
    int i,j,ij,ir,jj,l,kj,ig,jp1,jb,idiff=1,igg,is,k,ii,im1,nj,ier=0,jnt,threadnum, lasti, nthreads, verbose_nosim=0;
    double reps=1.45519152284e-11, ax=0.1, p1=0.1, f1, f2;
    extern TSingSpec SingSpec; nthreads=SingSpec.nthreads;
	pthread_t pth[nthreads]; 	
	lltype *llargs; llargs=new lltype[nthreads+1];
	for (i=0; i<=nthreads; i++) { 
		llargs[i].p=new double[n+1]; 
		llargs[i].dlta=0; 
		llargs[i].ii=llargs[i].idiff=1;
		llargs[i].npar=n; 
		llargs[i].pfunc=xxLike;
		llargs[i].rv=0; 
		for (j=0; j<n; j++) llargs[i].p[j]=x[j+1];
	}
	SingSpec.optmiz_done=0; SingSpec.ifn=0;  if (SingSpec.Verbose>1 && SingSpec.simulating==0) verbose_nosim=1;
	hh=sqrt(reps); eps=pow(10,(double)-nsig); ig=n; igg=n+n; is=igg; ir=n;
    w[1]=-1; w[2]=w[3]=0;

    for (ij=0, i=1; i<=n; i++) {  // SET H TO THE IDENTITY MATRIX
        for (j=1; j<=i; j++) h[++ij]=0;
        h[ij]=1;
    }
    f=xxLike(x+1,n); SingSpec.ifn=1; df=-1;  //    EVALUATE FUNCTION AT STARTING POINT
x110:
	for (i=0; i<nthreads; i++) for (j=0; j<n; j++) llargs[i].p[j]=x[j+1];
	//printf("eval gradient1...idiff=%d\n",idiff);	
	//for (i=1; i<=n; i++) printf("x(%d)=%20.15f\n",i,x[i]); printf("f=%20.15f hh=%20.15f ax=%20.15f\n",f,hh,ax);
        //  EVALUATE GRADIENT W(IG+I),I=1,...,N
	for (threadnum=i=lasti=0; i<n; i++) {
		z=hh*MAX(fabs(x[i+1]),ax); 
		llargs[threadnum].dlta=z; llargs[threadnum].ii=i; llargs[threadnum].idiff=idiff;
		pthread_create(&pth[threadnum],NULL,xxLike2,(void*)&llargs[threadnum]);  //f1=xxLike(x+1); 
		if ((++threadnum)>=nthreads) {
			for (k=0; k<nthreads; k++) {
				pthread_join(pth[k],NULL);	
				if (idiff == 2) w[ig+lasti+k+1]=llargs[k].rv/(llargs[k].dlta+llargs[k].dlta); // return value = f1-f2
				else w[ig+lasti+k+1]=(llargs[k].rv-f)/llargs[k].dlta; //  return value = f1
				//printf("%d %20.15f\n",SingSpec.ifn+k+1,llargs[k].rv);
			}
			threadnum=0; lasti+=nthreads; SingSpec.ifn+=nthreads*idiff;
			if (SingSpec.ifn>SingSpec.maxfn) { ier=131; goto x410; }
		}		
	}
	for (k=0; k<threadnum; k++) {
		pthread_join(pth[k],NULL); 
		if (idiff == 2) w[ig+lasti+k+1]=llargs[k].rv/(llargs[k].dlta+llargs[k].dlta); // return value = f1-f2
		else w[ig+lasti+k+1]=(llargs[k].rv-f)/llargs[k].dlta; //  return value = f1
		//printf("%d %20.15f\n",SingSpec.ifn+k+1,llargs[k].rv);
	}
	SingSpec.ifn+=threadnum*idiff;
	//printf("gradient done\n"); 	for (i=1; i<=n; i++) printf("%d %15.10f\n",i,w[ig+i]); exit(1);
		if(verbose_nosim>0)printf(">>>>%5d %22.16f            \r",SingSpec.ifn,llargs[k].rv); 
x120:       // BEGIN ITERATION LOOP
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++)  w[i]=-w[ig+i];

    if (ir>=n) {
      //  DETERMINE SEARCH DIRECTION W BY SOLVING H*W = -G WHERE H = L*D*L-TRANSPOSE
      g[1]=w[1];
      if (n>1) {
		for (ii=1, i=2; i<=n; i++) {  //  SOLVE L*W = -G
          ij=ii; ii+=i; v=w[i]; im1=i-1;
          for (j=1; j<=im1; j++) v-=h[++ij]*w[j];
          g[i]=w[i]=v;
        }
        w[n]/=h[ii];  // SOLVE (D*LT)*Z = W WHERE LT = L-TRANSPOSE
        jj=ii; 
        for (nj=1; nj<n; nj++) {
          j=n-nj; jp1=j+1; jj-=jp1; 
          v=w[j]/h[jj]; ij=jj;
          for (i=jp1; i<=n; i++) { ij+=i-1; v-=h[ij]*w[i]; }
          w[j]=v;
        }
	  } else  w[1]/=h[1];
	}
    relx=0; gs0=0;
    for (i=1; i<=n; i++) {
        w[is+i]=w[i]; 
        diff=fabs(w[i])/MAX(fabs(x[i]),ax);
        relx=MAX(relx,diff); gs0+=w[ig+i]*w[i];
    }
	//printf("relx=%f gs0=%f df=%f\n",relx,gs0,df);
    if (relx==0) goto x400;
    aeps=eps/relx; ier=130;
    if (gs0>=0 || df==0) goto x400;
    ier=0; alpha=(-df-df)/gs0;
    if (alpha<=0 || alpha>1) alpha=1;
    if (idiff==2) alpha=MAX(p1,alpha);
    ff=f; tot=0; jnt=0;
x210:
	//printf("x210 ifn=%d\n",SingSpec.ifn);
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++) w[i]=x[i]+alpha*w[is+i];
    f1=xxLike(w+1,n); SingSpec.ifn++; if(verbose_nosim>0)printf("%5d>%22.16f  %22.16f              \n",SingSpec.ifn,f1,f); 
    if (f1>=f) { 
      if (f==ff && idiff==2 && relx>eps) ier=130;
      if (alpha<aeps) goto x400;
      if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
      alpha=0.5*alpha;
      for (i=1; i<=n; i++) w[i]=x[i]+alpha*w[is+i];
      f2=xxLike(w+1,n); SingSpec.ifn++; 
      if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
      if (f2>=f) {
        z=p1; if ((f1+f)>(f2+f2)) z=1+0.5*(f-f1)/(f+f1-f2-f2);
        z=MAX(p1,z); alpha*=z; jnt=1;
        goto x210;
	  }
      tot+=alpha; ier=0; f=f2;
      for (i=1; i<=n; i++) x[i]=w[i];
      if (tot<aeps) goto x400; else goto x320;	
	}
    f2=f; tot+=alpha;
x230:
	//printf("x230 jnt=%d ifn=%d\n",jnt,SingSpec.ifn);
    ier=0; f=f1; for (i=1; i<=n; i++) x[i]=w[i];
    if ((jnt-1)==0) {
		if (tot<aeps) goto x400; else goto x320;
	}
    if ((jnt-1)>0) goto x320;
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++)  w[i]=x[i]+alpha*w[is+i];
    f1=xxLike(w+1,n); SingSpec.ifn++; if(verbose_nosim>0)printf("%5d>>%22.16f            \n",SingSpec.ifn,f1); 
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    if (f1>=f) goto x320;
    if ((f1+f2)>=(f+f) && (7*f1+5*f2)>(12*f)) jnt=2;
    tot+=alpha; alpha+=alpha;
    goto x230;


x320:
	//printf("eval gradient2...idiff=%d\n",idiff);
	//for (i=1; i<=n; i++) printf("x[%d]=%20.15f\n",i,x[i]); printf("f=%20.15f hh=%20.15f ax=%20.15f\n",f,hh,ax);
    alpha=tot; for (i=1; i<=n; i++) w[i]=w[ig+i]; //  SAVE OLD GRADIENT
	
	for (i=0; i<nthreads; i++) for (j=0; j<n; j++) llargs[i].p[j]=x[j+1];
        //  EVALUATE GRADIENT W(IG+I),I=1,...,N
	for (threadnum=i=lasti=0; i<n; i++) {
		z=hh*MAX(fabs(x[i+1]),ax); 
		llargs[threadnum].dlta=z; llargs[threadnum].ii=i; llargs[threadnum].idiff=idiff;
		pthread_create(&pth[threadnum],NULL,xxLike2,(void*)&llargs[threadnum]);  //f1=xxLike(x+1); 
		if ((++threadnum)>=nthreads) {
			for (k=0; k<nthreads; k++) {
				pthread_join(pth[k],NULL);	
				if (idiff == 2) w[ig+lasti+k+1]=llargs[k].rv/(llargs[k].dlta+llargs[k].dlta); // return value = f1-f2
				else w[ig+lasti+k+1]=(llargs[k].rv-f)/llargs[k].dlta; //  return value = f1
				//printf("%d %20.15f\n",SingSpec.ifn+k+1,llargs[k].rv);
			}
			threadnum=0; lasti+=nthreads; SingSpec.ifn+=nthreads*idiff;
			if (SingSpec.ifn>SingSpec.maxfn) { ier=131; goto x410; }
		}		
	}
	for (k=0; k<threadnum; k++) {
		pthread_join(pth[k],NULL); 
		if (idiff == 2) w[ig+lasti+k+1]=llargs[k].rv/(llargs[k].dlta+llargs[k].dlta); // return value = f1-f2
		else w[ig+lasti+k+1]=(llargs[k].rv-f)/llargs[k].dlta; //  return value = f1
		//printf("%d %20.15f\n",SingSpec.ifn+k+1,llargs[k].rv);
	}
	SingSpec.ifn+=threadnum*idiff;
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    gys=0;
    for (i=1; i<=n; i++) { gys+=w[ig+i]*w[is+i]; w[igg+i]=w[i]; }
    df=ff-f; dgs=gys-gs0;
    if (dgs>0) { 
      if ((dgs+alpha*gs0)>0) { //    UPDATE HESSIAN USING DFP FORMULA
        zz=alpha/(dgs-alpha*gs0);  sig=-zz;  optmz2(h,n,w,sig,g,&ir,0,reps);  z=dgs*zz-1;
        for (i=1; i<=n; i++) g[i]=w[ig+i]+z*w[igg+i];
        sig=1/(zz*dgs*dgs); optmz2(h,n,g,sig,w,&ir,0,0);	
	  } else {
        sig=1/gs0; ir=-ir;  optmz2(h,n,w,sig,g,&ir,0,0);  //  UPDATE HESSIAN H USING COMPLEMENTARY DFP FORMULA
        for (i=1; i<=n; i++) g[i]=w[ig+i]-w[igg+i];
        sig=1/(alpha*dgs); ir=-ir; optmz2(h,n,g,sig,w,&ir,0,0);
	  }
	}
    goto x120;
x400:
    if (idiff!=2) { idiff=2; goto x110;} //                CHANGE TO CENTRAL DIFFERENCES
x410:
    if (relx>eps && ier==0) goto x110;
    //             MOVE GRADIENT TO G AND RETURN
	//printf("move gradient to g...\n");
    for (gnrm=0, i=1; i<=n; i++) { g[i]=w[ig+i]; gnrm+=g[i]*g[i]; }
    w[1]=sqrt(gnrm); w[2]=SingSpec.ifn; w[3]=-log10(MAX(reps,relx));
    if(verbose_nosim>1)printf("norm(gradient)=%f No. function calls: %f  No sig digits:%f\n",w[1],w[2],w[3]);
    //             COMPUTE H = L*D*L-TRANSPOSE
    if (n!=1) {  
      jj=n*(n+1)/2;
      for (jb=1; jb<n; jb++) {
          jp1=n+1-jb; jj-=jp1; hjj=h[jj]; ij=jj; l=0;
          for (i=jp1; i<=n; i++) {
              l++; ij+=(i-1); v=h[ij]*hjj; kj=ij;
              for (k=i; k<=n; k++) { h[kj+l]+=h[kj]*v; kj+=k; }
              h[ij]=v;
          }
          hjj=h[jj];
      }
	}
	//printf("\n *** ERROR MESSAGE IER=%d FROM ROUTINE OPTMIZ\n",ier);
	//printf("500 #=%d %18.12f %18.12f %18.12f %18.12f\n",SingSpec.ifn,alpha,aeps,relx,eps);
    *rf=f; 
    if(verbose_nosim>0) printf("\noptmiz done... %1.0f iterations, %3.2f sig.digits.        \n",w[2],w[3]);
	SingSpec.optmiz_done=1;	
	f=xxLike(x+1,n); SingSpec.ifn++; //  call function after multi-threaded optmiz so exp values are right..
	for (i=0; i<n; i++) SingSpec.finalBetaEst[i]=x[i+1];
	for (i=0; i<=nthreads; i++) delete [] llargs[i].p; delete [] llargs;
    return(ier);
}

int optmiz(double xxLike(double *p, int N), int n, int nsig, int maxfn, int iopt, double *x, double *h, double *g, double *rf, double *w) {
    /*      SUBROUTINE OPTMIZ MAXIMIZES THE LOG LIKELIHOOD FUNCTION
C     SUBROUTINE OPTMZ2 IS CALLED.
C     FUNCTION            - A QUASI-NEWTON ALGORITHM FOR FINDING THE
C                           MINIMUM OF A FUNCTION OF N VARIABLES.
C     USAGE               - CALL OPTMIZ(FUNCT,N,NSIG,MAXFN,IOPT,X,H,G,F,W,IER)
C     PARAMETERS   FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           THE FUNCTION F FOR GIVEN PARAMETER VALUES
C                           X(1),X(2),...,X(N).
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
C                           CALL FUNCT(N,X,F)
C                           WHERE X IS A VECTOR OF LENGTH N.
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.
C                N      - THE NUMBER OF PARAMETERS (I.E., THE LENGTH
C                           OF X) (INPUT)
C                NSIG   - CONVERGENCE CRITERION. (INPUT). THE NUMBER
C                           OF DIGITS OF ACCURACY REQUIRED IN THE
C                           PARAMETER ESTIMATES.
C                           THIS CONVERGENCE CONDITION IS SATISIFIED IF
C                           ON TWO SUCCESSIVE ITERATIONS, THE PARAMETER
C                           ESTIMATES (I.E.,X(I), I=1,...,N) AGREE,
C                           COMPONENT BY COMPONENT, TO NSIG DIGITS.
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT)
C                IOPT   - INPUT OPTIONS SELECTOR.
C                         IOPT = 0 CAUSES OPTMIZ TO INITIALIZE THE
C                           HESSIAN MATRIX H TO THE IDENTITY MATRIX.
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
C                           VALUES.
C                         ON INPUT, X MUST CONTAIN THE INITIAL
C                           PARAMETER ESTIMATES.
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER
C                           ESTIMATES AS DETERMINED BY OPTMIZ.
C                H      - VECTOR OF LENGTH N*(N+1)/2 CONTAINING AN
C                           ESTIMATE OF THE HESSIAN MATRIX
C                           D**2F/(DX(I)DX(J)), I,J=1,...,N.
C                           H IS STORED IN SYMMETRIC STORAGE MODE.
C                         ON INPUT, IF IOPT = 0, OPTMIZ INITIALIZES H
C                           TO THE IDENTITY MATRIX. AN INITIAL SETTING
C                           OF H BY THE USER IS INDICATED BY IOPT=1.
C                           H MUST BE POSITIVE DEFINITE. IF IT IS NOT,
C                           A TERMINAL ERROR OCCURS.
C                         ON OUTPUT, H CONTAINS AN ESTIMATE OF THE
C                           HESSIAN AT THE FINAL PARAMETER ESTIMATES
C                           (I.E., AT X(1),X(2),...,X(N))
C                G      - A VECTOR OF LENGTH N CONTAINING AN ESTIMATE
C                           OF THE GRADIENT DF/DX(I),I=1,...,N AT THE
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
C                W      - A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
C                         ON OUTPUT, WORK(I), CONTAINS FOR
C                           I = 1, THE NORM OF THE GRADIENT (I.E.,
C                             SQRT(G(1)**2+G(2)**2+...+G(N)**2))
C                           I = 2, THE NUMBER OF FUNCTION EVALUATIONS
C                             PERFORMED.
C                           I = 3, AN ESTIMATE OF THE NUMBER OF
C                             SIGNIFICANT DIGITS IN THE FINAL
C                             PARAMETER ESTIMATES.
C                IER    - ERROR PARAMETER.
C                         IER = 0 IMPLIES THAT CONVERGENCE WAS
C                           ACHIEVED AND NO ERRORS OCCURRED.
C                         TERMINAL ERROR
C                           IER = 129 IMPLIES THAT THE INITIAL HESSIAN
C                             MATRIX IS NOT POSITIVE DEFINITE. THIS
C                             CAN OCCUR ONLY FOR IOPT = 1.
C                           IER = 130 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED DUE TO ROUNDING ERRORS
C                             BECOMING DOMINANT. THE PARAMETER
C                             ESTIMATES HAVE NOT BEEN DETERMINED TO
C                             NSIG DIGITS.
C                           IER = 131 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.
C     PRECISION           - SINGLE
C     REQD.  ROUTINES     - OPTMZ2
C*/

    double hh,eps,hjj,v,df,relx=0,gs0=0,diff,aeps=0,alpha=0,ff=0,tot,f=0,f1=0,f2=0,z,gys,dgs,sig,zz,gnrm;
    int i,j,ij,ir,nm1,np1,jj,l,kj,ig,jp1,jb,idiff=1,igg,is,k,ii,im1,nj,ier=0,jnt,verbose_nosim=0;
    double reps=1.45519152284e-11, ax=0.1, p1=0.1;
    extern TSingSpec SingSpec;
    SingSpec.optmiz_done=0; if (SingSpec.Verbose>1 && SingSpec.simulating==0) verbose_nosim=1;

    SingSpec.ifn=0; hh=sqrt(reps); eps=pow(10,(double)-nsig); ig=n; igg=n+n; is=igg; ir=n;
    w[1]=-1; w[2]=w[3]=0;

    for (ij=0, i=1; i<=n; i++) {  // SET H TO THE IDENTITY MATRIX
        for (j=1; j<=i; j++) h[++ij]=0;
        h[ij]=1;
    }

    f=xxLike(x+1,n); SingSpec.ifn=1; df=-1; //    EVALUATE FUNCTION AT STARTING POINT
x110:
        //  EVALUATE GRADIENT W(IG+I),I=1,...,N
	//printf("eval gradient1...idiff=%d\n",idiff); for (i=1; i<=n; i++) printf("x(%d)=%20.15f\n",i,x[i]);
	//printf("f=%20.15f hh=%20.15f ax=%20.15f\n",f,hh,ax);
	for (i=1; i<=n; i++) {
        z=hh*MAX(fabs(x[i]),ax); zz=x[i]; x[i]=zz+z; f1=xxLike(x+1,n); 
		if (idiff==2) {                    x[i]=zz-z; f2=xxLike(x+1,n); w[ig+i]=(f1-f2)/(z+z); }
		else                                                              w[ig+i]=(f1-f)/z; 
		x[i]=zz; SingSpec.ifn+=idiff; 
		if (verbose_nosim>0) printf("%d %20.15f\n",SingSpec.ifn,f1);
	}
	//printf("gradient done\n"); for (i=1; i<=n; i++) printf("%d %15.10f\n",i,w[ig+i]); exit(1);
	if(verbose_nosim>0)printf(">>>>%5d %22.16f            \r",SingSpec.ifn,f1); 	
x120:       // BEGIN ITERATION LOOP
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++)  w[i]=-w[ig+i]; 

    if (ir>=n) {
      //  DETERMINE SEARCH DIRECTION W BY SOLVING H*W = -G WHERE H = L*D*L-TRANSPOSE
      g[1]=w[1];
      if (n>1) {
		for (ii=1, i=2; i<=n; i++) {  //  SOLVE L*W = -G
          ij=ii; ii+=i; v=w[i]; im1=i-1;
          for (j=1; j<=im1; j++) v-=h[++ij]*w[j];
          g[i]=v; w[i]=v;
        }
        w[n]/=h[ii];  // SOLVE (D*LT)*Z = W WHERE LT = L-TRANSPOSE
        jj=ii; nm1=n-1;
        for (nj=1; nj<=nm1; nj++) {
          j=n-nj; jp1=j+1; jj-=jp1; 
          v=w[j]/h[jj]; ij=jj;
          for (i=jp1; i<=n; i++) { ij+=i-1; v-=h[ij]*w[i]; }
          w[j]=v;
        }
	  } else  w[1]/=h[1];
	}
    relx=0; gs0=0;
    for (i=1; i<=n; i++) {
        w[is+i]=w[i]; 
        diff=fabs(w[i])/MAX(fabs(x[i]),ax);
        relx=MAX(relx,diff); gs0+=w[ig+i]*w[i];
    }
	//printf("relx=%f gs0=%f df=%f\n",relx,gs0,df);
    if (relx==0) goto x400;
    aeps=eps/relx; ier=130;
    if (gs0>=0 || df==0) goto x400;
    ier=0; alpha=(-df-df)/gs0;
    if (alpha<=0 || alpha>1) alpha=1;
    if (idiff==2) alpha=MAX(p1,alpha);
    ff=f; tot=0; jnt=0;
x210:
	//printf("x210\n");
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++) w[i]=x[i]+alpha*w[is+i];
    f1=xxLike(w+1,n); SingSpec.ifn++; if(verbose_nosim>0)printf("%5d>%22.16f  %22.16f              \n",SingSpec.ifn,f1,f); 
    if (f1>=f) { 
      if (f==ff && idiff==2 && relx>eps) ier=130;
      if (alpha<aeps) goto x400;
      if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
      alpha=0.5*alpha;
      for (i=1; i<=n; i++) w[i]=x[i]+alpha*w[is+i];
      f2=xxLike(w+1,n); SingSpec.ifn++;
      if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
      if (f2>=f) {
        z=p1; if ((f1+f)>(f2+f2)) z=1+0.5*(f-f1)/(f+f1-f2-f2);
        z=MAX(p1,z); alpha*=z; jnt=1;
        goto x210;
	  }
      tot+=alpha; ier=0; f=f2;
      for (i=1; i<=n; i++) x[i]=w[i];
      if (tot<aeps) goto x400; else goto x320;	
	}
    f2=f; tot+=alpha;
x230:
    ier=0; f=f1; for (i=1; i<=n; i++) x[i]=w[i];
    if ((jnt-1)==0) {
		if (tot<aeps) goto x400; else goto x320;
	}
    if ((jnt-1)>0) goto x320;
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    for (i=1; i<=n; i++)  w[i]=x[i]+alpha*w[is+i];
    f1=xxLike(w+1,n); SingSpec.ifn++; if(verbose_nosim>0)printf("%5d>>%22.16f            \n",SingSpec.ifn,f1); 
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};
    if (f1>=f) goto x320;
    if ((f1+f2)>=(f+f) && (7*f1+5*f2)>(12*f)) jnt=2;
    tot+=alpha; alpha+=alpha;
    goto x230;

x320:
    alpha=tot; for (i=1; i<=n; i++) w[i]=w[ig+i]; //  SAVE OLD GRADIENT
       // EVALUATE GRADIENT W(IG+I), I=1,...,N
	//printf("eval gradient2...idiff=%d\n",idiff); for (i=1; i<=n; i++) printf("x[%d]=%20.15f\n",i,x[i]);
	//printf("f=%20.15f hh=%20.15f ax=%20.15f\n",f,hh,ax);
	for (i=1; i<=n; i++) {
        z=hh*MAX(fabs(x[i]),ax); zz=x[i]; x[i]=zz+z;  f1=xxLike(x+1,n); 
		if (idiff==2) {                    x[i]=zz-z;  f2=xxLike(x+1,n); w[ig+i]=(f1-f2)/(z+z); }
		else                                                               w[ig+i]=(f1-f)/z; 
		x[i]=zz; SingSpec.ifn+=idiff; 
		if (verbose_nosim>0) printf("%d %20.15f\n",SingSpec.ifn,f1);
	}
    if (SingSpec.ifn>=maxfn) {ier=131; goto x410;};	
    for (gys=0,i=1; i<=n; i++) { gys+=w[ig+i]*w[is+i]; w[igg+i]=w[i]; }
    df=ff-f; dgs=gys-gs0;
    if (dgs>0) { 
      if ((dgs+alpha*gs0)>0) { //    UPDATE HESSIAN USING DFP FORMULA
        zz=alpha/(dgs-alpha*gs0);  sig=-zz;  optmz2(h,n,w,sig,g,&ir,0,reps);  z=dgs*zz-1;
        for (i=1; i<=n; i++) g[i]=w[ig+i]+z*w[igg+i];
        sig=1/(zz*dgs*dgs); optmz2(h,n,g,sig,w,&ir,0,0);	
	  } else {
        sig=1/gs0; ir=-ir;  optmz2(h,n,w,sig,g,&ir,0,0);  //  UPDATE HESSIAN H USING COMPLEMENTARY DFP FORMULA
        for (i=1; i<=n; i++) g[i]=w[ig+i]-w[igg+i];
        sig=1/(alpha*dgs); ir=-ir; optmz2(h,n,g,sig,w,&ir,0,0);
	  }
	}
    goto x120;
x400:
    if (idiff!=2) { idiff=2; goto x110;} //                CHANGE TO CENTRAL DIFFERENCES
x410:
    if (relx>eps && ier==0) goto x110;
    //             MOVE GRADIENT TO G AND RETURN
    for (gnrm=0, i=1; i<=n; i++) { g[i]=w[ig+i]; gnrm+=g[i]*g[i]; }
    w[1]=sqrt(gnrm); w[2]=SingSpec.ifn; w[3]=-log10(MAX(reps,relx));
    if(verbose_nosim>1)printf("norm(gradient)=%f No. function calls: %f  No sig digits:%f\n",w[1],w[2],w[3]);
    //             COMPUTE H = L*D*L-TRANSPOSE
    if (n!=1) {  
      np1=n+1; nm1=n-1; jj=n*np1/2;
      for (jb=1; jb<=nm1; jb++) {
          jp1=np1-jb; jj-=jp1; hjj=h[jj]; ij=jj; l=0;
          for (i=jp1; i<=n; i++) {
              l++; ij+=(i-1); v=h[ij]*hjj; kj=ij;
              for (k=i; k<=n; k++) { h[kj+l]+=h[kj]*v; kj+=k; }
              h[ij]=v;
          }
          hjj=h[jj];
      }
	}
	//printf("\n *** ERROR MESSAGE IER=%d FROM ROUTINE OPTMIZ\n",ier);
	//printf("500 #=%d %18.12f %18.12f %18.12f %18.12f\n",SingSpec.ifn,alpha,aeps,relx,eps);
    *rf=f; 
    if(verbose_nosim>0) printf("\noptmiz done... %1.0f iterations, %3.2f sig.digits.        \n",w[2],w[3]);
	SingSpec.optmiz_done=1;	for (i=0; i<n; i++) SingSpec.finalBetaEst[i]=x[i+1];

    return(ier);
}

void try_optmz_randiv(double (*xxLinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *LL) {
	//         try optimization routine with user-specified number of random initial value vectors
	//         and save results from vector which gave the maximum log-likelihood
	double Random(void);
    int optmiz(double xxlike(double *Params, int NPar),int NPar, int Sig, int MaxFN, int iopt, double *Params, double *Hess, double *Grad, double *LL, double *Work);
	
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g; //	MSOpenLinkFn=&MSOpenLink; 
	
	double *Grad1, *iv, *bestIV, minLL=9e44, maxLL=-9e44, *LLsv, LLscale=1; time_t t1,t2;
	int i,j,kk,itry,mx=1, *LLhist, hist_size=20; 
	printf("try_optmz_randiv... NPar=%d\n",NPar);
	Grad1 = new double[NPar+1]; iv=new double[NPar]; bestIV=new double[NPar];  
	LLsv=new double[SingSpec.nrandIV+1]; LLhist=new int[hist_size]; for (i=0; i<hist_size; i++) LLhist[i]=0;
	
	for (i=0; i<NPar; i++) iv[i]=bestIV[i]=Params[i];
	for (itry=0; itry<=SingSpec.nrandIV; itry++) {
	  if (itry>0) for (i=0; i<NPar; i++) iv[i]=Params[i]=Random()*10-5; SingSpec.ifn=0;  time(&t1);
      if (SingSpec.nthreads<2) 
		  optmiz(xxLinkFn,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1, LL, Work); 
	  else {
		  printf("optimizing using %d threads...\n",SingSpec.nthreads);
		  fprintf(g,"optimizing using %d threads...\n",SingSpec.nthreads);
		  optmiz_mt(xxLinkFn,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1, LL, Work); 
	  }		  
	  *LL=(*LL); LLsv[itry]=(*LL); time(&t2);
	  if (SingSpec.Verbose>0) 
		  printf("attempt %d/%d (iv=%5.2f %5.2f %5.2f...), -2loglike=%9.4f,  optmiz time:%1.1f min\n",
				itry,SingSpec.nrandIV,iv[0],iv[1],iv[2],*LL*2,(t2-t1)/60.); 
	  fprintf(g,"attempt %d/%d (iv=%5.2f %5.2f %5.2f...), -2loglike=%9.4f,  optmiz time:%1.1f min.\n",
		itry,SingSpec.nrandIV,iv[0],iv[1],iv[2],*LL*2,(t2-t1)/60.); 
	  if (*LL < minLL) { minLL=*LL; for (i=0; i<NPar; i++) bestIV[i]=iv[i]; }
	  if (*LL > maxLL) maxLL=*LL;
	}
	if (SingSpec.nrandIV>0) {
	  LLscale=(maxLL-minLL)/(hist_size-1);
	  for (i=0; i<=SingSpec.nrandIV; i++) { j=floor((LLsv[i]-minLL)/LLscale); LLhist[j]++;}
	  for (i=0; i<hist_size; i++) if (LLhist[i]>mx) mx=LLhist[i];
	
      fprintf(g,"\nDistribution of -2logLike:\n");
      for (i=0; i<hist_size; i++) {
        fprintf(g,"%6.2f %5d:",minLL+i*LLscale,LLhist[i]);
        kk=LLhist[i]*50/mx; 
		if (kk>0) for (j=0; j<kk; j++) fprintf(g,"*"); fprintf(g,"\n");
        if (SingSpec.Verbose>0) printf("%6.2f %5d:",minLL+i*LLscale,LLhist[i]);
        if (kk>0) for (j=0; j<kk; j++) { if (SingSpec.Verbose>0) printf("*"); printf("\n");}
      }
	}
	if (SingSpec.nrandIV>0)  {
		for (i=0; i<NPar; i++) Params[i]=bestIV[i];
		optmiz(xxLinkFn,NPar, SingSpec.LikeNRSig, SingSpec.maxfn, 0, Params-1, Hess, Grad1, LL, Work); 
	}
	SingSpec.finalLL=*LL; SingSpec.finalBetaEst=Params; LLscale=xxLinkFn(Params,NPar); 
	
	if (SingSpec.Verbose>0) printf("%d random initial value vectors (+ default initial vector) attempted:\n",SingSpec.nrandIV);
	fprintf(g,"%d random initial value vectors attempted:\n",SingSpec.nrandIV);
	fprintf(g,"Initial values:\n"); 
	for (i=0; i<NPar; i++) { fprintf(g," %9.4f",bestIV[i]); if ((i%8)==7) fprintf(g,"\n"); }
	fprintf(g,"\n");

	delete[] iv; 
	delete[] bestIV; 
	delete[] Grad1; 
	delete[] LLsv; 
	delete [] LLhist;
}
