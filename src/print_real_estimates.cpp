#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

double xxLike(double *Params, int NPar);

int print_real_estimates(double xxLike(double *p, int N), int NPar, double *Params, double **covar) {
    extern TSingSpec SingSpec;
    int i,n,k,site,idm,iest,irow,betapar,N=SingSpec.N,	NrealPar=0; 
	double se,ci1,ci2,var,like,eps=.00001, **varEst, **dpsi_dth, *prd,*realpar;
	FILE *g=SingSpec.g; time_t t1,t2; float t3=0; time(&t1);

	if (SingSpec.novar>3) {
	  for (i=0; i<6; i++) NrealPar+=SingSpec.NrowsDM[i];
      prd=new double[NPar]; dpsi_dth=new double*[NPar]; for (k=0; k<NPar; k++) { prd[k]=0; dpsi_dth[k]=new double[N]; }
	  realpar = new double[N]; varEst = new double*[NrealPar]; for (i=0; i<NrealPar; i++) varEst[i] = new double[N];
	  for (i=0; i<NPar; i++) Params[i]=SingSpec.finalBetaEst[i];
	  SingSpec.optmiz_done=1; like=xxLike(Params,NPar); se=like;

	  fprintf(g,"\nReal parameter estimates\n=======================================================================\n");
	  for (idm=iest=0; idm<6; idm++) { 
	      for (irow=0; irow<SingSpec.NrowsDM[idm]; irow++,iest++) {   
	  		for (n=1,site=0; site<N; site++) {
	  			realpar[site]=SingSpec.realParmEst[iest][site];  
	  			if (site>0) if (realpar[site]!=realpar[site-1]) n=N;
	  		}
	  		for (betapar=0; betapar<NPar; betapar++) {
	  			Params[betapar]+=eps; like=xxLike(Params,NPar);	 // compute derivative of realParm wrt betaParm		
	  			for (site=0; site<n; site++) {
	  			    dpsi_dth[betapar][site]=(SingSpec.realParmEst[iest][site]-realpar[site])/eps;
	  			}
	  			Params[betapar]-=eps; like=xxLike(Params,NPar);
	  		}
	  		for (site=0; site<n; site++) {
	  			for (betapar=0; betapar<NPar; betapar++) { 
	  				for (k=0,prd[betapar]=0; k<NPar; k++) {
	  					prd[betapar]+=dpsi_dth[k][site]*covar[k][betapar];
	  				}
	  			}
	  			for (k=0,var=0; k<NPar; prd[k++]=0) var+=prd[k]*dpsi_dth[k][site];
	  			varEst[iest][site]=var;
	  		}
            if (irow==0) fprintf(g,"\n   parameter     site                    estimate    std.error       95%% conf. interval\n");
			if (n>SingSpec.lmt) n=SingSpec.lmt+1;
			for (site=0; site<n; site++) {
				if (SingSpec.fixed[iest]< -998) {
				    se=sqrt(varEst[iest][site]); ci1=realpar[site]-1.96*se; ci2=realpar[site]+1.96*se;
				    if (ci1<0) ci1=0; if (ci2>1) ci2=1;
                    fprintf(g,"%-12s %4d %-16s: %12.5f %12.5f  %12.5f - %8.5f\n",SingSpec.Realname[idm][irow],
					                site+1,SingSpec.sitename[site],realpar[site],se,ci1,ci2);
                }
                else {
                    if (site==0) fprintf(g,"%12s **** %12.5f                                      *** fixed\n",
					     SingSpec.Realname[idm][irow],SingSpec.fixed[iest]);   // realpar[site]);
				}
			}	//  end for (site=0...	
		}    //  end for(irow=0...
	  }     //  end for(idm=iest=0...
	  for (i=0; i<NrealPar; i++) free(varEst[i]); 
	  for (i=0; i<NPar; i++) free(dpsi_dth[i]);
	  free(realpar); free(varEst); free(dpsi_dth); free(prd);
	  time(&t2); t3=t2-t1; printf("cpu time to print real parameter estimates: %1.0f\n",t3);
	}  // end if (SingSpec.novar>3)
	return(0);
}
