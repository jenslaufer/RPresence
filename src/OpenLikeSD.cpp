#include <math.h>
#include "SingAnal.h"
#include <stdio.h>
#include <stdlib.h>
//#define DBG 1
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif
double OpenLikeSD(double ****Phi, double ***ClosedP, double ***Theta, double *th0pi, double pi) {
    int site, seasn, srvy, allmiss, no_det_yet, kk, nlog0, lastobs,missng,seen,PrmyPers; 
    double p,p2,prd,prd0a,prd0b,prd1a,prd1b,lprd0a,lprd0b,lprd1a,lprd1b,pstar,pstar2,qstar,qstar2,
          theta,theta_prime,sum=0,prdOcc,prdNotOcc,prdOcc_save,prdNotOcc_save,cellProb,
          prdOcc2,prdNotOcc2,prdOcc_save2,prdNotOcc_save2;
    int dbgsite=-15;
    extern TSingSpec SingSpec; PrmyPers=SingSpec.PrmyPeriods;
#ifdef DBG
    for (seasn=0; seasn<PrmyPers; seasn++) 
      printf("phi(%d):%f %f, %f %f\n",seasn,Phi[0][0][0][seasn],Phi[0][0][1][seasn],Phi[0][1][0][seasn],Phi[0][1][1][seasn]);
    printf("theta=%f %f  th0pi=%f  p=%f\n",Theta[0][0][0],Theta[0][0][1],th0pi[0],ClosedP[0][0][0]);
#endif
    for (site=nlog0=0; site<SingSpec.N; site++) {  cellProb=0;
        
        prdOcc = prdOcc2 = Phi[site][0][0][0]; prdNotOcc = prdNotOcc2 = 1-prdOcc;
        for (srvy=0, lastobs=-1; srvy<SingSpec.T; srvy++) 
            if (SingSpec.Data[site][srvy]!=-1) lastobs=srvy+1;
        
        for (seasn=srvy=0; seasn<PrmyPers; seasn++) {   //  prd0,prd1 = prod of p/q's within season         
            prd0a=prd0b=prd=pstar=pstar2=qstar=qstar2=1; no_det_yet=allmiss=1; 
            seen=missng=allmiss=0; prd1a=prd1b=0;
            if (site==dbgsite) printf("prdOcc:%f %f  prd:%f %f\n",prdOcc,prdNotOcc,prd0a,prd1a);
            for (kk=0; kk<SingSpec.SecPeriods[seasn]; kk++, srvy++) {   
                theta=Theta[site][srvy][0]; theta_prime=Theta[site][srvy][1];
                if (no_det_yet) {   //  if no detections yet...
                    if (th0pi[seasn]>1) {//  option for theta for starting in middle of trail
						if (theta<.1e-5) theta=.1e-5;
                        theta=theta/(theta+1-theta_prime);  //  using equalibrium assumption 
					}
                    else {
                        theta=th0pi[seasn]*theta_prime+(1-th0pi[seasn])*theta;  //  theta for 1st segment
                    }
                }
                p=p2=ClosedP[site][srvy][0]; if (pi>=0) p2=ClosedP[site][srvy][1];
                if (SingSpec.Data[site][srvy]!=-1 && p>=0) {  //  if not missing data in this survey...
                    qstar*=(1-p); qstar2*=(1-p2);
                    if (SingSpec.Data[site][srvy]>0) { pstar*=p; pstar2*=p2;}
                    else { pstar*=(1-p); pstar2*=(1-p2);}
                }
                lprd0a=prd0a; lprd1a=prd1a; lprd0b=prd0b; lprd1b=prd1b; 
                if (site==dbgsite) printf("%d,%d theta=%f %f p=%f p*=%f q*=%f\n",seasn,kk,theta,theta_prime,p,pstar,qstar);
                if (SingSpec.Data[site][srvy]>=0 && p>=0) {
                    if (SingSpec.Data[site][srvy]>0) { //   if detected in this survey...
                        prd0a= prd0b=0;
                        prd1a = (lprd0a*theta+lprd1a*theta_prime)*p;
                        prd1b = (lprd0b*theta+lprd1b*theta_prime)*p2;
                        seen=true;
                    } 
                    else {              //   not detected in this primary period
                        prd0a= lprd0a*(1-theta)+lprd1a*(1-theta_prime);
                        prd0b= lprd0b*(1-theta)+lprd1b*(1-theta_prime);
                        prd1a = (lprd0a*theta+lprd1a*theta_prime)*(1-p);
                        prd1b = (lprd0b*theta+lprd1b*theta_prime)*(1-p2);
                    }
                    allmiss=no_det_yet=0;
                    if (site==dbgsite) printf("prd=%15.12f %15.12f\n",prd0a,prd1a);
                }
            }   //  end of survey in season loop
            if (seen) prdNotOcc=prdNotOcc2=0;  
            if(missng<0)prd0a=prd0b=prd1a=prd1b=1;
            prdOcc  *= (prd0a+prd1a);   prdOcc2    *= (prd0b+prd1b);              
            prdOcc_save=prdOcc;  prdNotOcc_save=prdNotOcc; // save prds so they can be used below
            prdOcc_save2=prdOcc2;  prdNotOcc_save2=prdNotOcc2; // save prds so they can be used below
            if (seasn<(PrmyPers-1) && srvy<lastobs) {  //  Phi(0,0)=psi or 1-eps, Phi(0,1)=1-psi or eps
                prdOcc = prdOcc_save*Phi[site][0][0][seasn+1] + prdNotOcc_save*Phi[site][1][0][seasn+1];
                prdNotOcc = prdOcc_save*Phi[site][0][1][seasn+1] + prdNotOcc_save*Phi[site][1][1][seasn+1];
                prdOcc2 = prdOcc_save2*Phi[site][0][0][seasn+1] + prdNotOcc_save2*Phi[site][1][0][seasn+1];
                prdNotOcc2 = prdOcc_save2*Phi[site][0][1][seasn+1] + prdNotOcc_save2*Phi[site][1][1][seasn+1];
            }
            if (site==dbgsite) printf("prdOcc:%15.12f %15.12f  (save):%15.12f %15.12f\n",prdOcc,prdNotOcc,prdOcc_save,prdNotOcc_save);
        } // end PrmyPeriods loop
        cellProb += pi*(prdOcc+prdNotOcc)+(1-pi)*(prdOcc2+prdNotOcc2);
        if (site==dbgsite) printf("cellProb:%15.12f pi:%15.12f prdOcc:%15.12f %15.12f  (save):%15.12f %15.12f\n",cellProb,pi,prdOcc,prdNotOcc,prdOcc2,prdNotOcc2);
        //printf("cellProb(%d)=%15.12f %f %f\n",site,1500*cellProb,prdOcc,prdNotOcc);
        //printf("pi=%f prd0a=%f prd0b=%f prd1a=%f prd1b=%f\n",pi,prd0a,prd0b,prd1a,prd1b);
        //printf("pstar=%f pstar2=%f\n",pstar,pstar2);
        if (cellProb < SingSpec.nearzero) { nlog0++; cellProb=SingSpec.nearzero; }
        sum += SingSpec.det_hist_frq[site]*log(cellProb); 
        SingSpec.expval[site]=cellProb;
    } // end individual loop
    
    site=0;
#ifdef DBG
    printf("OpenLikeSD:2*sum=%f (%d) %f %f %f %f %f\n",-sum*2,nlog0,Phi[site][0][0][0],
                 Phi[site][0][1][1],Phi[site][1][0][1],Phi[site][1][1][1],ClosedP[0][0][0]);
#endif
	//printf("%f     %d \n",sum,SingSpec.ifn);
    return(-sum);
}
