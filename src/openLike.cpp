#include <math.h>
#include "SingAnal.h"
#include <stdio.h>
#include <stdlib.h>
//#define dbg 1
double OpenLike(double ****Phi, double ***ClosedP11, double ***ClosedP10, double ***ClosedB, double *pi) { 
    int site, seasn, srvy, kk, nlog0, lastobs,missng,seen,ngrps,PrmyPers,fpflag=0; 
    double OpenP0[2],OpenP2[2],sum=0,prdOcc,prdNotOcc,prdOcc_save,prdNotOcc_save,cellProb=0;
    extern TSingSpec SingSpec; PrmyPers=SingSpec.PrmyPeriods;
    ngrps=SingSpec.NrowsDM[4]+1; if (ngrps<1) ngrps=1;
    OpenP0[0]=OpenP0[1]=OpenP2[0]=OpenP2[1]=0;
    if (SingSpec.FalsePos==1) fpflag=1; //  see if (seen>fpflag) below... 
    for (site=nlog0=0; site<SingSpec.N; site++) {  cellProb=0;
        
        prdOcc = Phi[site][0][0][0]; prdNotOcc = 1-prdOcc;
        for (srvy=0, lastobs=-1; srvy<SingSpec.T; srvy++) 
            if (SingSpec.Data[site][srvy]!=-1) lastobs=srvy+1;
        for (seasn=srvy=0; seasn<PrmyPers; seasn++) {            
            missng=seen=0; OpenP0[0]=OpenP0[1]=OpenP2[0]=OpenP2[1]=1;
            for (kk=0; kk<(SingSpec.NMethods*SingSpec.SecPeriods[seasn]); kk++, srvy++) {   // calculate OpenP;
                if(SingSpec.Data[site][srvy]!=-1) {
                    if(SingSpec.Data[site][srvy]>1) {
                        OpenP0[0] *= ClosedP11[site][srvy][0]*ClosedB[site][srvy][0]; 
                        OpenP0[1] *= ClosedP11[site][srvy][1]*ClosedB[site][srvy][1]; 
                        OpenP2[0]=OpenP2[1]=0; seen=2; 
                    }
                    if(SingSpec.Data[site][srvy]==1) {
                        OpenP0[0] *= ClosedP11[site][srvy][0]*(1.-ClosedB[site][srvy][0]); 
                        OpenP0[1] *= ClosedP11[site][srvy][1]*(1.-ClosedB[site][srvy][1]); 
                        OpenP2[0] *= ClosedP10[site][srvy][0]; 
                        OpenP2[1] *= ClosedP10[site][srvy][1]; 
                        if (seen<1) seen=1;
                    }
                    if(SingSpec.Data[site][srvy]==0) {
                        OpenP0[0] *= 1.0-ClosedP11[site][srvy][0];
                        OpenP0[1] *= 1.0-ClosedP11[site][srvy][1];
                        OpenP2[0] *= (1.-ClosedP10[site][srvy][0]);
                        OpenP2[1] *= (1.-ClosedP10[site][srvy][1]);
                    }
                    if(ClosedP11[site][srvy][0]<0) missng=-1;
                }    //  end of  if not missing
            }   //  end of survey in season loop
            if(missng<0)OpenP0[0]=OpenP0[1]=1;
            prdOcc    *= (pi[seasn]*OpenP0[0]+(1-pi[seasn])*OpenP0[1]);  
            prdNotOcc *= (pi[seasn]*OpenP2[0]+(1-pi[seasn])*OpenP2[1]);
            if (seen>fpflag) prdNotOcc = 0;  //  seen > 1 if false-pos model, seen > 0 if not fp
            prdOcc_save=prdOcc;  prdNotOcc_save=prdNotOcc; // save prds so they can be used below
            if (seasn<(PrmyPers-1) && srvy<lastobs) {  //  Phi(0,0)=psi or 1-eps, Phi(0,1)=1-psi or eps
                prdOcc = prdOcc_save*Phi[site][0][0][seasn+1] + prdNotOcc_save*Phi[site][1][0][seasn+1];
                prdNotOcc = prdOcc_save*Phi[site][0][1][seasn+1] + prdNotOcc_save*Phi[site][1][1][seasn+1];
            }
        } // end PrmyPeriods loop
        cellProb += (prdOcc+prdNotOcc); 
        if (cellProb < SingSpec.nearzero) { nlog0++; cellProb=SingSpec.nearzero; }
        sum += SingSpec.det_hist_frq[site]*log(cellProb); 
		SingSpec.expval[site]=cellProb;
    } // end individual loop
    return(-sum);
}
