#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SingAnal.h"

double OpenModLinkMM(double *Params, int NPars) {
    
    double OpenLikeMM(double ****Phi, double **theta, double **closedp);
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int i,j,iie,iig,site, seasn,srvy, isrvy, iseasn, jj, kk,
	    N=SingSpec.N, T=SingSpec.T, PrmyPeriods=SingSpec.PrmyPeriods, pstart=SingSpec.NrowsDM[0];
    double temp=0.0, eps, psiseasn,  pnlty=0., ****Phi, **closedp, **theta; 
	
	Phi = new double***[N]; closedp = new double *[N]; theta = new double *[N];
    for (site=0; site<SingSpec.N; site++) {
        Phi[site] = new double**[2];  closedp[site] = new double[T]; theta[site] = new double[T];
        for (srvy=0; srvy<T; srvy++) closedp[site][srvy]=0;
        for (jj=0; jj<2; jj++) {    
			Phi[site][jj] = new double*[2];
            for (kk=0; kk<2; kk++) Phi[site][jj][kk] = new double[PrmyPeriods];
        }
    }
    for (site=0; site<N; site++) {
		if (SingSpec.NMethods>1)  
			for (iseasn=i=0; iseasn<SingSpec.PrmyPeriods; iseasn++) {
				for (isrvy=srvy=0; isrvy<(SingSpec.NMethods*SingSpec.SecPeriods[iseasn]); isrvy++,srvy++)
					if ((srvy % SingSpec.NMethods) == 0) { 
						theta[site][i]=getParam(i+1,0,i+1,site,srvy,Params,0,SingSpec.LinkFn); i++;
					}
				}
        for (srvy=isrvy=seasn=0; srvy<T; srvy++) {    //   isrvy = survey within season                                 
            if (SingSpec.Data[site][srvy]!=-1) { // if not missing data(site site, srvy srvy...
                closedp[site][srvy] = getParam(pstart+srvy,1,srvy,site,srvy,Params,0,SingSpec.LinkFn);
            }
     
        }
    } 
	iie=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]; iig=iie+SingSpec.NrowsDM[2];
    for (site=0; site<N; site++) { 
        // calculate psi
        Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);     //                Occ  [  psi     1-psi ]
        Phi[site][0][1][0] = 1.0-Phi[site][0][0][0];                     // initial PHI:   Unocc[  0         0   ]
		Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;                   //
        psiseasn=Phi[site][0][0][0];                                              //
        for (seasn=1; seasn<PrmyPeriods; seasn++) {
            // calculate epsilon
            eps = getParam(iie+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
            if (SingSpec.Alt_Param_Checked) eps=1-eps;  // alt-parm: use persistance instead of extinction
            Phi[site][0][0][seasn] = 1.0 - eps;                                         //                        Occ     Unocc
            Phi[site][0][1][seasn] = eps;                                               //               Occ  [  1-eps     eps   ]
            // calculate gamma                                                                              Phi:  Unocc[   gam      1-gam ]
            Phi[site][1][0][seasn] = getParam(iig+seasn-1,3,seasn-1,site,seasn-1,Params,0,1); // 
            Phi[site][1][1][seasn] = 1.0 - Phi[site][1][0][seasn];             //
            //   calc psi(t+1) and set 1000'th site covar to this (to use psi(t) as covar in eps or gam)
            psiseasn=psiseasn*Phi[site][0][0][seasn]+(1.-psiseasn)*Phi[site][1][0][seasn]; 
        }
    }  // end site loop

    if (temp == 0.0) temp = pnlty+OpenLikeMM(Phi, theta, closedp);
    
    for (site=0; site<N; site++) {
		for (i=0; i<2; i++) {
			for (j=0; j<2; j++) delete [] Phi[site][i][j]; 
			delete [] Phi[site][i]; 
		}
		delete [] Phi[site]; delete [] closedp[site]; delete [] theta[site];
    }
    delete [] Phi; delete [] closedp; delete [] theta;
    return (temp); 
}
