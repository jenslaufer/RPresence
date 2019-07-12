#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SingAnal.h"
//#define DBG 1
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif

double OpenModLinkSD(double *Params, int NPars) {   
    double OpenLikeSD(double ****Phi, double ***closedp, double ***theta, double *th0pi, double pi);
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int i,j,site=0, seasn,srvy=0, isrvy, pstart, th0pistrt, pistrt, othrsite, pDM, jj, kk, 
	    N=SingSpec.N, T=SingSpec.T, PrmyPeriods=SingSpec.PrmyPeriods;
    double temp=0.0, eps=0, pnlty=0.,pi=-1.,gam=0,  denom, x, qstar1, qstar2,*psi,**qstar,*th0pi,****Phi, ***closedp, ***theta;
#ifdef DBG
    printf("OpenModLinkSD...(model=%d)\n",SingSpec.Model); 
    printf("params:"); for (i=0; i<NPars; i++) printf(" %f",Params[i]); printf("\n");
#endif
	psi = new double[N]; qstar = new double*[N]; th0pi=new double[N];
	Phi = new double***[N]; closedp = new double **[N]; theta = new double **[N];
    for (site=0; site<SingSpec.N; site++) {
        Phi[site] = new double**[2];  theta[site] = new double*[T]; 
		qstar[site]= new double[SingSpec.PrmyPeriods]; closedp[site] = new double*[T]; 
		for (jj=0; jj<2; jj++) {
			Phi[site][jj] = new double*[2]; 
			for (kk=0; kk<2; kk++) Phi[site][jj][kk]=new double[T];
		}
        for (srvy=0; srvy<T; srvy++) { 
			closedp[site][srvy]=new double[2]; theta[site][srvy]=new double[2];
			for (jj=0; jj<2; jj++) closedp[site][srvy][jj]=theta[site][srvy][jj]=0;
        }
    }
    pstart=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]; pDM=3; 
	if (SingSpec.NrowsDM[2]>0) 
		if (SingSpec.Realname[2][0][0]=='t') { pstart=SingSpec.NrowsDM[0]; pDM=1;}
	th0pistrt=pstart+SingSpec.NrowsDM[pDM]; 
	pistrt=th0pistrt+SingSpec.NrowsDM[pDM+2]; 
    if (SingSpec.NrowsDM[pDM+2]>0) pi=getParam(pistrt,pDM+2,0,site,srvy,Params,0,SingSpec.LinkFn);
    for (site=0; site<N; site++) { 
        th0pi[site]=0;  //  prob(detect at least once in a year)
        if (SingSpec.NrowsDM[pDM+1]>0) th0pi[site]=getParam(th0pistrt,pDM+1,0,site,0,Params,0,SingSpec.LinkFn);
		qstar1=qstar2=1;
        for (srvy=isrvy=seasn=0; srvy<T; srvy++) {  //   isrvy = survey within season
            closedp[site][srvy][0] =  closedp[site][srvy][1] = 0;
            theta[site][srvy][0]=getParam(1+srvy,0,1+srvy,site,srvy,Params,0,SingSpec.LinkFn);
            theta[site][srvy][1]=getParam(T+1+srvy,0,T+1+srvy,site,srvy,Params,0,SingSpec.LinkFn);
            if (SingSpec.Data[site][srvy]!=-1) { // if not missing data(site site, srvy srvy...
                closedp[site][srvy][0] = getParam(pstart+srvy,pDM,srvy,site,srvy,Params,0,SingSpec.LinkFn);
                if (pi>=0)
                    closedp[site][srvy][1] = getParam(pstart+srvy+T,pDM,srvy+T,site,srvy,Params,0,SingSpec.LinkFn);
            }
			qstar1*=(1-theta[site][srvy][0]+theta[site][srvy][0]*(1-closedp[site][srvy][0]));
			qstar2*=(1-theta[site][srvy][srvy==0]+theta[site][srvy][srvy==0]*(1-closedp[site][srvy][0]));
            if ((++isrvy == SingSpec.SecPeriods[seasn]) && srvy<(T-1)) {
                qstar[site][seasn]=(1-theta[site][seasn][2])*qstar1+theta[site][seasn][2]*qstar2; // (1-th0pi,th0pi) * q*
				isrvy=0; theta[site][++seasn][2]=0;    //   this is really th0pi(site,season)
                if (SingSpec.NrowsDM[pDM+1]>0) {
                    theta[site][seasn][2]=getParam(th0pistrt+seasn,pDM+1,seasn,site,seasn,Params,0,SingSpec.LinkFn);
				}
            }
        }
    } 
#ifdef DBG
    printf("closedp=%f %f pi=%f\n",closedp[0][0][0],closedp[0][0][1],pi); 
    printf("theta=%f %f %f\n",theta[0][0][0],theta[0][0][1],theta[0][0][2]); 
#endif
        //  last subscript for psibar below is always zero in this model.
    for (seasn=0; seasn<=PrmyPeriods; seasn++) for (site=0; site<N; site++) SingSpec.psibar[seasn][site][0]=0;
    for (site=seasn=0; site<N; site++) {
        // transition matrix, Phi =  [  psi    1-psi ]  for t=0   [  1-eps   eps ]  for t > 0
        //                           [  0        0   ]            [    gam  1-gam] 
        Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
        Phi[site][0][1][0] = 1.0-Phi[site][0][0][0]; 
        Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;
        psi[site] = Phi[site][0][0][seasn];
        if (SingSpec.uncondpsi==0 && SingSpec.OpenData[site][seasn]>=0) {  //   use conditional psi (if at least one non-missing data for season)
            if (SingSpec.OpenData[site][seasn]>0) psi[site]=1.;  //  detect at least once in a year
            else {
                denom=1.0-psi[site]+psi[site]*qstar[site][seasn];
                if (denom!=0) psi[site]=psi[site]*qstar[site][seasn]/denom;
            }
        }
        if (SingSpec.UseNeighbor>0)
            for (othrsite=0; othrsite<N; othrsite++) {
#ifdef neighbor2					
				i=SingSpec.neighbor_lst[othrsite][site]; if (i<0) break;
				SingSpec.psibar[seasn][othrsite][0]+=psi[site]*SingSpec.neighbor_wgt[i];
#else				
                if (SingSpec.neighbor[othrsite][site]>1) break;  //  any digit > 1 indicates no more neighbors for site
                x=(double)SingSpec.neighbor[othrsite][site]*SingSpec.neighbor_wgt[othrsite];
                SingSpec.psibar[seasn][othrsite][0]+=psi[site]*x; 
#endif
            }
	}
    for (seasn=1; seasn<PrmyPeriods; seasn++) {
		for (site=0; site<N; site++) {
            eps = getParam(T*2+PrmyPeriods-1+seasn,2,seasn-1,site,seasn-1,Params,0,1);
            if (SingSpec.Alt_Param_Checked) eps=1-eps;  // alt-parm: use persistance instead of extinction
            Phi[site][0][0][seasn] = 1.0 - eps; Phi[site][0][1][seasn] = eps;
            gam=Phi[site][1][0][seasn] = getParam(T*2+seasn,1,seasn-1,site,seasn-1,Params,0,1);  //  gamma
            Phi[site][1][1][seasn] = 1.0 - gam                   ;             //  1-gamma
            //   calc unconditional psi(site,seasn) and sum to get avg psi(seasn)
            psi[site]=  psi[site]*Phi[site][0][0][seasn]+   //  psi(j-1)*(1-eps(j))+
                  (1.-psi[site])*Phi[site][1][0][seasn];  //  (1-psi(j-1))*gam(j)
            if (SingSpec.uncondpsi==0) {  //   use conditional psi
                if (SingSpec.OpenData[site][seasn]>0) psi[site]=1.;
                else if (qstar[site][seasn]>0.) psi[site]=psi[site]*qstar[site][seasn]/(1.-psi[site]+psi[site]*qstar[site][seasn]);
            }
            if (SingSpec.UseNeighbor>0)
                //  add psi(site) to each of it's (othrsite) neighbors
                for (othrsite=0; othrsite<N; othrsite++) {
#ifdef neighbor2					
					i=SingSpec.neighbor_lst[othrsite][site]; if (i<0) break;
					SingSpec.psibar[seasn][othrsite][0]+=psi[site]*SingSpec.neighbor_wgt[i];
#else  					
                    if (SingSpec.neighbor[othrsite][site]>1) break;       // last subscript for psibar always zero for 1 species model
                    x=(double)SingSpec.neighbor[othrsite][site]*SingSpec.neighbor_wgt[othrsite]; //  othrsite is focal site, site is neighbr site
                    SingSpec.psibar[seasn][othrsite][0]+=psi[site]*x;
#endif
                }
        }
    }  // end site loop
	 for (i=0; i<N; i++)  delete [] qstar[i]; 
	 delete [] qstar; delete [] psi;
#ifdef DBG
    printf("calling OpenLikeSD...\n");
#endif
    if (temp == 0.0) temp = pnlty+OpenLikeSD(Phi, closedp, theta, th0pi, pi);  
    for (site=0; site<N; site++) {
		for (i=0; i<2; i++) {
			for (j=0; j<2; j++) delete [] Phi[site][i][j]; 
			delete [] Phi[site][i]; 
		}
		for (i=0; i<T; i++) { delete [] closedp[site][i]; delete [] theta[site][i]; }
		delete [] Phi[site]; delete [] closedp[site]; delete [] theta[site];
    }
    delete [] Phi; delete [] closedp; delete [] theta; delete [] th0pi;
    return (temp); 
}
