#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SingAnal.h"
//#define DBG 1

double OpenModLink(double *Params, int NPars) {
    
    double OpenLike(double ****Phi,  double ***closedp, double ***closedp10, double ***closedb, double *pi);
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int i,site, seasn,srvy, isrvy, imeth,iseasn,kk, pstart, pistrt, gstrt, igrp, ngrps, jj,
	    N=SingSpec.N, T=SingSpec.T, PrmyPeriods=SingSpec.PrmyPeriods;
    
    double temp=0.0, temppsi, psit, eps, gam, *psi, psiseasn, prdq, **qstar, cp,denom, ***closedp10, ***closedb, *pi;
    double psit_lower_limit,psit_upper_limit, x, pnlty=0.;
#ifdef DBG    
    //printf("OpenModLink...(model=%d)\n",SingSpec.Model);
    //printf("params:"); for (kk=0; kk<NPars; kk++) printf(" %f",Params[kk]); printf("\n");
#endif
	double ****Phi, ***closedp; Phi = new double***[N]; closedp = new double **[N];
    for (site=0; site<SingSpec.N; site++) {
        Phi[site] = new double**[2];  closedp[site] = new double*[T];
        for (srvy=0; srvy<T; srvy++) {
            closedp[site][srvy]=new double[2];   //   limited to 2 groups
            closedp[site][srvy][0]=closedp[site][srvy][1]=0;
        }
        for (jj=0; jj<2; jj++) {    
			Phi[site][jj] = new double*[2];
            for (kk=0; kk<2; kk++) Phi[site][jj][kk] = new double[PrmyPeriods];
        }
    }
    pstart=1+PrmyPeriods*2-2; if (SingSpec.Model==4) pstart-=PrmyPeriods-1;
    closedb = new double**[N];    //   prob(detection is 'sure')
    closedp10 = new double**[N];  //   false-pos detection prob.
    qstar = new double*[N];       //   compute qstar = 1 - pstar, where pstar =
    ngrps=1; if (SingSpec.NrowsDM[4]>0) ngrps=2;

	for (site=0; site<N; site++) { 
        qstar[site]=new double[PrmyPeriods]; closedb[site]=new double*[T]; closedp10[site]=new double*[T];
		for (srvy=0; srvy<T; srvy++) {
            closedb[site][srvy]=new double[2]; closedp10[site][srvy]=new double[2];
            closedp10[site][srvy][0] = closedb[site][srvy][0] = closedp10[site][srvy][1] = closedb[site][srvy][1] = 0;
		}
	}
    pi = new double[T]; pistrt=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3];
    for (srvy=0; srvy<PrmyPeriods; srvy++) { //  only allow 2 groups : pi,1-pi
        pi[srvy]=1;
        if (SingSpec.NrowsDM[4]>0) pi[srvy]=getParam(pistrt+srvy,4,srvy,0,0,Params,0,SingSpec.LinkFn);
    }
    for (site=0; site<N; site++) { 
		if (SingSpec.NMethods>1)  
			for (iseasn=0; iseasn<PrmyPeriods; iseasn++)
				for (isrvy=0; isrvy<SingSpec.SecPeriods[iseasn]; isrvy++)
					for (srvy=imeth=0; imeth<SingSpec.NMethods; imeth++,srvy++)
						if (imeth==0) 
							SingSpec.theta[site][srvy][0]=getParam(srvy+1,0,srvy+1,site,srvy,Params,0,SingSpec.LinkFn);
        for (srvy=isrvy=seasn=0,prdq=1.; srvy<T; srvy++) {  //   isrvy = survey within season
           for (igrp=0; igrp<ngrps; igrp++) {
                closedp[site][srvy][igrp] = closedp10[site][srvy][igrp] = closedb[site][srvy][igrp] = cp = 0;
                if (SingSpec.Data[site][srvy]!=-1) { // if not missing data(site site, srvy srvy...
                    gstrt=T*(igrp); 
                    cp = getParam(pstart+gstrt+srvy,3,gstrt+srvy,site,srvy,Params,0,SingSpec.LinkFn);
                    closedp[site][srvy][igrp] = cp;
                    if (SingSpec.FalsePos) {
                        closedp10[site][srvy][igrp]=getParam(pstart+srvy+T,3,srvy+T,site,srvy,Params,0,SingSpec.LinkFn);
                        closedb[site][srvy][igrp]=getParam(pstart+srvy+T+T,3,srvy+T+T,site,srvy,Params,0,SingSpec.LinkFn);
                    }
                }
                prdq*=(1.-cp);    //   variables needed for autologistic model...
                if (igrp==0) 
                    if (++isrvy == SingSpec.SecPeriods[seasn]) { 
                        qstar[site][seasn]=prdq; isrvy=0; prdq=1; seasn++;  // qstar=prd of q's
                    }
            }
        }
    } 
    switch(SingSpec.Model) {
    case 2:    // psi and gamma
        for (site=0; site<N; site++) {
            // calculate psi[1]
            Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            Phi[site][0][1][0] = 1.0-Phi[site][0][0][0]; Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;
            temppsi = Phi[site][0][0][0];
            for (seasn=1; seasn<PrmyPeriods; seasn++) {
                // calculate gamma
                Phi[site][1][0][seasn] = getParam(PrmyPeriods+seasn-1,1,seasn-1,site,seasn-1,Params,0,1);
                Phi[site][1][1][seasn] = 1.0 - Phi[site][1][0][seasn];
                psit = getParam(seasn,0,seasn,site,seasn,Params,0,1);    // calculate psi[t]
                // calculate epsilon, phi[site][0][0][seasn] = epsilon
                Phi[site][0][0][seasn] = (psit-(1.0-temppsi)*Phi[site][1][0][seasn])/temppsi;
                temppsi = psit;
                // check 0<epsilon<1
                if (Phi[site][0][0][seasn]<0.0 ) 
                    temp = (1.0 + Phi[site][0][0][seasn]*Phi[site][0][0][seasn])*1.0e10;
                else {
                    if (Phi[site][0][0][seasn]>1.0) 
                        temp = (1.0 + (Phi[site][0][0][seasn]-1.0)*(Phi[site][0][0][seasn]-1.0))*1.0e10;
                }
                Phi[site][0][1][seasn] = 1.0 - Phi[site][0][0][seasn];
            }
        }  // end site loop
        break;
    case 3:    // psi and epsilon
        for (site=0; site<N; site++) {
            // calculate psi[1]
            Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            Phi[site][0][1][0] = 1.0-Phi[site][0][0][0]; Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;
            temppsi = Phi[site][0][0][0];
            for (seasn=1; seasn<PrmyPeriods; seasn++) {
                // calculate epsilon
                eps = getParam(PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
                if (SingSpec.Alt_Param_Checked) eps=1-eps;  // alt-parm: use persistance instead of extinction
                Phi[site][0][0][seasn] = 1.0 - eps;
                Phi[site][0][1][seasn] = eps;
                psit = getParam(seasn,0,seasn,site,seasn,Params,0,1);//SingSpec.LinkFn);       // calculate psi[t]
                if (temppsi==1.) temppsi=1-1.0e-20;
                psit_lower_limit=temppsi*(1-eps);
                psit_upper_limit=1-temppsi*eps;
                if (psit<psit_lower_limit) psit=psit_lower_limit+1.0e-20;
                if (psit>psit_upper_limit) psit=psit_upper_limit-1.0e-20;
                // calculate gamma using psi(i) and eps(i)
                gam = (psit-temppsi*(1.-eps))/(1.0-temppsi);
                // check 0<gamma<1
                if (gam<0){ pnlty=10000.; } //  temp=(1.0 + gam*gam)*1.0e10; gam=0;}
                if (gam>1){ pnlty=10000.; } //  temp=1.0e20*rand(); gam=1;}
                //printf("seasn=%d gam=%f psit=%f temppsi=%f eps=%f altparmchk=%d\n",seasn,gam,psit,temppsi,eps,SingSpec.Alt_Param_Checked);
                Phi[site][1][0][seasn] = gam;  Phi[site][1][1][seasn] = 1.0 - gam;
                temppsi = psit;
            }
        }  // end site loop
        break;
    case 4:    // psi(.), gamma(t), eps(t)=1-gamma(t)
        for (site=0; site<N; site++) { 
            // calculate psi
            Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            Phi[site][0][1][0] = 1.0-Phi[site][0][0][0]; Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;
            for (seasn=1; seasn<PrmyPeriods; seasn++) {
                // calculate gamma
                Phi[site][1][0][seasn] = getParam(seasn,0,seasn,site,seasn-1,Params,0,1);
                Phi[site][1][1][seasn] = 1.0 - Phi[site][1][0][seasn];
                // calculate epsilon
                Phi[site][0][0][seasn] = Phi[site][1][0][seasn];
                Phi[site][0][1][seasn] = 1.0 - Phi[site][0][0][seasn];
            }
        }  // end site loop
        break;
    case 5: // Standard Parameterization w/ gam (or eps) a function of psi
        //  actually, gam/eps function of total expected number of occupied sites in i-1
        //     which are neighbors of site site (ie., sum(psi(k)*ngbr(k,site), k=1..N)
        psi = new double[N];       //  last subscript for psibar below is always zero in this model.
        for (seasn=0; seasn<=PrmyPeriods; seasn++) for (site=0; site<N; site++) SingSpec.psibar[seasn][site][0]=0;
        seasn=0;
        for (site=0; site<N; site++) {
            // calculate psi
            Phi[site][0][0][seasn] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            Phi[site][0][1][seasn] = 1.0 - Phi[site][0][0][seasn]; 
            Phi[site][1][0][seasn] = Phi[site][1][1][seasn] = 0.0;
            psi[site] = Phi[site][0][0][seasn];
            if (SingSpec.uncondpsi==0 && SingSpec.OpenData[site][seasn]>=0) {  //   use conditional psi (if at least one non-missing data for season)
                if (SingSpec.OpenData[site][seasn]>0) psi[site]=1.;  //  detect at least once in a year
                else {
                    denom=1.0-psi[site]+psi[site]*qstar[site][seasn];
                    if (denom!=0) psi[site]=psi[site]*qstar[site][seasn]/denom;
                }
            }
			//  here's how the autologistic covar is computed:
			//   1)  compute psi(1)   {psi[site]}
			//   2)      psibar(2) = psi(1)*nei(2,1)*w(2)   {psibar[kk=2]=psi[site]*nei(kk,site)*wgt(kk)}
			//   3)      psibar(3) = psi(1)*nei(3,1)*w(3)
			//   4)      psibar(4) = psi(1)*nei(4,1)*w(4)
			//              :          :       :       :
			//   1)  compute psi(2)
			//   2)      psibar(1) = psi(2)*nei(1,2)*w(1)
			//   3)      psibar(3) = psi(2)*nei(3,2)*w(3)
			//   4)      psibar(4) = psi(2)*nei(4,2)*w(4)
			//              :          :       :       :
			//   1)  compute psi(3)
			//   2)      psibar(1) = psi(3)*nei(1,3)*w(1)
			//   3)      psibar(2) = psi(3)*nei(2,3)*w(2)
			//   4)      psibar(4) = psi(3)*nei(4,3)*w(4)
			//              :          :       :       :
            if (SingSpec.UseNeighbor>0) {
                for (kk=0; kk<N; kk++) {
#ifdef neighbor2					
					i=SingSpec.neighbor_lst[kk][site]; if (i<0) break;
					SingSpec.psibar[seasn][kk][0]+=psi[site]*SingSpec.neighbor_wgt[i];
#else
                    if (SingSpec.neighbor[kk][site]>1) break;  //  any digit > 1 indicates no more neighbors for site
                    x=(double)SingSpec.neighbor[kk][site]*SingSpec.neighbor_wgt[kk]; 
                    SingSpec.psibar[seasn][kk][0]+=psi[site]*x;
#endif
                }
			}
        }
        for (seasn=1; seasn<PrmyPeriods; seasn++) { 
            for (site=0; site<N; site++) {
                // calculate epsilon
                eps = getParam(PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1); 
                if (SingSpec.Alt_Param_Checked) eps = 1 - eps;  // alt-parm: use persistance instead of extinction
                Phi[site][0][0][seasn] = 1.0 - eps; 
                Phi[site][0][1][seasn] = eps;
                // calculate gamma
                gam = Phi[site][1][0][seasn] = getParam(seasn,1,seasn-1,site,seasn-1,Params,0,1);
                Phi[site][1][1][seasn] = 1.0 - gam;
                //   calc unconditional psi(site,seasn) and sum to get avg psi(seasn)
                psi[site]=  psi[site]*Phi[site][0][0][seasn]+   //  psi(j-1)*(1-eps(j))+
                      (1.-psi[site])*Phi[site][1][0][seasn];  //  (1-psi(j-1))*gam(j)
                if (SingSpec.uncondpsi==0) {  //   use conditional psi
                    if (SingSpec.OpenData[site][seasn]>0) psi[site]=1.;
                    else if (qstar[site][seasn]>0.) psi[site]=psi[site]*qstar[site][seasn]/(1.-psi[site]+psi[site]*qstar[site][seasn]);
                }
                if (SingSpec.UseNeighbor>0)
                    //  add psi(site) to each of it's (kk) neighbors
                    for (kk=0; kk<N; kk++) {
#ifdef neighbor2					
						i=SingSpec.neighbor_lst[kk][site]; if (i<0) break;
						SingSpec.psibar[seasn][kk][0]+=psi[site]*SingSpec.neighbor_wgt[i];
#else                        
						if (SingSpec.neighbor[kk][site]>1) break;       // last subscript for psibar always zero for 1 species model
                        x=(double)SingSpec.neighbor[kk][site]*SingSpec.neighbor_wgt[kk]; //  kk is focal site, site is neighbr site
                        SingSpec.psibar[seasn][kk][0]+=psi[site]*x;
#endif						
                    }
            }   // end site loop
        }  // end seasn loop
        delete [] psi;
        break;
    default: // Standard Parameterization
        for (site=0; site<N; site++) {
            // calculate psi
            Phi[site][0][0][0] = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            Phi[site][0][1][0] = 1.0-Phi[site][0][0][0]; Phi[site][1][0][0] = Phi[site][1][1][0] = 0.0;
            psiseasn=Phi[site][0][0][0];
            for (seasn=1; seasn<PrmyPeriods; seasn++) {
                // calculate epsilon
                eps = getParam(PrmyPeriods+seasn-1,2,seasn-1,site,seasn-1,Params,0,1);
                if (SingSpec.Alt_Param_Checked) eps=1-eps;  // alt-parm: use persistance instead of extinction
                Phi[site][0][0][seasn] = 1.0 - eps; 
                Phi[site][0][1][seasn] = eps;
                // calculate gamma
                Phi[site][1][0][seasn] = getParam(seasn,1,seasn-1,site,seasn-1,Params,0,1);
                Phi[site][1][1][seasn] = 1.0 - Phi[site][1][0][seasn];
                //   calc psi(t+1) and set 1000'th site covar to this (to use psi(t) as covar in eps or gam)
                psiseasn=psiseasn*Phi[site][0][0][seasn]+(1.-psiseasn)*Phi[site][1][0][seasn]; 
            }
        }  // end site loop
        break;
    } // end switch
    if (temp == 0.0) temp = pnlty+OpenLike(Phi, closedp, closedp10, closedb, pi);
    for (i=0; i<N; i++) {
        for (srvy=0; srvy<T; srvy++) {
            delete [] closedp10[i][srvy]; delete[] closedb[i][srvy]; delete [] closedp[i][srvy]; 
        }
        delete [] closedp10[i]; delete [] closedb[i]; delete [] qstar[i]; delete [] closedp[i];
		for (jj=0; jj<2; jj++) {
			for (kk=0; kk<2; kk++) delete [] Phi[i][jj][kk];
			delete [] Phi[i][jj];
		}
		delete [] Phi[i];
    }
    delete [] qstar; delete [] closedp10; delete [] closedb; delete [] closedp; delete [] pi; delete [] Phi;
    return (temp); 
}

double OpenModLink1(double *Params, int NPars) {   
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int site, srvy, lastobs, seasn, kk, missng, pstart, gstart=1, estart ;
    double sum=0, op, oq, cp, Psi, gam, eps, prdOcc, prdNotOcc, prdOcc_save, prdNotOcc_save, cellProb; 

	if (SingSpec.ifn<1) {
		printf("OpenModLink1 - primPer=%d secPer=",SingSpec.PrmyPeriods);
	    for (site=0; site<SingSpec.PrmyPeriods; site++) printf(" %d",SingSpec.SecPeriods[site]); 
		printf("\n"); 
	}
	estart=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1];
	pstart=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2];	
 // Standard Parameterization
 	//Psi=1/(1+exp(-Params[0])); cp=1/(1+exp(-Params[3])); eps=1/(1+exp(-Params[2])); gam=1/(1+exp(-Params[1]));
    for (site=0; site<SingSpec.N; site++) {  
		Psi = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); 
		prdOcc = Psi; prdNotOcc = 1-Psi;
		for (srvy=0, lastobs=-1; srvy<SingSpec.T; srvy++) if (SingSpec.Data[site][srvy]!=-1) lastobs=srvy+1;
		for (seasn=srvy=0; seasn<SingSpec.PrmyPeriods; seasn++) {            
			op=oq=1;                      // op = open p* = 1-(1-p1)(1-p2)...  oq={(1-p1)(1-p2)... or zero}
			for (kk=missng=0; kk<SingSpec.SecPeriods[seasn]; kk++, srvy++) {   // calculate OpenP;
				if(SingSpec.Data[site][srvy]!=-1) {
					cp = getParam(pstart+srvy,3,srvy,site,srvy,Params,0,SingSpec.LinkFn);
					if(SingSpec.Data[site][srvy]>=1) { op *= cp; oq = 0; }
					if(SingSpec.Data[site][srvy]==0) { op *= 1.0-cp;  }
					if(cp<0) missng=-1;
				}    //  end of  if not missing
			}   //  end of survey in season loop
			if(missng<0)op=1; 
			prdOcc    *= op;  prdNotOcc *= oq; 
			prdOcc_save=prdOcc;  prdNotOcc_save=prdNotOcc;  // save prds so they can be used below in matrix mult.
			if (seasn<(SingSpec.PrmyPeriods-1) && srvy<lastobs) {  //  Phi(0,0)=psi or 1-eps, Phi(0,1)=1-psi or eps
				eps = getParam(estart+seasn,2,seasn,site,seasn,Params,0,1); 
				gam = getParam(gstart+seasn,1,seasn,site,seasn,Params,0,1);  // calculate gamma
				prdOcc    = prdOcc_save*(1-eps) + prdNotOcc_save*gam;
				prdNotOcc = prdOcc_save*eps      + prdNotOcc_save*(1-gam); 
			}
		} // end PrmyPeriods loop
		cellProb = prdOcc+prdNotOcc; if (cellProb < SingSpec.nearzero) cellProb=SingSpec.nearzero;
		//SingSpec.expval[site]=cellProb;
		sum-=SingSpec.det_hist_frq[site]*log(cellProb);
		//printf("%d %25.20e %25.20e\n",site,cellProb,sum);
	} // end individual loop	
	//printf("-like:%25.20e  %d %25.20e %25.20e %25.20e %25.20e\n",sum,SingSpec.ifn,Params[0],Params[1],Params[2],Params[3]);
	//printf("%25.20e %25.20e %25.20e %25.20e\n",Psi,gam,eps,cp);
    return(sum); 
}
