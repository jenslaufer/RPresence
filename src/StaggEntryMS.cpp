#include <stdlib.h>
#include <math.h>
#include "SingAnal.h"
//#define DEBUG 1
double StaggEntryModLink(double *Params, int NPars) { //  Staggered entry (single and multiseason) model
    
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int i, site, srvy, j, T=SingSpec.T, T2, N=SingSpec.N, igam, ieps, seasn, isrvy, eptr, dptr, esrvy, dsrvy;
    bool seen,dbg=false;
    
    double psi,**D, **cp, **cq, sum=0, prod[3], newprod[3], rnewprod[3], like;
    double seasnvec[2], msprd[2], gam, eps;
    
    rnewprod[0]=rnewprod[1]=rnewprod[2]=msprd[0]=msprd[1]=0; 
    igam=1+T; ieps=igam+SingSpec.PrmyPeriods-1;

    D = new double*[3]; cp = new double*[3]; cq = new double*[3];
    for (i=0; i<3; i++) {
        D[i] = new double[3]; cp[i] = new double[3]; cq[i] = new double[3];
        for (j=0; j<3; j++) D[i][j]=cp[i][j]=cq[i][j]=0;
    }                                        //               not 
                                             //             entered
                                             //               yet  ent depart         
                                             //              [ 1-e  e   0 ]     e=prob(ent in i, | not ent yet)       
                                             //          D=  [  0  1-d  d ]     d=prob(depart in i)
                                             //              [  0   0   1 ]   
                                             //  e(K)=1 (must enter by last occ.), d(0)=0
    for (site=0, sum=0.; site<N; site++) {
        psi = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
        msprd[1]=1-psi;  //  msprd[0]=Pr(site occ) at end of season
		eptr=ieps+SingSpec.PrmyPeriods-1; dptr=eptr+T-1-SingSpec.PrmyPeriods;
        for (seasn=srvy=esrvy=dsrvy=0; seasn<SingSpec.PrmyPeriods; seasn++) { 
            prod[0]=psi; prod[1]=prod[2]=0; T2=SingSpec.SecPeriods[seasn]-1; seen=0;
            
            for (isrvy=0; isrvy<=T2; isrvy++,srvy++) {            
                
                D[0][1] = D[2][2] = 1.0; D[1][2] = 0.0;
                
                if (isrvy<T2) D[0][1] = getParam(eptr++,4,esrvy++,site,srvy,Params,0,SingSpec.LinkFn);// e(t)
                if (isrvy>0)  D[1][2] = getParam(dptr++,5,dsrvy++,site,srvy,Params,0,SingSpec.LinkFn);// d(t)
                
                D[0][0] = 1-D[0][1]; 
                D[1][1] = 1-D[1][2]; 
				cp[1][1] = 0;
                if (SingSpec.Data[site][srvy]!=-1) 
                    cp[1][1] = getParam(srvy+1,1,srvy,site,srvy-1,Params,0,SingSpec.LinkFn);  // p(t)
                cq[0][0] = cq[2][2] = 1; cq[1][1] = 1 - cp[1][1];
                
                         //newprod = prod * D  (matrix multiply vector * 3x3 matrix
                newprod[0]=prod[0]*D[0][0]+prod[1]*D[1][0]+prod[2]*D[2][0];
                newprod[1]=prod[0]*D[0][1]+prod[1]*D[1][1]+prod[2]*D[2][1];
                newprod[2]=prod[0]*D[0][2]+prod[1]*D[1][2]+prod[2]*D[2][2];
                if(SingSpec.Data[site][srvy]>0) { 
                    seen=1;                    //  rnewprod = newprod * cp
                    rnewprod[0]=newprod[0]*cp[0][0]+newprod[1]*cp[1][0]+newprod[2]*cp[2][0];
                    rnewprod[1]=newprod[0]*cp[0][1]+newprod[1]*cp[1][1]+newprod[2]*cp[2][1];
                    rnewprod[2]=newprod[0]*cp[0][2]+newprod[1]*cp[1][2]+newprod[2]*cp[2][2];
                } 
                else {                         //  rnewprod = newprod * (1-cp)
                    rnewprod[0]=newprod[0]*cq[0][0]+newprod[1]*cq[1][0]+newprod[2]*cq[2][0];
                    rnewprod[1]=newprod[0]*cq[0][1]+newprod[1]*cq[1][1]+newprod[2]*cq[2][1];
                    rnewprod[2]=newprod[0]*cq[0][2]+newprod[1]*cq[1][2]+newprod[2]*cq[2][2];
                }
                prod[0]=rnewprod[0]; prod[1]=rnewprod[1]; prod[2]=rnewprod[2];
            }   //  end season loop
            
            seasnvec[0]=rnewprod[0]+rnewprod[1]+rnewprod[2]; seasnvec[1]=msprd[1]; if (seen) seasnvec[1]=0;
            msprd[0]=seasnvec[0]; msprd[1]=seasnvec[1];
            
            if (seasn<(SingSpec.PrmyPeriods-1)) {
                gam = getParam(igam+seasn,2,seasn,site,srvy,Params,0,SingSpec.LinkFn);      //  gam(t)
                eps = getParam(ieps+seasn,3,seasn,site,srvy,Params,0,SingSpec.LinkFn);      //  eps(t)
                msprd[0]=seasnvec[0]*(1-eps)+seasnvec[1]*gam; msprd[1]=seasnvec[0]*eps+seasnvec[1]*(1-gam);
                psi=msprd[0]; 
            }          
        } // end PrmyPeriods loop
        like = msprd[0]+msprd[1];
        if (like < SingSpec.nearzero) like=SingSpec.nearzero;
        sum += SingSpec.det_hist_frq[site]*log(like); 
		SingSpec.expval[site]=like;
        if (dbg) {
            printf("site=%d ",site);
            for (seasn=srvy=0; seasn<SingSpec.PrmyPeriods; seasn++) 
                for (isrvy=0; isrvy<SingSpec.SecPeriods[0]; isrvy++,srvy++) printf("%d",SingSpec.Data[site][srvy]);
            printf(" %d %22.15f %22.15f\n",site,4000*like,sum);
        }
    } // end individual loop  
    
    if(dbg)printf("OpenLike:2*sum=%f %f %f %f %f\n",-sum*2,D[0][0],D[0][1],D[1][0],D[1][1]);
    
    for (j=0; j<3; j++) { delete[] D[j]; delete [] cp[j]; delete [] cq[j]; }
    delete[] D; delete [] cp; delete [] cq;
    if (SingSpec.Verbose>1) printf("\nloglike=%f\r",sum); if (dbg) abort();
    return (-sum);
}
