#include "SingAnal.h"
#include <math.h>
#define MISSNG -1
double CMLike(double *Params, int NPar) {
    //   CMLike - computes the log-likelihood value for the Custom model (single-season)
    //            for current parameter values (Params vector).
    //   
    //    This routine calls 'getParam' function to convert beta params to probabilities
    //    using the appropriate design matrix, and logit transformation.
    //
    //    Multi-method model likelihood can also be computed using this routine.
    //
    int ii, jj, kk, nmeth, inarea, seen, sure, allmiss, frstp,isurvey, T, N; 
    double psi=0, sum, sumll, theta, prd, like, p11, p11a, p10=0, b=0, prd2, prda, sumdiffpsi, condpsi, eps=.1e-10,pi=0;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec; T=SingSpec.T; N=SingSpec.N;
    nmeth=SingSpec.NMethods; frstp=1;
	if (nmeth>1) frstp=T/SingSpec.NMethods+1;
    for (ii=0; ii<N; ii++) { 
        SingSpec.psibar[0][ii][0]=0.5;   // initial values for spatial covariate
        SingSpec.psibar[0][ii][1]=0;
    }
    if (SingSpec.NrowsDM[2]>0) pi=getParam(T+frstp+T,2,0,0,0,Params,0,SingSpec.LnkFn[2][0]);      // pi
    for (kk=0; kk<100; kk++) { sumdiffpsi=0;   //  spatial dependence loop... 
        //  loops until psi's become stable
        for (ii=0,sumll=0; ii<N; ii++) {                                //    loop for each site (ii=0 to N)...
            psi=condpsi=getParam(0,0,0,ii,0,Params,0,SingSpec.LinkFn);                 //   get logit-transformed psi param
            like=psi; inarea=seen=sure=false; allmiss=true; theta=prd=prd2=prda=1;          
            for (jj=0; jj<T; jj++) {                                  //   loop for each survey (jj=0 to T)...
                p11=p11a=getParam(jj+frstp,1,jj,ii,jj,Params,0,SingSpec.LinkFn);      // p11=Pr(det | occupied)
                if (SingSpec.NrowsDM[2]>0) {
                    p11a=getParam(jj+frstp+T,1,jj+T,ii,jj,Params,0,SingSpec.LinkFn);  // p11a= det.prob. grp2
                }         
                if (SingSpec.FalsePos) {
                    p10=getParam(jj+frstp+T,1,jj+T,ii,jj,Params,0,SingSpec.LinkFn);   // p10=pr(det | unoccupied)
                    b=getParam(jj+frstp+T+T,1,jj+T+T,ii,jj,Params,0,SingSpec.LinkFn); //  b = pr(det is certain | detection)
                }
                if (SingSpec.Data[ii][jj]!=MISSNG  && p11>=0) {                 //   if not missing data...
                    allmiss=false;       
                    if (SingSpec.Data[ii][jj]>1)  { prd *= p11*b;      prd2=0;    prda*=p11a*b;    seen=inarea=sure=true;} 
                    if (SingSpec.Data[ii][jj]==1) { prd *= p11*(1.-b); prd2*=p10; prda*=p11a*(1.-b); seen=inarea=true;}  //  mult prd by p if detected,
                    if (SingSpec.Data[ii][jj]==0) { prd *= (1.0-p11);  prd2*=(1.-p10); prda*=(1-p11a);} //  mult prd by (1-p) if not detected,
                }
                if (nmeth>1) {
                    if ((jj+1) % nmeth == 0) {                                  //  if last method of survey...
                        if ( !allmiss ) {                                       //   if not all missing data in survey...
                            isurvey=(int)(jj/nmeth)+1;                          //   get logit-transformed theta param
                            theta=getParam(isurvey,0,isurvey,ii,jj,Params,0,SingSpec.LinkFn);  
                            prd*=theta; prda*=theta;                             //  mul prd by theta
                            if ( !seen ) { prd+=(1-theta); prda+=(1-theta);}     //  if never seen, add 1-theta to prd
                        }
                        like*=(1-pi)*prd+pi*prda; prda=prd=1; seen=false; // seen means seen by >= 1 method, inarea means seen >=1 in season
                    }
                }
            }  //  end jj loop
            if (nmeth<=1) like*=((1-pi)*prd+pi*prda); //  mult like by prd for not-multi-method (multi-method done above)
            if ( !inarea || !sure) like+=(1-psi)*prd2;     // if not in area (never seen), add 1-psi to likelihood
            //     
            if (seen) condpsi=1;
            else 
                if (!allmiss) condpsi=prd*psi/((1-psi)*(1-prd)+eps);
            if (SingSpec.UseNeighbor>0)
                for (jj=0; jj<N; jj++) {
                    if (SingSpec.neighbor[jj][ii]>1) break;  //  any digit > 1 indicates no more neighbors for site
                    SingSpec.psibar[0][jj][1]+=condpsi*(double)SingSpec.neighbor[jj][ii];
                }
            //
            SingSpec.expval[ii]=like; 
            if (like<SingSpec.nearzero) like=SingSpec.nearzero;  //  like <= 0, limit amt added to log-likelihood
            sumll -= SingSpec.det_hist_frq[ii]*log(like);        //  add -(# sites obs with this history)*log(like) to sum
        } // end ii loop
        
        if (SingSpec.UseNeighbor>0)         //   compute avg psi for neighbors of each site...
            for (ii=0,sumdiffpsi=0; ii<N; ii++) {
                for (jj=0,sum=0; jj<N; jj++) {
                    if (SingSpec.neighbor[ii][jj]>1) break;
                    sum+=SingSpec.neighbor[ii][jj];
                }
                SingSpec.psibar[0][ii][1]/=sum;
                sumdiffpsi+=fabs(SingSpec.psibar[0][ii][1]-SingSpec.psibar[0][ii][0]);  
                SingSpec.psibar[0][ii][0]=SingSpec.psibar[0][ii][1];
                SingSpec.psibar[0][ii][1]=0;
            }
        if (sumdiffpsi<.0001) break;
    }   //  end kk loop
    //fprintf(SingSpec.g,"-logLike=%f parms:%f %f %f %f pi=%f nmeth=%d npar=%d\n",sumll,Params[0],Params[1],Params[2],Params[3],pi,nmeth,NPar); 
    //for (ii=0; ii<N; ii++) fprintf(g,"%f\n",SingSpec.expval[ii]);
    //if (SingSpec.Verbose>1) printf("xll=%lf psi=%lf %30.22lf\n",sumll,psi,SingSpec.nearzero);
    return(sumll);
}
