#include "SingAnal.h"
//****  Likelihood for joint modelling of 2 species from presence/absence data  ****//
#include <math.h>
#include <stdio.h>

double TwoSpCorrDetLink(double *Params, int NPar) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    int i, j, survey, site, N=SingSpec.N, T=SingSpec.T, tA, tB, seenA, seenB, pstrt=3+6*T;
    double sum = 0, tmp=0, sum2=0, pA, pB, rA, psiA, likeA[3], likeB[3], likeAB[4], v[4], tr[4][4],  psiBA, psiBa;
    double rBA, rBa, thA, thBA, thBa, thA2, thBA2, thBa2, sum3=0, sum4=0, th0piA, th0piBA, th0piBa;
    for (i=0; i<4; i++) v[i]=0;
    
    for (site=0; site<N; site++) {  // loop through individual sites - calculate psi[A] and psi[B]      
        psiA=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
        psiBA=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
        psiBa=getParam(2,0,2,site,0,Params,0,SingSpec.LinkFn); 
        thA=getParam(3,0,3,site,0,Params,0,SingSpec.LinkFn);
        thA2=getParam(3+T,0,3+T,site,0,Params,0,SingSpec.LinkFn);
        thBA=getParam(3+2*T,0,3+2*T,site,0,Params,0,SingSpec.LinkFn);
        thBA2=getParam(3+3*T,0,3+3*T,site,0,Params,0,SingSpec.LinkFn);
        thBa=getParam(3+4*T,0,3+4*T,site,0,Params,0,SingSpec.LinkFn);
        thBa2=getParam(3+5*T,0,3+5*T,site,0,Params,0,SingSpec.LinkFn);
        th0piA=getParam(3+11*T,2,0,site,0,Params,0,SingSpec.LinkFn);
        th0piBA=getParam(4+11*T,2,1,site,0,Params,0,SingSpec.LinkFn);
        th0piBa=getParam(5+11*T,2,2,site,0,Params,0,SingSpec.LinkFn);
        thA=th0piA*thA2+(1-th0piA)*thA;
        thBA=th0piBA*thBA2+(1-th0piBA)*thBA;
        thBa=th0piBa*thBa2+(1-th0piBa)*thBa;
        likeAB[1]=psiA*psiBA*thA*(1-thBA);  // current local occ. state
        likeAB[2]=psiA*psiBA*(1-thA)*thBa; 
        likeAB[3]=psiA*psiBA*thA*thBA;
        likeAB[0]=psiA*psiBA*(1-thA)*(1-thBa);
        likeA[1]=psiA*(1-psiBA)*thA; likeA[0]=psiA*(1-psiBA)*(1-thA);
        likeB[1]=(1-psiA)*psiBa*thBa; likeB[0]=(1-psiA)*psiBa*(1-thBa);
        seenA=seenB=0; 
        //calculate pA's and pB's
        for (survey=0; survey<T; survey++) {   // set up 4 x 4 transition matrix - tr
            thA=getParam(survey+3,0,survey+3,site,survey,Params,0,SingSpec.LinkFn);
            thA2=getParam(survey+3+T,0,survey+3+T,site,survey,Params,0,SingSpec.LinkFn);
            thBA=getParam(survey+3+2*T,0,survey+3+2*T,site,survey,Params,0,SingSpec.LinkFn);
            thBA2=getParam(survey+3+3*T,0,survey+3+3*T,site,survey,Params,0,SingSpec.LinkFn);
            thBa=getParam(survey+3+4*T,0,survey+3+4*T,site,survey,Params,0,SingSpec.LinkFn);
            thBa2=getParam(survey+3+5*T,0,survey+3+5*T,site,survey,Params,0,SingSpec.LinkFn);
            if (SingSpec.Data[site][survey]!=-1) {
                pA=getParam(pstrt+survey,1,survey,site,survey,Params,0,SingSpec.LinkFn);
                pB=getParam(pstrt+survey+T,1,survey+T,site,survey,Params,0,SingSpec.LinkFn);
                rA=getParam(pstrt+survey+2*T,1,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
                rBA=getParam(pstrt+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn); 
                rBa=getParam(pstrt+survey+4*T,1,survey+4*T,site,survey,Params,0,SingSpec.LinkFn);
                tA=tB=0;
                if (SingSpec.Data[site][survey]==1 || SingSpec.Data[site][survey]==3) tA = 1;
                if (SingSpec.Data[site][survey]==2 || SingSpec.Data[site][survey]==3) tB = 1;
                if (SingSpec.Nstates==2) if (SingSpec.Data[site+N][survey]==1) tB = 1;
                switch (tA+2*tB) {
                case 0:     // neither detected
                    likeAB[3] *= (1-rA)*(1-rBa);    // [3] -> both sp locally present
                    likeAB[2] *= (1-pB);             // [2] -> sp B only locally present
                    likeAB[1] *= (1-pA);               // [1] -> sp A only locally present
                    likeA[1] *= (1-pA);               // [1] -> sp A only locally present
                    likeB[1] *= (1-pB);               // [1] -> sp B only locally present
                    break;
                case 1: 
                    likeAB[3] *= rA*(1-rBA);
                    likeAB[2] = likeAB[0] = 0; 
                    likeAB[1] *= pA;    
                    likeA[1] *= pA; likeA[0]=0;              // [1] -> sp A only locally present
                    likeB[1] = likeB[0]=0;               // [1] -> sp B only locally present
                    seenA = 1; break;   // A detected, not B
                case 2:
                    likeAB[3] *= (1-rA)*rBa;
                    likeAB[2] *= pB;    
                    likeAB[1] = likeAB[0] = 0;    
                    likeA[1] = likeA[0]=0;               // [1] -> sp A only locally present
                    likeB[1] *= pB; likeB[0]=0;              // [1] -> sp B only locally present
                    seenB = 1; break;   // B detected, not A
                case 3: 
                    likeAB[3] *= rA*rBA;    
                    likeAB[2] = likeAB[1] = likeAB[0] = 0;    
                    likeA[1]=likeB[1]=likeA[0]=likeB[0]=0;
                    seenA = seenB = 1;   // both A and B detected
                    break;    // A and B detected
                }   // end switch
            } // end if missing
            if (survey<(T-1)) {
                tr[0][0]=(1-thA)*(1-thBa);   tr[0][1]=thA*(1-thBA);   tr[0][2]=(1-thA)*thBa;   tr[0][3]=thA*thBA;
                tr[1][0]=(1-thA2)*(1-thBa);  tr[1][1]=thA2*(1-thBA);  tr[1][2]=(1-thA2)*thBa;  tr[1][3]=thA2*thBA;
                tr[2][0]=(1-thA)*(1-thBa2);  tr[2][1]=thA*(1-thBA2);  tr[2][2]=(1-thA)*thBa2;  tr[2][3]=thA*thBA2;
                tr[3][0]=(1-thA2)*(1-thBa2); tr[3][1]=thA2*(1-thBA2); tr[3][2]=(1-thA2)*thBa2; tr[3][3]=thA2*thBA2;
                //   matrix multply current state (likeAB) * transition matrix.
                for (i=0; i<4; i++) for (j=0; j<4; j++) v[i]+=likeAB[j]*tr[j][i];
                //  store result in likeAB
                for (i=0; i<4; i++) { likeAB[i]=v[i]; v[i]=0; }
                //  matrix mult current state (likeA) * trans matrix
                tmp=likeA[0]*(1-thA)+likeA[1]*(1-thA2); likeA[1]=likeA[0]*thA+likeA[1]*thA2; likeA[0]=tmp;
                //  matrix mult current state (likeB) * trans matrix
                tmp=likeB[0]*(1-thBa)+likeB[1]*(1-thBa2); likeB[1]=likeB[0]*thBa+likeB[1]*thBa2; likeB[0]=tmp;
            }
        } // end survey=0..T
        if (seenA>0) likeB[1]=likeB[0]=0;
        if (seenB>0) likeA[1]=likeA[0]=0;
        tmp = likeAB[3] + likeAB[2] + likeAB[1] + likeAB[0] +
              likeA[1] + likeA[0] + likeB[1] + likeB[0];
        if (seenA==0 && seenB==0) tmp+= (1-psiA)*(1-psiBa);
        if (tmp < SingSpec.nearzero) tmp = SingSpec.nearzero;
        sum-=SingSpec.det_hist_frq[site]*log(tmp);
        SingSpec.expval[site]=tmp; sum2+=SingSpec.det_hist_frq[site];
    } // end individual site loop
    for (site=0; site<N; site++) {
        sum4+=SingSpec.expval[site];
        SingSpec.expval[site]*=sum2; 
        sum3+=(SingSpec.det_hist_frq[site]-SingSpec.expval[site])*(SingSpec.det_hist_frq[site]-SingSpec.expval[site]);
    }
    //printf("loglik=%f chisq=%f sumprb=%f\n",sum,sum3,sum4);
    return(sum);
}
