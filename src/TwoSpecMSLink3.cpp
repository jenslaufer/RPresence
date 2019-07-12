#include "SingAnal.h"
//****  Likelihood for joint modelling of 2 species from presence/absence data  ****//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double TwoSpeciesMSLink(double *Params, int NPar) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    double invlogit(double x);
    int i, j, kk, tA, tB, seenA, seenB, seasn, srvy, survey, site, PHIPARMOPT, GAMBAAdifferent, iroweps1, irow; 
    int PrmyPeriods=SingSpec.PrmyPeriods, N=SingSpec.N, T=SingSpec.T, ip, Intervals;
    double sum = 0.0, tmp=0, pA=0, pB=0, rA=0, rB=0, rAB=0, like0, likeA, likeB, likeAB, psiBA=0, psiBa=0, rBA=0, rBa=0;
    double lower, upper, phi, delta, alf, gam, prVec[4]={0,0,0,0}, xphi[4][4], tmpVec[4];
    double **psiA,**psiB,**psiAB,**gamAb,**gamAB,**gamBaA,**gamBaa,**gamBAA,**gamBAa,
                                 **epsAb,**epsAB,**epsBaA,**epsBaa,**epsBAA,**epsBAa;
    double epsA,epsB,gamA,gamB,x,xsum,pnlty=0., rAwoB,rBwoA,a,b,c,K,psiAwoB,psiBwoA;
    
    Intervals=PrmyPeriods-1;
    //printf("TwoSpeciesMSLink altparm=%d\n",SingSpec.Alt_Param_Checked);
    psiA=new double*[N]; psiB=new double*[N]; psiAB=new double*[N]; 
    gamAB=new double*[N]; gamAb=new double*[N]; 
    gamBAA=new double*[N]; gamBAa=new double*[N]; gamBaA=new double*[N]; gamBaa=new double*[N];
    epsAB=new double*[N]; epsAb=new double*[N]; 
    epsBAA=new double*[N]; epsBAa=new double*[N]; epsBaA=new double*[N]; epsBaa=new double*[N];
    for (seasn=0; seasn<PrmyPeriods; seasn++) 
        for (kk=0; kk<N; kk++) SingSpec.psibar[seasn][kk][0]=SingSpec.psibar[seasn][kk][1]=0;
        //  if param name is 'phi' instead of 'eps'...
    PHIPARMOPT=(SingSpec.Realname[2][0][0]=='p');
    if (SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]<2)
       PHIPARMOPT=(SingSpec.Realname[0][3+6*Intervals][0]=='p' && SingSpec.Realname[0][3+6*Intervals][1]=='h');
    // 13[AB] = psi[A]*psi[B]*phi; phi is parameter being estimated
    // alt parameterization: Psi[2] = psi[Ba] = Pr(occ by B|not occ by A)
    for (site=0; site<N; site++) {  // loop through individual sites
        psiA[site]=new double[PrmyPeriods];  psiB[site]=new double[PrmyPeriods]; 
        psiAB[site]=new double[PrmyPeriods];
        gamAB[site]=new double[PrmyPeriods];  gamAb[site]=new double[PrmyPeriods];
        gamBAA[site]=new double[PrmyPeriods];  gamBaA[site]=new double[PrmyPeriods];
        gamBAa[site]=new double[PrmyPeriods];  gamBaa[site]=new double[PrmyPeriods];
        epsAB[site]=new double[PrmyPeriods];  epsAb[site]=new double[PrmyPeriods];
        epsBAA[site]=new double[PrmyPeriods];  epsBaA[site]=new double[PrmyPeriods];
        epsBAa[site]=new double[PrmyPeriods];  epsBaa[site]=new double[PrmyPeriods];
        // calculate psi[A] and psi[B]
        psiA[site][0]=getParam(0,0,0,site,0,Params,0,0);
        if (SingSpec.Alt_Param_Checked==PSIBA_ODDS_RATIO) {  // newer parameterzation - psiA, psiBa, odds-ratio(psiBa:psiBA) 
            alf=getParam(1,0,1,site,0,Params,0,IDLnk);                 // (Alt_Param_Checked=2)
            gam=getParam(2,0,2,site,0,Params,0,IDLnk);
            psiBA=invlogit(alf); 
            psiBa=invlogit(alf+gam); 
            psiB[site][0]=psiA[site][0]*psiBA+(1-psiA[site][0])*psiBa; 
            psiAB[site][0]=psiA[site][0]*psiBA; 
        }
        if (SingSpec.Alt_Param_Checked==PSIAB_ODDS_RATIO) {  // newest parameterzation - psiAb, psiBa, odds-ratio(psiBa:psiBA) 
            a=fmax(-50,fmin(50,getParam(0,0,0,site,0,Params,0,IDLnk)));  //  a=logit(psiA|notB)
            b=fmax(-50,fmin(50,getParam(1,0,1,site,0,Params,0,IDLnk)));  //  b=logit(psiB|notA)
            c=fmax(-50,fmin(50,getParam(2,0,2,site,0,Params,0,IDLnk)));  //  c=log(tau)
            K=1+exp(a+b+c)+exp(a)+exp(b);
            psiAB[site][0]=exp(a+b+c)/K; psiAwoB=exp(a)/K; psiBwoA=exp(b)/K;   
            psiA[site][0]=psiAB[site][0]+psiAwoB; psiB[site][0]=psiAB[site][0]+psiBwoA;
        }
        if (SingSpec.Alt_Param_Checked==PSIBA_PSIBa) {  // new parameterzation - psiA, psiBa, psiBA(Alt_Param_Checked=1)
            psiBA=getParam(1,0,1,site,0,Params,0,0);
            psiBa=getParam(2,0,2,site,0,Params,0,0); 
            psiAB[site][0]=psiA[site][0]*psiBA;                //  Pr(occ by B) = Pr(occ by A)*pr(Occ by B | A present)
            psiB[site][0]=psiA[site][0]*psiBA+(1-psiA[site][0])*psiBa;  //    + Pr(not occ by A)*Pr(occ by B | A not present)
        }
        if (SingSpec.Alt_Param_Checked==PSIB_PHI) {  //  orig parameterization - psiA,psiB,phi (Alt_Param_Checked=0)
            psiB[site][0]=getParam(1,0,1,site,0,Params,0,0); 
            psiAB[site][0] = psiA[site][0]*psiB[site][0];
            phi=getParam(2,0,2,site,0,Params,1,expLnk);
            psiAB[site][0] = psiA[site][0]*psiB[site][0]*phi;
            // check whether value is in allowable bounds
            if (psiAB[site][0]<0.) { psiAB[site][0]=0.; pnlty=99999.;}
            if (psiAB[site][0]>1.) { psiAB[site][0]=1.; pnlty=99999.;}
        }
        if (SingSpec.UseNeighbor>0)
            for (kk=0; kk<N; kk++) {
#ifdef neighbor2					
				i=SingSpec.neighbor_lst[kk][site]; if (i<0) break;
				SingSpec.psibar[0][kk][0]+=psiA[site][0]*SingSpec.neighbor_wgt[i];
				SingSpec.psibar[0][kk][1]+=psiB[site][0]*SingSpec.neighbor_wgt[i];
#else				
                if (SingSpec.neighbor[kk][site]>1) break;
                x=(double)SingSpec.neighbor[kk][site]*SingSpec.neighbor_wgt[site];
                SingSpec.psibar[0][kk][0]+=psiA[site][0]*x;
                SingSpec.psibar[0][kk][1]+=psiB[site][0]*x;
#endif
            }
    }
    //  variables for option of gamBAA and gamBAa different in design matrix
    GAMBAAdifferent=(SingSpec.NrowsDM[1]>(Intervals*4+1));
    iroweps1=Intervals*4+3;
    if (GAMBAAdifferent) iroweps1=Intervals*6+3;
	//printf("altparmchecked=%d, gambaadifferent=%d\n",SingSpec.Alt_Param_Checked,GAMBAAdifferent);

    for (seasn=0; seasn<Intervals; seasn++) {   //  loop by season so we can get avg psi over sites
        for (site=0; site<N; site++) {  // loop through individual sites
            gamAB[site][seasn] =getParam(3+seasn,                 1,seasn,                site,seasn,Params,0,0);
            gamAb[site][seasn] =getParam(3+seasn+Intervals      , 1,seasn+Intervals,site,seasn,Params,0,0);
            epsAB[site][seasn] =getParam(iroweps1+seasn         , 2,seasn,                site,seasn,Params,0,0);
            epsAb[site][seasn] =getParam(iroweps1+seasn+Intervals, 2,seasn+Intervals,site,seasn,Params,0,0);
            if (PHIPARMOPT) {   //  if param name is 'phi' instead of 'eps'...
                epsAB[site][seasn]=1-epsAB[site][seasn];
                epsAb[site][seasn]=1-epsAb[site][seasn];
            }
            if (SingSpec.Alt_Param_Checked==PSIBA_ODDS_RATIO) {   //  newer parameterzation - psiA, psiBA, odds-ratio(psiBa:psiBA) 
                alf=getParam(3+seasn+2*Intervals, 1,seasn+2*Intervals,site,seasn,Params,0,IDLnk);
                gam=getParam(3+seasn+3*Intervals, 1,seasn+3*Intervals,site,seasn,Params,0,IDLnk);
                gamBAA[site][seasn]=invlogit(alf);
                gamBAa[site][seasn]=invlogit(alf+gam);
                alf=getParam(3+seasn+4*Intervals, 1,seasn+4*Intervals,site,seasn,Params,0,IDLnk);
                gam=getParam(3+seasn+5*Intervals, 1,seasn+5*Intervals,site,seasn,Params,0,IDLnk);
                gamBaA[site][seasn]=invlogit(alf);
                gamBaa[site][seasn]=invlogit(alf+gam);
                alf=getParam(3+seasn+8*Intervals, 2,seasn+2*Intervals,site,seasn,Params,0,IDLnk);
                gam=getParam(3+seasn+9*Intervals, 2,seasn+3*Intervals,site,seasn,Params,0,IDLnk);
                epsBAA[site][seasn]=invlogit(alf);
                epsBAa[site][seasn]=invlogit(alf+gam);
                alf=getParam(3+seasn+10*Intervals, 2,seasn+4*Intervals,site,seasn,Params,0,IDLnk);
                gam=getParam(3+seasn+11*Intervals, 2,seasn+5*Intervals,site,seasn,Params,0,IDLnk);
                epsBaA[site][seasn]=invlogit(alf);
                epsBaa[site][seasn]=invlogit(alf+gam);
            }
            if (SingSpec.Alt_Param_Checked==PSIBA_PSIBa || SingSpec.Alt_Param_Checked==PSIB_PHI) {
                irow=seasn+2*Intervals;
                gamBAA[site][seasn]=gamBAa[site][seasn]=getParam(3+irow,        1,irow,site,seasn,Params,0,0);
				//if (site==0) printf(">>gamBAA=%f par=%f %f %f cov=%f\n",gamBAA[site][seasn],Params[2],Params[3],Params[4],
				 //           SingSpec.SampCov[0][0][0]);
                epsBAA[site][seasn]=epsBAa[site][seasn]=getParam(iroweps1+irow, 2,irow,site,seasn,Params,0,0);
                if (GAMBAAdifferent) {
					irow+=Intervals;
					gamBAa[site][seasn]=getParam(3+irow,        1,irow,site,seasn,Params,0,0);
					epsBAa[site][seasn]=getParam(iroweps1+irow, 2,irow,site,seasn,Params,0,0);
				}
                irow+=Intervals;
                gamBaA[site][seasn]=gamBaa[site][seasn]=getParam(3+irow,        1,irow,site,seasn,Params,0,0);
                epsBaA[site][seasn]=epsBaa[site][seasn]=getParam(iroweps1+irow, 2,irow,site,seasn,Params,0,0);
                if (GAMBAAdifferent) {
					irow+=Intervals;
					gamBaa[site][seasn]=getParam(3+irow,        1,irow,site,seasn,Params,0,0);
					epsBaa[site][seasn]=getParam(iroweps1+irow, 2,irow,site,seasn,Params,0,0);
				}
            }
            if (PHIPARMOPT) {   //  if param name is 'phi' instead of 'eps'...
                epsBAA[site][seasn]=1-epsBAA[site][seasn];
                epsBAa[site][seasn]=1-epsBAa[site][seasn];
                epsBaA[site][seasn]=1-epsBaA[site][seasn];
                epsBaa[site][seasn]=1-epsBaa[site][seasn];
            }
            epsA=psiB[site][seasn]*epsAB[site][seasn]+(1-psiB[site][seasn])*epsAb[site][seasn];
            gamA=psiB[site][seasn]*gamAB[site][seasn]+(1-psiB[site][seasn])*gamAb[site][seasn];
            psiA[site][seasn+1]=psiA[site][seasn]*(1-epsA)+(1-psiA[site][seasn])*gamA;
            
            epsB=    psiA[site][seasn] *((1-epsA)*epsBAA[site][seasn]+epsA*epsBAa[site][seasn])+
                  (1-psiA[site][seasn])*(gamA*epsBaA[site][seasn]+(1-gamA)*epsBaa[site][seasn]);
            gamB=    psiA[site][seasn] *((1-epsA)*gamBAA[site][seasn]+epsA*gamBAa[site][seasn])+
                  (1-psiA[site][seasn])*(gamA*gamBaA[site][seasn]+(1-gamA)*gamBaa[site][seasn]);
            psiB[site][seasn+1]=psiB[site][seasn]*(1-epsB)+(1-psiB[site][seasn])*gamB;
            if (SingSpec.UseNeighbor>0)
                for (kk=0; kk<N; kk++) {
#ifdef neighbor2					
					i=SingSpec.neighbor_lst[kk][site]; if (i<0) break;
					SingSpec.psibar[seasn+1][kk][0]+=psiA[site][seasn+1]*SingSpec.neighbor_wgt[i];
					SingSpec.psibar[seasn+1][kk][1]+=psiB[site][seasn+1]*SingSpec.neighbor_wgt[i];
#else						
                    if (SingSpec.neighbor[kk][site]>1) break;
                    x=(double)SingSpec.neighbor[kk][site]*SingSpec.neighbor_wgt[site];
                    SingSpec.psibar[seasn+1][kk][0]+=psiA[site][seasn+1]*x;
                    SingSpec.psibar[seasn+1][kk][1]+=psiB[site][seasn+1]*x;
#endif
                }
            if(site<0)
                printf("site=%d psiA(%d)=%f psiA(%d)=%f epsA=%f gamA=%f\npsiB=%f epsAB=%f epsAb=%f gamAB=%f gamAb=%f\n",
                       site,seasn+1,psiA[site][seasn+1],seasn,psiA[site][seasn],
                       epsA,gamA,psiB[site][seasn],epsAB[site][seasn],epsAb[site][seasn],
                       gamAB[site][seasn],gamAb[site][seasn]);
        }
    }
    //printf("\n%f %f %f\n%f %f %f %f %f %f\n%f %f %f %f %f %f\n",psiA[0][0],psiB[0][0],psiAB[0][0],
	  //gamAB[0][0],gamAb[0][0],gamBAA[0][0],gamBAa[0][0],gamBaA[0][0],gamBaa[0][0],
	  //epsAB[0][0],epsAb[0][0],epsBAA[0][0],epsBAa[0][0],epsBaA[0][0],epsBaa[0][0]);
    for (site=0; site<N; site++) {  // loop through individual sites
        likeA = psiA[site][0]-psiAB[site][0]; likeB = psiB[site][0]-psiAB[site][0]; 
        likeAB = psiAB[site][0]; like0 = 1-(likeA+likeB+likeAB);
        
        // P[0]=pA  pB=pB  P[2]=pAB  P[3]=pAB'  P[4]=pA'B
        // new definition of detection probabilities
        // P[2] = rA = Pr(detect A | A & B present)
        // P[3] = rB = Pr(detect B | A & B present)
        // P[4] = rAB = Pr(detect A&B | A&B present) = pA x pB if independent
        
        //calculate pA's and pB's
        survey=-1; ip=Intervals*8+3;
        if (GAMBAAdifferent) ip=Intervals*12+3;
        for (seasn=0; seasn<PrmyPeriods; seasn++) { seenA=seenB=0;
            for (srvy=0; srvy<SingSpec.SecPeriods[seasn]; srvy++) { survey++;
                if (SingSpec.Data[site][survey]!=-1) {  //   if not missing data for survey ...
                    pA=rA=getParam(ip+survey,3,survey,site,survey,Params,0,0);
                    pB=rB=getParam(ip+survey+T,3,survey+T,site,survey,Params,0,0);
                    //    rBa = Pr(detect B | not detected A)         case is important!
                    //    rBA = Pr(detect B | detected A)
                    //    rA = Pr(detect A)   rB = Pr(detect B)   rAB = Pr(detect A and B)
                    rA=getParam(ip+survey+2*T,3,survey+2*T,site,survey,Params,0,0);
                    rB=getParam(ip+survey+3*T,3,survey+3*T,site,survey,Params,0,0); 
                    //if(site==0)printf("pA=getparm(%d)=%f pB=%f rA=%f rB=%f\n",ip+survey,pA,pB,rA,rB);
                    rAB = rA*rB;
                    if (SingSpec.Alt_Param_Checked==0) {
                        delta=getParam(ip+survey+4*T,3,survey+4*T,site,survey,Params,1,0);
                        // check whether value is in allowable bounds
                        lower = (rA+rB-1.0)/(rA*rB); if (lower<0.0) lower=0.0;
                        upper = 1.0/rB; if (upper>(1.0/rA)) upper=1.0/rA;
                        if(delta < lower) return((fabs(delta-lower)+1.0)*1.0e10);
                        if(delta > upper) return((fabs(delta-upper)+1.0)*1.0e10);
                        rAB = rA*rB*delta;
                    }
                    if (SingSpec.Alt_Param_Checked==1) {
                        rBA=rB; rBa=getParam(ip+survey+4*T,3,survey+4*T,site,survey,Params,0,0);
                        rAB=rA*rBA; rB=rA*rBA+(1-rA)*rBa;
                    }
                    if (SingSpec.Alt_Param_Checked==2) {
                        rBa=getParam(ip+survey+3*T,3,survey+3*T,site,survey,Params,0,0);
                        alf=getParam(ip+survey+3*T,3,survey+3*T,site,survey,Params,0,IDLnk);
                        gam=getParam(ip+survey+4*T,3,survey+4*T,site,survey,Params,0,IDLnk);  //  logit(rBA) = alpha + gamma 
                        rBA=invlogit(alf+gam); rB=rA*rBA+(1-rA)*rBa; rAB=rA*rBA; 
                    }
                    if (SingSpec.Alt_Param_Checked==3) {   //  newest parameterzation - psiAb, psiBa, odds-ratio(psiAb:psiA|B)
                        a=getParam(ip+survey+2*T,3,survey+2*T,site,0,Params,0,IDLnk);  //  untransformed = logit(rAb)
                        b=getParam(ip+survey+3*T,3,survey+3*T,site,0,Params,0,IDLnk);  //  untransformed = logit(rBa)
                        c=getParam(ip+survey+4*T,3,survey+4*T,site,0,Params,0,IDLnk);  
                        K=1+exp(a+b+c)+exp(a)+exp(b);
                        rAB=exp(a+b+c)/K; rAwoB=exp(a)/K; rBwoA=exp(b)/K;   
                        rA=rAB+rAwoB; rB=rAB+rBwoA;
                    }
                    tA=tB=0;
                    if (SingSpec.Data[site][survey]==1 || SingSpec.Data[site][survey]==3) tA = 1;
                    if (SingSpec.Data[site][survey]==2 || SingSpec.Data[site][survey]==3) tB = 1;
                    if (SingSpec.Nstates==2) if (SingSpec.Data[site+N][survey]==1) tB = 1;
                    switch (tA + 2*tB) {
                    case 0:     // neither detected
                        likeAB *= (1.0-rA-rB+rAB); likeA *= (1.0-pA); likeB *= (1.0-pB); break;
                    case 1: likeAB *= rA-rAB; likeA *= pA; seenA = 1; break;   // A detected, not B
                    case 2: likeAB *= rB-rAB; likeB *= pB; seenB = 1; break;   // not A, B detected
                    case 3: likeAB *= rAB;  seenA = seenB = 1;  break;    // A and B detected
                    }   // end switch
                    //if(site==0) printf("%d %d tA=%d tB=%d rA=%f rB=%f rAB=%f likeA=%f likeB=%f likeAB=%f\n",seasn,srvy,tA,tB,rA,rB,rAB,likeA,likeB,likeAB);
                } // end if not missing..
            } // end survey=0..nsrvys/season
            
            prVec[1]=likeA; prVec[2]=likeB; prVec[3]=likeAB; prVec[0]=like0;
            
            switch(1*seenA + 2*seenB) {
            case 0:                      break;    //neither species detected  
            case 1: prVec[0]=prVec[2]=0; break;    // A seen, not B
            case 2: prVec[0]=prVec[1]=0; break;    // not A, B seen
            case 3: prVec[0]=prVec[1]=prVec[2]=0;          break;    // both seen
            }
			//              transition matrix  
//                  0                     A                    B                    AB
//      0      (1-gamAb)*(1-gamBaa) gamAb*(1-gamBaA)    (1-gamAb)*gamBaa      gamAb*gamBaA
//      A      epsAb*(1-gamBAa)    (1-epsAb)*(1-gamBAA)   epsAb*gamBAa      (1-epsAb)*gamBAA
//      B      (1-gamAB)*epsBaa      gamAB*epsBaA     (1-gamAB)*(1-epsBaa)  gamAB*(1-epsBaA)
//     AB      epsAB*epsBAa        (1-epsAB)*epsBAA     epsAB*(1-epsBAa)   (1-epsAB)*(1-epsBAA)
             
            if (seasn<Intervals) {
                
                xphi[0][0]=(1-gamAb[site][seasn])*(1-gamBaa[site][seasn]);  // see transition matrix above  
                xphi[0][1]=gamAb[site][seasn]*(1-gamBaA[site][seasn]);        
                xphi[0][2]=(1-gamAb[site][seasn])*gamBaa[site][seasn];        
                xphi[0][3]=gamAb[site][seasn]*gamBaA[site][seasn];            
                xphi[1][0]=epsAb[site][seasn]*(1-gamBAa[site][seasn]);        
                xphi[1][1]=(1-epsAb[site][seasn])*(1-gamBAA[site][seasn]);    
                xphi[1][2]=epsAb[site][seasn]*gamBAa[site][seasn];            
                xphi[1][3]=(1-epsAb[site][seasn])*gamBAA[site][seasn];        
                xphi[2][0]=(1-gamAB[site][seasn])*epsBaa[site][seasn];        
                xphi[2][1]=gamAB[site][seasn]*epsBaA[site][seasn];            
                xphi[2][2]=(1-gamAB[site][seasn])*(1-epsBaa[site][seasn]);    
                xphi[2][3]=gamAB[site][seasn]*(1-epsBaA[site][seasn]);        
                xphi[3][0]=epsAB[site][seasn]*epsBAa[site][seasn];            
                xphi[3][1]=(1-epsAB[site][seasn])*epsBAA[site][seasn];        
                xphi[3][2]=epsAB[site][seasn]*(1-epsBAa[site][seasn]);        
                xphi[3][3]=(1-epsAB[site][seasn])*(1-epsBAA[site][seasn]);    
                
                for (i=0; i<4; i++) { //  matrix mult 1x4 vector by 4x4 matrix
                    xsum=0; for (j=0; j<4; j++) xsum+=prVec[j]*xphi[j][i]; tmpVec[i]=xsum;
                }
                like0=tmpVec[0]; likeA=tmpVec[1]; likeB=tmpVec[2]; likeAB=tmpVec[3];
            }   // end if seasn<Intervals
        } // end for seasn=1..prmyPeriods
        tmp=prVec[0]+prVec[1]+prVec[2]+prVec[3]; 
        if(site<0)printf(">>%d %f %f %f %f %f %f %f\n",site+1,prVec[0],prVec[1],prVec[2],prVec[3],
                          tmp*10000-SingSpec.det_hist_frq[site],tmp*10000,SingSpec.det_hist_frq[site]);
        if (tmp < SingSpec.nearzero) tmp = SingSpec.nearzero; 
        if (tmp>1.0000000000001) { 
            printf("tmp>1(%22.20e) site=%d prVec=%f %f %f %f\n",tmp,site+1,prVec[0],prVec[1],prVec[2],prVec[3]); 
            exit(1);
        }
        sum-=SingSpec.det_hist_frq[site]*log(tmp);
        if (site<0) {
            printf("\npsiA,psiAB,psiBA,psiBa,psiB:%f %f %f %f %f\n",psiA[0][0],psiAB[0][0],psiBA,psiBa,psiB[0][0]);
            printf("pA,pB,rA,rBA,rBa:%f %f %f %f %f\n",pA,pB,rA,rAB,rB);
            printf("gam AB,Ab,BAA,BAa,BaA,Baa:%f %f %f %f %f %f\n",
		       gamAB[0][0],gamAb[0][0],gamBAA[0][0],gamBAa[0][0],gamBaA[0][0],gamBaa[0][0]);
            printf("eps AB,Ab,BAA,BAa,BaA,Baa:%f %f %f %f %f %f\n",
			   epsAB[0][0],epsAb[0][0],epsBAA[0][0],epsBAa[0][0],epsBaA[0][0],epsBaa[0][0]);
			printf("parms:%f %f %f %f %f\n%f %f %f %f %f\n",Params[0],Params[1],Params[2],Params[3],Params[4],Params[5],Params[6],Params[7],Params[8],Params[9]);
       
        }
        delete [] psiA[site]; delete [] psiB[site]; delete [] psiAB[site];
        delete [] gamAB[site]; delete [] gamAb[site]; 
        delete [] gamBAA[site]; delete [] gamBAa[site]; delete [] gamBaA[site]; delete [] gamBaa[site];
        delete [] epsAB[site]; delete [] epsAb[site]; 
        delete [] epsBAA[site]; delete [] epsBAa[site]; delete [] epsBaA[site]; delete [] epsBaa[site];
    } // end individual site loop
    delete [] psiA; delete [] psiB; delete [] psiAB; 
    delete [] gamAB; delete [] gamAb; delete [] gamBAA; delete [] gamBAa; delete [] gamBaA; delete [] gamBaa;
    delete [] epsAB; delete [] epsAb; delete [] epsBAA; delete [] epsBAa; delete [] epsBaA; delete [] epsBaa;
    if (SingSpec.Verbose>0) printf("%f %f %f\n",sum,pnlty,sum+pnlty);
    //printf("%6d %16.4f\r",SingSpec.ifn,sum);
    return(sum+pnlty);
}
