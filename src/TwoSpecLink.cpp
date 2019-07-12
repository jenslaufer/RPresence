#include "SingAnal.h"
//****  Likelihood for joint modelling of 2 species from presence/absence data  ****//
#include <math.h>
#include <stdio.h>

double TwoSpeciesLink(double *Params, int NPar) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    double invlogit(double x);
    int survey, site, seenA, seenB, N=SingSpec.N, T=SingSpec.T;
    double sum = 0.0, tmp=0, pA, pB, rA=0, rB=0, rAB=0, rBA=0, rBa=0,
	             psiA, psiB=0, psiAB=0, likeA, likeB, likeAB, psiBA=0, psiBa=0;
    double lower, upper, phi, delta, alf, gam, xnu, denom,
	   logitPsiAb,logitPsiBa,logtau,logitPsiAgB,logitPsiBgA,psiAgB,psiBgA,
       logitrAb,logitrBa,logitrAgB,logitrBgA,rAgB,rBgA;
    // psi[AB] = psi[A]*psi[B]*phi; phi is parameter being estimated
    // alt parameterization: Psi[2] = psi[Ba] = Pr(occ by B|not occ by A)
    for (site=0; site<N; site++) {  // loop through individual sites
        // calculate psi[A] and psi[B]
        psiA=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
        if (SingSpec.Alt_Param_Checked==0)  {
            psiB=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn); psiAB = psiA*psiB;
            phi=getParam(2,0,2,site,0,Params,1,expLnk);
            // check whether value is in allowable bounds
            lower = (psiA+psiB-1.0)/(psiA*psiB); if (lower<0.0) lower=0.0;
            upper = 1.0/psiA; if (upper>(1.0/psiB)) upper=1.0/psiB;
            if(phi < lower) return((fabs(phi-lower)+1.0)*1.0e10);
            if(phi > upper) return((fabs(phi-upper)+1.0)*1.0e10);
            psiAB = psiA*psiB*phi;
        }
        if (SingSpec.Alt_Param_Checked==1) {            //  new parmeterization - psiA, psiBA, psiBa
            psiBA=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
            psiBa=getParam(2,0,2,site,0,Params,0,SingSpec.LinkFn); 
            psiAB=psiA*psiBA;                //  Pr(occ by B) = Pr(occ by A)*pr(Occ by B | A present)
            psiB=psiA*psiBA+(1-psiA)*psiBa;  //    + Pr(not occ by A)*Pr(occ by B | A not present)
        }
        if (SingSpec.Alt_Param_Checked==2) {   //  newer parameterzation - psiA, psiBa, odds-ratio(psiBa:psiBA)
            psiBa=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
            alf=50; if (psiBa<1) alf=log(psiBa/(1-psiBa));
            //alf=getParam(1,0,1,site,0,Params,0,IDLnk);    //  logit(psiBa) = alpha, 
            xnu=getParam(2,0,2,site,0,Params,0,expLnk);    //  logit(psiBA) = alpha + gamma 
            gam=-50; if (xnu>0) gam=log(xnu);
            psiBA=invlogit(alf+gam); psiB=psiA*psiBA+(1-psiA)*psiBa; psiAB=psiA*psiBA; 
        }
        if (SingSpec.Alt_Param_Checked==3) {   //  newest parameterzation - psiAb, psiBa, odds-ratio(psiAb:psiA|B)
            logitPsiAb=getParam(0,0,0,site,0,Params,0,IDLnk);  //  untransformed = logit(psiAb)
            logitPsiBa=getParam(1,0,1,site,0,Params,0,IDLnk);  //  untransformed = logit(psiBa)
            logtau=getParam(2,0,2,site,0,Params,0,IDLnk);  
            logitPsiAgB=logitPsiAb+logtau;  logitPsiBgA=logitPsiBa+logtau;  
			psiAgB=invlogit(logitPsiAgB);	psiBgA=invlogit(logitPsiBgA);
		    
			tmp=exp(logitPsiAb+logitPsiBa+logtau); 
			denom=1+tmp+exp(logitPsiAb)+exp(logitPsiBa);
			psiAB=tmp/denom; 
			//tmp=psiAgB/psiBgA;                 //  do some algebra to get psiAB
			//psiAB=(psiAb-tmp*psiBa)/(tmp-1);
			psiA=psiAB/psiBgA; psiB=psiAB/psiAgB;
        }
        likeA = psiA-psiAB; likeB = psiB-psiAB; likeAB = psiAB;  seenA=0; seenB=0; 
        //printf("%f %f %f %f %f\n",psiA,psiB,psiAB,psiBA,psiBa);
        
        // P[0]=pA  pB=pB  P[2]=pAB  P[3]=pAB'  P[4]=pA'B
        // new definition of detection probabilities
        // P[2] = rA = Pr(detect A | A & B present)
        // P[3] = rB = Pr(detect B | A & B present)
        // P[4] = rAB = Pr(detect A&B | A&B present) = pA x pB if independent
        
        //calculate pA's and pB's
        for (survey=0; survey<T; survey++) {
            if (SingSpec.Data[site][survey]!=-1) {
                pA=rA=getParam(3+survey,1,survey,site,survey,Params,0,SingSpec.LinkFn);
                pB=rB=getParam(3+survey+T,1,survey+T,site,survey,Params,0,SingSpec.LinkFn);
                rA=getParam(3+survey+2*T,1,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
                rB=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn);                
                rAB = rA*rB;
				switch(SingSpec.Alt_Param_Checked) {
                case 0:
                    delta=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,1,SingSpec.LinkFn);
                    // check whether value is in allowable bounds
                    lower = (rA+rB-1.0)/(rA*rB); if (lower<0.0) lower=0.0;
                    upper = 1.0/rB; if (upper>(1.0/rA)) upper=1.0/rA;
                    if(delta < lower) return((fabs(delta-lower)+1.0)*1.0e10);
                    if(delta > upper) return((fabs(delta-upper)+1.0)*1.0e10);
                    rAB = rA*rB*delta;
					break;
                case 1:
                    rBA=rB; rBa=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,0,SingSpec.LinkFn);
                    rAB=rA*rBA; rB=rA*rBA+(1-rA)*rBa;
					break;
                case 2:
                    rBa=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn);
                    alf=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,IDLnk);
                    gam=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,0,IDLnk);  //  logit(rBA) = alpha + gamma
					if (gam>0) gam=log(gam); else gam=-50;
                    rBA=invlogit(alf+gam); rB=rA*rBA+(1-rA)*rBa; rAB=rA*rBA; 
					break;
                case 3:                    //  newest parameterzation - psiAb, psiBa, odds-ratio(psiAb:psiA|B)
                    logitrAb=getParam(0,0,0,site,0,Params,0,IDLnk);  //  untransformed = logit(rAb)
                    logitrBa=getParam(1,0,1,site,0,Params,0,IDLnk);  //  untransformed = logit(rBa)
                    logtau=getParam(2,0,2,site,0,Params,0,IDLnk);  
                    logitrAgB=logitrAb+logtau;  logitrBgA=logitrBa+logtau;  
			        rAgB=invlogit(logitrAgB);	rBgA=invlogit(logitrBgA);
					
					tmp=exp(logitrAb+logitrBa+logtau); 
					denom=1+tmp+exp(logitrAb)+exp(logitrBa);
					rAB=tmp/denom; 
			        rA=rAB/rBgA; rB=rAB/rAgB;
                }
                //    rBa = Pr(detect B | not detected A)         case is important!
                //    rBA = Pr(detect B | detected A)
                //    rAB = Pr(detect A and B)
                //    rA = Pr(detect A)
                //    rB = Pr(detect B)
                switch (SingSpec.Data[site][survey]) {
                case 0:     // neither detected
                    likeAB *= (1.0-rA-rB+rAB); likeA *= (1.0-pA); likeB *= (1.0-pB); break;
                case 1: likeAB *= rA-rAB; likeA *= pA; seenA = 1; break;   // A detected, not B
                case 2: likeAB *= rB-rAB; likeB *= pB; seenB = 1; break;   // not A, B detected
                case 3: likeAB *= rAB;  seenA = seenB = 1;  break;         // A and B detected
				
				case 4: likeAB *= (1-rA); likeA *= (1-pA);  break;   // A not detected, B missing
				case 5: likeAB *= rA; likeA *= pA; seenA=1; break;             // A detected, B missing
				
				case 6: likeAB *= (1-rB); likeB *= (1-pB); break; // B not detected, A missing
				case 7: likeAB *= rB; likeB *= pB; seenB = 1;  break;        // B detected, A missing			
                }   // end switch
            } // end if missing
        } // end survey=0..T
        switch(1*seenA + 2*seenB) {
        case 0: tmp = likeAB + likeA + likeB + (1.0-psiA-psiB+psiAB); break;    //neither species detected  
        case 1: tmp = likeAB + likeA; break;    // A seen, not B
        case 2: tmp = likeAB + likeB; break;    // not A, B seen
        case 3: tmp = likeAB;  break;    // both seen
        } // end switch 
        if (tmp < SingSpec.nearzero) tmp = SingSpec.nearzero;
        sum-=SingSpec.det_hist_frq[site]*log(tmp);
    } // end individual site loop
    return(sum);
}
