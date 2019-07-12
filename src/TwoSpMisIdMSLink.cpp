#include "SingAnal.h"
//****  Likelihood for joint modelling of 2 species w/ mis-identification as per Chambert (2018)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MISSNG -1
double TwoSpMisIdMSLink(double *Params, int NPar) {
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    extern TSingSpec SingSpec;
    double invlogit(double x);
    int survey, site, N=SingSpec.N, T=SingSpec.T;
    int i,j,k,l,*w,*y,g_ij,seasn,Intervals,epsptr,irow,GAMBAAdifferent,oaptr,gamptr;
    double sum = 0.0,  pA, pB, rA=0, rB=0, rAB=0, rBA=0, rBa=0,
                 psiA, psiB=0, psiAB=0, psiBA=0, psiBa=0,
                 gamphi,oA,oa,oB,ob,cA,cB,cAB,den,lower,upper,po30,po31,po32,po33,
				 trmat00,trmat01,trmat02,trmat03,trmat10,trmat11,trmat12,trmat13,
				 trmat20,trmat21,trmat22,trmat23,trmat30,trmat31,trmat32,trmat33,
				 gamAB,gamAb,gamBAA,gamBAa,gamBaA,gamBaa,epsAB,epsAb,epsBAA,epsBAa,epsBaA,epsBaa,
				 newphi0,newphi1,newphi2,newphi3;

    double po[4][4]={1,0,0,0,  0,1,0,0,   0,0,1,0,  0,0,0,0};
    double pcg1[4][4]={1,0,0,0,  0,1,0,0,   0,0,1,0,  0,0,0,0};
    double mu[4],alpha[4],pcg0[4][4],phi[4],pc[4][4];

    y=(int*)malloc(T*sizeof(int)); w=(int*)malloc(T*sizeof(int));
    pcg0[0][0]=1; pcg0[0][1]=pcg0[0][2]=pcg0[0][3]=0;
    pcg0[1][0]=1; pcg0[1][1]=pcg0[1][2]=pcg0[1][3]=0;
    pcg0[2][0]=1; pcg0[2][1]=pcg0[2][2]=pcg0[2][3]=0;
    pcg0[3][0]=1; pcg0[3][1]=pcg0[3][2]=pcg0[3][3]=0;

    pcg1[0][0]=pcg1[1][1]=pcg1[2][2]=1;
    pcg1[0][1]=pcg1[0][2]=pcg1[0][3]=pcg1[1][0]=pcg1[1][2]=pcg1[1][3]=pcg1[2][0]=pcg1[2][1]=pcg1[2][3]=pcg1[3][0]=0;
    Intervals=SingSpec.PrmyPeriods-1;
	epsptr=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2];
	GAMBAAdifferent=(SingSpec.NrowsDM[2]>(4*SingSpec.T));
	oaptr=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3];
	gamptr=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1];
    // psi[AB] = psi[A]*psi[B]*phi; phi is parameter being estimated
    // alt parameterization: Psi[2] = psi[Ba] = Pr(occ by B|not occ by A)
    for (site=0; site<N; site++) {  // loop through individual sites
        i=site;
        psiA=       getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); // get psi[A] and psi[B]
        psiBA=psiB= getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
        if (SingSpec.Alt_Param_Checked!=0)  {// interaction = dominance ...  psiA,psiBA,psiBa parameterization
            psiBa=getParam(2,0,2,site,0,Params,0,SingSpec.LinkFn);
            phi[0]=(1-psiA)*(1-psiBa); phi[1]=psiA*(1-psiBA); phi[2]=(1-psiA)*psiBa; phi[3]=psiA*psiBA;
        }
        else {
            gamphi=getParam(2,0,2,site,0,Params,1,SingSpec.LinkFn);
            // check whether value is in allowable bounds
            lower = (psiA+psiB-1.0)/(psiA*psiB); if (lower<0.0) lower=0.0;
            upper = 1.0/psiA; if (upper>(1.0/psiB)) upper=1.0/psiB;
            if(gamphi < lower) {sum=9.9999E22; printf("err=%f\n",sum); return(sum);}
            if(gamphi > upper) {sum=8.8888E22; printf("Err=%f\n",sum); return(sum);}
            psiAB=psiA*psiB*gamphi;   // interaction = general
            //     phi = state vector... 0=unocc, 1=spA only, 2=spB only, 3=both sp.
            phi[0]=1-psiA-psiB+psiAB; phi[1]=psiA-psiAB; phi[2]=psiB-psiAB; phi[3]=psiAB;
        }       
        for (seasn=survey=j=0; seasn<SingSpec.PrmyPeriods; seasn++) {
            mu[0]=mu[1]=mu[2]=mu[3]=alpha[0]=alpha[1]=alpha[2]=alpha[3]=1;			
			//calculate pA's and pB's
			for (survey=0; survey<SingSpec.SecPeriods[seasn]; survey++) {
				j=survey;
				if (SingSpec.Data[i][j]!=MISSNG) {
					y[j]=(SingSpec.Data[i][j] & 3);    //  y is the observed data (in rightmost 2 bits)
					w[j]=(SingSpec.Data[i][j] >>2);      // w is the confirmed data. (in 2nd pair of bits from right)

									//### DETECTION
					pA=    getParam(3+survey,    1,survey,    site,survey,Params,0,SingSpec.LinkFn);
					pB=    getParam(3+survey+T,  1,survey+T,  site,survey,Params,0,SingSpec.LinkFn);
					rA=    getParam(3+survey+2*T,1,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
					rBA=rB=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn);
					if (SingSpec.Alt_Param_Checked!=0)  {    // interaction = dominance ...  rA,rBA,rBa parameterization
						rBa=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,0,SingSpec.LinkFn);
						po30=(1-rA)*(1-rBa); po31=rA*(1-rBA); po32=(1-rA)*rBa; po33=rA*rBA;
					}
					else { // interaction = general ...  rA,rB,delta parameterization
						gamphi=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,1,SingSpec.LinkFn);
						// check whether value is in allowable bounds
						lower = (rA+rB-1.0)/(rA*rB); if (lower<0.0) lower=0.0;
						upper = 1.0/rA; if (upper>(1.0/rB)) upper=1.0/rB;
						if(gamphi < lower) {sum=9.9999E22; printf("err=%f\n",sum); return(sum);}
						if(gamphi > upper) {sum=8.8888E22; printf("Err=%f\n",sum); return(sum);}
						rAB=rA*rB*gamphi; po30=1-rA-rB+rAB; po31=rA-rAB; po32=rB-rAB; po33=rAB;
					}
								//### MISIDENTIFICATION
					oA=getParam(oaptr+survey,4,survey,site,survey,Params,0,SingSpec.LinkFn);
					oa=getParam(oaptr+T+survey,4,survey+T,site,survey,Params,0,SingSpec.LinkFn);
					oB=getParam(oaptr+2*T+survey,4,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
					ob=getParam(oaptr+3*T+survey,4,survey+3*T,site,survey,Params,0,SingSpec.LinkFn);
									//### CONFIRMATION        (use MLOGIT here)
					cA=getParam(oaptr+4*T+survey,5,survey,site,survey,Params,0,expLnk);
					cB=getParam(oaptr+5*T+survey,5,survey+T,site,survey,Params,0,expLnk);
					den=1+cA+cB; cA=cA/den; cB=cB/den; cAB=1-cA-cB;
					SingSpec.realParmEst[3+survey+9*T][site]=cA;
					SingSpec.realParmEst[3+survey+10*T][site]=cB;
					// ## DETECTION MATRIX 'po' (Non-Confirmed Observations)
					//po[0][0]=1; po[0][1]=po[0][2]=po[0][3]=0;
					po[1][0]=(1-pA)*(1-oa); po[1][1]=pA*(1-oA); po[1][2]=(1-pA)*oa; po[1][3]=pA*oA;
					po[2][0]=(1-pB)*(1-ob); po[2][1]=(1-pB)*ob; po[2][2]=pB*(1-oB); po[2][3]=pB*oB;
					po[3][0]=po30; po[3][1]=po31; po[3][2]=po32; po[3][3]=po33;
					//## SPECIES CONFIRMATION  MATRIX 'pc' (Confirmed Observations)
					pcg1[3][1]=cA; pcg1[3][2]=cB; pcg1[3][3]=cAB;

					g_ij=(w[j]>0);
					for (k=0; k<4; k++) {
						for (l=0; l<4; l++) pc[k][l]=(1-g_ij)*pcg0[k][l]+g_ij*pcg1[k][l];
						mu[k]*=po[k][y[j]]; // y(j)=SingSpec.Data(i,j)// ## Probability MATRIX for the OBSERVATION PIECE
					}
					for (k=0; k<4; k++) alpha[k]*=pc[k][w[j]];
				}       // end if Data[i][j] not missing...
			}   //   end for (survey=0 to SecPeriods[seasn]...
			for (k=0; k<4; k++) phi[k]=phi[k]*mu[k]*alpha[k];
			if (seasn<SingSpec.PrmyPeriods) {
				gamAB =getParam(gamptr+seasn,                 2,seasn,                site,seasn,Params,0,0);
				gamAb =getParam(gamptr+seasn+Intervals      , 2,seasn+Intervals,site,seasn,Params,0,0);
				epsAB =getParam(epsptr+seasn          , 3,seasn,                site,seasn,Params,0,0);
				epsAb =getParam(epsptr+seasn+Intervals, 3,seasn+Intervals,site,seasn,Params,0,0); 
				irow=seasn+2*Intervals;
				gamBAA=gamBAa=getParam(gamptr+irow, 2,irow,site,seasn,Params,0,0);
				epsBAA=epsBAa=getParam(epsptr+irow, 3,irow,site,seasn,Params,0,0);
				if (GAMBAAdifferent) {
					irow+=Intervals;
					gamBAa=getParam(gamptr+irow, 2,irow,site,seasn,Params,0,0);
					epsBAa=getParam(epsptr+irow, 3,irow,site,seasn,Params,0,0);
				}
				irow+=Intervals;
				gamBaA=gamBaa=getParam(gamptr+irow, 2,irow,site,seasn,Params,0,0);
				epsBaA=epsBaa=getParam(epsptr+irow, 3,irow,site,seasn,Params,0,0);
				if (GAMBAAdifferent) {
					irow+=Intervals;
					gamBaa=getParam(gamptr+irow, 2,irow,site,seasn,Params,0,0);
					epsBaa=getParam(epsptr+irow, 3,irow,site,seasn,Params,0,0);
				}    //  compute transition matrix...
				trmat00=(1-gamAb)*(1-gamBaa); trmat01=gamAb*(1-gamBaA);     trmat02=(1-gamAb)*gamBaa;      trmat03=gamAb*gamBaA;
				trmat10=epsAb*(1-gamBAa);     trmat11=(1-epsAb)*(1-gamBAA); trmat12=epsAb*gamBAa;          trmat13=(1-epsAb)*gamBAA;
				trmat20=(1-gamAB)*epsBaa;     trmat21=gamAB*epsBaA;          trmat22=(1-gamAB)*(1-epsBaa); trmat23=gamAB*(1-epsBaA);
				trmat30=epsAB*epsBAa;          trmat31=(1-epsAB)*epsBAA;      trmat32=epsAB*(1-epsBAa);     trmat33=(1-epsAB)*(1-epsBAA);
				//         matrix multiply phi * transition matrix , then store in phi...
				newphi0=phi[0]*trmat00+phi[1]*trmat10+phi[2]*trmat20+phi[3]*trmat30;
				newphi1=phi[0]*trmat01+phi[1]*trmat11+phi[2]*trmat21+phi[3]*trmat31;
				newphi2=phi[0]*trmat02+phi[1]*trmat12+phi[2]*trmat22+phi[3]*trmat32;
				newphi3=phi[0]*trmat03+phi[1]*trmat13+phi[2]*trmat23+phi[3]*trmat33;
				phi[0]=newphi0; phi[1]=newphi1; phi[2]=newphi2; phi[3]=newphi3;
			}   //  end of if(seasn<PrmyPeriods)...
		}    // end of seasn=0..PrmyPeriods...
		sum-=SingSpec.det_hist_frq[i]*log(phi[0]+phi[1]+phi[2]+phi[3]);
	} // end individual site loop
    //  free memory and return
    free(w); free(y);
    //printf("%f %f %f %f %f %f %f %f %f %f %f\n",psiA,psiB,gamphi,pA,pB,rA,rB,oa,ob,cA,cB);
    if (!isnormal(sum)) sum=9.999999e22;
    return(sum);
}
