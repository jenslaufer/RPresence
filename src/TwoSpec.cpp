#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SingAnal.h"

void TwoSpec(int NBoot, char *s, int NPar, double *Params) { 
  FILE *g;
  void TwoSpCorrDet(int NBoot, char *s, int NPar, double *Params);
  double TwoSpeciesLink(double *Params, int NPar);
  void print_individual_estimates2(FILE *g, int parmnum, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn, int fsite, int nsites);
  void prntBetas(double xxLike(double *p, int N), int NPar, double *Params, double *Work, double MaxLL, double **covar);
  void DoAmoeba( double *Params, int NPar, double (*funk)(double [], int), double xincr);
  double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
  void try_optmz_randiv(double (*xxLinkFn)(double Params[], int npar), int NPar, double *Params, double *Work, double *Hess, double *MaxLL);

  extern TSingSpec SingSpec;
  int i, ii, jj, T, N, ndec=4, site, survey, naiveA=0, naiveB=0, dA, dB;
  double *Hess, *Work, MaxLL, **covar;
  
  if (SingSpec.NrowsDM[0]>3) 
	TwoSpCorrDet(NBoot,s,NPar,Params);
  else {
    g=SingSpec.g; T = SingSpec.T;  N = SingSpec.N;
    fprintf(g,"Two Species, Single Season Model:\n---------------------------------\n");
    SingSpec.Alt_Param_Checked=0;
    if (strstr(SingSpec.Realname[0][2],"psiBa") != (void*)NULL || strstr(SingSpec.Realname[0][2],"psiB2") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=1;
        fprintf(g,"-- alternate parameterization requested: psiA,psiBA,psiBa,pA,pB,rA,rBA,rBa\n");
    }
    if (strstr(SingSpec.Realname[0][2],"nu") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=2;
        fprintf(g,"-- alternate parameterization requested: psiA,psiBa,nu,pA,pB,rA,rBa,rho\n");
    }
    if (strstr(SingSpec.Realname[0][2],"tau") != (void*)NULL) {
        SingSpec.Alt_Param_Checked=3;
        fprintf(g,"-- alternate parameterization requested: psiAb,psiBa,tau,pA,pB,rAb,rBa,ups\n");
    }
    if (SingSpec.Alt_Param_Checked==0) {
            fprintf(g,"-- standard parameterization requested: psiA,psiB,phi,pA,pB,rA,rB,delta\n");
    }
    if (SingSpec.Nstates==2) {   //  if data only contains 0's and 1's... then stacked input!
        SingSpec.N/=2; N/=2; SingSpec.NMiss/=2; 
        fprintf(g,"-- stacked input: 1st %d records=species A, last %d records=species B\n",N,N);
		for (site=0; site<N; site++) {
			for (survey=dA=dB=0; survey<=T; survey++) {                              //  combined codes:  1=det of A only
				if (SingSpec.Data[site][survey]>0) dA=1;
				if (SingSpec.Data[N+site][survey]>0) dB=1;
				if (SingSpec.Data[site][survey]>=0) {                           //                  2=det of B only
					if (SingSpec.Data[N+site][survey]>=0)                        //                 3=det of  A & B  
						SingSpec.Data[site][survey]=SingSpec.Data[site][survey]+2*SingSpec.Data[N+site][survey];
					else                                                              //              4=no det of A, B not sampled(missing)
 						SingSpec.Data[site][survey]+=4; //  A not det or det, B missing              5=det of A, B not sampled(missing)
				}                                                                      //             6=no det of B, A not sampled     
				else {   //  species A missing...                                                     7=det of B, A not sampled
					if (SingSpec.Data[N+site][survey]>=0) 
						SingSpec.Data[site][survey]=SingSpec.Data[N+site][survey]+6;       //        -1=Neither A or B sampled
					else
						SingSpec.Data[site][survey]=-1; //  both A and B missing					
				}
			}
			naiveA+=dA; naiveB+=dB;
		}
    }
    else  {
		fprintf(g,"-- compressed input detected: 1=species A, 2=species B, 3=both\n");
		for (site=0; site<N; site++) {        //  compute naive occupancy for each species...
			for (survey=dA=dB=0; survey<=T; survey++) {   
				i=SingSpec.Data[site][survey];
				if (i==1 || i==3 || i==5) dA=1;
				if (i==2 || i==3 || i==7) dB=1;
			}
			naiveA+=dA; naiveB+=dB;
		}
	}
	fprintf(g,"Naive occupancy for species A:%8.4f\n",(double)naiveA/(double)N);
	fprintf(g,"Naive occupancy for species B:%8.4f\n",naiveB/(double)N);
	
			// try amoeba routine to get starting values
    
    fprintf(g,"\nModel:%s  NPar=%d  s=%s\n",SingSpec.modname,NPar,s);

    if (SingSpec.UseAmoeba) DoAmoeba(Params, NPar, TwoSpeciesLink, 1.0);
    Hess = new double[NPar*(NPar+1)/2+2]; Work = new double[3*NPar+2]; 
    covar = new double*[NPar]; for (ii=0; ii<NPar; ii++) covar[ii]=new double[NPar];

	try_optmz_randiv(TwoSpeciesLink, NPar, Params, Work, Hess, &MaxLL);

    // output results
    prntBetas(TwoSpeciesLink,NPar, Params, Work, MaxLL, covar);
    
    //   call link function to print re-parameterized estimates
    //SingSpec.ifn=-2; fx=TwoSpeciesLink(Params,NPar);
    
    //  print_individual_estimates(FILE *g, char *lbl, int pgrp, int istrt, int nrealparms, double **covar, double *Params, int linkfn);
    fprintf(g,"\n============================================================\n");
    //strcpy(labl,"PsiA"); 

	print_individual_estimates2(g,0,0,0,1,covar,Params,SingSpec.LinkFn,0,0); // PsiA
    print_individual_estimates2(g,1,0,1,1,covar,Params,SingSpec.LinkFn,0,0);  // psiBA
    //strcpy(tmpstr,SingSpec.Realname[0][2]);
    if (SingSpec.Alt_Param_Checked!=1) print_individual_estimates2(g,2,0,2,1,covar,Params,expLnk,0,0); // psiBa
    else                               print_individual_estimates2(g,2,0,2,1,covar,Params,SingSpec.LinkFn,0,0); // phi
    ii=0;
    //strcpy(labl,"pA"); 
	print_individual_estimates2(g,3,1,0,T,covar,Params,SingSpec.LinkFn,0,0);  // pA
    //strcpy(labl,"pB"); 
	print_individual_estimates2(g,3+T,1,T,T,covar,Params,SingSpec.LinkFn,0,0);   // pB
    //strcpy(labl,"rA"); 
	print_individual_estimates2(g,3+T*2,1,2*T,T,covar,Params,SingSpec.LinkFn,0,0);  // rA
    //strcpy(tmpstr,SingSpec.Realname[1][3*T]); ptr=strstr(tmpstr,"["); if (ptr>0) *ptr='\0';
    print_individual_estimates2(g,3+T*3,1,3*T,T,covar,Params,SingSpec.LinkFn,0,0);
    //if (SingSpec.Alt_Param_Checked==0) {strcpy(labl,"delta"); print_individual_estimates(g,labl,1,-4*T,T,covar,Params,expLnk,0,0);}
    //if (SingSpec.Alt_Param_Checked==1) {strcpy(labl,"rBa");   print_individual_estimates(g,labl,1,-4*T,T,covar,Params,SingSpec.LinkFn,0,0);}
    //if (SingSpec.Alt_Param_Checked==2) {strcpy(labl,"rho");   print_individual_estimates(g,labl,1,-4*T,T,covar,Params,expLnk,0,0);}
    if (SingSpec.Alt_Param_Checked!=1) print_individual_estimates2(g,3+T*4,1,4*T,T,covar,Params,expLnk,0,0);  //  delta
    else                               print_individual_estimates2(g,3+T*4,1,4*T,T,covar,Params,SingSpec.LinkFn,0,0); // rBa
	//
	if (SingSpec.Alt_Param_Checked==1) {     //  if conditional parameterization... psiA,psiBA,psiBa
	                                       //            print derived parameters... phi, dlta
		char *lbl,*fmt2; lbl=new char[80]; fmt2=new char[80]; int site,survey,kk;
		double psiA,psiBA,psiBa,psiAB,psiB,*prd,var,phise,ci1,ci2,phi,phi2,*grad=new double[NPar]; 
		double rA,rBA,rBa,rAB,rB,eps=.1e-10,lastphi; prd=new double[NPar];
		ii=sprintf(lbl,"%%%d.%df",ndec+4,ndec);
		ii=sprintf(fmt2,"%%-12s %%4d %%-16s:%s %s   %s -%s\n",lbl,lbl,lbl,lbl);
		strcpy(lbl,"phi"); phi=-1;
		fprintf(g,"\n    Derived parameter: phi (Species Interaction Factor, SIF = psiAB/(psiA*psiB))\n");
		fprintf(g,"                Site               estimate  Std.err   95%% conf. interval\n");
		for (site=0; site<SingSpec.N; site++) {
			psiA=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
            psiBA=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
            psiBa=getParam(2,0,2,site,0,Params,0,SingSpec.LinkFn); 
            psiAB=psiA*psiBA;                //  Pr(occ by B) = Pr(occ by A)*pr(Occ by B | A present)
            psiB=psiA*psiBA+(1-psiA)*psiBa;  //    + Pr(not occ by A)*Pr(occ by B | A not present)
			lastphi=phi; phi=psiAB/psiA/psiB;
			for (jj=0; jj<NPar; jj++) prd[jj]=0;
			for (jj=0; jj<NPar; jj++) {
				Params[jj]+=eps;
				psiA=getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn);
				psiBA=getParam(1,0,1,site,0,Params,0,SingSpec.LinkFn);
				psiBa=getParam(2,0,2,site,0,Params,0,SingSpec.LinkFn); 
				psiAB=psiA*psiBA;                //  Pr(occ by B) = Pr(occ by A)*pr(Occ by B | A present)
				psiB=psiA*psiBA+(1-psiA)*psiBa;  //    + Pr(not occ by A)*Pr(occ by B | A not present)
				phi2=psiAB/psiA/psiB;
				grad[jj]=(phi2-phi)/eps;
				for (kk=0; kk<NPar; kk++) prd[kk]+=grad[jj]*covar[jj][kk];
				Params[jj]-=eps;
			}
			var=0; for (jj=0; jj<NPar; jj++) var+=prd[jj]*grad[jj];
			phise=sqrt(var); ci1=phi-1.96*phise; ci2=phi+1.96*phise;
	        if (phi!=lastphi) fprintf(g,fmt2,lbl,site+1,SingSpec.sitename[site],
                                phi,phise,ci1,ci2);
			
		}
		fprintf(g,"\n    Derived parameter: dlta (Detection Species Interaction Factor, = rAB/(rA*rB))\n");
		fprintf(g,"                Site               estimate  Std.err   95%% conf. interval\n");
		for (survey=0; survey<T; survey++) {
			ii=sprintf(lbl,"dlta[%d]",survey+1); phi=-1;
			for (site=0; site<SingSpec.N; site++) {
				rA =getParam(3+survey+2*T,1,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
                rBA=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn); 
				rBa=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,0,SingSpec.LinkFn);
                rAB=rA*rBA; rB=rA*rBA+(1-rA)*rBa; lastphi=phi; phi=rAB/rA/rB;
				for (jj=0; jj<NPar; jj++) prd[jj]=0;
				for (jj=0; jj<NPar; jj++) {
					Params[jj]+=eps;
					rA =getParam(3+survey+2*T,1,survey+2*T,site,survey,Params,0,SingSpec.LinkFn);
					rBA=getParam(3+survey+3*T,1,survey+3*T,site,survey,Params,0,SingSpec.LinkFn); 
					rBa=getParam(3+survey+4*T,1,survey+4*T,site,survey,Params,0,SingSpec.LinkFn);
					rAB=rA*rBA; rB=rA*rBA+(1-rA)*rBa; phi2=rAB/rA/rB;
					grad[jj]=(phi2-phi)/eps;
					for (kk=0; kk<NPar; kk++) prd[kk]+=grad[jj]*covar[jj][kk];
					Params[jj]-=eps;
				}
				var=0; for (jj=0; jj<NPar; jj++) var+=prd[jj]*grad[jj];
				phise=sqrt(var); ci1=phi-1.96*phise; ci2=phi+1.96*phise;
				if (phi!=lastphi) fprintf(g,fmt2,lbl,site+1,SingSpec.sitename[site],
                                phi,phise,ci1,ci2);
			}	
		}
		delete[] lbl; delete[] fmt2; delete[] grad; delete[] prd;
    }
    // delete dynamic variables
    delete[] Hess;    delete[] Work; 
    for (ii=0; ii<NPar; ii++) delete[] covar[ii];    delete[] covar;
  }  // end if (SingSpec.NrowsDM[0]>3) .. else
}
