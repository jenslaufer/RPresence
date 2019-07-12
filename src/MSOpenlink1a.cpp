#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>

double MSOpenLink1a(double Params[], int npar) {  // special case of multi-state psi/R
                                                  // parameterization w/ 1 survey per season
    extern TSingSpec SingSpec;
    double getParam(int param_num, int param_grp, int row, int site, int seasn, double *Params, int transf, int LinkFn);
    int i=0, site, FROMstate, seasn, ns1,obsstate,
          nseasns = SingSpec.PrmyPeriods, Nstates=SingSpec.Nstates, linkFn, dstart;
    double sum=0, psi, R, p1, p2, de, Phi[3][3], ClosedP[3][3], 
	     s0,s1,s2,t0,t1,t2,like,OpenP0,OpenP1,OpenP2,just_below_one=1-1e-12,verysmallnumber=1e-28;
    linkFn=SingSpec.LinkFn; nseasns=SingSpec.T; ns1=nseasns-1; 
	ClosedP[1][0]=ClosedP[2][1]=ClosedP[2][0]=0; ClosedP[0][0] = 1;   //  if true state=0(unocc), then obs state must be zero
    for (site=0; site<SingSpec.N; site++) {
        //     For psi-R parameterization...    1st season, only compute init psi, and init R
        //        other seasons, compute conditional psi for each state and cond. R.
		// temp[t0,t1,t2] used to cary forward the probability of observing the detection history 
        //     while going through the matrix multiplication
        t0=1; t1=t2=0; i=0;
        psi=getParam(i,  0,i,  site,0,Params,0,linkFn);       //  initial psi
        R  =getParam(i+1,0,i+1,site,0,Params,0,linkFn);  //  initial R
        Phi[0][0]=1-psi; Phi[0][1]=psi*(1-R); Phi[0][2]=psi*R; Phi[1][0]=Phi[1][1]=Phi[1][2]=Phi[2][0]=Phi[2][1]=Phi[2][2]=0;
        for (seasn=0; seasn<nseasns; seasn++) { 
		  if (seasn>0)
            for (FROMstate=0; FROMstate<Nstates; FROMstate++) {
                i=2 + FROMstate*ns1 + seasn-1;
                psi = getParam(i,0,i,site,seasn-1,Params,0,linkFn);         //  conditional psi's
                i=2 + (3+FROMstate)*ns1 + seasn-1;
                R   = getParam(i,0,i,site,seasn-1,Params,0,linkFn); // conditional R's
                Phi[FROMstate][0] = 1-psi; Phi[FROMstate][1] = psi*(1-R); Phi[FROMstate][2] = psi*R; 
            }
		    // set up parameters for detection,    ClosedP[ObsState][TrueState]
		  dstart=2+2*Nstates*ns1; i=seasn; 
          de=getParam(dstart          +i, 4,           i, site,seasn,Params,0,linkFn); // delta
          p1=getParam(dstart+  nseasns+i, 4,   nseasns+i, site,seasn,Params,0,linkFn);
          p2=getParam(dstart+2*nseasns+i, 4, 2*nseasns+i, site,seasn,Params,0,linkFn);
          ClosedP[0][1]=1-p1; ClosedP[1][1]=p1; ClosedP[0][2]=1-p2; ClosedP[1][2]=p2*(1-de); ClosedP[2][2]=p2*de; 
		  OpenP0=OpenP1=OpenP2=1; obsstate=SingSpec.Data[site][seasn]; 
		  s0=t0; s1=t1; s2=t2;         //  matrix mult. state vector * phi
		  t0=s0*Phi[0][0]*(0>=obsstate)+ s1*Phi[1][0]*(0>=obsstate)+ s2*Phi[2][0]*(0>=obsstate);
		  t1=s0*Phi[0][1]*(1>=obsstate)+ s1*Phi[1][1]*(1>=obsstate)+ s2*Phi[2][1]*(1>=obsstate);
		  t2=s0*Phi[0][2]*(2>=obsstate)+ s1*Phi[1][2]*(2>=obsstate)+ s2*Phi[2][2]*(2>=obsstate);
          if(obsstate>(-1)) { // if not missing...
             OpenP0=ClosedP[obsstate][0]*(0>=obsstate);  OpenP1=ClosedP[obsstate][1]*(1>=obsstate);
             OpenP2=ClosedP[obsstate][2]*(2>=obsstate);
          }
          t0 *= OpenP0; t1 *= OpenP1; t2 *= OpenP2;  // mult. state vector * P
		  //if (site==165) printf("%d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f | %g %g %g\n",grp,seasn,i,dstart,p1start,p2start,obsstate,
		  //Phi[0][0],Phi[0][1],Phi[0][2],Phi[1][0],Phi[1][1],Phi[1][2],OpenP0,OpenP1,OpenP2,t0,t1,t2);        
	    } // end PrmyPeriods loop
        like=t0+t1+t2;  //  sum state vector to get cell prob (like)
	    //printf("%d %d%d%d%d %g %g %g %f %f %12.8f\n",site,SingSpec.Data[site][0],SingSpec.Data[site][1],SingSpec.Data[site][2],
		//     SingSpec.Data[site][3],t0,t1,t2,like,SingSpec.det_hist_frq[site]/1000.,like-SingSpec.det_hist_frq[site]/1000.);
		if (like > -0.0000000001 && like < 0.0) like=0.0;  // just do some checking that we haven't got some funky answers
        if (like <0.0 || like > 1.01) { 
            printf("error in likelihoood for site %d: like=%g t0=%g t1=%g t2=%g\n",site,like,t0,t1,t2); exit(1); 
        }
        // the log funciton will fail to evaluate if the likelihood value is too small
        if (like > just_below_one) like=just_below_one; 
		if (like < verysmallnumber) like=verysmallnumber;
        sum += SingSpec.det_hist_frq[site]*log(like); 
		SingSpec.expval[site]=like;
		//sum2 += (SingSpec.det_hist_frq[site]/1000.-like)*(SingSpec.det_hist_frq[site]/1000.-like);
    } // end individual (site) loop
    if (SingSpec.Verbose>1) printf("XLL=%f    %d   %f %f %f %f %f %f\n",-sum,SingSpec.ifn,Params[0],Params[1],Params[2],Params[3],Params[4],Params[5]);
	//exit(1);
    return(-sum);  // return the negative log-likelihood
}
