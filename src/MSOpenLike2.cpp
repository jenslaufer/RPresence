#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

//#include "SPO.h"

double MSOpenLike(double ****Phi, double ****ClosedP) {
/*  This routine returns the negative log-likelihood for Multi-State occupancy model
    according to the parameterization used in my notes which differes slightly from
    that in Chap 10 of the book. Passed are the transition probability matrix (Phi)
   and detection probabtility for each survey (ClosedP).
   
   Phi[site][state @ t][state @ t+1][season]
   ClosedP[site][ObsState][TrueState][survey number]
   
   Indexing of states is:
   0 = Unoccupied/Nondetection
   1 = Occupancy 1 (eg species present but not breeding)/ditto detection
   2 = Occupancy 2 (eg species present with breeding)/ditto detection
   
   Global variables declared externally
   SingSpec.N = # sites
   SingSpec.PrmyPeriods = # of seasons
   SingSpec.SecPeriods[seasn] = vector with # of surveys conducted within each season seasn
   SingSpec.Data[site][srvy] integer matrix with 2,1,0 values for occupacny 2, occupancy 1 and
   nondetection for each survey. and srvy keeps track of survey number
   SingSpec.Missing[site][srvy] binary matrix indicating whether survey was NOT
   conducted (1), 0 otherwise.
 */
    extern TSingSpec SingSpec;
	
    int site, t, fstate, tstate, Nstates, srvy, seasn, obsstate, trustate, N=SingSpec.N, maxstate;
    double like,sum=0.,verysmallnumber=SingSpec.nearzero, just_below_one; just_below_one=1-verysmallnumber;
    double *OpenP,*OldTemp,*temp; 

    Nstates=SingSpec.Nstates; OldTemp=new double[Nstates]; temp=new double[Nstates];
    
    // OpenP is equivalent to the diagonalized p vector from MacKenzie et al. (2003)
    //               doesn't need to be a square matrix though
    OpenP = new double[Nstates]; // 3 'states', 2 occupied, and 1 unoccupied


    for (site=0; site<N; site++) { //printf("\nSite=%d\n",site);
        // temp[] used to cary forward the probability of observing the detection history 
        //     while going through the matrix multiplication
        for (tstate=0; tstate<Nstates; tstate++) temp[tstate]=1;
        
        for (srvy=seasn=0; seasn<SingSpec.PrmyPeriods; seasn++) {
            for (fstate=0; fstate<Nstates; fstate++) {     
                OldTemp[fstate] = temp[fstate];	temp[fstate] = 0; OpenP[fstate]=1.0;
            }
            
            for (t=0,maxstate=-1; t<SingSpec.SecPeriods[seasn]; t++) {    // find max state in season
                obsstate=SingSpec.Data[site][srvy+t]; 
				if (obsstate>maxstate) maxstate=obsstate; 
				if (obsstate>=Nstates) { printf("\nerror : obs state(%d) > Nstates(%d) (site %d,srvy %d,t=%d)\n",obsstate,Nstates,site+1,srvy+t+1,t); exit(1); }
			}
			//if (SingSpec.maxfn == 123456. && maxstate==2) maxstate=3; //  A. Hitch data special case!!
			
			for (fstate=0; fstate<Nstates; fstate++)       //  matrix mult. state vector * phi
                for (tstate=0; tstate<Nstates; tstate++) 
                    temp[tstate]+=OldTemp[fstate]*Phi[site][fstate][tstate][seasn]*(tstate>=maxstate);
            
            for (t=0; t<SingSpec.SecPeriods[seasn]; t++, srvy++) {    // calculate OpenP;
                obsstate=SingSpec.Data[site][srvy]; 
                if(obsstate>-1) // if not missing...
                    for (trustate=0; trustate<Nstates; trustate++) { 
                        OpenP[trustate] *= ClosedP[site][obsstate][trustate][srvy]*(trustate>=maxstate);
						if (OpenP[trustate]>1) 
						printf("error:%d %d %d %d %f %f\n",site,srvy,trustate,obsstate,OpenP[trustate],ClosedP[site][obsstate][trustate][srvy]);
					}
            }
            for (fstate=0; fstate<Nstates; fstate++) temp[fstate] *= OpenP[fstate];  // mult. state vector * P
            
        } // end PrmyPeriods loop
        for (fstate=0,like=0; fstate<Nstates; fstate++) like+=temp[fstate];  //  sum state vector to get cell prob (like)
        
        // just do some checking that we haven't got some funky answers
		if (like > -0.0000000001 && like < 0.0) like=0.0;
        if (like <0.0 || like > 1.01) { 
            printf("error in likelihoood for site %d: like=%g\n",site,like);
            printf("p=\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n",
                   ClosedP[0][0][0][0],ClosedP[0][0][1][0],ClosedP[0][0][2][0],ClosedP[0][0][3][0],
                   ClosedP[0][1][0][0],ClosedP[0][1][1][0],ClosedP[0][1][2][0],ClosedP[0][1][3][0],
                   ClosedP[0][2][0][0],ClosedP[0][2][1][0],ClosedP[0][2][2][0],ClosedP[0][2][3][0],
                   ClosedP[0][3][0][0],ClosedP[0][3][1][0],ClosedP[0][3][2][0],ClosedP[0][3][3][0]);
            exit(1); 
        }
        // the log funciton will fail to evaluate if the likelihood value is too small
        if (like > just_below_one) like=just_below_one;
        if (like < verysmallnumber) like=verysmallnumber;
        sum += SingSpec.det_hist_frq[site]*log(like);
		SingSpec.expval[site]=like;
    } // end individual (site) loop
    
    // delete dynamic variables
    delete[] OpenP;	delete temp; delete[] OldTemp;
    
    // return the negative log-likelihood
    if (SingSpec.Verbose>1) printf("XLL=%f\r",-sum);
    return(-sum);   
}
