#include <math.h>
#include "SingAnal.h"
#include <stdio.h>
#include <stdlib.h>

double OpenLikeMM(double ****Phi, double **theta, double **ClosedP) { 
    int site, seasn, srvy, nlog0, missng,seenInSeasn,seenInSrvy,PrmyPers,tsrvy,msrvy,imeth; 
    double prdp[4],prdLocOcc[3],sum=0,prdOcc,prdNotOcc,prdOcc_save,prdNotOcc_save,cellProb=0,pstar;
    extern TSingSpec SingSpec; PrmyPers=SingSpec.PrmyPeriods;

    for (site=nlog0=0; site<SingSpec.N; site++) {  cellProb=0;
        
        prdOcc = Phi[site][0][0][0]; prdNotOcc = 1-prdOcc;
        //for (srvy=0, lastobs=-1; srvy<SingSpec.T; srvy++) if (SingSpec.Data[site][srvy]!=-1) lastobs=srvy+1;
		
        for (seasn=tsrvy=msrvy=0; seasn<PrmyPers; seasn++) {            // prb=psi  *  theta  *  p  * z  * phi  *  theta ...
                                                                        //       1x2       2x3   3x3  3x2    2x2       2x3
            for (srvy=missng=seenInSeasn=0; srvy<SingSpec.SecPeriods[seasn]; srvy++,tsrvy++) {
                pstar=1;                                                        //         Occ  UnOcc
				for (imeth=seenInSrvy=0; imeth<SingSpec.NMethods; imeth++,msrvy++) {  //  psi = [psi  1-psi]
                                                                                      //  
					if(SingSpec.Data[site][msrvy]!=-1) {  // calculate OpenP;
						                                               //              h* = history for srvy j
						if(SingSpec.Data[site][msrvy]==1) {                   //        h*>0    h*=0  unocc
							                                                 //     LO[ p*       q*   0    ] 
							pstar *= ClosedP[site][msrvy];              // p = LU[ 0        1    0    ] 
							seenInSeasn=seenInSrvy=1;                        //     Un[ 0        0    1    ]
						}
						if(SingSpec.Data[site][msrvy]==0) pstar *= (1-ClosedP[site][msrvy]);  //       Occ UnOcc
					}    //  end of  if not missing                                               z = h*>0[ 1    0   ] 
				}   //  end of imeth loop                                                             h*=0[ 1    0   ]
				                                                                     //              unocc[ 0    1   ]
				prdLocOcc[0]=prdOcc*theta[site][tsrvy];                // psi * theta           LO    LU     UnOcc
				prdLocOcc[1]=prdOcc*(1-theta[site][tsrvy])*(1-seenInSrvy);   //  theta = Occ  [theta 1-theta  0]
				prdLocOcc[2]=prdNotOcc;                      //                             UnOcc[   0      0    1] 
				
				prdp[0]=prdLocOcc[0]*pstar;                      //  psi * theta * p
				prdp[1]=prdLocOcc[1];
				prdp[2]=prdLocOcc[2];

				prdOcc    = prdp[0]+prdp[1];                   //  psi * theta * p * z
				prdNotOcc = prdp[2];
				
            }   //  end of survey in season loop
			prdNotOcc *= (1-seenInSeasn);
            prdOcc_save=prdOcc;  prdNotOcc_save=prdNotOcc; // save prds so they can be used below
            if (seasn<(PrmyPers-1)) {  // psi*theta*p*z*phi                                                             Occ  UnOcc
                prdOcc = prdOcc_save*Phi[site][0][0][seasn+1] + prdNotOcc_save*Phi[site][1][0][seasn+1];  // phi=  Occ[1-eps  eps ]
                prdNotOcc = prdOcc_save*Phi[site][0][1][seasn+1] + prdNotOcc_save*Phi[site][1][1][seasn+1];//    UnOcc[ gam  1-gam]
            }
			
        } // end PrmyPeriods loop
        cellProb += (prdOcc+prdNotOcc); 
        if (cellProb < SingSpec.nearzero) { nlog0++; cellProb=SingSpec.nearzero; }
        sum += SingSpec.det_hist_frq[site]*log(cellProb); 
		SingSpec.expval[site]=cellProb; 
    } // end individual loop
    return(-sum);
}
