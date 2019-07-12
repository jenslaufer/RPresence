#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>

double MSOpenLink(double Params[], int npar) {
    
    extern TSingSpec SingSpec;
    double MSOpenLike(double ****fi, double ****CP);
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int i, site, FROMstate, TOstate, OBSstate, TRUstate, srvy, seasn, pstart, nsrvys, rstart, row,j,k,ns1,
          nseasns = SingSpec.PrmyPeriods, Nstates=SingSpec.Nstates, linkFn,*fixed;
    double sum, sum2, psi, R, p1, p2, de, ****Phi, ****ClosedP, xll, eta[2], gameps[4], pi, psiA, psiB, sumfix;
    
	if (SingSpec.Alt_Param_Checked==1) Nstates=3;
    fixed=new int[Nstates];
    linkFn=SingSpec.LinkFn;
    nsrvys=SingSpec.T; 
    Phi=new double***[SingSpec.N]; ClosedP=new double***[SingSpec.N];
    for (site=0; site<SingSpec.N; site++) {
        Phi[site]=new double**[Nstates]; ClosedP[site]=new double**[Nstates];
        for (FROMstate=0; FROMstate<Nstates; FROMstate++) {
            Phi[site][FROMstate]=new double*[Nstates]; 
            ClosedP[site][FROMstate]=new double*[Nstates];
            for (TOstate=0; TOstate<Nstates; TOstate++) {
                Phi[site][FROMstate][TOstate]=new double[nsrvys]; 
                ClosedP[site][FROMstate][TOstate]=new double[nsrvys];
                for (srvy=0; srvy<nsrvys; srvy++) {
                    Phi[site][FROMstate][TOstate][srvy] = 0; 
                    ClosedP[site][FROMstate][TOstate][srvy] = 0;
                }
            }
        }
    }
    ns1=nseasns-1; rstart=1+ns1*3;
    for (site=0; site<SingSpec.N; site++) {
        //     For psi-R parameterization...
        //        1st season, only compute init psi, and init R
        //        other seasons, compute conditional psi for each state and cond. R.
        switch (SingSpec.Alt_Param_Checked) {
        case 1:                                 //  psi,R parameterization
            psi=getParam(0,0,0,site,0,Params,0,linkFn);
            R=getParam(rstart,0,rstart,site,0,Params,0,linkFn);
            Phi[site][0][1][0] = psi*(1-R);
            Phi[site][0][2][0] = psi*R; sum=psi;
            break;
        case 2:                                 //   integrated habitat model parameterization
            pi=getParam(0,0,0,site,0,Params,0,linkFn);
            psiA=getParam(1,0,1,site,0,Params,0,linkFn);
            psiB=getParam(2,0,2,site,0,Params,0,linkFn);
            Phi[site][0][1][0] = pi*psiA;                  // state=1 -> hab=A,occupied
            Phi[site][0][2][0] = (1-pi)*(1-psiB);          // state=2 -> hab=B,unoccupied
            Phi[site][0][3][0] = (1-pi)*psiB;              // state=3 -> hab=B,occupied
            sum=Phi[site][0][1][0]+Phi[site][0][2][0]+Phi[site][0][3][0]; // state=0 -> hab=A,unoccupied
            break;
        case 3:                                 //   general parmeterization with conditional probs
            for (TOstate=1,sum=0; TOstate<Nstates; TOstate++) {
                psi=getParam(TOstate-1,0,TOstate-1,site,0,Params,0,linkFn)*(1-sum);
                Phi[site][0][TOstate][0]=psi;    //  psi = prob(unocc at start to state 'TOstate')
                sum+=psi;
            }
            break;
        default:                                //   general parameterization
            for (TOstate=1,sum=sum2=0; TOstate<Nstates; TOstate++) {
                psi=getParam(TOstate-1,0,TOstate-1,site,0,Params,0,expLnk);
                Phi[site][0][TOstate][0]=psi;    //  psi = prob(unocc at start to state 'TOstate')
                sum2+=psi;
            }
            for (TOstate=1,sum=0; TOstate<Nstates; TOstate++) {
                Phi[site][0][TOstate][0]=Phi[site][0][TOstate][0]/(1+sum2);
                sum+=Phi[site][0][TOstate][0];
				SingSpec.realParmEst[TOstate-1][site]=Phi[site][0][TOstate][0];
            }
        }
        Phi[site][0][0][0]=1-sum; 
        for (seasn=1; seasn<nseasns; seasn++) { 
            for (FROMstate=0; FROMstate<Nstates; FROMstate++) {
                switch (SingSpec.Alt_Param_Checked) {
                case 1:
                    i=1+(seasn-1)*3+FROMstate; 
                    psi = getParam(i,0,i,site,seasn-1,Params,0,linkFn);
                    R   = getParam(i+rstart,0,i+rstart,site,seasn-1,Params,0,linkFn);
                    Phi[site][FROMstate][1][seasn] = psi*(1-R);
                    Phi[site][FROMstate][2][seasn] = psi*R; sum=psi;
                    break;
                case 2:    //  subscripting...   assume 4 states (AaBb) for now...
                    
                    /*  transition matrix  
                       
                          state=0    |   state=1     |    state=2        |   state=3   
                       not occ in A  |   occ in A    |  not occ in B     |   occ in B     
                       |-------------|---------------|-------------------|-----------------|
           not occ in A| N0AA*(1-gAA)|  N0AA*gAA     |  (1-N0AA)*(1-gAB) | (1-N0AA)*gAB    |
             occ in A  | N1AA*eAA    |  N1AA*(1-eAA) |  (1-N1AA)*eAB     | (1-N1AA)*(1-eAB)|
           not occ in B| N0BA*(1-gBA)|  N0BA*gBA     |  (1-N0BA)*(1-gBB) | (1-N0BA)*gBB    |
             occ in B  | N1BA*eBA    |  N1BA*(1-eBA) |  (1-N1BA)*eBB     | (1-N1BA)*(1-eBB)|
                       |-------------|---------------|-------------------|-----------------|
                       Where N0AA=eta(0,AA)=Pr(transiton A-A | not occupied ) */
                    
                    j=FROMstate%2; k=FROMstate/2; 
                    row=seasn-1+(nseasns-1)*FROMstate;  i=SingSpec.NrowsDM[0]+row; 
                    eta[0]=getParam(i,1,row,site,seasn-1,Params,0,linkFn); eta[1]=1-eta[0];
                    if (j==0) {   //   even rows (rows 0,2,4,...), use 1-gamma, gamma
                        row=k*2*ns1+seasn-1; i=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+row;
                        gameps[1]=getParam(i,2,row,site,seasn-1,Params,0,linkFn);   //  gamma(xA)
                        gameps[0]=1-gameps[1];
                        gameps[3]=getParam(i+ns1,2,row+ns1,site,seasn-1,Params,0,linkFn); //  gamma(xB)
                        gameps[2]=1-gameps[3];
                    }
                    else  {              //    odd rows, use epsilon, 1-epsilon
                        row=k*2*ns1+seasn-1; i=SingSpec.NrowsDM[0]+SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+row;
                        gameps[0]=getParam(i,3,row,site,seasn-1,Params,0,linkFn);   //  eps(xA)
                        gameps[1]=1-gameps[0];
                        gameps[2]=getParam(i+ns1,3,row+ns1,site,seasn-1,Params,0,linkFn); //  eps(xB)
                        gameps[3]=1-gameps[2];
                    }
                    for (TOstate=1,sum=0; TOstate<Nstates; TOstate++) {
                        Phi[site][FROMstate][TOstate][seasn]=eta[TOstate/2]*gameps[TOstate];	
                        sum+=Phi[site][FROMstate][TOstate][seasn];
                    }
                    break;
                case 3:    //  general transitions  with conditional probs
                    for (TOstate=1,sum=0; TOstate<Nstates; TOstate++) {
                        i=Nstates+FROMstate*(Nstates-1)*ns1+(TOstate-1)*ns1+seasn-2;
                        psi = getParam(i,0,i,site,seasn-1,Params,0,linkFn)*(1.-sum);
                        Phi[site][FROMstate][TOstate][seasn]=psi;
                        sum+=psi;
                    }
                    break;
                default:
                    for (TOstate=1,sum=sum2=0; TOstate<Nstates; TOstate++) {
                        i=Nstates+FROMstate*(Nstates-1)*ns1+(TOstate-1)*ns1+seasn-2;
                        psi = getParam(i,0,i,site,seasn-1,Params,0,expLnk); 
                        Phi[site][FROMstate][TOstate][seasn]=psi;
                        sum2+=psi;
                    }
                    for (TOstate=1,sum=0; TOstate<Nstates; TOstate++) {
                        Phi[site][FROMstate][TOstate][seasn]=Phi[site][FROMstate][TOstate][seasn]/(1+sum2);
                        sum+=Phi[site][FROMstate][TOstate][seasn];
                        i=Nstates+FROMstate*(Nstates-1)*ns1+(TOstate-1)*ns1+seasn-2;						
						SingSpec.realParmEst[i][site]=Phi[site][FROMstate][TOstate][seasn];
                    }
                }
                Phi[site][FROMstate][0][seasn] = 1 - sum;
            }
        }
    }
    // set up parameters for detection,    ClosedP[site][ObsState][TrueState][survey]
    pstart=Nstates-1+Nstates*(Nstates-1)*ns1;
    if (SingSpec.Alt_Param_Checked) pstart=2+6*ns1;
    for (site=0; site<SingSpec.N; site++) {
        for (srvy=0; srvy<nsrvys; srvy++) { 
            ClosedP[site][0][0][srvy] = 1;   //  if true state=0(unocc), then obs state must be zero
                                             //  special parameterization w/ 3 states, using p1,p2,delta
            switch (SingSpec.Alt_Param_Checked) {
            case 1:
                de=getParam(pstart+srvy,4,srvy,site,srvy,Params,0,linkFn);     // delta
                i=pstart+srvy+nsrvys;
                p1=getParam(i,       4,i-pstart,       site,srvy,Params,0,linkFn);
                p2=getParam(i+nsrvys,4,i-pstart+nsrvys,site,srvy,Params,0,linkFn);
                ClosedP[site][0][1][srvy] = 1 - p1;      
                ClosedP[site][1][1][srvy] = p1;          
                //ClosedP[site][2][1][srvy] = 0;           // can't observe a 2 if true state=1
                ClosedP[site][0][2][srvy] = 1-p2;        
                ClosedP[site][1][2][srvy] = p2*(1-de); 
                ClosedP[site][2][2][srvy] = p2*de; 
                break;
            case 2:                               //  integrated habitat/occupancy model
                row=srvy; i=row+ns1*12+3;
                p1=getParam(i,       4,row,       site,srvy,Params,0,linkFn);
                p2=getParam(i+nsrvys,4,row+nsrvys,site,srvy,Params,0,linkFn);
                ClosedP[site][0][0][srvy] = 1;      // obs state=0, tru state=0 -> hab=A, not detected
                ClosedP[site][0][1][srvy] = 1-p1;      // obs state=A, tru state=a -> hab=A, 
                //ClosedP[site][0][2][srvy] = 0;      // can't observe a 0 if true state=2
                //ClosedP[site][0][3][srvy] = 0;      // can't observe a 0 if true state=3
                ClosedP[site][1][0][srvy] = 0;        // obs state=a, tru state=A
                ClosedP[site][1][1][srvy] = p1;        // obs state=a, tru state=a          
                //ClosedP[site][2][0][srvy] = 0;           // can't observe...  whole ClosedP matrix is
                //ClosedP[site][2][1][srvy] = 0;           // can't observe...     initialized to zero above
                ClosedP[site][2][2][srvy] = 1;           //  obs state=B, tru state=B
                ClosedP[site][2][3][srvy] = 1-p2;           //  obs state=B, tru state=b
                ClosedP[site][3][2][srvy] = 0;            //  obs state=b, tru state=B
                ClosedP[site][3][3][srvy] = p2;            //  obs state=b, tru state=b
                if (SingSpec.Data[site][srvy]>3) {    //   if det.hist>3 then missing habitat type
                    ClosedP[site][0][2][srvy] = 1;     
                    ClosedP[site][0][3][srvy] = 1-p2;  
                    ClosedP[site][2][0][srvy] = 1;     
                    ClosedP[site][2][1][srvy] = 1-p1;  
                }
                break;
            case 3:
                for (TRUstate=0; TRUstate<Nstates; TRUstate++) {  //parm. order- p01(srvy),p02(srvy),...p0k(srvy),
                    //             p11,p12,...p1k,
                    for (OBSstate=1,sum=0; OBSstate<Nstates; OBSstate++) {//     p21,p22,...p2k,
                        i=pstart+srvy+(OBSstate-1)*nsrvys+TRUstate*(Nstates-1)*nsrvys; //     conditional probs,
                        fixed[OBSstate]=0;                                                   //  use logit transf
                        if (SingSpec.fixed[i]>-998) { fixed[OBSstate]=1; p1=SingSpec.fixed[i];}
                        else { p1=getParam(i,4,i-pstart,site,srvy,Params,0,linkFn)*(1.-sum); }
                        ClosedP[site][OBSstate][TRUstate][srvy]=p1; 
                        sum+=p1;  
                    }
                    ClosedP[site][0][TRUstate][srvy]=1-sum;
                }
            default:
                for (TRUstate=1; TRUstate<Nstates; TRUstate++) {  //parm. order- p11(srvy),p12(srvy),...p1k(srvy),
                    //             p11,p12,...p1k,
                    for (OBSstate=1,sum2=sumfix=1; OBSstate<Nstates; OBSstate++) {//     p21,p22,...p2k,
                        i=pstart+srvy+(OBSstate-1)*nsrvys+(TRUstate-1)*(Nstates-1)*nsrvys; //       :
                        fixed[OBSstate]=0;                                                   //  use mlogit transformation
                        if (SingSpec.fixed[i]>-998) { fixed[OBSstate]=1; p1=SingSpec.fixed[i]; sumfix-=p1;}
                        else { p1=getParam(i,4,i-pstart,site,srvy,Params,0,expLnk); sum2+=p1; }
                        ClosedP[site][OBSstate][TRUstate][srvy]=p1;					
                    }
					ClosedP[site][0][TRUstate][srvy]=1; 
                    for (OBSstate=1,sum=0; OBSstate<Nstates; OBSstate++) {  //  sum of cols must be <= (1.0-sum of fixed values)
                        i=pstart+srvy+(OBSstate-1)*nsrvys+(TRUstate-1)*(Nstates-1)*nsrvys; 
						if (fixed[OBSstate]==0)   //  if not fixed
                            ClosedP[site][OBSstate][TRUstate][srvy]=sumfix*ClosedP[site][OBSstate][TRUstate][srvy]/sum2;
						SingSpec.realParmEst[i][site]=ClosedP[site][OBSstate][TRUstate][srvy];
						ClosedP[site][0][TRUstate][srvy]-=ClosedP[site][OBSstate][TRUstate][srvy];
						if (ClosedP[site][0][TRUstate][srvy]<0) ClosedP[site][0][TRUstate][srvy]=0;				
                    }                    
                }
            }
        } // end of survey (srvy) loop
    }    	
    xll=MSOpenLike(Phi,ClosedP);
    
    for (site=0; site<SingSpec.N; site++) {
        for (FROMstate=0; FROMstate<Nstates; FROMstate++) {
            for (TOstate=0; TOstate<Nstates; TOstate++) { 
                delete [] Phi[site][FROMstate][TOstate]; 
                delete [] ClosedP[site][FROMstate][TOstate];
            }
            delete [] Phi[site][FROMstate]; delete [] ClosedP[site][FROMstate];
        }
        delete [] Phi[site]; delete [] ClosedP[site];
    }
    delete [] Phi; delete [] ClosedP; delete[] fixed;
    
    return(xll);  // call OpenLike
}
