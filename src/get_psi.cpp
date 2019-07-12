#include "SingAnal.h"
#include <math.h>
#include <stdlib.h>
#define PLIMIT 20.0

double invlogit(double x) { return( (x > -60. ? (x < 60. ? 1./(1.+exp(-x)) : 1.) : 0.)  );}
double logit(double x) {    return( (x > 0. ? (x < 1. ? log(x/(1-x)) : 40.) : -40.)  );}

double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params,int transf, int LinkFn) {
    
    //   get 'real' parameter using beta's and design matrices
    
    extern TSingSpec SingSpec;  
	int i,j,dmptr, kk,missng=0,indexPar=0; double rslt,sum=0,dm,xn=(double)SingSpec.N;
    rslt = SingSpec.realParmEst[param_num][site] = SingSpec.fixed[param_num]; xn=1.;
    //if (param_num==10 && site<10)  
	//	printf("getParam: param_num=%d grp=%d row=%d site=%d srvy=%d %d %d\n",param_num,param_grp,row,site,srvy,transf,LinkFn);
   
    if (rslt>-998) return(rslt);
	
    //    count number of parms used in preceeding design matrices and store in 'IndexPar'
    for (j=0; j<param_grp; j++) 
        for (kk=0; kk<SingSpec.NParKK[j]; kk++) {
            if (SingSpec.BetaFixed[j][kk]>1.0e44) indexPar++;  //  if beta(j) not fixed, inc indexpar
        }
    if (SingSpec.NParKK[1]+SingSpec.NParKK[2]+SingSpec.NParKK[3]+SingSpec.NParKK[4]==0) { // if only 1 des.mat.
        param_grp=0; row=param_num; indexPar=0;                                             //  set index to 0.
    }
    if (row<0) 	for (i=0,row=param_num; i<param_grp; i++) row-=SingSpec.NrowsDM[i];
    if (row<0) { fprintf(SingSpec.g,"\n\n***** error in getParam (row<0) ***\n\n"); exit(1);}
    
    if (SingSpec.LnkFn[param_grp][row]>0) LinkFn=SingSpec.LnkFn[param_grp][row];
    	
    if (SingSpec.fixed[param_num]<-998.) {
        for (kk=0; kk<SingSpec.NParKK[param_grp]; kk++) { 
            if (row>=SingSpec.NrowsDM[param_grp]) { 
                fprintf(SingSpec.g,"**** error in getParam (parnum=%d mat=%d,row=%d,nrows=%d) ***\n",
					param_num,param_grp,row+1,SingSpec.NrowsDM[param_grp]); 
                fclose(SingSpec.g); 
                exit(1); 
            }
            dm=SingSpec.DMat[param_grp][row][kk]; dmptr=SingSpec.DMat_ptr[param_grp][row][kk];
            if(dmptr>3000) { 
			  dm=SingSpec.psibar[srvy][site][dmptr-3001]/xn; 
			  //if (site==1) printf(">>getParam(param_num=%d param_grp=%d dmptr=%d psibar=%f xn=%f dm=%f\n",
			  //param_num,param_grp,dmptr,SingSpec.psibar[srvy][site][dmptr-3001],xn,dm);
			} //  3001=>psiA, 3002=>psiB
            else if(dmptr<0) { dm=SingSpec.SampCov[-dmptr-1][site][srvy];	}
            else if(dmptr>0) { dm=SingSpec.SiteCov[site][dmptr-1];
			}
            if(dmptr<-9990)missng=-1;
			//if (site==1 && param_num==45) printf("sum=%f dm=%f par=%f\n",sum,dm,Params[indexPar]);
	
            if (SingSpec.BetaFixed[param_grp][kk]>1.0e44) sum += dm*Params[indexPar++];
            else                                          sum += dm*SingSpec.BetaFixed[param_grp][kk];
        }   //  end for (kk=0 ...	
		
        if (transf==0) {
            switch (LinkFn) {
            case 0: 
            case logitLnk: rslt = invlogit(sum); break;
            case logitLnk2: rslt = invlogit(sum)/2.; break;
            case loglogLnk: rslt = (sum>-60 ? (sum<60 ? 1-exp(-exp(sum)) : 1-1.0e-50) : 1.0e-50); break;
            case expLnk: rslt = (sum<50. ? exp(sum) : exp(50.)); break;
            case sinLnk: rslt = (sin(sum)+1.)/2.; break;
            case sinLnk2: rslt = (sin(sum)+1.)/4.; break;
            case IDLnk: rslt = sum; break;
            default: rslt = (sum>-60 ? (sum<60 ? 1/(1+exp(-sum)) : 1) : 0);
            }
            if (rslt == 0.0) rslt = 1.0e-50;  // asdf make value something other than 0 or 1
            if (rslt == 1.0) rslt -= 1.0e-50;
            if (LinkFn==IDLnk) rslt=sum;
        }
        else rslt=exp(sum);
    }
	//if (site==1 && param_num==45) printf("sum=%f rslt=%f\n",sum,rslt);
    if(missng<0) rslt=-1;
	if (SingSpec.optmiz_done==1) SingSpec.realParmEst[param_num][site]=rslt;
	if (SingSpec.Verbose>2) printf("site=%d getParam(%d) = %f\n",site,param_num,rslt);
    return rslt;
}
