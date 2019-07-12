#include "SingAnal.h"
#include <math.h>
#include <stdlib.h>
#define II -85
double MS1CorrLike(double *Params, int NPar) {

    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int ii, jj, seen, k1, kp, jp, T, N, i, j;
    double psi1,  // =Pr(patch i is occupied in state 1)
	       psi2,  //=Pr(patch i is occupied in state 2)
		   psi3,  //=Pr(patch i is occupied in state 3)
		   pi11,  // =Pr(site is initially in local occ state 1 | global occ state=1)
		   pi12,  // =Pr(site is initially in local occ state 1 | global occ state=2)
		   pi22,  // =Pr(site is initially in local occ state 2 | global occ state=2)
		   pi13,  // =Pr(site is initially in local occ state 1 | global occ state=3)
		   pi23,  // =Pr(site is initially in local occ state 2 | global occ state=3)
		   pi33,  // =Pr(site is initially in local occ state 3 | global occ state=3)
		   th1[2][2],  // =Pr(transition from local occ state r to s | global occ state=1)
		   th2[3][3],  // =Pr(transition from local occ state r to s | global occ state=2)
		   th3[4][4],  // =Pr(transition from local occ state r to s | global occ state=3)
		   p11,        // =Pr(det in state 1 | true local occ state=1)
		   p12,        // =Pr(det in state 1 | true local occ state=2)
		   p13,        // =Pr(det in state 1 | true local occ state=3)
		   p22,        // =Pr(det in state 2 | true local occ state=2)
		   p23,        // =Pr(det in state 2 | true local occ state=3)
		   p33,        // =Pr(det in state 3 | true local occ state=3)
	 sum, sumll=0, sumchi=0, sumfrq=0, like=0,like1[2]={0,0},like2[3]={0,0,0}, 
	 like3[4]={0,0,0,0}, v[4]={0,0,0,0};
    extern TSingSpec SingSpec; T=SingSpec.T; N=SingSpec.N;
/*      * computes estimates under multi-strata-correlated-detections PRESENCE model with 3 strata 
	 
	 initial state:  X = [ X0   X1   X2   X3 ]    X0=1-psi1-psi2-psi3
	 
	        X1: Occupancy state 1                    X2: Occupancy state 2
	           0             1                   0                   1            2
		  ---                    ---        ---                                       ---
	     | psi1*(1-pi11)  psi1*pi11 |      | psi2*(1-pi12-pi22)   psi2*pi12   psi2*pi22  |
	      ---                    ---        ---                                       ---  
	         X3: Occupancy state 3
	                     0                1           2              3
		  ---                                                            ---
	     |  psi3*(1-pi13-pi23-pi33)   psi3*pi13   psi3*pi23      psi3*pi33  |
	      ---                    ---        ---                          ---  
	 
	 possible transitions: 
	 
	 Occpancy state 1:          Occupancy state 2:           Occupancy state 3:
	        0     1                  0     1      2           0      1      2      3
		---         ---      ---                 ---      ---                        ---
	 0 | th1.00 th1.01 |  0 | th2.00 th2.01 th2.02  |  0 | th3.00 th3.01 th3.02 th3.03  |
	 1 | th1.10 th1.11 |  1 | th2.10 th2.11 th2.12  |  1 | th3.00 th3.01 th3.02 th3.13  | 
        ---         ---   2 | th2.20 th2.21 th2.22  |  2 | th3.00 th3.01 th3.02 th3.23  |
	                         ---                 ---   3 | th3.00 th3.01 th3.02 th3.33  |
	                                                      ---                        ---   
     *  note:  th1.01 = theta, th1.11 = theta' in single-state model
	 *  also:  thX.r0 = 1 - sum(thX.ri)  i=1..max state
            
          Detection probs:
		     True Local Occupancy state
	Obs        0           1      2     3
   state  ---                          ---
	 0   |     1           0      0     0 |  
	 1   |   1-p11       p11      0     0 |  
	 2   | 1-p12-p22     p12    p22     0 |  
	 3   | 1-p13-p23-p33 p13    p23   p33 |  
	      ---                          ---   
	 */
    // go through data and evaluate likelihood
    for (ii=0, sum=0.; ii<N; ii++) {
        psi1=getParam(0,0,0,ii,0,Params,0,expLnk);
        psi2=getParam(1,0,1,ii,0,Params,0,expLnk);
        psi3=getParam(2,0,2,ii,0,Params,0,expLnk);
		sum=1+psi1+psi2+psi3; psi1/=sum; psi2/=sum; psi3/=sum; 
		
		kp=26*T+3;
        pi11=getParam(kp,  2,0,ii,0,Params,0,SingSpec.LinkFn);
        pi12=getParam(kp+1,2,1,ii,0,Params,0,expLnk);
        pi22=getParam(kp+2,2,2,ii,0,Params,0,expLnk);
		sum=1+pi12+pi22; pi12/=sum; pi22/=sum; 
        pi13=getParam(kp+3,2,3,ii,0,Params,0,expLnk);
        pi23=getParam(kp+4,2,4,ii,0,Params,0,expLnk);
        pi33=getParam(kp+5,2,5,ii,0,Params,0,expLnk);
		sum=1+pi13+pi23+pi33; pi13/=sum; pi23/=sum; pi33/=sum; 
		if(ii==II) printf("psi:%f %f %f pi:%f %f %f %f %f %f\n",psi1,psi2,psi3,pi11,pi12,pi22,pi13,pi23,pi33);
		
        like1[1]=psi1*pi11; like1[0]=psi1*(1-pi11);
		like2[2]=psi2*pi22; like2[1]=psi2*pi12; like2[0]=psi2*(1-pi22-pi12);
		like3[3]=psi3*pi33; like3[2]=psi3*pi23; like3[1]=psi3*pi13; like3[0]=psi3*(1-pi33-pi23-pi13);
		if(ii==II)printf("like3=%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",like3[0],like3[1],like3[2],like3[3],like2[0],like2[1],like2[0],like1[1],like1[0]);
		seen=0;
        for (jj=0; jj<T; jj++) {
		    k1=3+jj; kp=20*T+k1; jp=jj;
            th1[0][1]=getParam(k1,0,k1,ii,jj,Params,0,SingSpec.LinkFn); th1[0][0]=1-th1[0][1]; k1+=T;
            th1[1][1]=getParam(k1,0,k1,ii,jj,Params,0,SingSpec.LinkFn); th1[1][0]=1-th1[1][1]; k1+=T;
            th2[0][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th2[0][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th2[0][1]+th2[0][2]; th2[0][1]/=sum; th2[0][2]/=sum; th2[0][0]=1-th2[0][1]-th2[0][2];
            th2[1][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th2[1][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th2[1][1]+th2[1][2]; th2[1][1]/=sum; th2[1][2]/=sum; th2[1][0]=1-th2[1][1]-th2[1][2];
            th2[2][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th2[2][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th2[2][1]+th2[2][2]; th2[2][1]/=sum; th2[2][2]/=sum; th2[2][0]=1-th2[2][1]-th2[2][2];
            th3[0][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[0][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[0][3]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th3[0][1]+th3[0][2]+th3[0][3]; th3[0][1]/=sum; th3[0][2]/=sum; th3[0][3]/=sum;
			th3[0][0]=1-th3[0][1]-th3[0][2]-th3[0][3];
            th3[1][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[1][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[1][3]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th3[1][1]+th3[1][2]+th3[1][3]; th3[1][1]/=sum; th3[1][2]/=sum; th3[1][3]/=sum;
			th3[1][0]=1-th3[1][1]-th3[1][2]-th3[1][3];
            th3[2][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[2][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[2][3]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
			sum=1+th3[2][1]+th3[2][2]+th3[2][3]; th3[2][1]/=sum; th3[2][2]/=sum; th3[2][3]/=sum;
			th3[2][0]=1-th3[2][1]-th3[2][2]-th3[2][3];
            th3[3][1]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[3][2]=getParam(k1,0,k1,ii,jj,Params,0,expLnk); k1+=T;
            th3[3][3]=getParam(k1,0,k1,ii,jj,Params,0,expLnk);
			sum=1+th3[3][1]+th3[3][2]+th3[3][3]; th3[3][1]/=sum; th3[3][2]/=sum; th3[3][3]/=sum;
			th3[3][0]=1-th3[3][1]-th3[3][2]-th3[3][3];
            p11=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
            p12=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
            p22=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
            p13=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
            p23=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
            p33=getParam(kp,1,jp,ii,jj,Params,0,SingSpec.LinkFn); kp+=T; jp+=T;
		
                //   matrix multply current state (like3) * transition matrix.
            for (i=0; i<4; i++) for (j=0; j<4; j++) v[i]+=like3[j]*th3[j][i];   
            for (i=0; i<4; i++) { like3[i]=v[i]; v[i]=0; }//  store result in likeAB
			if(ii==II)printf("j=%d like3=%15.9f %15.9f %15.9f %15.9f (%d)\n",
			  jj,like3[0],like3[1],like3[2],like3[3],SingSpec.Data[ii][jj]);

			//   matrix multply current state (like2) * transition matrix.
            for (i=0; i<3; i++) for (j=0; j<3; j++) v[i]+=like2[j]*th2[j][i];
            for (i=0; i<3; i++) { like2[i]=v[i]; v[i]=0; }//  store result in like2
			
			//   matrix multply current state (like1) * transition matrix.
            for (i=0; i<2; i++) for (j=0; j<2; j++) v[i]+=like1[j]*th1[j][i];
            for (i=0; i<2; i++) { like1[i]=v[i]; v[i]=0; } //  store result in like1
			
            if (SingSpec.Data[ii][jj]==3) {
				like3[3]*=p33;             // site must be in global and local state 3, obs in state 3
				like3[2]=like3[1]=like3[0]=0;   // so states 0,1,2 impossible
				like2[2]=like2[1]=like2[0]=like1[1]=like1[0]=0; seen=3; 
			}
            if (SingSpec.Data[ii][jj]==2) {      // obs in state 2
                like3[3]*=p23;    //  so, could be in global state 3, local state 3, but obs in 2
				like3[2]*=p22;    //  or, could be in global state 3, local state 2, obs in 2
				like3[1]=like3[0]=0;   //  can't be in local states 0 or 1
				like2[2]*=p22;    //  or, could be in global state 2, local state 2, obs in 2
				like2[1]=like2[0]=like1[1]=like1[0]=0;  if (seen<2) seen=2;
            } 
            if (SingSpec.Data[ii][jj]==1) {     // obs in state 1
                like3[3]*=p13;   //  so, could be in global state 3, local state 3, obs in 1
				like3[2]*=p12;   //  or, could be in global state 3, local state 2, obs in 1
				like3[1]*=p11;   //  or, could be in global state 3, local state 2, obs in 1 
				like3[0]=0;
				like2[2]*=p12;   //  or, could be in global state 2, local state 2, obs in 1 
				like2[1]*=p11;   //  or, could be in global state 2, local state 1, obs in 1 
				like2[0]=0;
				like1[1]*=p11;    //  or, could be in global state 1, local state 1, obs in 1
				like1[0]=0; if (seen<1) seen=1;
            } 
            if (SingSpec.Data[ii][jj]==0) {  // obs in state 0 = unobserved
                like3[3]*= (1.0-p33-p23-p13);    //  so, could be in global state 3, local state 3, unobs
				like3[2]*=(1-p22-p12);          //  or, could be in global state 3, local state 2, unobs
				like3[1]*=(1-p11);             //  or, could be in global state 3, local state 1, unobs 
                like2[2]*=(1-p22-p12);        //  or, could be in global state 2, local state 2, unobs
				like2[1]*=(1-p11);            //  or, could be in global state 2, local state 1, unobs
                like1[1]*=(1-p11);            //  or, could be in global state 1, local state 1, unobs
            }       //  end if Data==0
			if(ii==II)printf("j:%d like3=%15.9f %15.9f %15.9f %15.9f\n",jj,like3[0],like3[1],like3[2],like3[3]);
        }         //  end for jj=0 to T-1
        like=like3[0]+like3[1]+like3[2]+like3[3]+like2[0]+like2[1]+like2[2]+like1[0]+like1[1];
        if (seen==0) like += (1.0-psi1-psi2-psi3); 
		if (ii==II) printf("like=%15.9f  %15.9f  %15.9f\n",like, 20000*like,SingSpec.det_hist_frq[ii]);
        if (like<SingSpec.nearzero) like=SingSpec.nearzero;
        sumll -= SingSpec.det_hist_frq[ii]*log(like); SingSpec.expval[ii]=like; 
		sumfrq += SingSpec.det_hist_frq[ii];
    }    //  end for ii=0 to N-1
	for (ii=0; ii<N; ii++)
		sumchi += (SingSpec.det_hist_frq[ii]-SingSpec.expval[ii]*sumfrq)*
		          (SingSpec.det_hist_frq[ii]-SingSpec.expval[ii]*sumfrq)/(SingSpec.expval[ii]*sumfrq);
				  
    printf("ifn=%d like=%f chi=%f  \n",SingSpec.ifn,sumll,sumchi); return(sumll);
}
