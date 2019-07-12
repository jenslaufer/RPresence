#include <stdio.h>
#include <math.h>
#include "SingAnal.h"

double PDLike(double *Params, int NPar)
{
    int ii, jj, kk,k; double psi, *mix, **p, sum, twopi=6.283185307179586476925286766559;
    extern TSingSpec SingSpec;
    double like,temp1,xlmt=twopi+.1e-9;  bool seen;
    // get current parameter values
    //fprintf(g,"pdlike:%f %f %f: ",Params[0],Params[1],floor(Params[0]/twopi));
    //if (Params[0]>twopi) Params[0]-=(floor(Params[0]/twopi)*twopi);
    //if (Params[0]<-twopi) Params[0]-=(floor(Params[0]/twopi)*twopi);
	//printf("pdlike NPar=%d xlmt=%f params=%f %f\n",NPar,xlmt,Params[0],Params[1]);
    if (SingSpec.MissClass) {psi=1; jj=0;}
    else {
        jj=1; 
        for (k=0; k<NPar; k++) 
            if (k==0 || k>=SingSpec.Groups) {
                while (Params[k]>xlmt) Params[k]-=xlmt; 
                while (Params[k]<-xlmt) Params[k]+=xlmt; 
            }
        psi = (sin(Params[0])+1.0)/2.0;
    }
    mix = new double [SingSpec.Groups];	mix[0]=1.0;
    if (SingSpec.Groups==2) { mix[0]=(sin(Params[1])+1.)/2.; mix[1]=1-mix[0]; }
    if (SingSpec.Groups>2) {
        for (ii=0,sum=1.0; ii<SingSpec.Groups-1; ii++) sum += exp(Params[ii+jj]); 
        mix[0] = 1.0/sum; for (ii=1; ii<SingSpec.Groups; ii++) mix[ii] = exp(Params[ii+jj-1])/sum; 
    }
    p = new double*[SingSpec.Groups]; kk=SingSpec.Groups-1;
    for (ii=0; ii<SingSpec.Groups; ii++) {
        p[ii] = new double[SingSpec.T]; kk++;
        for (jj=0; jj<SingSpec.T; jj++) {
            if (SingSpec.TSpecificP && jj>0) kk++;
            p[ii][jj] = (sin(Params[kk])+1.0)/2.;
        }
    }
    // go through data and evaluate likelihood
    for (ii=0,sum=0.0; ii<SingSpec.N; ii++) {
        for (kk=0,like=0.0,seen=false; kk<SingSpec.Groups; kk++) {
            temp1 = mix[kk];
            for (jj=0; jj<SingSpec.T; jj++) {
                if (SingSpec.Data[ii][jj]!=-1) {
                    if (SingSpec.Data[ii][jj]>=1) { temp1 *= p[kk][jj]; seen = true; } 
                    else { temp1 *= (1.0-p[kk][jj]); }
                }
            }
            like += temp1;
        }
        like *= psi; if (!seen) like += (1.0-psi);  SingSpec.expval[ii]=like;
        sum -= SingSpec.det_hist_frq[ii]*log(like);
    }
    //printf("ll=%f parms:%f %f %f %f %f %f\n",2*sum,psi,mix[0],mix[1],p[0][1],p[1][1],like);
    delete[] mix;  // delete dynamic variables
    for (ii=0; ii<SingSpec.Groups; ii++) { delete[] p[ii]; } delete[] p;
    return(sum);
}
