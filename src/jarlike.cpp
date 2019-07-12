#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SingAnal.h"

double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);

double JARLike(double *Params, int NPar) {            //  Royle Point-count model w/ poisson distribution
    int i, k, x, t, lmt1, lmt2, site, errcnt=0, zip=0;
    double *p, lam, prd, *f, sumf=0, like, loglike=0, eps=1.e-18,pnlty=0, psi=1;
    extern TSingSpec SingSpec; static int nfcall=0;  int Abund_Index_Model=SingSpec.TSpecificP;

    f=new double[SingSpec.rlmt+10];  p=new double[SingSpec.T];
    //printf("npar=%d lmt=%d parms:%f %f\n",NPar,SingSpec.rlmt,Params[0],Params[1]);
    for (i=0; i<NPar; i++) { if(Params[i]<-33) Params[i]=-33; if(Params[i]>33) Params[i]=33;}
    lmt1=SingSpec.rlmt; if (SingSpec.NrowsDM[0]>1) zip=1;

    for (site=0,loglike=0.0; site<SingSpec.N; site++) {      // go through data and evaluate likelihood
		for (t=0; t<SingSpec.T; t++) p[t]=getParam(1+t+zip,1,t,site,t,Params,0,1);
		if (zip>0) psi=getParam(1,0,1,site,0,Params,0,1)*SingSpec.area[site]; //printf("lmt=%d site=%d psi=%f area=%f\n",lmt1,site,psi,SingSpec.area[site]);
		lam=getParam(0,0,0,site,0,Params,1,1)*SingSpec.area[site];	//printf("lmt=%d site=%d lam=%f\n",lmt1,site,lam);
        if (lam<.1e-8)lam=.1e-8;
        if(!Abund_Index_Model && lam>800) lam=800;  // avoid numerical errors if lam too large
        f[0]=exp(-lam); lmt2=lmt1;
        for (k=1; k<=lmt1; k++) {
			f[k]=f[k-1]*lam/(double)k;  f[k-1]*=psi; //f=dpois(j,lam)=lam^k*exp(-lam)/k!
			if (f[k]<eps && lmt2==lmt1) lmt2=k;
		}
		f[0]=f[0]+1-psi; f[lmt1]*=psi; for (k=0; k<=lmt1; k++) sumf+=f[k];
		//printf("f(%d)=%22.18e f(%d)=%22.18e maxx(site)=%d sumf=%f\n",lmt2,f[lmt2],lmt1,f[lmt1],SingSpec.maxx[site],sumf);
        for (k=0,like=0.0; k<=lmt2; k++) {
            if (k>=SingSpec.maxx[site]) {
                for (t=0, prd=1; t<SingSpec.T; t++) {
                    x=SingSpec.Data[site][t];		//printf("t=%d k=%d x=%f\n",t,k,x);
                    if(x>=0)prd*=SingSpec.choose[k][x]*pow(p[t],x)*pow((1.-p[t]),(k-x));//a=prd(dbinom(x(site,),k,p));
                }
                if (sumf!=0) like +=(prd*f[k]/sumf); 
                else {
                    pnlty=99999;
                    if (++errcnt<1) 
                        printf("k,maxx(site),prd,f(k),sumf,like:%d %d %f %f %f %f\n",k,SingSpec.maxx[site],prd,f[k],sumf,like);
                }
            }   //  end if k>=SingSpec.maxx[site]
        }    //  end for(k=0 to lmt1)
        if (like==0) like=SingSpec.nearzero;
        loglike -= SingSpec.det_hist_frq[site]*log(like); SingSpec.expval[site]=like;
    }    //  end for(site=0 to N)
    nfcall++; delete [] f; delete [] p;
    return(loglike+pnlty);
}

double JARLikeNB(double *Params, int NPar) {            //  Royle Point-count model w/ neg. binomial prior dist.
    int i, k, x, t, site, lmt1, lmt2, zip=0; 
	double *p, prd, *f, sumf=0, like, loglike, u=0, a=0, term2, term3, eps=1.e-18, psi=1; 
    extern TSingSpec SingSpec; static int nfcall=0;
    
    f=new double[SingSpec.rlmt+1];  p=new double[SingSpec.T]; 
	lmt1=SingSpec.rlmt; if (SingSpec.NrowsDM[0]>2) zip=1;
    for (i=0; i<NPar; i++) { if(Params[i]<-22) Params[i]=-22; if(Params[i]>22) Params[i]=22;}
	
    for (site=0,loglike=0.0; site<SingSpec.N; site++) {      // go through data and evaluate likelihood
 	    for (t=0; t<SingSpec.T; t++) p[t]=getParam(zip+2+t,1,t,site,t,Params,0,1);
        u=getParam(0,0,0,site,0,Params,1,1);  a=u*u/(getParam(1,0,1,site,0,Params,1,1));
        if (zip>0) psi=getParam(2,0,2,site,0,Params,0,1)*SingSpec.area[site];      
		//r=a/(u+a); f[k]=gamma(k+a)/gamma(a)/gamma(k-1)*pow(r,a)*pow(1-r,k);
		if (u>100.) u=100;
        if (u<SingSpec.nearzero)u=SingSpec.nearzero;
        term2=u/(a+u+eps); term3=pow(a/(a+u+eps),a); //printf("a=%f u=%f term3=%f\n",a,u,term3); // 1/(1+u/a) = a/(a+u)
        f[0]=term3*psi; prd=1; lmt2=SingSpec.rlmt;
        for (k=1; k<=lmt1; k++) {// recursive func for neg.bin.
            prd*=((a+k-1)/(double)(k)) * term2; f[k]=term3*prd*psi; 
			if (f[k]<eps && lmt2==SingSpec.rlmt) lmt2=k;
        }
		f[0]=f[0]+1-psi; for (k=0; k<=lmt1; k++) sumf+=f[k];
        for (k=0,like=0.0; k<=lmt2; k++) {
            if (k>=SingSpec.maxx[site]) {
                for (t=0, prd=1; t<SingSpec.T; t++) {
                    x=SingSpec.Data[site][t];
                    if(x>=0)prd*=SingSpec.choose[k][x]*pow(p[t],x)*pow((1.-p[t]),(k-x));
                }
                like +=prd*f[k]/(sumf+eps);
            }
        }  //  end for (k=0 to SingSpec.rlmt)
        loglike -= SingSpec.det_hist_frq[site]*log(like+eps); SingSpec.expval[site]=like;
    }  //  end for (site=0 to N)
    if(SingSpec.Verbose>1) printf("lik3([%f %f %f %f %f %f]=%f (%d)\n ",Params[0],Params[1],Params[2],p[0],u,a,2*loglike,SingSpec.rlmt); 
    nfcall++; delete f; delete p;
    return(loglike);
}

double NHetLike(double *Params, int NPar) {           //  Royle-Nichols heterogeneity model
    int site, k, t, lmt1, lmt2,zip=0;  double p, *r, lam, prd, *f, loglike, exlam, sm, tfac, eps=1.e-18, psi=1;
    extern TSingSpec SingSpec;    r=new double[SingSpec.T]; f=new double[SingSpec.rlmt+10];
    lmt1=SingSpec.rlmt; if (SingSpec.NrowsDM[0]>1) zip=1;
    for (site=0,loglike=0.; site<SingSpec.N; site++) {  // go through data and evaluate likelihood
        for (t=0; t<SingSpec.T; t++) r[t]=getParam(zip+1+t,1,t,site,t,Params,0,1);
        lam=getParam(0,0,0,site,0,Params,1,1); exlam=exp(-lam);
		if (zip>0) psi=getParam(1,0,1,site,0,Params,0,1)*SingSpec.area[site];   
		for (t=0,sm=0.,tfac=1; sm<1; tfac=tfac*t++) sm=pow(lam,t)*exlam/tfac; 
        exlam=f[0]=exp(-lam); lmt2=SingSpec.rlmt;
		for (k=1; k<=SingSpec.rlmt; k++) {
			f[k]=f[k-1]*lam/(double)k;  f[k-1]*=psi;//f=dpois(j,lam)=lam^k*exp(-lam)/k!
			if (f[k]<eps && lmt2==SingSpec.rlmt) lmt2=k;
		}
		f[0]=f[0]+1-psi; f[lmt1]*=psi;
        for (k=0,sm=0.; k<=lmt2; k++) {
            for (t=0,prd=1; t<SingSpec.T; t++) {
                if (SingSpec.Data[site][t]>=0) {
                    p=1-pow((1-r[t]),k); 
                    if (SingSpec.Data[site][t]==0) prd*=(1.-p);
                    else prd*=p;
                }
            }    // end for (t=0 to T)
            sm+=prd*f[k];
        }   //  end for (k=0 to SingSpec.lmt)
        if(sm<SingSpec.nearzero)sm=SingSpec.nearzero;
        loglike -=SingSpec.det_hist_frq[site]*log(sm); SingSpec.expval[site]=sm;
    }   //  end for (site=0 to N)
    delete[] r;  delete [] f;
	if(SingSpec.Verbose>1) printf("%f\n",loglike);
	return(loglike);
}

double NHetLikeNB(double *Params, int NPar) {           //  Royle-Nichols heterogeneity model w/ neg-binomial
    int site, k, t, lmt2, zip=0;  
	double p, u, a, *r, prd, *f, loglike, sm, term2, term3, eps=1.e-18, psi=1;
    extern TSingSpec SingSpec;    r=new double[SingSpec.T]; f=new double[SingSpec.rlmt+10];
    if (SingSpec.NrowsDM[0]>1) zip=1;
    for (site=0,loglike=0.; site<SingSpec.N; site++) {  // go through data and evaluate likelihood
        for (t=0; t<SingSpec.T; t++) r[t]=getParam(zip+2+t,1,t,site,t,Params,0,1);
        u=getParam(0,0,0,site,0,Params,1,1);  a=u*u/getParam(1,0,1,site,0,Params,1,1);
		if (zip>0) psi=getParam(2,0,2,site,0,Params,0,1)*SingSpec.area[site];   
        if (u<SingSpec.nearzero)u=SingSpec.nearzero; if (u>100.) u=100.;
        term2=u/(a+u); term3=pow(a/(a+u),a);   // 1/(1+u/a) = a/(a+u)
        f[0]=term3*psi; prd=1; lmt2=SingSpec.rlmt;
        for (k=1; k<=SingSpec.rlmt; k++) {// recursive func for neg.bin.
            prd*=((a+k-1)/(double)(k)) * term2; f[k]=term3*prd*psi; 
			if (f[k]<eps && lmt2==SingSpec.rlmt) lmt2=k;
        }
		f[0]=f[0]+1-psi; 
        for (k=0,sm=0.; k<=lmt2; k++) {
            for (t=0,prd=1; t<SingSpec.T; t++) {
                if (SingSpec.Data[site][t]>=0) {
                    p=1-pow((1-r[t]),k); 
                    if (SingSpec.Data[site][t]==0) prd*=(1.-p);
                    else prd*=p;
                }
            }    // end for (t=0 to T)
            sm+=prd*f[k];
        }   //  end for (k=0 to SingSpec.lmt)
        if(sm<SingSpec.nearzero)sm=SingSpec.nearzero;
        loglike -=SingSpec.det_hist_frq[site]*log(sm); SingSpec.expval[site]=sm;
    }   //  end for (site=0 to N)
    delete[] r;  delete f;
	if(SingSpec.Verbose>1) printf("%f\n",loglike);
	return(loglike);
}
