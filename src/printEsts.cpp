#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
// #include <afx.h>
#include "SingAnal.h"

void print_individual_estimates2(FILE *g, int realParmNum, int pgrp, int istrt, int nrealparms, double **covar, 
                                 double *Params, int iLinkFn, int firstSite, int nsites) {
    extern TSingSpec SingSpec; 
    double invlogit(double x);
    int i,site,srvy,current_real_param,kk,ll,kkk,lll,j,ipar,k1=0,k2=0,N,pmiss=0,sitelmt=1,nbetas,prntflag,ndec=4;
    double xpsi,xpsise,sum,sum2,beta,betaSE,dm,dm2,ci1,ci2,xpar,avgpsi=0,avgse=0,nout=0,rtmp,xn=(double)SingSpec.N;
    N=SingSpec.N; char lbl[80],fmt1[80],fmt2[80],fmt3[80];
    if (nrealparms<1) nrealparms=SingSpec.NrowsDM[pgrp];
    if (nrealparms<1) return;
    if (nsites==0) nsites=N; xn=1.;
    i=sprintf(lbl,"%%%d.%df",ndec+4,ndec);
    i=sprintf(fmt1,"%%-12s all-sites             :%s %s   %s -%s fixed\n",lbl,lbl,lbl,lbl);
    i=sprintf(fmt2,"%%-12s %%4d %%-16s:%s %s   %s -%s\n",lbl,lbl,lbl,lbl);
    i=sprintf(fmt3,"\n   avg. %%s (n=%%1.0f)   :%s %s   %s -%s\n",lbl,lbl,lbl,lbl);
    if(SingSpec.NrowsDM[1]+SingSpec.NrowsDM[2]+SingSpec.NrowsDM[3]+SingSpec.NrowsDM[4]==0 || pgrp==0) {
	  istrt=realParmNum; pgrp=0;
	}
    
    for (i=0; i<pgrp; i++) { // k1= # beta parms in list before current parm group
        for (j=0; j<SingSpec.NParKK[i]; j++) if (SingSpec.BetaFixed[i][j]>1.0e44) k1++;                                          
        k2+=SingSpec.NrowsDM[i];   // k2=# real parms in list before current parm group
    }
    strcpy(lbl,SingSpec.Realname[pgrp][istrt]);
    fprintf(g,"\n   Individual Site estimates of <%s>",lbl);
    if (firstSite>0) { fprintf(g," (non-surveyed sites)"); N=nsites-1; }
    fprintf(g,"\n                Site               estimate  Std.err   95%% conf. interval\n");
    current_real_param=istrt;
    for (ipar=0; ipar<nrealparms; ipar++) { prntflag=0; srvy=ipar;
        xpsi=SingSpec.fixed[k2+current_real_param]; xpsise=ci1=ci2=nout=avgpsi=avgse=0;
        strcpy(lbl,SingSpec.Realname[pgrp][istrt+ipar]);
        if (SingSpec.fixed[k2+current_real_param]>=-998) {     //  fixed parameter
            fprintf(g,fmt1,lbl,xpsi,xpsise,ci1,ci2); 
            nout=1; avgpsi=xpsi;
        }
        else {
            sitelmt=1;              //   if covariates in design matrix, set lmt to num of sites,  else set to 1.
            nbetas=SingSpec.NParKK[pgrp];
            if (SingSpec.lmt!=1)
                for (kk=0; kk<SingSpec.NParKK[pgrp]; kk++) {
                    if (SingSpec.DMat_ptr[pgrp][current_real_param][kk]!=0) sitelmt=SingSpec.N+SingSpec.nphantom;
                }
            if (SingSpec.LnkFn[pgrp][current_real_param]>0) iLinkFn=SingSpec.LnkFn[pgrp][current_real_param];
			if (sitelmt>SingSpec.lmt) sitelmt=SingSpec.lmt;
            for (site=0; site<sitelmt; site++) {                            //    For each site...
                for (kk=kkk=pmiss=0, sum=sum2=0.; kk<nbetas; kk++) {  //       For each beta column...
                    dm=SingSpec.DMat[pgrp][current_real_param][kk]; 
					j=SingSpec.DMat_ptr[pgrp][current_real_param][kk];
                    if(j>3000) dm=SingSpec.psibar[ipar][site][j-3001]/xn;
					else {
                        if(j<0) dm=SingSpec.SampCov[-j-1][site+firstSite][srvy];
                        if(j>0) dm=SingSpec.SiteCov[site+firstSite][j-1];
					}
                    xpar=SingSpec.BetaFixed[pgrp][kk];
                    if (xpar>1.0e44) xpar=Params[k1+kkk];
                    if (dm<-9998) pmiss=1; else sum += dm*xpar;
                    if (SingSpec.BetaFixed[pgrp][kk]>1.0e44) {
                        for (ll=lll=0; ll<nbetas; ll++) {
                            if (SingSpec.BetaFixed[pgrp][ll]>1.0e44) { 
                                dm2=SingSpec.DMat[pgrp][current_real_param][ll]; 
								j=SingSpec.DMat_ptr[pgrp][current_real_param][ll];
                                if(j>3000) dm2=SingSpec.psibar[ipar][site][j-3001]/xn;
								else {
                                    if(j<0) dm2=SingSpec.SampCov[-j-1][site+firstSite][srvy];
                                    if(j>0) dm2=SingSpec.SiteCov[site+firstSite][j-1];
								}
                                if (dm2>-9999) sum2 += dm*dm2*covar[k1+kkk][k1+lll];
                                lll++;
                            }
                        }
                        kkk++;
                    }
                }
                beta=sum; betaSE=sqrt(sum2);
                if(pmiss==0) {          //  not missing data
                    if (iLinkFn==sinLnk) {      //    sin link
                        xpsi=(sin(beta)+1.)/2.; xpsise=fabs(cos(beta))*betaSE/2.;
                        beta-=1.96*betaSE; ci1=(sin(beta)+1.)/2.; beta+=3.92*betaSE; ci2=(sin(beta)+1)/2.; 
                    }
                    if (iLinkFn==sinLnk2) {      //    sin link / 2
                        xpsi=(sin(beta)+1.)/4.; xpsise=fabs(cos(beta))*betaSE/4.;
                        beta-=1.96*betaSE; ci1=(sin(beta)+1.)/4.; beta+=3.92*betaSE; ci2=(sin(beta)+1)/4.; 
                    }
                    if(beta < -52.0) beta =  -52.0; if(beta > 52.0) beta =  52.0;
                    if (iLinkFn<1 || iLinkFn>8 || iLinkFn==logitLnk) {    //   logit link = default
                        xpsi=invlogit(beta); xpsise=exp(-beta)*xpsi*xpsi*betaSE;
                        beta-=1.96*betaSE; ci1=1./(1.+exp(-beta)); beta+=3.92*betaSE; ci2=1./(1.+exp(-beta)); 
                    }
                    if (iLinkFn==logitLnk2) {    //   logit link2 = logit link / 2
                        xpsi=invlogit(beta)/2.; xpsise=exp(-beta)*xpsi*xpsi*2*betaSE;
                        beta-=1.96*betaSE; ci1=1./(1.+exp(-beta))/2.; beta+=3.92*betaSE; ci2=1./(1.+exp(-beta))/2.; 
                    }
                    if (iLinkFn==loglogLnk) {     //   cloglog link
                        xpsi=1.-exp(-exp(beta)); xpsise=exp(beta)*(1.-xpsi)*betaSE;
                        beta-=1.96*betaSE; ci1=1.-exp(-beta); beta+=3.92*betaSE; ci2=1.-exp(-beta); 
                    }
                    if (iLinkFn==expLnk) {      //    exp link
                        xpsi=exp(beta); xpsise=xpsi*betaSE;
                        beta-=1.96*betaSE; ci1=exp(beta); beta+=3.92*betaSE; ci2=exp(beta); 
                    }
                }
                else {    // missing data 
                    xpsi=-999;
                }
                if (SingSpec.NParKK[pgrp]>0) {
                    if (pmiss==0) {
                        fprintf(g,fmt2,lbl,site+firstSite+1,SingSpec.sitename[site+firstSite],
                                xpsi,xpsise,ci1,ci2);
                        rtmp=1; if (firstSite==0) rtmp=SingSpec.det_hist_frq[site];
                        nout+=rtmp; avgpsi+=xpsi*rtmp; avgse+=xpsise*rtmp; prntflag=1;
                    }
                    else {
                        fprintf(g,"%-12s %4d %-16s:   .          .            .         -   .\n",
                                lbl,site+1,SingSpec.sitename[site]);
                    }
                }
                else fprintf(g,"\n");
            }
        }
        if (SingSpec.lmt==12345 && nout>1 && prntflag>0) {
            avgse/=nout; avgpsi/=nout; ci1=avgpsi-1.96*avgse; ci2=avgpsi+1.96*avgse;
            fprintf(g,fmt3,lbl,nout,avgpsi,avgse,ci1,ci2);
        }
        current_real_param++; if (strstr(lbl,"Cond")!=NULL) current_real_param+=2;
    }
}
