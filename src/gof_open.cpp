#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

double gof_open(double *Params, int prnt) {
    
    //   gof test for multi-season (open) model 
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    int i,ii,j,jj,k,kk,j2,N=SingSpec.N,T=SingSpec.T,debug=0; 
    double expect=0,*cohort,temp,TestStat=0,Nx=0,SeenSum,sumexp=0,df,lastexp=0,lasttemp=0;
    double lastsumexp=0, lastobs=0,sumobs,lastsumobs=0,ts1,ts2,eps=0,gam,cp; char fmt[120];
    
    int site,srvy, isrvy, seasn, PrmyPeriods=SingSpec.PrmyPeriods;
    double psi=0, prd1=1, prd0, lprd0, lprd1;
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    
    int *uniqdata; uniqdata=new int[N];
    int *uniqloc;  uniqloc=new int[N];
    int *uniqcht; uniqcht=new int[N];
    int *chtptr; chtptr=new int[N];
    int nuniqdata=0,nuniqcht=0,misseq,chtfound,found,equal; double *obs; char **hist;
    obs = new double[N]; cohort = new double[N];  hist = new char*[N];               prnt=2; 
    for (i=0; i<N; i++) { 
        obs[i]=cohort[i]=0; hist[i]=new char[T+1]; hist[i][T]='\0';
        for (j=0; j<T; j++) { 
            if (SingSpec.Data[i][j]>=0) hist[i][j]='0'+SingSpec.Data[i][j]; 
            else hist[i][j]='.';
        }
    }
    for (i=0,TestStat=0.;i<N; i++) {    // for each history...
        Nx+=SingSpec.det_hist_frq[i];
        for (found=-1,j=0; j<nuniqdata; j++) { // search table of unique histories
            jj=uniqdata[j];
            for (k=0,equal=1; k<T; k++) {
                if (SingSpec.Data[i][k]!=SingSpec.Data[jj][k]) equal=0;
            }
            if (equal==1 && found<0) found=j;
        }                
        for (chtfound=-1,j=0; j<nuniqcht; j++) { // search table of unique histories
            j2=uniqcht[j];
            for (k=0,misseq=1; k<T; k++) {
                if ((SingSpec.Data[i][k]<0 && SingSpec.Data[j2][k]>=0)||
                    (SingSpec.Data[i][k]>=0 && SingSpec.Data[j2][k]<0))misseq=0;
            }
            if (misseq==1 && chtfound<0) chtfound=j;
        }                
        if (found<0) { found=nuniqdata++; uniqdata[found]=i;}
        if (chtfound<0) {chtfound=nuniqcht++; uniqcht[chtfound]=i;}
        obs[found]+=SingSpec.det_hist_frq[found]; 
        cohort[chtfound]+=SingSpec.det_hist_frq[found]; 
        chtptr[found]=chtfound; uniqloc[i]=found;
    }
    if (debug) {
        fprintf(g,"\n i   uniqdata(i) chtptr(i) hist(uniqdata(i))\n");
        for (j=0; j<nuniqdata; j++) { 
            fprintf(g,"%d  %d %d : %s\n",j,uniqdata[j],chtptr[j],hist[uniqdata[j]]);
        }
        fprintf(g,"\n i   uniqcht(i)  hist(uniqcht(i))\n");
        for (j=0; j<nuniqcht; j++) { 
            fprintf(g,"%d  %d : %s\n",j,uniqcht[j],hist[uniqcht[j]]);
        }
        fprintf(g,"\n i uniqloc(i)\n");
        for (i=0; i<N; i++) {
            fprintf(g,"%d  %d : %s\n",i,uniqloc[i],hist[i]);
        }
    }
    
    // calculate test stat
    if(prnt>1) { 
        fprintf(g,"\nHistory(cohort)");
        strcpy(fmt,"  ** %s **        ");
        for (kk=15; kk<T; kk++) { fprintf(g," "); strcat(fmt," ");}
        fprintf(g,"		Observed    Expected    Chi-square\n");
        strcat(fmt,"%12.4f %16.4f %12.2f\n");
    }
    sumexp=sumobs=TestStat=SeenSum=df=0; int pstart=(SingSpec.PrmyPeriods-1)*2+1;
    for(ii=0; ii<nuniqdata; ii++) {    
        j=uniqdata[ii]; k=chtptr[ii];      //   sum prob of ii'th history over all histories in cohort
        for (site=0,expect=0; site<N; site++) {
			psi = getParam(0,0,0,site,0,Params,0,SingSpec.LinkFn); 
			prd0=1-psi; prd1=psi; if (SingSpec.OpenData[j][0]>0) prd0=0;
			for (srvy=isrvy=seasn=0; srvy<T; srvy++) {
				cp = getParam(pstart+srvy,3,srvy,site,srvy,Params,0,SingSpec.LinkFn);
                if (SingSpec.Data[j][srvy]!=-1 && cp>=0) {  //  if not missing data in this survey...
                    if (SingSpec.Data[j][srvy]>0) prd1*=cp;
                    else prd1*=(1-cp);
                }
				if (++isrvy==SingSpec.SecPeriods[seasn]) {
					if (srvy<(T-1)) {
						gam = getParam(seasn+1,1,seasn,site,seasn,Params,0,1);	eps=1-gam;
						if (SingSpec.Model<2) {
							eps = getParam(PrmyPeriods+seasn,2,seasn,site,seasn,Params,0,1);
							if (SingSpec.Alt_Param_Checked) eps=1-eps;
                        }
						lprd0=prd0; lprd1=prd1; prd0=lprd0*(1-gam);
						///prd0=lprd0*(1-gam)+lprd1*(SingSpec.OpenData[j][seasn+1]>0 ? 0 : eps);
						if (SingSpec.OpenData[j][seasn+1]==0) prd0+=lprd1*eps;
						prd1=lprd0*gam+lprd1*(1-eps);
					}
					seasn++; isrvy=0;
				}
            }
            temp=(prd0+prd1); expect+=temp*SingSpec.det_hist_frq[site];
        } // end site loop
        
        temp = (obs[ii]-expect)*(obs[ii]-expect)/(expect>0 ? expect : eps);
        TestStat += temp; SeenSum+=expect;
        if(prnt>1)fprintf(g,"%s(%2d %d) %15.4f %16.9f %12.2f\n",hist[j],k,j,obs[ii],expect,temp);
        sumexp=0; sumobs=0; df++; lastexp=expect; lastobs=obs[ii]; lasttemp=temp;
    }

    if (sumobs>0) {      //  if small exp/obs values remaining in pool at end
        if (lastexp==0) { lastobs=lastsumobs; lastexp=lastsumexp;}
        sumobs+=lastobs; sumexp+=lastexp;     // of histories, pool with last 'good' history
        temp = (sumobs-sumexp)*(sumobs-sumexp)/(sumexp>0 ? sumexp : eps); ++df;   //   and subtract obs/exp/chi 
        TestStat += temp; SeenSum+=expect;                     //   from TestStat
        if(prnt>1)fprintf(g,fmt,"pooled ",sumobs,sumexp,temp);
        TestStat-=lasttemp; --df; SeenSum-=lastexp;
        if(prnt>1)fprintf(g,fmt,"deleted",lastobs,lastexp,lasttemp);
    }
    
    ts1=TestStat; TestStat += Nx - SeenSum;  ts2=TestStat;
    if (prnt>0) fprintf(g,">sum(chisq)=%15.4f  sum+N-seensum=%15.4f,  Test Statistic = %11.4f\n",ts1,ts2,TestStat);
    for (i=0; i<N; i++) delete [] hist[i]; delete [] hist;
    delete [] uniqdata; delete [] uniqloc; delete [] uniqcht; delete [] chtptr; 
    delete [] obs; delete [] cohort;
    return(TestStat);
}
