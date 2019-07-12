#include "SingAnal.h"
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

int cmphist(int *h1, int *h2) {
	extern TSingSpec SingSpec;  int i,rv=0;
	for (i=0; (i<SingSpec.T && h1[i]==h2[i]); i++);
	if (i<SingSpec.T) rv=h1[i]-h2[i];
	return(rv);
}
int minhist(int *used) {
	extern TSingSpec SingSpec;   int i,rv=-1;
	for (i=0; i<SingSpec.N && rv<0; i++) if (used[i]==0) rv=i;  //  find 1st unused history...
	if (rv>=0) {
		for (i=0; i<SingSpec.N; i++) if (cmphist(SingSpec.Data[i],SingSpec.Data[rv])<0 && used[i]==0) rv=i;
		for (i=0; i<SingSpec.N; i++) if (cmphist(SingSpec.Data[i],SingSpec.Data[rv])==0) used[i]=1;
	}
	return(rv);
}
double sumfrq(int new1) {
	extern TSingSpec SingSpec;   int i; double sum=0;
	for (i=0; i<SingSpec.N; i++) 
		if (cmphist(SingSpec.Data[i],SingSpec.Data[new1])==0) 	sum+=SingSpec.det_hist_frq[i];
	return(sum);
}
void chisq_test(void) {
    extern TSingSpec SingSpec;     //  compute expected vals for GOF test for Single-season model
	double obs,exp,chi,sum=0,sumchi=0; 
	int site,k,N=SingSpec.N,T=SingSpec.T,new1,nout=0,*used; 
	FILE *g=SingSpec.g; time_t t1, t2; float t3=0; time(&t1);
	
	if (SingSpec.do_chisq>0) {
		used=new int[N]; for (site=0; site<N; site++) used[site]=0;
		printf("\nChi-square test...(NO POOLING of small expected values)\n           "); 
		fprintf(g,"\nChi-square test...(NO POOLING of small expected values)\n           "); 
		for (k=0; k<(T-4); k++) fprintf(g," ");
		fprintf(g,"hist         obs         exp         chi-sq\n"); 
		for (site=0; site<N; site++) sum+=SingSpec.det_hist_frq[site];
		if (sum==0) sum=N; 
		for (nout=0,new1=minhist(used); new1>=0; new1=minhist(used)) {
			fprintf(g,"%5d %5d ",++nout,new1+1);
			for (k=0; k<T; k++) 
				if (SingSpec.Data[new1][k]>=0) fprintf(g,"%d",SingSpec.Data[new1][k]);
				else fprintf(g,".");
			obs=sumfrq(new1);
			exp=SingSpec.expval[new1]*sum;
			chi=(obs-exp)*(obs-exp)/(exp+1e-15); sumchi+=chi; 
			fprintf(g," %12.5f %12.5f %12.3f\n",obs,exp,chi); 
		}
		fprintf(g,"\nTotal chi-square=%f\n",sumchi);
		time(&t2); t3=t2-t1; printf("chisq test CPU time: %f min.\n",t3/60.);
		fprintf(g,"chisq test CPU time: %f min.\n",t3/60.);
	}
}

/*double gofnew(int prnt, double poolcutoff) {    
    //   gof test for custom model 
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    int site,i,ii,j,jj,k,kk,j2,N=SingSpec.N,T=SingSpec.T,debug=1; char fmt[120];
    double expect=0,*cohort,temp,TestStat=0,Nx=0,SeenSum,df;
    double eps=1.0e-10,minexp=9.99e99;
    
    int *uniqdata; uniqdata=new int[N];    int *uniqloc;  uniqloc=new int[N];
    int *uniqcht; uniqcht=new int[N];    int *chtptr; chtptr=new int[N];
    int nuniqdata=0,nuniqcht=0,misseq,chtfound,found,equal; double *obs; 
    obs = new double[N]; cohort = new double[N];
    for (i=0; i<nboot; i++) {
		for (site=0; site<N; site++) {
			for (j=0,x=0; j<SingSpec.det_hist_frq[site]; j++) x+=(rand()<SingSpec.expval[site])
		    obs[site]=x;
     // calculate test stat
    if(prnt>1) { 
        chsq[i]=chisq_test(obs);
        strcpy(fmt,"  ** %s **        ");
        for (kk=5; kk<T; kk++) { fprintf(g," "); strcat(fmt," ");}
        fprintf(g,"		Observed        Expected        Chi-square\n");
        strcat(fmt,"%12.4f %16.4f %12.2f\n");
    }
    TestStat=SeenSum=df=0; 
    for(ii=0; ii<nuniqdata; ii++) {
        j=uniqdata[ii]; k=chtptr[ii];      //   sum prob of ii'th history over all histories in cohort
        for (site=0,expect=0; site<N; site++) {
            if (chtptr[uniqloc[site]]==k) {       //   if current hist is in same cohort as hist ii
                //  sum cell probs for each site in current cohort
                //  using psi,p,theta for site 'site', and data from site 'j'
                if (SingSpec.NMethods>1) expect+=comp_expval_mm(site,j,psi_hat,p_hat,th_hat);
                else  
                    if (SingSpec.Model==7) expect+=comp_expval_sd(site,j,psi_hat,p_hat,th_hat,pi_hat);
                    else 
                        if (SingSpec.NrowsDM[1]<1) expect+=comp_expval_het(site,j,psi_hat,th_hat,pi_hat);
                        else expect+=comp_expval(site,j,psi_hat,p_hat);
            }
        }   //  end of site loop
        if (prnt>1) 
            for(kk=0; kk<T; kk++) fprintf(g,"%c",(SingSpec.Data[j][kk]<0?'-':SingSpec.Data[j][kk]+'0'));				
        temp = (obs[ii]-expect)*(obs[ii]-expect)/(expect<eps?eps:expect);
        TestStat += temp; SeenSum+=expect;
        if(prnt>1)fprintf(g,"(%2d %d) %15.4f %16.9f %12.2f\n",k,j,obs[ii],expect,temp);
        df++;  
        if (expect<minexp) minexp=expect;
    }
    TestStat += Nx - SeenSum;  if (poolcutoff>.5) TestStat/=df;
    if (prnt>1) fprintf(g,"Test Statistic = %11.4f min(expect)=%e\n",TestStat,minexp);
    delete [] uniqdata; delete [] uniqloc; delete [] uniqcht; delete [] chtptr; 
    delete [] obs; delete [] cohort;
    return(TestStat);
}
*/
double comp_expval(int site, int j, double *psi_hat, double **p_hat) {
    extern TSingSpec SingSpec;     //  compute expected vals for GOF test for Single-season model
    double prd,expect=0; int t; bool seen; 
    for (t=0,prd=psi_hat[site],seen=false; t<SingSpec.T; t++) {  // compute exp val of unique hist ii
        if(SingSpec.Data[j][t]>=0) {
            if(SingSpec.Data[j][t]>0) {seen=true; prd*=p_hat[site][t];}
            else prd*=(1.0-p_hat[site][t]);
        }
    }
    if (!seen) prd+=(1.0-psi_hat[site]); 
    expect+=prd*SingSpec.det_hist_frq[site];  //  so, expval(hist ii) = num in cohort * prob of hist ii
    return(expect);
}
double comp_expval_fp(int site, int j, double *psi_hat, double **p_hat, double **p10_hat, double **b_hat) {
    extern TSingSpec SingSpec;     //  compute expected vals for GOF test for Single-season model
    double prd=psi_hat[site],prd2=(1-psi_hat[site]),expect=0; int t; bool seen=false;
    for (t=0,seen=false; t<SingSpec.T; t++) {  // compute exp val of unique hist ii
        if(SingSpec.Data[j][t]>=0) {
            if(SingSpec.Data[j][t]>0) {
				seen=true; 
			if (SingSpec.Data[j][t]>1) { prd*=p_hat[site][t]*b_hat[site][t]; prd2=0; }
				else { prd*=p_hat[site][t]*(1-b_hat[site][t]); prd2*=(1-p10_hat[site][t]); }
			}
            else prd*=(1-p_hat[site][t]);
        }
    }
    if (seen) prd2=0;
    expect=(prd+prd2)*SingSpec.det_hist_frq[site];  //  so, expval(hist ii) = num in cohort * prob of hist ii
    return(expect);
}
double comp_expval_het(int site, int j, double *psi_hat, double ***p_hat, double *theta) {
    extern TSingSpec SingSpec;     //  compute expected vals for GOF test for Single-season predefined model
    double prd1,prd2,expect=0; int t; bool seen;

    for (t=0,prd1=prd2=1,seen=false; t<SingSpec.T; t++) {  // compute exp val of unique hist ii
        if(SingSpec.Data[j][t]>=0) {
            if(SingSpec.Data[j][t]>0) {seen=true; prd1*=p_hat[0][t][0]; prd2*=p_hat[0][t][1];}
            else { prd1*=(1.0-p_hat[0][t][0]); prd2*=(1-p_hat[0][t][1]);}
        }
    }
	prd1=psi_hat[0]*(theta[0]*prd1+(1-theta[0])*prd2);
    if (!seen) prd1+=(1.0-psi_hat[site]); 
    expect+=prd1*SingSpec.det_hist_frq[site];  //  so, expval(hist ii) = num in cohort * prob of hist ii
    return(expect);
}
double comp_expval_mm(int site, int j, double *psi_hat, double **p_hat, double ***th_hat) {
    extern TSingSpec SingSpec; double prb,prd1=1,prd0=0,lprd0,lprd1,expect=0; 
    int t,srvy,seen_in_srvy=0,meth,seen=0;    //  compute expected vals for GOF test for Multi-method model
    for (t=srvy=meth=0; t<SingSpec.T; t++) {  // compute exp val of unique hist ii
        if (t % SingSpec.NMethods == 0) {
            lprd0=prd0; lprd1=prd1; seen_in_srvy=(SingSpec.Data[j][t]>0);
            prd0=(lprd0+lprd1)*(1-th_hat[site][srvy][0]);
            prd1=(lprd0+lprd1)*th_hat[site][srvy][0];
            srvy++; meth=0;
        }
        if(SingSpec.Data[j][t]>=0) {
            if(SingSpec.Data[j][t]>0) { seen=seen_in_srvy=1; prd1*=p_hat[site][t];}
            else prd1*=(1.0-p_hat[site][t]);
        }
		if (++meth==SingSpec.NMethods && seen_in_srvy>0) prd0=0;
    }
    prb=psi_hat[site]*(prd0+prd1);
    if (seen==0) prb+=(1.0-psi_hat[site]); 
    expect+=prb*SingSpec.det_hist_frq[site];  //  so, expval(hist ii) = num in cohort * prob of hist ii
    return(expect);
}

double comp_expval_sd(int site, int j, double *psi_hat, double **p_hat, double ***th_hat, double *pi_hat) {
    extern TSingSpec SingSpec;    //  compute expected vals for GOF test for Spatial Dependence model
    double prd1=0,prd0=1,lprd0,lprd1,expect=0,pstar=1,theta,theta_prime,pi;
    int srvy; 
    bool inarea=false,seen=false,no_det_yet=true;  //  theta is now 'season' specific.  So, .  
    //  order of parameters is: theta(1),theta(2),..theta(T),theta'(1),theta'(2),...theta'(T)
    for (srvy=0; srvy<SingSpec.T; srvy++) {        
        theta      =th_hat[site][srvy][0];   //  theta(not prev occupied)
        theta_prime=th_hat[site][srvy][1];   //  theta(prev occupied)
        if (no_det_yet && th_hat[0][0][1]>=0) {
            if (SingSpec.Alt_Param_Checked==1) theta=theta/(theta+1-theta_prime);  //  theta starting in middle of trail
            else {
                pi=0;   //  pi = prob that site is occupied before 1st segment
                if (SingSpec.NrowsDM[2]>0) pi=pi_hat[site];
                theta=pi*theta_prime+(1-pi)*theta;  //  theta for starting in middle of trail
            }
        }
        if (SingSpec.Data[j][srvy]!=-1 && p_hat[site][srvy]>=0) {  //  if not missing data in this survey...
            pstar=(SingSpec.Data[j][srvy]>0?p_hat[site][srvy]:1-p_hat[site][srvy]);
        }
        lprd0=prd0; lprd1=prd1; 
        if (SingSpec.Data[j][srvy]>=0 && pstar>=0) {
            if (SingSpec.Data[j][srvy]>0) { //   if detected in this survey...
                prd0=0; prd1=(lprd0*theta+lprd1*theta_prime)*pstar; inarea=seen=true;
            } 
            else {              //   not detected in this primary period
                prd0=lprd0*(1-theta)+lprd1*(1-theta_prime);
                prd1=(lprd0*theta+lprd1*theta_prime)*pstar;
            }
            no_det_yet=false;
        }
    }
    expect=(prd0+prd1)*psi_hat[site] + (1-inarea)*(1-psi_hat[site]);
    return(expect);
}
double gof(double *psi_hat, double **p_hat, double ***th_hat, double *pi_hat, double **p10_hat, double **b_hat, int prnt, double poolcutoff, int mtype) {    
    //   gof test for custom model 
    extern TSingSpec SingSpec; FILE *g; g=SingSpec.g;
    int site,i,ii,j,jj,k,kk,j2,N=SingSpec.N,T=SingSpec.T; char fmt[120];
    double expect=0,*cohort,temp,TestStat=0,Nx=0,SeenSum,df;
    double eps=1.0e-10,minexp=9.99e99;
    
    int *uniqdata; uniqdata=new int[N];    int *uniqloc;  uniqloc=new int[N];
    int *uniqcht; uniqcht=new int[N];    int *chtptr; chtptr=new int[N];
    int nuniqdata=0,nuniqcht=0,misseq,chtfound,found,equal; double *obs; 
    obs = new double[N]; cohort = new double[N];
    for (i=0; i<N; i++) obs[i]=cohort[i]=0;
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
    // calculate test stat
    if(prnt>1) { 
        fprintf(g,"\nHistory(cohort)"); strcpy(fmt,"  ** %s **        ");
        for (kk=5; kk<T; kk++) { fprintf(g," "); strcat(fmt," ");}
        fprintf(g,"		Observed        Expected        Chi-square\n");
        strcat(fmt,"%12.4f %16.4f %12.2f\n");
    }
    TestStat=SeenSum=df=0; 
    for(ii=0; ii<nuniqdata; ii++) {
        j=uniqdata[ii]; k=chtptr[ii];      //   sum prob of ii'th history over all histories in cohort
        for (site=0,expect=0; site<N; site++) {
            if (chtptr[uniqloc[site]]==k) {       //   if current hist is in same cohort as hist ii
                //  sum cell probs for each site in current cohort
                //  using psi,p,theta for site 'site', and data from site 'j'
				switch(mtype) {
                    case 1: expect+=comp_expval_mm(site,j,psi_hat,p_hat,th_hat); break;
					case 2: expect+=comp_expval_sd(site,j,psi_hat,p_hat,th_hat,pi_hat); break;
                    case 3: expect+=comp_expval_het(site,j,psi_hat,th_hat,pi_hat); break;
					case 4: expect+=comp_expval_fp(site,j,psi_hat,p_hat,p10_hat,b_hat); break;
                    default: expect+=comp_expval(site,j,psi_hat,p_hat);
				}
            }
        }   //  end of site loop
        if (prnt>1) for(kk=0; kk<T; kk++) fprintf(g,"%c",(SingSpec.Data[j][kk]<0?'-':SingSpec.Data[j][kk]+'0'));				
        temp = (obs[ii]-expect)*(obs[ii]-expect)/(expect<eps?eps:expect);
        TestStat += temp; SeenSum+=expect;
        if(prnt>1)fprintf(g,"(%2d %d) %15.4f %16.9f %12.2f\n",k,j,obs[ii],expect,temp);
        df++;  
        if (expect<minexp) minexp=expect;
    }
    TestStat += Nx - SeenSum;  if (poolcutoff>.5) TestStat/=df;
    if (prnt>1) fprintf(g,"Test Statistic = %11.4f min(expect)=%10.4f\n",TestStat,minexp);
    delete [] uniqdata; delete [] uniqloc; delete [] uniqcht; delete [] chtptr; 
    delete [] obs; delete [] cohort;
    return(TestStat);
}
