#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>

#include "SingAnal.h"
#include "rngs.h"
#define FALSE (1==0)

//TSingSpec SingSpec;

int rpresence(int *NN, int *kk, int *nseasons, int *nmeth,                           // nsites, nsurveys, nseasons, nmethods
         int *dd, double *frq, int *nunitcov, double *unitcov, char **unitcovnames, // det.data, frequencies, nunitcovs, unitcov matrix (as vector)
         int *nsurvcov, double *survcov, char **survcovnames,                       // nsurvey covs, survcov array (as vector), survcovnames
         int *nsurveyseason, char **title,                                         // vector of nsurveys per season, title
         char **unitnames, char **surveynames,                                    // vector of unitnames, surveynames
         int *dmdims,                                                             // vector of... dim(dm0), dim(dm1), dim(dm2),...
         char **dm0, char **rnames0,                                             // vectorized dm0,dm1,dm2,dm3,dm4,dm5,  vector of dm rownames
         double *iv, double *fv,                                                // vector of init values,  fixed values
		 int *neighbor,                                                        // neighbor matrix (NxN) as vector
		 double *neighbor_wgt,                                                // neighbor weights (Nx1) vector
         char **gname,                                                       //  output filename (should include random chars for parallel runs)
         char **argv,                                                       // char vector of other run options
		 double *rvals) {

	extern TSingSpec SingSpec;
    int i,ii,j=0,jj,dmnum,irows,icols,k,nbootstraps=0,nboot2=0,prnt=0,modtype=1,nsites,csum,nn,NPar;
    int nsitecov,nsampCov,CustomModel,notallmiss,seen,gof_err=0,maxstate=0,N,T,NrealPar=0;
    time_t t1; time_t t2; FILE *g; float t3;
    char s[65536], cwd[512]=".",tmps[1024];
    double sim_psi=0.7, sim_p=0.4, poolcutoff=0.0,n=0,tsites=0, *Params,totsites=0, eps;
    int sim_nsites=200, sim_nintsites=100, sim_nintvisits=9, sim_nother=5, sim_removal=0, sim_nsims=100, do_sims=0;
    char sim_fname[256]=""; const char *autocovs[4]={"psi1","psiA","upsi","psiB"};

    double atof_comma(char *s);
    void Simulation(double PSI,double P,int NTot,int NInt,int TInt,int TOther,int Removal, int NSims, char *fname);
    void PreDefMod(int NBoot, int LoFBoot);
    void CustomMod(int NBoot, int nboot2, int NPar, double *Params,double poolcutoff);
    void OpenPopn(int NBoot, char *s, int NPar, double *Params, int nboot2);
    void OpenPopnMM(int NBoot, char *s, int NPar, double *Params, int nboot2);
    void TwoSpec(int NBoot, char *s, int NPar, double *Params);
    void TwoSpecMultiSeason(int NBoot, char *s, int NPar, double *Params);
    void RoyleMod(int NBoot, int goftest, int NPar, double *Params);
    void MSMod(int NBoot, int LoFBoot, int NPar, double *Params);
    void MSOpenPop(int NBoot, char *s, int NPar, double *Params);
    void StagEntModl(int nbootstraps, char *s, int NPar, double *Params);  //  single-season staggered entry model
    void OpenPopSD(int NPar, double *Params, int nboot2);
	void TwoSpecMisId(int NPar, double *Params);
    int chk_real_labels(char *lbl, int *modtype, int i, int dmnum);
	int check_desmat_rowlabels(int *modtype);
    int get_modl_num(char *s);
	void get_opendata(int modtype);

    SingSpec.DMat = new double**[9]; SingSpec.Betaname=new char **[9];  SingSpec.DMat_ptr = new int**[9];
    SingSpec.Realname=new char **[9]; SingSpec.BetaFixed=new double*[9];
    SingSpec.LnkFn=new int *[9];
    char timstr[40]; time(&t1); strcpy(timstr,ctime(&t1)); timstr[strlen(timstr)-1]='\0';
    CustomModel=true;
    SingSpec.seed=(unsigned long) -time(NULL);  modtype=1;

	SingSpec.NSiteCov=0;        //  Number of site covariates
	SingSpec.NSampCov=0;        //  Number of sample(survey) covariates
	SingSpec.PrmyPeriods=1;     //  Number of Primary periods
	SingSpec.Nstates=1;                    //  Number of states(strata) in detection-history records
	SingSpec.NMethods=1;                   //  Number of methods per survey
	SingSpec.lmt=999999;           //  max number of sites to print for each real parameter
	SingSpec.rlmt=200;
	SingSpec.ifn=0;
	SingSpec.maxfn=50000;
	SingSpec.Verbose=1;
	SingSpec.novar=4;
    SingSpec.Groups=1;
    SingSpec.UseNeighbor=SingSpec.Alt_Param_Checked=0;
	SingSpec.uncondpsi=SingSpec.do_chisq=SingSpec.nrandIV=SingSpec.nphantom=SingSpec.NMiss=0;
	SingSpec.Model=SingSpec.LinkFn=1;
	SingSpec.LikeNRSig=6;
	SingSpec.BootNRSig=4;             // numerical covergence options
    SingSpec.TSpecificP=SingSpec.MissClass=SingSpec.FalsePos=SingSpec.optmiz_done=SingSpec.UseAmoeba=SingSpec.multiseason=false;
	SingSpec.prnt_derived=SingSpec.dups=true;
    SingSpec.nearzero=1e-200;
	SingSpec.VCeps=.01;
	SingSpec.nthreads=1;
	SingSpec.simulating=0;

    for (i=0; i<6; i++) SingSpec.NParKK[i]=SingSpec.NrowsDM[i]=0;

    if ((g=SingSpec.g=fopen(gname[0],"wt"))==NULL)
	    if (strstr(gname[0],"/dev/null")==NULL) { printf("\nError opening output file (%s)\n",gname[0]); return(1); }
    fprintf(g,"\n\nPRESENCE - Presence/Absence-Site Occupancy data analysis\n");
    fprintf(g,"%s,       Version %s.\n",timstr,VERSN);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");


    for (i=0; strlen(argv[i])>0; i++) {
		strcpy(s,argv[i]);
        //if (strstr(s,"quiet")  !=NULL || strstr(s,"QUIET")!=NULL) SingSpec.Verbose=0;
		if (strstr(s,"verbose")!=NULL) SingSpec.Verbose=2;
        if (strstr(s,"debug")  !=NULL) SingSpec.Verbose=3;
    }
    if (SingSpec.Verbose>0) printf("PRESENCE Version %s.\n",VERSN);

	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		if (SingSpec.Verbose>0) printf("\nCurrent working dir: %s\n", cwd);
	}
	else
		perror("getcwd() error");

    for (i=0; strlen(argv[i])>0; i++) {
      strcpy(s,argv[i]); j=(int)strlen(s); fprintf(g,"==>%s\n",s); printf(">>%s\n",s); if (SingSpec.Verbose) printf("==>%s\n",s);
      for (k=0; k<j && s[k]!='='; k++) s[k]=tolower(s[k]);
      if (strstr(s,"echo")    !=NULL) prnt=1;
      if (strstr(s,"name=")   !=NULL) strcpy(SingSpec.modname,s+5);
      if (strstr(s,"maxfn=")  !=NULL)  SingSpec.maxfn=atoi(s+6);
      if (strstr(s,"lmt=")    !=NULL && strstr(s,"rlmt")==NULL) SingSpec.lmt=atoi(s+4);
      if (strstr(s,"rlmt=")   !=NULL) SingSpec.rlmt=atoi(s+5);
      if (strstr(s,"link=c")  !=NULL) SingSpec.LinkFn=2;
      if (strstr(s,"link=s")  !=NULL) SingSpec.LinkFn=4;
      if (strstr(s,"nsig=")   !=NULL) SingSpec.LikeNRSig=atoi(s+5);
      if (strstr(s,"zero=")   !=NULL) SingSpec.nearzero=pow(10.,(double)(-atoi(s+5)));
      if (strstr(s,"nboot=")  !=NULL) nbootstraps=atoi(s+6);
      if (strstr(s,"boot2=")  !=NULL) nboot2=atoi(s+6);
      if (strstr(s,"pool=")   !=NULL) poolcutoff=atof(s+5);
      if (strstr(s,"nstates=")!=NULL) SingSpec.Nstates=atoi(s+8);
      if (strstr(s,"amoeba")  !=NULL) SingSpec.UseAmoeba=1;
      if (strstr(s,"seed")    !=NULL) SingSpec.seed = atol(s+5);
      if (strstr(s,"vceps")   !=NULL) SingSpec.VCeps = atof(s+6);
      if (strstr(s,"neighb")  !=NULL) SingSpec.UseNeighbor = 2;
      if (strstr(s,"nrand=")  !=NULL) SingSpec.nrandIV=atoi(s+6);
	  if (strstr(s,"nthreads=")  !=NULL) SingSpec.nthreads=atoi(s+9);

//        output options for SE's and VC matrix
//	  Option 1:  only compute likelihood and beta estimates, no beta var-cov matrix or beta SE's (useful for faster model selection of big models or simulations)
//    Option 2: compute but don't print beta var-cov matrix so we get beta SE's, don't compute real params
//    Option 3: compute and print beta beta var-cov matrix, don't compute real params (I think this is what you want for RPresence)
//    Option 4: (default) compute real params and real var-cov matrix, don't print either var-cov matrix.
//    Option 5: compute real params and real var-cov matrix, print only beta var-cov matrix.
//    Option 6: print both var-cov matrices.

      if (strstr(s,"nose")      !=NULL) SingSpec.novar = 1;
      if (strstr(s,"betavc")    !=NULL) SingSpec.novar = 2;
      if (strstr(s,"noreal")    !=NULL) SingSpec.novar = 3;
      if (strstr(s,"realbetavc")!=NULL) SingSpec.novar = 5;
      if (strstr(s,"bothvc")    !=NULL) SingSpec.novar = 6;
      if (strstr(s,"chsq")           !=NULL) SingSpec.do_chisq = 1;
      if (strstr(s,"noderived")      !=NULL) SingSpec.prnt_derived = false;
      if (strstr(s,"nodups")         !=NULL) SingSpec.dups = false;
      if (strstr(s,"model=")         !=NULL) { modtype=get_modl_num(s); if (SingSpec.Verbose) printf("modtype=%d\n",modtype);}
      if (strstr(s,"sim_psi=")       !=NULL) { sim_psi=atof(s+8); do_sims=1;}
      if (strstr(s,"sim_p=")         !=NULL) { sim_p=atof(s+6); do_sims=1;}
      if (strstr(s,"sim_nsites=")    !=NULL) { sim_nsites=atoi(s+11); do_sims=1;}
      if (strstr(s,"sim_nintsites=") !=NULL) { sim_nintsites=atoi(s+14); do_sims=1;}
      if (strstr(s,"sim_nintvisits=")!=NULL) { sim_nintvisits=atoi(s+15); do_sims=1;}
      if (strstr(s,"sim_nother=")    !=NULL) { sim_nother=atoi(s+11); do_sims=1;}
      if (strstr(s,"sim_removal")    !=NULL) { sim_removal=1; do_sims=1;}
      if (strstr(s,"sim_nsims=")     !=NULL) { sim_nsims=atoi(s+10); do_sims=1;}
      if (strstr(s,"sim_fname=")     !=NULL) { strcpy(sim_fname,s+10); do_sims=1;}
    }
    SelectStream(0); PutSeed(SingSpec.seed);

    eps=pow(pow(10.0, -SingSpec.LikeNRSig), 1.0/3.0);
    fprintf(g,"varcov: nsig=%d eps=%e\n",SingSpec.LikeNRSig,eps);

    if (do_sims>0) {
        if (SingSpec.Verbose) printf("calling Simulation(psi=%f p=%f N=%d n=%d k=%d nother=%d removal=%d nsims=%d fname=%s)\n",
               sim_psi,sim_p,sim_nsites,sim_nintsites,sim_nintvisits,sim_nother,sim_removal,sim_nsims,sim_fname);
        Simulation(sim_psi,sim_p,sim_nsites,sim_nintsites,sim_nintvisits,sim_nother,sim_removal,sim_nsims,sim_fname);
        return(0);
    }
    N=SingSpec.N=*NN; T=SingSpec.T=*kk;
    if (SingSpec.Verbose>0)printf("N,T-->%d,%d\n",N,T);
	if (N<1) {
		printf("\n******* Input error - number of sites = zero ****\n");
		fprintf(g,"\n******* Input error - number of sites = zero ****\n");
		return(1);
	}
    if (SingSpec.Verbose>1)printf("nearzero=%g\n",SingSpec.nearzero);
    SingSpec.Data         = new int*[N];    SingSpec.OpenData     = new int*[N];
    SingSpec.area         = new double[N];  SingSpec.det_hist_frq =new double[N];
    SingSpec.expval       = new double[N+1];
    SingSpec.neighbor     = new char*[N];   SingSpec.neighbor_wgt = new double[N];
    SingSpec.neighbor_lst = new int*[N];
	for (i=ii=0; i<N; i++) {
		SingSpec.neighbor[i]=new char[N]; SingSpec.neighbor_lst[i]=new int[N];
		for (j=0; j<N; j++) {SingSpec.neighbor[i][j]=neighbor[ii++]; SingSpec.neighbor_lst[i][j]=-1;}
	}
    k=N; if (SingSpec.UseNeighbor!=2) k=1;
    for (i=ii=0; i<N; i++) SingSpec.neighbor_wgt[i]=neighbor_wgt[ii++];  // wgt[i]= 1/((double)N);
    i=sizeof(s); gof_err=0;
    //              ****************    get detection history data *********************
    fprintf(g,"\n\n");
    for (ii=i=0,csum=0; i<N; i++) {
        SingSpec.area[i]=1; 	SingSpec.Data[i] = new int[T+1]; SingSpec.OpenData[i] = new int[T+1];
        for (j=0; j<=T; j++) SingSpec.Data[i][j]=SingSpec.OpenData[i][j]=-1;
		seen=notallmiss=0;
        for (j=0; j<T; j++,ii++) {
            if (dd[ii]<0) SingSpec.NMiss++;
            else { SingSpec.Data[i][j]=dd[ii]; notallmiss=1;}
            if(SingSpec.Data[i][j]>0) seen=1;
            if(SingSpec.Data[i][j]>maxstate) maxstate=SingSpec.Data[i][j];
            csum+=(j+1)*('2'+SingSpec.Data[i][j]); if (csum>=65536) csum-=65536;
        }
        SingSpec.det_hist_frq[i]=frq[i];
        totsites+=SingSpec.det_hist_frq[i];
        if (seen>0) n+=SingSpec.det_hist_frq[i];
        if (notallmiss==1) tsites+=SingSpec.det_hist_frq[i];
        if (SingSpec.det_hist_frq[i]>0 && SingSpec.det_hist_frq[i]!=1) gof_err=1;
    }
    if (prnt>=1) {
        for (i=0; i<N; i++) { fprintf(g,"%4d:",i+1);
            for (j=0; j<T; j++) fprintf(g,"%c",(SingSpec.Data[i][j]>=0?SingSpec.Data[i][j]+'0':'.'));
            if (SingSpec.det_hist_frq[0]!=1) fprintf(g," %f",SingSpec.det_hist_frq[i]);
            fprintf(g,"\n");
        }
    }
    SingSpec.Nstates=maxstate+1;
    printf("Data checksum = %d\n\n",csum);
    printf("\n\n********* Input Data summary *******\n");
    printf("Number of sites                = %d\n", (int)(totsites+.5));
    printf("Number of sampling occasions   = %d\n", T);
    printf("Number of states               = %d\n", SingSpec.Nstates);
    printf("Number of missing observations = %d\n", SingSpec.NMiss);
    printf("Data checksum = %d\n\n",csum);		
    fprintf(g,"\n\n********* Input Data summary *******\n");
    fprintf(g,"Number of sites                = %d\n", (int)(totsites+.5));
    fprintf(g,"Number of sampling occasions   = %d\n", T);
    fprintf(g,"Number of states               = %d\n", SingSpec.Nstates);
    fprintf(g,"Number of missing observations = %d\n", SingSpec.NMiss);
    fprintf(g,"Data checksum = %d\n\n",csum);

    SingSpec.NSiteCov=nsitecov=*nunitcov; SingSpec.NSampCov=nsampCov=*nsurvcov; nsites=N;  nn=nsites;
    fprintf(g,"NSiteCovs-->%d\n",SingSpec.NSiteCov);
    if(SingSpec.Verbose>0) printf("NSiteCovs-->%d\n",SingSpec.NSiteCov);
    SingSpec.CovNames=new char*[nsitecov];
	SingSpec.nphantom=0; nn=N+SingSpec.nphantom; SingSpec.SiteCov=new double*[nn];
    for (i=0; i<nsitecov; i++) {
        SingSpec.CovNames[i]=new char[64];
		strcpy(SingSpec.CovNames[i],unitcovnames[i]);
        if (SingSpec.Verbose>0) printf("site_covname[%d]=%s\n",i,SingSpec.CovNames[i]);
        fprintf(g,"site_covname[%d]=%s\n",i,SingSpec.CovNames[i]);
    }
    for (ii=i=0; i<nn; i++) {
        SingSpec.SiteCov[i]=new double[nsitecov];
        for (j=0; j<*nunitcov; j++) SingSpec.SiteCov[i][j]=unitcov[ii++];
    }
    //           ************** get sample covars ********

    SingSpec.NSampCov=nsampCov=*nsurvcov;
	if(SingSpec.Verbose>0) printf("NSampCovs-->%d\n",nsampCov);
    fprintf(g,"NSampCovs-->%d\n",nsampCov);
    SingSpec.CovNames2=new char*[nsampCov]; SingSpec.SampCov=new double**[nsampCov];
    for (ii=k=0; k<nsampCov; k++) {
        SingSpec.CovNames2[k]=new char[64];
        strcpy(SingSpec.CovNames2[k],survcovnames[k]);
        if (SingSpec.Verbose>0) printf("samp_covname[%d]=%s\n",k,SingSpec.CovNames2[k]);
        fprintf(g,"samp_covname[%d]=%s\n",k,SingSpec.CovNames2[k]);
        SingSpec.SampCov[k]=new double*[N];
        for (j=0; j<nsites; j++) SingSpec.SampCov[k][j]=new double[T];
	}
	for (i=0; i<T; i++)
		for (j=0; j<nsites; j++)
			for (k=0; k<nsampCov; k++) {
				SingSpec.SampCov[k][j][i]=survcov[ii++];
				if (strcmp(survcovnames[k],"confirm")==0)
					SingSpec.Data[j][i]= ((int)SingSpec.SampCov[k][j][i])<<2 | SingSpec.Data[j][i];
			}

    if (prnt) fprintf(g,"\nSite Covariates:\n");
    if (SingSpec.NSiteCov>0) {
        if (prnt>0) {
			for (i=0; i<SingSpec.NSiteCov; i++) fprintf(g,"%s ",SingSpec.CovNames[i]);
            fprintf(g,"\n");
            for (j=0; j<nsites; j++) { fprintf(g,"%d:",j);
                for (i=0; i<SingSpec.NSiteCov; i++) fprintf(g,"%f ",SingSpec.SiteCov[j][i]);
                fprintf(g,"\n");
            }
        }
	}
    if (prnt) fprintf(g,"\n Sample Covariates:");
    for (k=0; k<SingSpec.NSampCov; k++) {
        if (prnt) fprintf(g,"\n  Covariate %d: %s\n",k,SingSpec.CovNames2[k]);
        for (j=0; j<nsites; j++) {
            for (i=0; i<T; i++) {
                if (prnt) fprintf(g,"%1.1f ",SingSpec.SampCov[k][j][i]);
                if (SingSpec.SampCov[k][j][i]<-9998 && SingSpec.Data[j][i]>=0) {
                    fprintf(g,"\n***\nWARNING: missing covariate (%s) with non-missing data\n***         site:%d survey:%d\n***\n",SingSpec.CovNames2[k],j,i);
                }
            }
            if (prnt) fprintf(g,"\n");
        }
    }
	//  *********    get secondary periods (and number of methods)  ******************

    SingSpec.SecPeriods = new int[T+1]; SingSpec.PrmyPeriods=*nseasons;
	for (i=0; i<=*nseasons; i++) SingSpec.SecPeriods[i]=nsurveyseason[i];
	SingSpec.NMethods=*nmeth;
    if (SingSpec.Verbose>0) printf("Primary periods=%d Secondary periods:",SingSpec.PrmyPeriods);
    fprintf(g,"Primary periods=%d Secondary periods:",SingSpec.PrmyPeriods);
    for (i=0;i<SingSpec.PrmyPeriods;i++) {
	  if (SingSpec.Verbose>0) printf(" %d",SingSpec.SecPeriods[i]);
	  fprintf(g," %d",SingSpec.SecPeriods[i]);
	}
    fprintf(g,"\n");
	if (SingSpec.Verbose>0) printf("\nNumber of methods/survey = %d\n",SingSpec.NMethods); 
	//       **************  get sitenames *****************8
    SingSpec.sitename=new char*[N]; SingSpec.surveyname=new char*[T]; printf("N=%d T=%d\n",N,T);
    for (i=0; i<N; i++) { SingSpec.sitename[i]=new char[1024]; strcpy(SingSpec.sitename[i],unitnames[i]); }
	//   **********************   get survey names *****************
    for (i=0; i<T; i++) { 
		SingSpec.surveyname[i]=new char[1024]; 
		strcpy(SingSpec.surveyname[i],surveynames[i]); 
	} 

    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    fprintf(g," %s\n",title[0]);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    time(&t2); if (SingSpec.Verbose>0) printf("N=%d T=%d Groups=%d bootstraps=%d\n",N,T,SingSpec.Groups,nbootstraps);
    fprintf(g,"N=%d T=%d Groups=%d bootstraps=%d\n",N,T,SingSpec.Groups,nbootstraps);
    //printf("reading design matrix file...\n");
	NrealPar=NPar=0;
	for (ii=jj=dmnum=0; dmnum<6; dmnum++) {
		irows=dmdims[dmnum*2]; icols=dmdims[dmnum*2+1];
		if (SingSpec.Verbose>0) printf("\nMatrix %d: rows=%d, cols=%d modtype=%d\n",dmnum+1,irows,icols,modtype);
		fprintf(g,"\nMatrix %d: rows=%d, cols=%d\n",dmnum+1,irows+(irows>0),icols+(icols>0));
		if (icols>0) {
		  if (SingSpec.Verbose) printf("            -");
		  for (j=0; j<icols; j++) if (SingSpec.Verbose) printf(",%c%d",'a'+dmnum,j+1);
		  if (SingSpec.Verbose) printf("\n");
		  fprintf(g,"            -"); for (j=0; j<icols; j++) fprintf(g,",%c%d",'a'+dmnum,j+1); fprintf(g,"\n");
		}
		SingSpec.NrowsDM[dmnum]=irows; NrealPar+=irows;
		SingSpec.NParKK[dmnum]=icols; NPar+=icols;
        SingSpec.DMat[dmnum] = new double*[irows+1];   //  SingSpec.DMat, SingSpec.DMat_ptr used to obtain real value in design matrix
        SingSpec.DMat_ptr[dmnum] = new int*[irows+1];  //  If DMat_ptr==0, use value in SingSpec.DMat
        SingSpec.Realname[dmnum]=new char*[irows+1];   //  If DMat_ptr==3001, use "psi1" or "PsiA" in des.mat.
        SingSpec.Betaname[dmnum]=new char*[icols+1];   //  If DMat_ptr==3002, use "psiB" in des.mat.
        SingSpec.BetaFixed[dmnum]=new double[icols+1]; //  if Dmat_ptr==i (where i>0), use i'th site covar.
        SingSpec.LnkFn[dmnum]=new int[irows+1];        //  if Dmat_ptr==i (i<0), use -i'th sample covar.
        for (j=0; j<icols; j++) {
		    SingSpec.Betaname[dmnum][j]=new char[1024];
			strcpy(SingSpec.Betaname[dmnum][j],"x");
            SingSpec.BetaFixed[dmnum][j]=1.1e44;
		}
		if (irows<1) jj++;   //  des.mat from R cannot be NULL, so each des.mat has at least a "0"
		for (i=0; i<irows; i++) {
            SingSpec.DMat[dmnum][i]=new double[icols+1]; SingSpec.DMat_ptr[dmnum][i]=new int[icols+1];
            SingSpec.Realname[dmnum][i]=new char[64];   strcpy(SingSpec.Realname[dmnum][i],rnames0[ii++]);
			fprintf(g,"%-13s",SingSpec.Realname[dmnum][i]); if (SingSpec.Verbose) printf("%-13s",SingSpec.Realname[dmnum][i]);
			j=chk_real_labels(SingSpec.Realname[dmnum][i],&modtype, i, dmnum);
            for (j=0; j<icols; j++,jj++) {
				strcpy(tmps,""); SingSpec.DMat[dmnum][i][j]=atof(dm0[jj]);
				SingSpec.DMat_ptr[dmnum][i][j]=0;      //  check list of site covars
                for (k=0; k<SingSpec.NSiteCov; k++) if (strcmp(dm0[jj],SingSpec.CovNames[k])==0) break;
				if (k<SingSpec.NSiteCov) { SingSpec.DMat_ptr[dmnum][i][j]=k+1; strcpy(tmps,SingSpec.CovNames[k]); }
                else {				 //  check list of sample covars
                    for (k=0; k<SingSpec.NSampCov; k++) if (strcmp(dm0[jj],SingSpec.CovNames2[k])==0) break;
					if (k<SingSpec.NSampCov) { SingSpec.DMat_ptr[dmnum][i][j]=-(k+1); strcpy(tmps,SingSpec.CovNames2[k]); }
                    else {				//  check list of autologistic covars
				        if (strncmp(dm0[jj],"psi1",4)==0) { SingSpec.DMat_ptr[dmnum][i][j]=3001; SingSpec.Model=5;}
                        if (strncmp(dm0[jj],"psi1",4)==0 || strncmp(dm0[jj],"psiA",4)==0) { SingSpec.UseNeighbor=2;}
                        if (strncmp(dm0[jj],"psiB",4)==0) { SingSpec.DMat_ptr[dmnum][i][j]=3002;}
                        if (strncmp(dm0[jj],"upsi1",5)==0) { SingSpec.Model=5; SingSpec.uncondpsi=1; }
				        for (k=0; k<4; k++)  if (strcmp(dm0[jj],autocovs[k])==0) break;
						if (k<4) {
				        	SingSpec.DMat_ptr[dmnum][i][j]=3001+(k>2);
				        	SingSpec.UseNeighbor=2;
				        	if (k<3) SingSpec.Model=5;
				        	if (k==2) SingSpec.uncondpsi=1;
				        	strcpy(tmps,autocovs[k]);
				        }
					}
				}
				if ((SingSpec.DMat_ptr[dmnum][i][j]!=0) | (SingSpec.DMat[dmnum][i][j]!=0))
					if (strlen(SingSpec.Betaname[dmnum][j])<2) {
						if (dm0[jj][0]!='1') sprintf(SingSpec.Betaname[dmnum][j],"%s.%s",SingSpec.Realname[dmnum][i],dm0[jj]);
					    else sprintf(SingSpec.Betaname[dmnum][j],"%s",SingSpec.Realname[dmnum][i]);
					}
				fprintf(g,"  %s",dm0[jj]); if (SingSpec.Verbose) printf(" %s",dm0[jj]);
			}
            fprintf(g,"\n"); if (SingSpec.Verbose) printf("\n");
            SingSpec.LnkFn[dmnum][i]=0;
		}
		for (j=0; j<icols; j++) if (strlen(SingSpec.Betaname[dmnum][j])<2) sprintf(SingSpec.Betaname[dmnum][j],"%c%d",'a'+dmnum,j+1);
	}
	SingSpec.fixed = new double[NrealPar];  for (i=0; i<NrealPar; i++) SingSpec.fixed[i]=-999;
	j=check_desmat_rowlabels(&modtype);
    Params=new double[NPar+1]; for (i=0; i<=NPar; i++) Params[i]=0;
    //readdesmat(&NPar,&Params,f,&modtype); fclose(f);
	//                  allocate space for real param estimates for each site
	SingSpec.realParmEst = new double*[NrealPar+1];
	SingSpec.finalBetaEst = new double[NPar+1];
	for (i=0; i<NrealPar; i++) SingSpec.realParmEst[i] = new double[N]; SingSpec.realParmEst[NrealPar]=NULL;
	SingSpec.finalVCmat = new double*[NPar]; for (i=0; i<NPar; i++) SingSpec.finalVCmat[i]=new double[NPar];
	if (modtype==7 && SingSpec.NrowsDM[4]>0) modtype=10;

	for (i=k=0; i<6; i++) {
	    for (j=0; j<SingSpec.NrowsDM[i]; j++,k++) {
			SingSpec.fixed[k]=fv[k];
			if (fv[k]>-990) fprintf(g,"fix(%d) %s=%f\n",k,SingSpec.Realname[i][j],fv[k]);
		}
	}
	for (i=0; i<NPar; i++) Params[i]=iv[i];

	if (SingSpec.Verbose>0) printf("modtype=%d\n",modtype);
    if (SingSpec.PrmyPeriods<2 && modtype!=3) fprintf(g,"Naive occupancy estimate       = %6.4f\n\n",n / tsites);
	get_opendata(modtype); fprintf(g,"\nmodtype=%d\n",modtype);
    if (modtype==1 && nbootstraps>0 && gof_err>0) {
        fprintf(g,"\n============================================================\n");
        fprintf(g,"Note: GOF only works with unsummarized data\n");
        fprintf(g,"============================================================\n");
        nbootstraps=0;
    }
	if (SingSpec.UseNeighbor>0) {
		for (i=0; i<N; i++)
			for (j=k=0; j<N; j++)
				if (SingSpec.neighbor[i][j]>0) SingSpec.neighbor_lst[i][k++]=j;
	}
    switch (modtype) {
	  case 9:
    case 1:    // single-season model  (possibly multi-method
		if (CustomModel) CustomMod(nbootstraps,nboot2,NPar,Params,poolcutoff);
        else PreDefMod(nbootstraps,nboot2);
        break;
    case 2:            //  multi-season model
		if (SingSpec.NrowsDM[0]>1) {
			if (strncmp(SingSpec.Realname[0][1],"theta",5)==0) OpenPopnMM(nbootstraps,s,NPar,Params,nboot2);
			else                                                 OpenPopn(nbootstraps,s,NPar,Params,nboot2);
		} else                                                 OpenPopn(nbootstraps,s,NPar,Params,nboot2);
		break;
    case 3:
		if (!SingSpec.multiseason) TwoSpec(nbootstraps,s,NPar,Params);    //   two-species model
	    else TwoSpecMultiSeason(nbootstraps,s,NPar,Params);
		break;
    case 4: RoyleMod(nbootstraps,nboot2,NPar,Params); break;  //  Royle model

    case 5:
		if (!SingSpec.multiseason) MSMod(nbootstraps,nboot2,NPar,Params);  //  single-season,multi-strata model
    	else                       MSOpenPop(nbootstraps,s,NPar,Params);
		break;
    case 6: MSOpenPop(nbootstraps,s,NPar,Params); break;      //  multi-season,multi-strata model
    //case 7: SpaDepMod(nbootstraps,nboot2,NPar,Params,poolcutoff); break; //  single-season spatial dependence model
    case 8: StagEntModl(nbootstraps,s,NPar,Params); break; //  single-season staggered entry model
    //case 9: FalsePosMod(NPar,Params); break; //  single-season false-positive detection model
    case 7:
    case 10: OpenPopSD(NPar,Params,nboot2);  break; //  multi-season spatial-dependence model
	case 11: TwoSpecMisId(NPar,Params);
    }

    time(&t2); t3=t2-t1; if (SingSpec.Verbose>0) printf("\nCPU time= %1.0f seconds (%1.2f min)\n",t3,t3/60.);
    fprintf(g,"\nCPU time= %1.0f seconds (%1.2f min)\n",t3,t3/60.);
	fprintf(g,"\n\nPRESENCE - Presence/Absence-Site Occupancy data analysis\n");
    fprintf(g,"%s,       Version %s.\n",timstr,VERSN);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
	fclose(g); 
	rvals[0]=2*SingSpec.finalLL;      //   return values to R...
	for (i=0; i<NPar; i++) rvals[i+1]=SingSpec.finalBetaEst[i];
	for (i=0,ii=NPar+1; i<NPar; i++) for (j=0; j<NPar; j++) rvals[ii++]=SingSpec.finalVCmat[i][j];

	delete[] Params;
	for (i=0; i<N; i++) {
		delete [] SingSpec.Data[i];
		delete [] SingSpec.neighbor_lst[i];
		delete [] SingSpec.neighbor[i];
	}
	delete [] SingSpec.Data;
	delete [] SingSpec.neighbor; delete [] SingSpec.neighbor_lst;
    delete [] SingSpec.area; delete [] SingSpec.det_hist_frq;  delete [] SingSpec.expval;
    delete [] SingSpec.neighbor_wgt;
    printf("exit rpresence (main2b)\n"); return(0);
}
