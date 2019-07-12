#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>

#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif

#include "SingAnal.h"
#include "rngs.h"
#define FALSE (1==0)

TSingSpec SingSpec;

int presence(int *iargc, char **argv) {

    // Data[][] contains detection/nondetection info
    // Missing[][] contains whether observations were missing
    // N is number of sites
    // T is number of sampling occasions
    // Groups is number of groups for mixture models
    // TSpecificP indicates whether detection probabilities are time specific
	extern TSingSpec SingSpec;
    int argc=*iargc,i,j=0,k,l,nbootstraps=0,nboot2=0,prnt=0,dprnt=0,modtype=1,TT,nsites,csum,nn,NPar,phantom=0;
    int nsampCov,CustomModel,notallmiss,seen,gof_err=0,site,maxstate=0,N=SingSpec.N,T=SingSpec.T,NrealPar=0,iopt=0;
    time_t t1, t2; FILE *f,*g,*fcov; float t3;
    char s[65536], *s2, *s3, titstr[512]="no title", secper[2560]="", modsave[64]="no model name",
          cmtchar='|', fname[256]="genpres.pao", gname[256]="presence.out", jname[256]="",
		  cname[256]="",phname[256]="", *p2;
    double sim_psi=0.7, sim_p=0.4, tmp,poolcutoff=0.0,n=0,tsites=0, *Params,totsites=0, eps;
    int sim_nsites=200, sim_nintsites=100, sim_nintvisits=9, sim_nother=5, sim_removal=0, sim_nsims=100, do_sims=0;
    char sim_fname[256]="",abc[]="_abc_AFH_DBI_EGC";

    double atof_comma(char *s);
    void readdesmat(int *NPar, double **Params, FILE *f, int *modtype);
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
  
    int get_modl_num(char *s);
	void get_opendata(int modtype);

    SingSpec.DMat = new double**[9]; SingSpec.Betaname=new char **[9];  SingSpec.DMat_ptr = new int**[9];
    SingSpec.Realname=new char **[9]; SingSpec.BetaFixed=new double*[9];
    SingSpec.LnkFn=new int *[9];
    char timstr[40]; time(&t1); strcpy(timstr,ctime(&t1)); timstr[strlen(timstr)-1]='\0';
	
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

    printf(">>timstr=%s\n",timstr);
    CustomModel=true;
    SingSpec.seed=(unsigned long) -time(NULL);
	char cwd[1024];
	
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(stdout, "Current working dir: %s\n", cwd);
	else
		perror("getcwd() error");

    for (i=0; i<6; i++) SingSpec.NParKK[i]=SingSpec.NrowsDM[i]=0;
    for (i=0; i<argc; i++) {
			strcpy(s,argv[i]);
        if (strstr(s,"quiet")!=NULL || strstr(s,"QUIET")!=NULL) SingSpec.Verbose=0;
		if (strstr(s,"verbose")!=NULL) SingSpec.Verbose=2;
        if (strstr(s,"debug")!=NULL) SingSpec.Verbose=3;
		if (SingSpec.Verbose>0) printf("%d:%s\n",i,argv[i]); 
    }
	if (SingSpec.Verbose>0) printf("PRESENCE Version %s.\n",VERSN);	
    for (i=0; i<argc; i++) {
      strcpy(s,argv[i]); j=(int)strlen(s);
      if (SingSpec.Verbose>0) printf("==>%s\n",s);
      for (k=0; k<j && s[k]!='='; k++) s[k]=tolower(s[k]);
      if (strstr(s,"echo")!=NULL) dprnt=1;
		  if (s[1]=='=') {
			  if (s[0]=='i') strcpy(fname,s+2);   //  input pao filename
			  if (s[0]=='j') strcpy(jname,s+2);   //  input design matrix filename
			  if (s[0]=='l' || s[0]=='o') strcpy(gname,s+2);   //    output filename
			  if (s[0]=='t') { strcpy(secper,s+2); printf("s2=%s\n",secper);}
			  if (s[0]=='c') { strcpy(cname,s+2); } //   input separate site covariate filename
			  if (s[0]=='f') { strcpy(phname,s+2); phantom=1; } //   input extra phantom detection file (site covar data only)
		  }
      if (strstr(s,"name=")!=NULL) strcpy(SingSpec.modname,s+5);
      if (strstr(s,"maxfn=")!=NULL)  SingSpec.maxfn=atoi(s+6);
      if (strstr(s,"lmt=")!=NULL && strstr(s,"rlmt")==NULL) SingSpec.lmt=atoi(s+4);
      if (strstr(s,"rlmt=")!=NULL) SingSpec.rlmt=atoi(s+5);
      if (strstr(s,"link=c")!=NULL) SingSpec.LinkFn=2;
      if (strstr(s,"link=s")!=NULL) SingSpec.LinkFn=4;
      if (strstr(s,"nsig=")!=NULL) SingSpec.LikeNRSig=atoi(s+5);
      if (strstr(s,"zero=")!=NULL) SingSpec.nearzero=pow(10.,(double)(-atoi(s+5)));
      if (strstr(s,"nboot=")!=NULL) nbootstraps=atoi(s+6);
      if (strstr(s,"boot2=")!=NULL) nboot2=atoi(s+6);
      if (strstr(s,"pool=")!=NULL) poolcutoff=atof(s+5);
      if (strstr(s,"nstates=")!=NULL) SingSpec.Nstates=atoi(s+8);
      if (strstr(s,"amoeba")!=NULL) SingSpec.UseAmoeba=1;
      if (strstr(s,"seed")!=NULL) SingSpec.seed = atol(s+5);
      if (strstr(s,"vceps")!=NULL) SingSpec.VCeps = atof(s+6);
      if (strstr(s,"neighb")!=NULL) SingSpec.UseNeighbor = 2;
      if (strstr(s,"nrand=")!=NULL) SingSpec.nrandIV=atoi(s+6);
	  if (strstr(s,"nthreads=")  !=NULL) SingSpec.nthreads=atoi(s+9);
	
//        output options for SE's and VC matrix
//	  Option 1:  only compute likelihood and beta estimates, no beta var-cov matrix or beta SE's (useful for faster model selection of big models or simulations)
//    Option 2: compute but don't print beta var-cov matrix so we get beta SE's, don't compute real params
//    Option 3: compute and print beta beta var-cov matrix, don't compute real params (I think this is what you want for RPresence)
//    Option 4: (default) compute real params and real var-cov matrix, don't print either var-cov matrix.
//    Option 5: compute real params and real var-cov matrix, print only beta var-cov matrix.
//    Option 6: print both var-cov matrices.

      if (strstr(s,"nose")!=NULL) SingSpec.novar = 1;
      if (strstr(s,"betavc")!=NULL) SingSpec.novar = 2;
      if (strstr(s,"noreal")!=NULL) SingSpec.novar = 3;
      if (strstr(s,"realbetavc")!=NULL) SingSpec.novar = 5;
      if (strstr(s,"bothvc")!=NULL) SingSpec.novar = 6;


      if (strstr(s,"chsq")!=NULL) SingSpec.do_chisq = 1;
      if (strstr(s,"noderived")!=NULL) SingSpec.prnt_derived = false;
      if (strstr(s,"nodups")!=NULL) SingSpec.dups = false;
      if (strstr(s,"model=")!=NULL) { strcpy(modsave,s); modtype=get_modl_num(s); if (modtype==1) CustomModel=(modsave[7]=='0'); }
      if (strstr(s,"sim_psi=")!=NULL) { sim_psi=atof(s+8); do_sims=1;}
      if (strstr(s,"sim_p=")!=NULL) { sim_p=atof(s+6); do_sims=1;}
      if (strstr(s,"sim_nsites=")!=NULL) { sim_nsites=atoi(s+11); do_sims=1;}
      if (strstr(s,"sim_nintsites=")!=NULL) { sim_nintsites=atoi(s+14); do_sims=1;}
      if (strstr(s,"sim_nintvisits=")!=NULL) { sim_nintvisits=atoi(s+15); do_sims=1;}
      if (strstr(s,"sim_nother=")!=NULL) { sim_nother=atoi(s+11); do_sims=1;}
      if (strstr(s,"sim_removal")!=NULL) { sim_removal=1; do_sims=1;}
      if (strstr(s,"sim_nsims=")!=NULL) { sim_nsims=atoi(s+10); do_sims=1;}
      if (strstr(s,"sim_fname=")!=NULL) { strcpy(sim_fname,s+10); do_sims=1;}
      if (strstr(s,"help")!=NULL) {
          printf("usage: presence i=infile j=desmatfile l=outfile model=xxx name=model_name\n");
          printf("                maxfn=? lmt=? link=[c|s] nsig=? zero=? nboot=? boot2=?\n");
          printf("                amoebo hess verbose debug seed=? vc help\n");
          printf("                sim_psi=.7 sim_p=.4 sim_nsites=200 sim_nintsites=100 sim_nintvisits=9\n");
          printf("                sim_nother=5 sim_removal sim_nsims=100 sim_fname=?\n");
		  printf("\n       PRESENCE command-line arguments:\n");
          printf("i=infile   : name of input pao file created by PRESENCE GUI\n");
          printf("j=desmatfile : name of design matrix file \n");
          printf("l=outfile  : name of output file\n");
          printf("model=xxx  : indicate pre-defined model (by number) to run, not used for other models\n");
          printf("name=model_name : name of model to appear in output (user-defined)\n");
          printf("maxfn=?    : max number of iterations before aborting max likelihood search\n");
          printf("lmt=?      : max number of sites to print in output (to save paper) \n");
          printf("link=[c|s] : change link function to (c)loglog or (s)sin\n");
          printf("nsig=?  : number of significant digits in likelihood when maximum is assumed\n");
          printf("zero=?  : highest value above zero assumed to equal zero for exp func.(default=10^(-200))\n");
          printf("nboot=? : number of bootstrap replicates for computation of var-cov matrix\n");
          printf("nboot2=?: number of bootstrap replicates for GOF \n");
          printf("amoeba  : use simplex algorithm for starting values\n");
          printf("hess    : *	no longer used *\n");
          printf("verbose : output input info to cmd console\n");
          printf("debug   :  output debugging info to cmd console\n");
          printf("seed=?  : set random number generator seed (for simulations and/or bootstraps)\n");
          printf("vc      : print var-cov matrix of real parameters in output\n");
          printf("nose    : dont waste time computing std. errs or real params for the model (betas only)\n");
          printf("betavc  : compute but don't print beta var-cov matrix so we get beta SE's, don't compute real params\n");
          printf("noreal  : compute and print beta beta var-cov matrix, don't compute real params\n");
          printf("realbetavc: compute real params and real var-cov matrix, print only beta var-cov matrix.\n");
          printf("bothvc    : print both var-cov matrices.\n");
          printf("sim_psi, sim_p, sim_nsites,... : parameters for legacy simulation routine (in Run menu)\n");
          return(1);
      }
    }
    if ((g=fopen(gname,"wt"))==NULL) {
        printf("\nError opening output file (%s)\n",gname); return(1);
    } SingSpec.g=g;
    SelectStream(0); PutSeed(SingSpec.seed);
    fprintf(g,"\n\nPRESENCE - Presence/Absence-Site Occupancy data analysis\n");
    fprintf(g,"%s,       Version %s.\n",timstr,VERSN);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    for (i=0; i<argc; i++) fprintf(g,"==>%s\n",argv[i]);
    eps=pow(pow(10.0, -SingSpec.LikeNRSig), 1.0/3.0);
    fprintf(g,"varcov: nsig=%d eps=%e\n",SingSpec.LikeNRSig,eps);

    if (do_sims>0) {
        printf("calling Simulation(psi=%f p=%f N=%d n=%d k=%d nother=%d removal=%d nsims=%d fname=%s)\n",
               sim_psi,sim_p,sim_nsites,sim_nintsites,sim_nintvisits,sim_nother,sim_removal,sim_nsims,sim_fname);
        Simulation(sim_psi,sim_p,sim_nsites,sim_nintsites,sim_nintvisits,sim_nother,sim_removal,sim_nsims,sim_fname);
		return(0);
    }
    if ((f=fopen(fname,"rt"))==NULL) {
		printf("\nError opening input file (%s)\n",fname); return(1);
	}
    for (p2=fgets(s,sizeof(s),f); s[0]==cmtchar && p2!=NULL; p2=fgets(s,sizeof(s),f)) 
		if (SingSpec.Verbose>2) printf("..>%s",s);
    i=sscanf(s,"%d %d %d",&N,&T,&iopt); if (SingSpec.Verbose>2) printf(">>s=%s\n",s);
    fprintf(g,"%s N,T-->%d,%d\n",modsave,N,T); SingSpec.N=N; SingSpec.T=T;
    if (SingSpec.Verbose>0)printf("N,T-->%d,%d\n",N,T);
	if (N<1) {
		printf("\n******* Input error - number of sites = zero ****\n");
		fprintf(g,"\n******* Input error - number of sites = zero ****\n");
		return(1);
	}

    if (SingSpec.Verbose>1)printf("nearzero=%g\n",SingSpec.nearzero);
    printf("allocating memory: N=%d\n",N);
    SingSpec.Data = new int*[N];  SingSpec.OpenData = new int*[N];  char tab[4]="x ,"; tab[0]=9;   // delims = tab, comma, space
    SingSpec.area = new double[N]; SingSpec.det_hist_frq=new double[N];  SingSpec.expval = new double[N+1];
    SingSpec.neighbor = new char*[N]; SingSpec.neighbor_wgt = new double[N]; SingSpec.neighbor_lst = new int*[N];
	for (i=0; i<SingSpec.N; i++) {
		SingSpec.neighbor[i]=new char[SingSpec.N]; SingSpec.neighbor_lst[i]=new int[SingSpec.N];
		for (j=0; j<SingSpec.N; j++) {SingSpec.neighbor[i][j]=1; SingSpec.neighbor_lst[i][j]=-1;}
	}
    for (i=0; i<N; i++) SingSpec.neighbor_wgt[i]=1/((double)N);
    //              ****************    get detection history data *********************
    fprintf(g,"\n\n");
    for (i=0,csum=0; i<N; i++) {
        SingSpec.area[i]=1;
        p2=fgets(s,sizeof(s),f); l=strlen(s); if (s[l-1]<' ') s[l-1]='\0';
        if (SingSpec.Verbose>1) printf("p/a data(%d):%s\n",i,s);
        s2=strtok(s,tab);
		SingSpec.Data[i] = new int[T+1];	SingSpec.OpenData[i] = new int[T+1];
        for (j=0; j<=T; j++) SingSpec.Data[i][j]=SingSpec.OpenData[i][j]=-1;
		seen=notallmiss=0;
        for (j=0; j<T; j++) {
            if(s2!=NULL) {
				if (s2[0]>='A' && s2[0]<='c') {
					k=SingSpec.Data[i][j]=(int)(strchr(abc,s2[0])-abc);
				}
				else 
                if (s2[0]>='-' && s2[0]<='9') {
                    if (s2[0]=='-' || s2[0]=='.') SingSpec.NMiss++;
                    else { SingSpec.Data[i][j]=atoi(s2); notallmiss=1;}
                    if(SingSpec.Data[i][j]>0) seen=1;
                    if(SingSpec.Data[i][j]>maxstate) maxstate=SingSpec.Data[i][j];
                    if(SingSpec.Data[i][j]>maxstate) maxstate=SingSpec.Data[i][j];
                }
                else
                    fprintf(g,"\n Error: invalid data, site:%d - %s\n",i,s);
            }
            s2=strtok(NULL,tab); csum+=(j+1)*('2'+SingSpec.Data[i][j]); if (csum>=65536) csum-=65536;
        }
        SingSpec.det_hist_frq[i]=1;
        if (s2!=NULL) if (s2[0]>='0' && s2[0]<='9') {
            s3=strchr(s2,','); if (s3!=NULL) s3[0]='.';   //  change "," to "." for spanish decimal point
            SingSpec.det_hist_frq[i]=atof(s2);
        }
        totsites+=SingSpec.det_hist_frq[i];
        if (seen>0) n+=SingSpec.det_hist_frq[i];
        if (notallmiss==1) tsites+=SingSpec.det_hist_frq[i];
        if (SingSpec.det_hist_frq[i]>0 && SingSpec.det_hist_frq[i]!=1) gof_err=1;
    }
    if (dprnt>=1) {
        for (i=0; i<N; i++) { fprintf(g,"%4d:",i+1);
            for (j=0; j<T; j++) fprintf(g,"%c",(SingSpec.Data[i][j]>=0?SingSpec.Data[i][j]+'0':'.'));
            if (SingSpec.det_hist_frq[0]!=1) fprintf(g," %f",SingSpec.det_hist_frq[i]);
            fprintf(g,"\n");
        }
    }
    SingSpec.Nstates=maxstate+1;
    fprintf(g,"\n\n********* Input Data summary *******\n");
    fprintf(g,"Number of sites                = %d\n", (int)(totsites+.5));
    fprintf(g,"Number of sampling occasions   = %d\n", T);
    fprintf(g,"Number of states               = %d\n", SingSpec.Nstates);
    fprintf(g,"Number of missing observations = %d\n", SingSpec.NMiss);

    for(s2=fgets(s,sizeof(s),f); s[0]==cmtchar && strlen(s)>0; s2=fgets(s,sizeof(s),f)) 
		if (SingSpec.Verbose>1) printf("..>%s",s);
    SingSpec.NSiteCov=0; SingSpec.NSampCov=0; nsampCov=0; nsites=N;  nn=nsites;
    if (s2!=NULL) SingSpec.NSiteCov=atoi(s); //i=sscanf(s,"%d",&SingSpec.NSiteCov);
    fprintf(g,"NSiteCovs-->%d\n",SingSpec.NSiteCov);
    if(SingSpec.Verbose>0) printf("NSiteCovs-->%d\n",SingSpec.NSiteCov);
	i=j=0;
    if (strlen(cname)>1) {       // count # site covars from separate file if specified
        fcov=fopen(cname,"rt"); p2=fgets(s,sizeof(s),fcov); j=strchr(s,'\n')-s; s[j]='\0';
        s3=strtok(s,tab);
        for (i=0; s3; i++) s3=strtok(NULL,tab);
		for (j=0; (p2=fgets(s,sizeof(s),fcov)); j++)
			;
		fclose(fcov);  printf("site cov file: nsitecovs=%d sites=%d\n",i,j);
    }
    int nsitecov=SingSpec.NSiteCov+i; SingSpec.CovNames=new char*[nsitecov];
	i=j=0;
    if (strlen(phname)>1) {       // count # sites from phantom data file if specified
        fcov=fopen(phname,"rt"); p2=fgets(s,sizeof(s),fcov); j=strchr(s,'\n')-s; s[j]='\0';
        s3=strtok(s,tab);
        for (i=0; s3; i++) s3=strtok(NULL,tab);
		for (j=0; (p2=fgets(s,sizeof(s),fcov)); j++)
			;
		fclose(fcov); phantom=SingSpec.nphantom=j;
        printf("phantom data file: nsitecovs=%d sites=%d\n",i,j);
		SingSpec.phantom = new int*[phantom];
		for (i=0; i<phantom; i++) SingSpec.phantom[i] = new int[T];
    }
	nn=N+SingSpec.nphantom; SingSpec.SiteCov=new double*[nn];
    for (i=0; i<nsitecov; i++) {
        SingSpec.CovNames[i]=new char[64]; strcpy(SingSpec.CovNames[i],"");
    }
    for (i=0; i<nn; i++) {
        SingSpec.SiteCov[i]=new double[nsitecov]; 
        for (j=0; j<nsitecov; j++) SingSpec.SiteCov[i][j]=-9999;
    }

	//              ****************    get site covariate data from pao file *********************

    for (i=0; i<SingSpec.NSiteCov; i++) {
        p2=fgets(s,sizeof(s),f);  if (SingSpec.Verbose>1) printf("..>%s",s);
        s[strlen(s)-1]=0;  l=strlen(s)-1; if (l>0 && s[l]<' ') s[l]='\0';
        strcpy(SingSpec.CovNames[i],s);
        if(SingSpec.Verbose>0) printf("site_covname[%d]=%s\n",i,SingSpec.CovNames[i]);
        fprintf(g,"site_covname[%d]=%s\n",i,SingSpec.CovNames[i]); nsitecov++;
        for (site=0; site<nsites; site++) {
            p2=fgets(s,sizeof(s),f);
            if (strpbrk(s,"0123456789")!=NULL) {
                s3=strchr(s,','); if (s3!=NULL) s3[0]='.';   //  change "," to "." for spanish decimal point
                SingSpec.SiteCov[site][i]=atof(s);
            }
            if (strcmp(SingSpec.CovNames[i],"area")==0) SingSpec.area[site]=SingSpec.SiteCov[site][i];
        }
    }
    if (strlen(cname)>1) {
	//           ************** get site covars from separate file if specified ********
        fcov=fopen(cname,"rt"); p2=fgets(s,sizeof(s),fcov); j=strchr(s,'\n')-s; s[j]='\0';
        s3=strtok(s,tab);
        for (i=0; s3; i++,SingSpec.NSiteCov++) {
            strcpy(SingSpec.CovNames[SingSpec.NSiteCov],s3);
            if(SingSpec.Verbose>0)
                printf("site_covname[%d]=%s\n",SingSpec.NSiteCov,SingSpec.CovNames[SingSpec.NSiteCov]);
            fprintf(g,"site_covname[%d]=%s\n",SingSpec.NSiteCov,SingSpec.CovNames[SingSpec.NSiteCov]);
            s3=strtok(NULL,tab);
        }

        for (site=0; (p2=fgets(s,sizeof(s),fcov)); site++) {   //  for each line (site) ...
            j=strchr(s,'\n')-s; s[j]='\0'; s3=strtok(s,tab);
            for (j=0; s3; j++) {                    //   for each field (covariate) ...
                tmp=atof(s3); if(SingSpec.Verbose>0) { printf("%f ",tmp); fprintf(g,"%f ",tmp); }
				if (site<N && j<SingSpec.NSiteCov) SingSpec.SiteCov[site][j]=tmp;
				else {fprintf(g,"\nERROR in reading site covariate file\n"); return(1);}
                s3=strtok(NULL,tab);
            }
            if(SingSpec.Verbose>0) { printf("\n"); fprintf(g,"\n"); }
        }
        fclose(fcov);
    }
    if (strlen(phname)>1) {
	//           ************** get site covars from separate file if specified ********
        fcov=fopen(phname,"rt"); p2=fgets(s,sizeof(s),fcov); j=strchr(s,'\n')-s; s[j]='\0';
        s3=strtok(s,tab);
        for (site=0; (p2=fgets(s,sizeof(s),fcov)); site++) {   //  for each line (site) ...
            j=strchr(s,'\n')-s; s[j]='\0'; s3=strtok(s,tab);
            for (j=0; s3; j++) {                    //   for each field (covariate) ...
                tmp=atof(s3); if(SingSpec.Verbose>0) { printf("%f ",tmp); fprintf(g,"%f ",tmp); }
				if ((site+N)<nn && j<SingSpec.NSiteCov) SingSpec.SiteCov[site+N][j]=tmp;
				else {fprintf(g,"\nERROR in reading site covariate file\n"); return(1);}
                s3=strtok(NULL,tab);
            }
            if(SingSpec.Verbose>0) { printf("\n"); fprintf(g,"\n"); }
        }
        fclose(fcov);
    }
    //           ************** get sample covars ********
    for(s2=fgets(s,sizeof(s),f); s[0]==cmtchar && strlen(s)>0; p2=fgets(s,sizeof(s),f)) 
		if (SingSpec.Verbose>1) printf("..>%s",s);
    if (s2!=NULL) nsampCov=atoi(s); //i=sscanf(s,"%d",&nsampCov);
    SingSpec.NSampCov=nsampCov; if(SingSpec.Verbose>0) printf("NSampCovs-->%d\n",nsampCov);
    fprintf(g,"NSampCovs-->%d\n",nsampCov);
    SingSpec.CovNames2=new char*[nsampCov]; SingSpec.SampCov=new double**[nsampCov];
    for (k=0; k<nsampCov; k++) {
        SingSpec.CovNames2[k]=new char[64];
        for(p2=fgets(s,sizeof(s),f); s[0]==cmtchar && strlen(s)>0; p2=fgets(s,sizeof(s),f))
			printf("..>%s",s);
        s[strlen(s)-1]=0; l=strlen(s)-1; if (l>0 && s[l]<' ') s[l]='\0';
        strcpy(SingSpec.CovNames2[k],s);
        if(SingSpec.Verbose>0) printf("samp_covname[%d]=%s.\n",k,SingSpec.CovNames2[k]);
        fprintf(g,"samp_covname[%d]=%s\n",k,SingSpec.CovNames2[k]);
        SingSpec.SampCov[k]=new double*[N];
        for (j=0; j<nsites; j++) {
            SingSpec.SampCov[k][j]=new double[T]; for (i=0; i<T; i++) SingSpec.SampCov[k][j][i]=-9999;
            p2=fgets(s,sizeof(s),f); s2=strtok(s,tab);
            for (i=0; i<T; i++) {
                if (s2!=NULL) {
                    s3=strchr(s2,','); if (s3!=NULL) s3[0]='.';   //  change "," to "." for spanish decimal point
                    tmp=atof(s2);
                    if ((s2[0]=='-' || s2[0]=='.') && s2[1]<=' ') SingSpec.SampCov[k][j][i]=-9999;
                    else {
						SingSpec.SampCov[k][j][i]=tmp;
						if (strcmp(SingSpec.CovNames2[k],"confirm")==0) {
							SingSpec.Data[j][i]= atoi(s2)<<2 | SingSpec.Data[j][i];
						}
					}
                }
                if (s2[strlen(s2)-1]!='\n')	s2=strtok(NULL,tab);
            }
        }
    }
    if (prnt>0) {
        fprintf(g,"\nSite Covariates:\n");
        if (SingSpec.NSiteCov>0) {
            for (i=0; i<SingSpec.NSiteCov; i++) fprintf(g,"%s ",SingSpec.CovNames[i]);
            fprintf(g,"\n");
            for (j=0; j<nsites; j++) { 
				fprintf(g,"%d:",j);
                for (i=0; i<SingSpec.NSiteCov; i++) fprintf(g,"%f ",SingSpec.SiteCov[j][i]);
                fprintf(g,"\n");
            }
        }
        else fprintf(g,"\nNo Site covariates\n");
        fprintf(g,"\n Sample Covariates:");
        for (k=0; k<SingSpec.NSampCov; k++) {
            fprintf(g,"\n  Covariate %d: %s\n",k,SingSpec.CovNames2[k]);
            for (j=0; j<nsites; j++) {
                for (i=0; i<T; i++) {
                    fprintf(g,"%1.1f ",SingSpec.SampCov[k][j][i]);
                    if (SingSpec.SampCov[k][j][i]<-9998 && SingSpec.Data[j][i]>=0) {
                        fprintf(g,"\n***\nWARNING: missing covariate with non-missing data\n");
                        fprintf(g,"***         site:%d survey:%d\n***\n",j,i);
                    }
                }
                fprintf(g,"\n");
            }
        }
    }
    if (!feof(f)) {p2=fgets(secper,sizeof(secper),f); j=strchr(secper,'\n')-secper; if(j>0)secper[j]='\0';}
    else strcpy(secper,"1");
    if (!feof(f)) {p2=fgets(titstr,sizeof(titstr),f);j=strchr(titstr,'\n')-titstr; if(j>0)titstr[j]='\0';}
    else strcpy(titstr,"no title specified");

	//  *********    get secondary periods (and number of methods)  ******************

    SingSpec.SecPeriods = new int[T+1]; for (i=0; i<=T; i++) SingSpec.SecPeriods[i]=1;
	if (strchr(secper,';')!=NULL) {
		strcpy(s,secper); s2=strtok(s,";"); s2=strtok(NULL,";"); if (strlen(s2)>0) SingSpec.NMethods=atoi(s2);
	}
	if (SingSpec.NMethods<2) SingSpec.NMethods=1;
    s2=strtok(secper,",");
	if (strlen(secper)>0) {
		j=atoi(s2);
		for(TT=T,i=1;s2;i++) {
			if(atoi(s2)>0) j=atoi(s2); if (j<1) j=1;
			if (i<=T) SingSpec.SecPeriods[i-1]=j; TT-=j;
			s2=strtok(NULL,",");
		}
		if (i<=T && j>1) SingSpec.SecPeriods[i-1]=j;
		for (;TT>0;TT-=j) if (i<=T) SingSpec.SecPeriods[i++]=j;
		SingSpec.PrmyPeriods=--i; if (i>T) i=T;
	}
    printf("Primary periods=%d Secondary periods:",SingSpec.PrmyPeriods);
    fprintf(g,"Primary periods=%d Secondary periods:",SingSpec.PrmyPeriods);
    for (i=0;i<SingSpec.PrmyPeriods;i++) {
		printf(" %d",SingSpec.SecPeriods[i]);
		fprintf(g," %d",SingSpec.SecPeriods[i]);
	}
    printf("\n"); fprintf(g,"\n");
	printf("Number of methods/survey = %d\n",SingSpec.NMethods);
	//       **************  get sitenames *****************8
    SingSpec.sitename=new char*[N]; SingSpec.surveyname=new char*[T];
    strcpy(s,"");
    for (i=0; i<N; i++) {
        SingSpec.sitename[i]=new char[32]; sprintf(SingSpec.sitename[i],"site_%d",i+1);
        if (!feof(f)) {
            p2=fgets(s,sizeof(s),f); l=strlen(s)-1;
            if (l>0) {
                j=strchr(s,'\n')-s; if (j>=0) s[j]='\0';
                if (s[l]<' ') s[l]='\0';
            	strncpy(SingSpec.sitename[i],s,31); s[31]='\0';
            }
        }
    }

	//   **********************   get survey names *****************
    for (i=0; i<T; i++) {
        SingSpec.surveyname[i]=new char[32];
        sprintf(SingSpec.surveyname[i],"survey_%d",i+1);
        if (!feof(f)) {
            p2=fgets(s,sizeof(s),f); l=strlen(s)-1;
            if (l>0) {
            	j=strchr(s,'\n')-s;
				if (j>=0 && j<32) s[j]='\0';
            	if (s[l]<' ') s[l]='\0';
            	if (j>0 && j<32) strcpy(SingSpec.surveyname[i],s);
			}
        }
    }
    fclose(f); time(&t2); 
	fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n %s\n",titstr);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    fprintf(g,"N=%d T=%d Groups=%d bootstraps=%d\n\n-->%s",N,T,SingSpec.Groups,nbootstraps,s);

    if (modtype!=1 || CustomModel) {
        if (strlen(jname)<1) {           //  if no design matrix file specified...
            strcpy(jname,fname); i=strlen(jname);   //   use pao filename, replacing 'pao' with 'dm'
            jname[i-3]='d'; jname[i-2]='m'; jname[i-1]='\0';
        }
        if (strlen(jname)>1) {
            printf("opening %s...\n",jname);
            if ((f=fopen(jname,"rt"))==NULL) {
                printf("\nError opening input dm file (%s)\n",jname);
                fprintf(g,"\nError opening input dm file (%s)\n",jname);
                return(1);
            }
            printf("reading design matrix file...\n");
            readdesmat(&NPar,&Params,f,&modtype); fclose(f); printf("done\n");
        }
    }
    if(SingSpec.Verbose>0) printf("==>dmfile=[%s] modtype=%d\n",jname,modtype);
    if (strlen(jname)>0 && SingSpec.NrowsDM[0]>1) {
		printf("nrowsdm=%d %d %d %d %d %d\n",SingSpec.NrowsDM[0],SingSpec.NrowsDM[1],SingSpec.NrowsDM[2],SingSpec.NrowsDM[3],SingSpec.NrowsDM[4],SingSpec.NrowsDM[5]);
	    if (strncmp(SingSpec.Realname[0][1],"thet",4)==0 && strncmp(SingSpec.Realname[1][0],"p1",2)==0)	 {
			SingSpec.NMethods=T/(SingSpec.NrowsDM[0]-1);
			printf("\nNumber of Methods per survey = %d\n",SingSpec.NMethods);
			fprintf(g,"\nNumber of Methods per survey = %d\n",SingSpec.NMethods);
		}
	}
	//                  allocate space for real param estimates for each site
	
    if (modtype==1 && !CustomModel) {
	    NPar =  1 /*psi*/ + (SingSpec.Groups-1) + SingSpec.Groups*(1 + (SingSpec.TSpecificP*(SingSpec.T-1)));
		SingSpec.NrowsDM[0]=SingSpec.NrowsDM[1]=1; 
		SingSpec.DMat[0]=new double*[1]; SingSpec.DMat[0][0]=new double[1]; SingSpec.DMat[0][0][0]=1;
		SingSpec.DMat[1]=new double*[SingSpec.T];
		for (i=0; i<SingSpec.T; i++) SingSpec.DMat[1][i]=new double[1]; SingSpec.DMat[1][i][0]=1;
    }
		  
	for (i=0; i<6; i++) NrealPar+=SingSpec.NrowsDM[i];
	SingSpec.realParmEst = new double*[NrealPar+1];
	SingSpec.finalBetaEst = new double[NPar+1];
	for (i=0; i<NrealPar; i++) SingSpec.realParmEst[i] = new double[N]; SingSpec.realParmEst[NrealPar]=NULL;
	SingSpec.finalVCmat = new double*[NPar]; for (i=0; i<NPar; i++) SingSpec.finalVCmat[i]=new double[NPar];
	printf("modtype=%d\n",modtype);
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
		if (SingSpec.NrowsDM[0]>2) {
			if (strncmp(SingSpec.Realname[0][1],"theta",5)==0) OpenPopnMM(nbootstraps,s,NPar,Params,nboot2);
			else                                                 OpenPopn(nbootstraps,s,NPar,Params,nboot2); 
		} else  OpenPopn(nbootstraps,s,NPar,Params,nboot2); 
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
    time(&t2); t3=t2-t1; printf("\nCPU time= %1.0f seconds (%1.2f min)\n",t3,t3/60.);
    fprintf(g,"\nCPU time= %1.0f seconds (%1.2f min)\n",t3,t3/60.);
    fprintf(g,"\n\nPRESENCE - Presence/Absence-Site Occupancy data analysis\n");
    fprintf(g,"%s,       Version %s.\n",timstr,VERSN);
    fprintf(g,"- - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");	
	fclose(g); if (SingSpec.Verbose>1) printf("free memory...\n");

	delete[] Params;
	for (i=0; i<SingSpec.N; i++) delete [] SingSpec.Data[i];
	for (i=0; i<N; i++) {
		delete [] SingSpec.neighbor_lst[i];
		delete [] SingSpec.neighbor[i];
	}
	delete [] SingSpec.Data;
	delete [] SingSpec.neighbor; delete [] SingSpec.neighbor_lst;
    delete [] SingSpec.area; delete [] SingSpec.det_hist_frq;  delete [] SingSpec.expval;
    delete [] SingSpec.neighbor_wgt; if (SingSpec.Verbose>0) printf("done\n");
    return(0);
}

int covIndex(char **covlist, char *srchstr) {
    int i,j=-1;
    for (i=0; strlen(covlist[i])>0; i++) {
        printf("i=%d covlist=%s\n",i,covlist[i]);
        if (strcmp(srchstr,covlist[i])==0) j=i;
    }
    if (j<0) { strcpy(covlist[i],srchstr); j=i; strcpy(covlist[i+1],"");}
    return(j);
}
int covlist_sort(char **covlist) {
    int i,j,n,lo; char tmp[512];
    for (n=0; strlen(covlist[n])>0; n++);
    for (i=0; i<(n-1); i++) { lo=i;
        for (j=i+1; j<n; j++) if (strcmp(covlist[j],covlist[lo])<0) lo=j;
        strcpy(tmp,covlist[i]); strcpy(covlist[i],covlist[lo]); strcpy(covlist[lo],tmp);
    }
    return(n);
}

int get_modl_num(char *s) { int modtype;
    // get which parameterization; 1=psi,gamma,phi,p; 2=psi,gamma,p; 3=psi,phi,p
    modtype=s[6]-'0';
    switch (modtype) {
    case 1:      //  single-season model (and possibly multi-method
        if (s[7]>'0' && s[7]<'4') SingSpec.Groups=s[7]-'0';
        if (s[7]=='4') modtype=4;
        if (s[8]>='2' && s[8]<='9') SingSpec.NMethods=s[8]-'0';
        if (s[8]=='t') SingSpec.TSpecificP=true;
        if (s[8]=='m') SingSpec.MissClass=true;
        if (s[8]=='f') SingSpec.FalsePos=true;
        break;
    case 2:     //   multi-season models
	    if (s[7]=='6') {SingSpec.Model=6;}
        if (s[7]=='5') {SingSpec.Model=5;}
        if (s[7]=='4') {SingSpec.Model=4;}
        if (s[7]=='3') {SingSpec.Model=3;}
        if (s[7]=='2') { SingSpec.Model=2;}
        if (s[8]=='1') SingSpec.Alt_Param_Checked=1;
        break;
    case 3:  // 2-species model
        if (s[10]>'0') SingSpec.Alt_Param_Checked=s[10]-'0';
        break;
    case 4:   //  Royle model
        if (s[7]=='0') SingSpec.Alt_Param_Checked=0;
        if (s[7]=='1') SingSpec.Alt_Param_Checked=1;
        if (s[7]=='2') { SingSpec.TSpecificP=true; SingSpec.Alt_Param_Checked=2;}
        if (s[7]=='3') SingSpec.Alt_Param_Checked=3;
        break;
    case 5:    //  single-season multi-state model
        if (s[8]=='1') SingSpec.Alt_Param_Checked=1;
        break;
    case 6: SingSpec.Alt_Param_Checked=s[7]-'0'; break;  //  multi-season,multi-strata model
    case 7:
        if (s[7]=='1') SingSpec.Alt_Param_Checked=1;   //  single-season, spatial dependence model
        if (s[8]>='2' && s[8]<='9') SingSpec.NMethods=s[8]-'0';
        break;
    }
    return(modtype);
}

void get_opendata(int modtype) {
	extern TSingSpec SingSpec;  FILE *g=SingSpec.g;
	int site, srvy, seasn, dsum=0, dsum2=0, srvy1, nseasns=1, N=SingSpec.N;
	if (SingSpec.PrmyPeriods>1) nseasns=SingSpec.PrmyPeriods;
	//if (SingSpec.NMethods>1) nseasns=T / SingSpec.NMethods;
    for (site=0; site<N; site++) {
        for (srvy=seasn=0; seasn<nseasns; seasn++) {
            for (srvy1=dsum=dsum2=0; srvy1<SingSpec.SecPeriods[seasn]; srvy++,srvy1++) {
                if (SingSpec.Data[site][srvy]>0 && !SingSpec.FalsePos &&
					(modtype<3 || modtype==7 || modtype==10)) SingSpec.Data[site][srvy]=1;
                if (SingSpec.Data[site][srvy]>MISSNG) dsum+=SingSpec.Data[site][srvy];
                else  dsum2++;
            }
            SingSpec.OpenData[site][seasn] = MISSNG;
            if (dsum2<SingSpec.SecPeriods[seasn]) SingSpec.OpenData[site][seasn] = (dsum>0);
        }
    }
    if (SingSpec.Verbose>1) {
        printf("opendata:\n"); fprintf(g,"opendata:\n");
        for (site=0; site<N; site++) { 
			printf("%d:",site+1); fprintf(g,"%d:",site+1);
            for (seasn=0; seasn<nseasns; seasn++) {
				printf(" %d",SingSpec.OpenData[site][seasn]);
				fprintf(g," %d",SingSpec.OpenData[site][seasn]);
			}
            printf("\n"); fprintf(g,"\n");
        }
    }
}
#ifdef noRPRES
int main(int argc, char **argv) {
	int i; i=presence(&argc, argv);
	printf("exiting main...(%d)\n",i);
	exit(0);
}
#endif
