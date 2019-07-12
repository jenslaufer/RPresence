//   SingAnal.h

#include <stdio.h>
#define VERSN "2.12.31"
#ifndef SingAnalH

#define SingAnalH
#define logitLnk 1
#define logitLnk2 8
#define loglogLnk 2
#define expLnk 3
#define sinLnk 4
#define sinLnk2 7
#define IDLnk 5
#define mlogitLnk 6
#define PSIBA_ODDS_RATIO 2
#define PSIAB_ODDS_RATIO 3
#define PSIBA_PSIBa 1
#define PSIB_PHI 0
#define MISSNG -1

struct TSingSpec {
    //  Stuff in PAO file...
    int N, T;              //  Number of sites, Number of surveys
	int **Data;            //  Input detection data
	double *det_hist_frq;  //  Number of sites each det. hist. represents
	int NSiteCov;        //  Number of site covariates
	char **CovNames;       //  vector[NSiteCov] of site-covariate names
	double **SiteCov;      //  Matrix[N,NSiteCov] of site covariates
	int NSampCov;        //  Number of sample(survey) covariates
	char **CovNames2;      //  vector[NSampCov] of sample-covariate names
	double ***SampCov;     //  Matrix[NSampCov,N,T] of sample covariates
	int PrmyPeriods, *SecPeriods;   //  Number of Primary periods, Vector of number of surveys in each primary period
	char **sitename, **surveyname;   // Vector[N] of sitenames, Vector[T] of surveynames
	int Nstates;                    //  Number of states(strata) in detection-history records
	int NMethods;                   //  Number of methods per survey
	int NMiss;                      //  Number of missing values in data
	int **OpenData;                   //  Pooled version of **Data (pooled by primary period)
	int mthreads;                   //  Number of threads for openmodmt.

	//  Stuff in design-matrix file...
	int NrowsDM[9], NParKK[9];  // Number of rows, cols in each input design matrix
	double ***DMat;             //  Design matrices[i, NrowsDM[i], NParKK[i]], i:6
	int ***DMat_ptr;            //  Design matrices pointers[i, NrowsDM[i], NParKK[i]]
	double *fixed;              //   Vector of fixed parameter values
	char ***Realname;           //   Matrices[i,NrowsDM[i]] of real parameter (eg., psi,gam,...) names
	char ***Betaname;           //   matrices[i,NParKK[i]] of beta parameter (eg., A1,A2,B1,B2,...) names

	int *maxx,
	     lmt,           //  max number of sites to print for each real parameter
		 rlmt,           //  upper limit for N in Royal model
		 ifn,            //  current number of times likelihood function has bee called
		 maxfn,          //  Max number of likelihood function calls
		 Verbose,        //  switch to print on console (0=not much printing, 1=some printing, 2=lots of printint)
		 novar,          //  switch to determine if variances are computed/printed
         Groups,         //  number of groups for pre-defined heterogeneity models
		  UseNeighbor,
          Alt_Param_Checked,
		  uncondpsi,
		  do_chisq,
		  nrandIV,
		  nphantom, **phantom,
		  Model,
		  LinkFn,
		  **LnkFn,
		  **neighbor_lst,
		  LikeNRSig,
		  BootNRSig,             // numerical covergence options
		  nthreads,
		  simulating;
    bool TSpecificP,
	      MissClass,
		  FalsePos,
		  prnt_derived,
		  optmiz_done,
		  UseAmoeba,
		  dups,
		  multiseason;

    double **choose,
	        *area,
			**BetaFixed,    //  fixed value for beta's (beta not fixed if betafixed>1e44)
		   **realParmEst,
		   ****OpenP,
		   *expval,
		   ***closedp,
		   ***theta,
		   ***psibar,
		   *neighbor_wgt,
		   nearzero,
		   VCeps,
		   *finalBetaEst,
		   finalLL,
		   **finalVCmat;
	char modname[256],
		 **neighbor;
    long seed;     FILE *g;
};
#endif

#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif
