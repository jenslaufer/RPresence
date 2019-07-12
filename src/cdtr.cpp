#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#define dabs(x) (x<0 ? -x : x)
double ndtr(double x) {
/*C        SUBROUTINE NDTR
C
C        PURPOSE
C           COMPUTES Y = P(X) = PROBABILITY THAT THE RANDOM VARIABLE  U,
C           DISTRIBUTED NORMALLY(0,1), IS LESS THAN OR EQUAL TO X.
C           F(X), THE ORDINATE OF THE NORMAL DENSITY AT X, IS ALSO
C           COMPUTED.
C
C        USAGE
C           CALL NDTR(X,P)
C
C        DESCRIPTION OF PARAMETERS
C           X--INPUT SCALAR FOR WHICH P(X) IS COMPUTED.
C           P--OUTPUT PROBABILITY.
C
C        REMARKS
C           MAXIMUM ERROR IS 0.0000007.
C
C        SUBROUTINES AND SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           BASED ON APPROXIMATIONS IN C. HASTINGS, APPROXIMATIONS FOR
C           DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON, N.J.,
C           1955.  SEE EQUATION 26.2.17, HANDBOOK OF MATHEMATICAL
C           FUNCTIONS, ABRAMOWITZ AND STEGUN, DOVER PUBLICATIONS, INC.,
C           NEW YORK.
C
   SUBROUTINE NDTR (X,P) */
    double p,t,d;
    t=1./(1.+.2316419*dabs(x));
    if (-x*x*0.5<-700) d=0.3989423*exp(-700.);
    else d=0.3989423*exp(-x*x*0.5);
    p=1.0-d*t*((((1.330274*t-1.821256)*t+1.781478)*t-0.3565638)*t+0.3193815);
    if (x<0.) p=1.-p;
    return(p);
}

double cdtr(double chisq, int idf, int *ier) {
/*C        SUBROUTINE CDTR
C        PURPOSE
C           COMPUTES P(CHISQ) = PROBABILITY THAT THE RANDOM VARIABLE U,
C           DISTRIBUTED ACCORDING TO THE CHI-SQUARE DISTRIBUTION WITH G
C           DEGREES OF FREEDOM, IS LESS THAN OR EQUAL TO CHISQ.  F(G,CHI
C           ORDINATE OF THE CHI-SQUARE DENSITY AT CHISQ, IS ALSO COMPUTE
C
C        USAGE
C           CALL CDTR(CHISQ,IDF,SIGCHI,IER)
C
C        DESCRIPTION OF PARAMETERS
C           CHISQ   - INPUT SCALAR FOR WHICH P(CHISQ) IS COMPUTED.
C           IDF     - NUMBER OF DEGREES OF FREEDOM OF THE CHI-SQUARE
C                     DISTRIBUTION.  IDF IS AN INTEGER.
C           SIGCHI  - OUTPUT PROBABILITY.
C           IER     - RESULTANT ERROR CODE WHERE
C               IER= 0 --- NO ERROR
C               IER=-1 --- AN INPUT PARAMETER IS INVALID.  X IS LESS
C                          THAN 0.0, OR IDF IS LESS THAN 0.0 OR GREATER
C                          THAN 2*10**(+5).  P AND D ARE SET TO -1.E25.
C               IER=-2 --- IDF = 0, P SET TO 0.
C               IER=+1 --- INVALID OUTPUT.  P IS LESS THAN ZERO OR
C                          GREATER THAN ONE, OR SERIES FOR T1 (SEE
C                          MATHEMATICAL DESCRIPTION) HAS FAILED TO
C                          CONVERGE.  P IS SET TO 1.E25.
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NDTR
C
C        METHOD
C
C     THIS FUNCTION WAS LIFTED FROM "NWAY" WRITTEN BY THE STANFORD
C     COMPUTATION CENTER. THE SCC BORROWED IT FROM JOHN MORRIS OF THE
C     COMPUTER INSTITUTE FOR SOCIAL SCIENCE RESEARCH AT MICHIGAN STATE
C     UNIVERSITY.
*/
    double anschi, ff, p, xi, pp, term, argnrm, gam, txx, sigchi=0;
    long i;

    *ier=0; anschi=1.;
//    printf("cdtr:chisq=%f, idf=%d)\n",chisq,idf);
    if (chisq<0. || idf<0){*ier=-1; sigchi=-1000000000000000000000000.; }
    else if (idf==0) { }
    else if (chisq+100==100) {}
    else {
        //    c     gt 60 d.f. or chi gt 100 - use a normal approximation;
        if (idf>60 || chisq>100){ argnrm=sqrt(2*chisq)-sqrt(2.*idf-1.); sigchi=ndtr(argnrm); }
        else {
            ff=idf-2;  p=0.5*ff;  xi=0.5*chisq; pp=p+2.0;  term=xi/pp;
           // printf("ff=%f p=%f xi=%f pp=%f term=%f\n",ff,p,xi,pp,term);
            for (i=1; i<=100; i++){
                anschi+=term;// printf("i=%d term=%f abs(term)=%f\n",i,term,fabs(term));
                if (fabs(term)<.000001) break;
                pp=pp+1.0; term=term*xi/pp;
              //  printf("%d: %f %f %f\n",i,anschi,term,pp);
            }
            gam=1.0; //printf("pp=%f, term=%f\n",pp,term);
            if (floor(idf/2.)*2==idf)
                for (i=1; i<=idf/2; i++) gam*=i; // factorial(idf/2)
            else {
                for (i=1; i<=idf; i+=2) gam*=i*0.5;
                gam*=1.7724588509;
            }
            txx=(p+1.0)*log(xi)- xi + log(anschi/gam);
            if (txx<=88.0)
                if (txx>-88.0) anschi=exp(txx);
                else  anschi=0.;
            else anschi=1.;
          //  printf("gam=%f, txx=%f\n",gam,txx);
            sigchi=anschi;
        }
        if (sigchi<0. || sigchi>1.) { *ier=1; sigchi=1000000000000000000000.; }
    }
    //   printf("cdtr:result=%f\n",*sigchi);
    return(sigchi);
}
